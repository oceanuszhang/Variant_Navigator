"""Extract frozen ESM-2 mutation-site delta embeddings on Modal."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import modal

MODEL_NAME = "facebook/esm2_t12_35M_UR50D"

image = modal.Image.debian_slim(python_version="3.11").pip_install(
    "numpy>=1.26,<3",
    "safetensors>=0.4",
    "torch>=2.4,<3",
    "transformers>=4.45,<5",
)
app = modal.App("variant-navigator-esm2-35m")


def reconstruct_sequence(reference: str, variant: dict[str, object]) -> str:
    position = int(variant["position"])
    sequence_index = position - 1
    return (
        reference[:sequence_index]
        + str(variant["alt_aa"])
        + reference[sequence_index + 1 :]
    )


def artifact_stem(row_index: int, mutant: str) -> str:
    safe_mutant = re.sub(r"[^A-Za-z0-9_.-]+", "_", mutant)
    return f"{row_index:04d}_{safe_mutant}"


def read_reference(path: str | Path) -> str:
    sequence = "".join(
        line.strip()
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.startswith(">")
    ).upper()
    if not sequence or not sequence.isalpha():
        raise ValueError(f"invalid protein sequence in {path}")
    return sequence


def read_variants(path: str | Path, reference: str) -> list[dict[str, object]]:
    variants: list[dict[str, object]] = []
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "mutant",
            "protein_start",
            "ref_aa",
            "alt_aa",
            "netgrowth_mean_corrected",
            "netgr_corr_rep1.high",
            "netgr_corr_rep2.high",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"CSV is missing columns: {sorted(missing)}")

        for row_number, row in enumerate(reader, start=2):
            position = int(row["protein_start"])
            ref_aa = row["ref_aa"].strip().upper()
            alt_aa = row["alt_aa"].strip().upper()
            if not 1 <= position <= len(reference):
                raise ValueError(f"row {row_number}: position {position} is outside reference")
            observed = reference[position - 1]
            if observed != ref_aa:
                raise ValueError(
                    f"row {row_number}: {row['mutant']} expects {ref_aa} at {position}, "
                    f"but reference has {observed}"
                )
            if len(ref_aa) != 1 or len(alt_aa) != 1 or ref_aa == alt_aa:
                raise ValueError(f"row {row_number}: invalid substitution {row['mutant']}")
            variants.append(
                {
                    "mutant": row["mutant"],
                    "position": position,
                    "ref_aa": ref_aa,
                    "alt_aa": alt_aa,
                    "target": float(row["netgrowth_mean_corrected"]),
                    "replicate_1": float(row["netgr_corr_rep1.high"]),
                    "replicate_2": float(row["netgr_corr_rep2.high"]),
                }
            )
    if not variants:
        raise ValueError(f"no variants found in {path}")
    return variants


@app.function(image=image, gpu="T4", timeout=3600)
def embed_variants(
    reference: str,
    variants: list[dict[str, object]],
    batch_size: int = 32,
):
    import numpy as np
    import torch
    from transformers import AutoModel, AutoTokenizer

    if batch_size < 1:
        raise ValueError("batch_size must be positive")

    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = AutoModel.from_pretrained(MODEL_NAME).to("cuda").eval()

    reference_tokens = tokenizer(reference, return_tensors="pt").to("cuda")
    with torch.inference_mode():
        reference_hidden = model(**reference_tokens).last_hidden_state[0]

    embeddings = np.empty((len(variants), model.config.hidden_size), dtype=np.float32)

    for start in range(0, len(variants), batch_size):
        batch = variants[start : start + batch_size]
        mutant_sequences: list[str] = []
        token_positions: list[int] = []

        for variant in batch:
            position = int(variant["position"])
            mutant_sequences.append(reconstruct_sequence(reference, variant))
            # ESM prepends a BOS token, so the one-based residue position is
            # also its token index in last_hidden_state.
            token_positions.append(position)

        tokens = tokenizer(mutant_sequences, return_tensors="pt", padding=True).to("cuda")
        with torch.inference_mode():
            mutant_hidden = model(**tokens).last_hidden_state

        row_indices = torch.arange(len(batch), device="cuda")
        position_indices = torch.tensor(token_positions, device="cuda")
        mutant_vectors = mutant_hidden[row_indices, position_indices]
        reference_vectors = reference_hidden[position_indices]
        delta = (mutant_vectors - reference_vectors).float().cpu().numpy()
        embeddings[start : start + len(batch)] = delta

    return embeddings


@app.local_entrypoint()
def main(
    sequence_path: str,
    csv_path: str,
    output_path: str,
    batch_size: int = 32,
):
    import numpy as np

    reference = read_reference(sequence_path)
    variants = read_variants(csv_path, reference)
    print(f"Validated {len(variants)} variants against a {len(reference)}-aa reference")

    embeddings = embed_variants.remote(reference, variants, batch_size)
    if embeddings.shape != (len(variants), 480):
        raise RuntimeError(f"unexpected embedding shape: {embeddings.shape}")

    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    sequences = [reconstruct_sequence(reference, variant) for variant in variants]
    np.savez_compressed(
        destination,
        embeddings=embeddings,
        sequence=np.asarray(sequences),
        mutant=np.asarray([v["mutant"] for v in variants]),
        protein_start=np.asarray([v["position"] for v in variants], dtype=np.int32),
        ref_aa=np.asarray([v["ref_aa"] for v in variants]),
        alt_aa=np.asarray([v["alt_aa"] for v in variants]),
        target=np.asarray([v["target"] for v in variants], dtype=np.float32),
        replicate_1=np.asarray([v["replicate_1"] for v in variants], dtype=np.float32),
        replicate_2=np.asarray([v["replicate_2"] for v in variants], dtype=np.float32),
        model=np.asarray(MODEL_NAME),
        feature=np.asarray("mutant_minus_reference_embedding_at_mutated_residue"),
    )
    print(f"Saved {embeddings.shape} embeddings to {destination}")

    fasta_path = destination.parent / "variant_sequences.fasta"
    fasta_lines: list[str] = []
    for row_index, (variant, sequence) in enumerate(zip(variants, sequences, strict=True)):
        fasta_lines.extend(
            [
                (
                    f">{artifact_stem(row_index, str(variant['mutant']))} "
                    f"position={variant['position']} ref={variant['ref_aa']} "
                    f"alt={variant['alt_aa']}"
                ),
                sequence,
            ]
        )
    fasta_path.write_text("\n".join(fasta_lines) + "\n")

    embedding_dir = destination.parent / "variant_embeddings"
    embedding_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = destination.parent / "variant_manifest.csv"
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "row_index",
                "mutant",
                "protein_start",
                "ref_aa",
                "alt_aa",
                "sequence_length",
                "embedding_file",
            ],
        )
        writer.writeheader()
        for row_index, (variant, sequence, embedding) in enumerate(
            zip(variants, sequences, embeddings, strict=True)
        ):
            filename = artifact_stem(row_index, str(variant["mutant"])) + ".npy"
            relative_embedding_path = Path("variant_embeddings") / filename
            np.save(destination.parent / relative_embedding_path, embedding, allow_pickle=False)
            writer.writerow(
                {
                    "row_index": row_index,
                    "mutant": variant["mutant"],
                    "protein_start": variant["position"],
                    "ref_aa": variant["ref_aa"],
                    "alt_aa": variant["alt_aa"],
                    "sequence_length": len(sequence),
                    "embedding_file": relative_embedding_path.as_posix(),
                }
            )

    print(f"Saved {len(sequences)} full mutant sequences to {fasta_path}")
    print(f"Saved {len(embeddings)} per-variant embedding files to {embedding_dir}")
    print(f"Saved artifact manifest to {manifest_path}")
