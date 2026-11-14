"""Boltz-2 structure + binding prediction for the mutant protein.

Runs the `boltz-api` CLI (v0.44.0, OAuth session in the system keyring) as a
subprocess, so this module stays importable on the stdlib interpreter like the
rest of the pipeline. No torch, no Modal.

WHY WE CROP TO THE KINASE DOMAIN. Boltz bills by crop size, and the price for
full-length P00519 (1130 aa) is 4x the price of the kinase domain alone:

    full-length 1130 aa + imatinib, 1 sample     $0.20
    kinase crop 242-512 + imatinib, 1 sample     $0.05

Cost is only half the argument. ABL1's N-terminal half carries the SH3/SH2
domains and long disordered linkers; folding them costs money and adds nothing
to a pocket measurement, while the extra chains give the sampler more ways to
produce a low-confidence model. The crop is the better prediction *and* the
cheaper one.

WHY 242-512 AND NOT THE STRUCTURE'S OWN RANGE. The obvious crop is the range
1IEP actually resolves, 225-498. That is a trap: the BCR-ABL DMS screen measures
positions 242-512, so cropping to 1IEP would silently drop the 14 positions
above 498. It would also save nothing -- see the $0.05 floor below. The crop is
therefore driven by the *data*, not by the reference structure.

THE $0.05 FLOOR. The published price table implies a cheaper 0-256 aa tier at
$0.025, but it is unreachable in practice: a 256-aa crop and even a 200-aa crop
both estimate at $0.05, presumably MSA-generation overhead. So there is no
cost argument for shaving the crop, and shaving it is exactly the change that
would corrupt the position coverage. Do not shave the crop.

BINDING METRICS ARE FREE. Structure-only and structure+affinity estimate at the
identical $0.05, so we always request `binding`. drug_screen.py wants those
metrics and they cost nothing.

FIVE SAMPLES, NOT ONE. Sample pricing is strongly sublinear and plateaus:

    1 sample  $0.05        5 samples  $0.10   <- best value, $0.020/unit
    3 samples $0.10       10 samples  $0.20

5 samples cost the same as 3 and reach the per-unit floor. One sample is a
single draw from a diffusion model with no spread to report; 5 give a usable
standard deviation, so `Structure.confidence_sd` is an honest number rather
than a missing one.

NUMBERING. The crop introduces a FIFTH frame on top of the four already in
play, and the Boltz API indexes residues 0-based within the submitted chain:

    UniProt P00519 position   p          (1-based, == DMS, == 1IEP, offset 0)
    crop-local 1-based        p - 241
    Boltz residue_index       p - 242    (0-based)

`_to_crop_index()` is the only place that conversion happens and it validates
the residue identity it lands on, raising rather than returning a plausible
wrong index. Every cached result records the crop bounds, so a confidence can
never be silently compared across two different crops.

EVERY CALL COSTS MONEY, so `predict_structure()` runs the free `estimate-cost`
endpoint first and refuses to submit if the estimate exceeds `max_usd`
(default $0.25, override with VN_BOLTZ_MAX_USD). Results are cached on disk
keyed by the full request, so re-running a variant never spends twice.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

MODEL = "boltz-2.1"
CACHE_DIR = Path("output/structures")
COMPOUNDS = Path("data/compounds/compounds.json")

# Kinase-domain crop, in UniProt P00519 numbering. Driven by the DMS position
# range (242-512), not by 1IEP's resolved range (225-498). See module docstring.
CROP_START = 242
CROP_END = 512

NUM_SAMPLES = 5           # same price as 3, reaches the per-unit floor
DEFAULT_MAX_USD = 0.25    # per call; a 5-sample kinase-crop run estimates $0.10
TIMEOUT = 1800            # a 271-aa complex finishes well inside this

# Boltz reports confidence on 0-1; the Structure contract is pLDDT on 0-100.
PLDDT_SCALE = 100.0


class BoltzUnavailable(RuntimeError):
    """boltz-api is missing, or its OAuth session is not usable."""


class BoltzTooExpensive(RuntimeError):
    """The free cost estimate came back above the caller's budget cap."""


@dataclass
class Structure:
    path: Path                              # predicted mutant structure (mmCIF/PDB)
    sequence: str
    confidence: float                       # mean pLDDT, 0-100
    pocket_residues: list[int] = field(default_factory=list)
    per_residue_confidence: dict[int, float] = field(default_factory=dict)
    source: str = "boltz"
    # --- added for Boltz-2; all optional so existing callers keep working ---
    binding_confidence: float | None = None
    ptm: float | None = None
    iptm: float | None = None
    confidence_sd: float | None = None      # spread across the 5 samples
    num_samples: int = 1
    crop: tuple[int, int] | None = None     # UniProt bounds of the submitted crop
    cost_usd: float | None = None
    prediction_id: str | None = None

    @property
    def reliable(self) -> bool:
        return self.confidence >= 70.0


# --------------------------------------------------------------------------
# availability
# --------------------------------------------------------------------------

def _cli() -> str:
    """Absolute path to boltz-api, or raise."""
    # Installed by its own installer to ~/.local/bin, which is not always on
    # PATH for a non-login shell, so look there explicitly before giving up.
    found = shutil.which("boltz-api") or shutil.which(
        "boltz-api", path=str(Path.home() / ".local/bin")
    )
    if not found:
        raise BoltzUnavailable(
            "boltz-api not found on PATH or in ~/.local/bin; "
            "install the CLI and run `boltz-api auth login`"
        )
    return found


def available() -> bool:
    """True when boltz-api exists and reports a usable auth state.

    NOTE: this deliberately does NOT look at MODAL_TOKEN_ID. An earlier stub
    gated this stage on the Modal credentials, but Boltz is a separate service
    with its own OAuth session; the Modal check would have reported the wrong
    answer in both directions.
    """
    try:
        cli = _cli()
    except BoltzUnavailable:
        return False
    try:
        out = subprocess.run(
            [cli, "auth", "status", "--format", "json"],
            capture_output=True, text=True, timeout=30,
        )
    except (OSError, subprocess.SubprocessError):
        return False
    if out.returncode != 0:
        return False
    try:
        return bool(json.loads(out.stdout).get("authenticated"))
    except (json.JSONDecodeError, AttributeError):
        return False


# --------------------------------------------------------------------------
# numbering
# --------------------------------------------------------------------------

def _to_crop_index(
    position: int, sequence: str, *, expect: str | None = None
) -> int:
    """UniProt 1-based position -> Boltz 0-based residue_index in the crop.

    Validates the residue it lands on when `expect` is given. Raises rather
        than returning a plausible wrong index -- a silent off-by-one here would
    place a pocket constraint on the wrong residue and quietly degrade every
    prediction.
    """
    if not CROP_START <= position <= CROP_END:
        raise ValueError(
            f"position {position} is outside the submitted crop "
            f"{CROP_START}-{CROP_END}; widen the crop rather than clamping"
        )
    idx = position - CROP_START
    if idx >= len(sequence):
        raise ValueError(
            f"position {position} maps to crop index {idx} but the crop is only "
            f"{len(sequence)} residues -- crop bounds and sequence disagree"
        )
    if expect is not None and sequence[idx] != expect:
        raise ValueError(
            f"numbering check failed: UniProt {position} should be {expect!r} "
            f"but crop index {idx} holds {sequence[idx]!r}. Do not adjust the "
            f"offset to make this pass -- find out which frame the caller used."
        )
    return idx


# Residue identities at six spread-out positions in P00519 numbering. These are
# what actually catch a wrongly-numbered reference: if a caller hands us a
# kinase-domain construct numbered from 1 (the "518-aa ABL1 reference" that a
# plan once proposed, against P00519's real 1130 aa), the crop still slices
# cleanly and still has the right *length* -- only the residue identities give
# it away. A length check alone sails straight past it.
CROP_ANCHORS = {244: "M", 271: "K", 286: "E", 381: "D", 393: "Y", 468: "V"}
MAX_ANCHOR_MISMATCHES = 2


def crop_sequence(full_sequence: str, *, strict: bool = True) -> str:
    """Slice the kinase domain out of a full-length P00519-numbered sequence.

    Validates residue identity at `CROP_ANCHORS`, tolerating up to
    `MAX_ANCHOR_MISMATCHES` because the caller may legitimately be passing a
    *mutant* sequence whose substitution lands on an anchor. Anchors therefore
    deliberately avoid the famous resistance positions (315, 250, 253, 396).
    """
    if len(full_sequence) < CROP_END:
        raise ValueError(
            f"sequence is {len(full_sequence)} aa, shorter than the crop end "
            f"{CROP_END}; cannot slice residues {CROP_START}-{CROP_END}"
        )
    bad = {
        p: (aa, full_sequence[p - 1])
        for p, aa in CROP_ANCHORS.items()
        if full_sequence[p - 1] != aa
    }
    if strict and len(bad) > MAX_ANCHOR_MISMATCHES:
        raise ValueError(
            f"{len(bad)}/{len(CROP_ANCHORS)} anchor residues do not match "
            f"P00519 numbering: "
            + ", ".join(f"{p} expected {e!r} got {g!r}" for p, (e, g) in
                        sorted(bad.items()))
            + ". This sequence is almost certainly in a different numbering "
              "frame (a kinase-only construct numbered from 1?). Do NOT add an "
              "offset to make this pass -- establish which frame it is first."
        )
    return full_sequence[CROP_START - 1:CROP_END]


# --------------------------------------------------------------------------
# request assembly
# --------------------------------------------------------------------------

def compound_smiles(name: str) -> str:
    """SMILES for a compound in data/compounds/compounds.json."""
    for entry in json.loads(COMPOUNDS.read_text()):
        if entry.get("name", "").lower() == name.lower():
            return entry["smiles"]
    raise KeyError(f"no SMILES for compound {name!r} in {COMPOUNDS}")


def build_input(
    sequence: str,
    *,
    compound: str | None = None,
    pocket_hint: list[int] | None = None,
    num_samples: int = NUM_SAMPLES,
) -> dict[str, Any]:
    """Assemble the prediction body for an already-cropped sequence."""
    entities: list[dict[str, Any]] = [
        {"type": "protein", "value": sequence, "chain_ids": ["A"]}
    ]
    body: dict[str, Any] = {"entities": entities, "num_samples": num_samples}

    if compound:
        entities.append({
            "type": "ligand_smiles",
            "value": compound_smiles(compound),
            "chain_ids": ["B"],
        })
        # Free, and drug_screen.py consumes it. Always on when there's a ligand.
        body["binding"] = {
            "type": "ligand_protein_binding", "binder_chain_id": "B",
        }

        # A pocket constraint tells the sampler where the ATP site is. That is
        # wild-type structural annotation from UniProt, identical for every
        # variant, so it does not leak anything about THIS variant's
        # resistance -- unlike literature, which evidence.py is kept away from
        # prediction for exactly that reason. force=False biases rather than
        # enforces, so a variant that genuinely reshapes the pocket can still
        # move the ligand.
        contacts = [
            _to_crop_index(p, sequence)
            for p in sorted(pocket_hint or [])
            if CROP_START <= p <= CROP_END
        ]
        if contacts:
            body["constraints"] = [{
                "type": "pocket",
                "binder_chain_id": "B",
                "contact_residues": {"A": contacts},
                "max_distance_angstrom": 6.0,
                "force": False,
            }]
    return body


# --------------------------------------------------------------------------
# cost preflight
# --------------------------------------------------------------------------

def estimate_cost(body: dict[str, Any], *, model: str = MODEL) -> dict[str, Any]:
    """Free cost estimate. Creates no resource and consumes no GPU."""
    out = subprocess.run(
        [_cli(), "predictions:structure-and-binding", "estimate-cost",
         "--input", json.dumps(body), "--model", model, "--format", "json"],
        capture_output=True, text=True, timeout=120,
    )
    if out.returncode != 0:
        raise BoltzUnavailable(
            f"estimate-cost failed: {(out.stderr or out.stdout).strip()[:300]}"
        )
    return json.loads(out.stdout)


# --------------------------------------------------------------------------
# result parsing
# --------------------------------------------------------------------------

def _mean_plddt_per_residue(cif_path: Path, sequence: str) -> dict[int, float]:
    """Per-residue confidence from the mmCIF B-factor column, keyed by UniProt.

    Boltz writes per-atom pLDDT into B_iso_or_equiv. We average over each
    residue's atoms and translate label_seq_id (1-based within the submitted
    chain) back to UniProt numbering so callers never see crop-local indices.
    """
    cols: list[str] = []
    sums: dict[int, list[float]] = {}
    in_loop = False
    for raw in cif_path.read_text().splitlines():
        line = raw.strip()
        if line.startswith("_atom_site."):
            cols.append(line.split(".", 1)[1])
            in_loop = True
            continue
        if in_loop:
            if not line or line.startswith(("#", "loop_", "_", "data_")):
                break
            parts = line.split()
            if len(parts) != len(cols):
                continue
            row = dict(zip(cols, parts))
            try:
                seq_id = int(row["label_seq_id"])
                b = float(row["B_iso_or_equiv"])
            except (KeyError, ValueError):
                continue
            if row.get("label_asym_id", "A") != "A":
                continue  # protein chain only; ligand has no residue numbering
            sums.setdefault(seq_id, []).append(b)

    scale = PLDDT_SCALE if sums and max(
        max(v) for v in sums.values()
    ) <= 1.0 else 1.0
    out: dict[int, float] = {}
    for seq_id, vals in sums.items():
        if 1 <= seq_id <= len(sequence):
            out[seq_id + CROP_START - 1] = sum(vals) / len(vals) * scale
    return out


def _confidence(run_json: dict[str, Any]) -> tuple[float, float | None, dict]:
    """Mean pLDDT on 0-100, its spread across samples, and the best sample."""
    samples = (
        run_json.get("results", {}).get("samples")
        or run_json.get("samples")
        or []
    )
    scored: list[tuple[float, dict]] = []
    for s in samples:
        c = s.get("confidence", s)
        val = c.get("complex_plddt", c.get("structure_confidence"))
        if val is None:
            continue
        val = float(val)
        scored.append(((val * PLDDT_SCALE if val <= 1.0 else val), s))
    if not scored:
        raise BoltzUnavailable(
            "prediction returned no per-sample confidence; cannot report pLDDT"
        )
    vals = [v for v, _ in scored]
    mean = sum(vals) / len(vals)
    sd = None
    if len(vals) > 1:
        var = sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)
        sd = var ** 0.5
    best = max(scored, key=lambda t: t[0])[1]
    return mean, sd, best


# --------------------------------------------------------------------------
# the stage
# --------------------------------------------------------------------------

def _cache_key(body: dict[str, Any], model: str) -> str:
    blob = json.dumps(
        {"body": body, "model": model, "crop": [CROP_START, CROP_END]},
        sort_keys=True,
    )
    return hashlib.sha256(blob.encode()).hexdigest()[:16]


def predict_structure(
    sequence: str,
    *,
    name: str = "mutant",
    out_dir: Path | str = CACHE_DIR,
    pocket_hint: list[int] | None = None,
    compound: str | None = "imatinib",
    num_samples: int = NUM_SAMPLES,
    max_usd: float | None = None,
    model: str = MODEL,
    already_cropped: bool = False,
    dry_run: bool = False,
) -> Structure:
    """Predict the mutant structure, with a mandatory free cost preflight.

    `sequence` is full-length in P00519 numbering unless `already_cropped`.
    Set `dry_run=True` to stop after the estimate and spend nothing; the
    returned Structure then carries only `cost_usd`.
    """
    if not available():
        raise BoltzUnavailable(
            "boltz-api is not authenticated; run `boltz-api auth login`"
        )

    cropped = sequence if already_cropped else crop_sequence(sequence)
    body = build_input(
        cropped, compound=compound, pocket_hint=pocket_hint,
        num_samples=num_samples,
    )

    out_dir = Path(out_dir)
    run_dir = out_dir / f"{name}_{_cache_key(body, model)}"
    cached = run_dir / "run.json"

    cap = max_usd if max_usd is not None else float(
        os.getenv("VN_BOLTZ_MAX_USD", DEFAULT_MAX_USD)
    )

    # Preflight runs even on a cache hit -- it is free, and it keeps the
    # reported cost honest if the cached run was made under a different price.
    est = estimate_cost(body, model=model)
    cost = float(est["estimated_cost_usd"])
    if cost > cap:
        raise BoltzTooExpensive(
            f"estimate ${cost:.4f} exceeds the ${cap:.2f} cap for {name}. "
            f"Raise max_usd / VN_BOLTZ_MAX_USD deliberately, or shrink the "
            f"request ({est.get('breakdown')})."
        )
    if dry_run:
        return Structure(
            path=run_dir, sequence=cropped, confidence=0.0,
            pocket_residues=list(pocket_hint or []), num_samples=num_samples,
            crop=(CROP_START, CROP_END), cost_usd=cost, source="boltz-estimate",
        )

    if not cached.exists():
        run_dir.parent.mkdir(parents=True, exist_ok=True)
        # @json:// keeps a 271-aa sequence plus SMILES off the argv length limit.
        payload = run_dir.parent / f".{run_dir.name}.input.json"
        payload.write_text(json.dumps(body))
        try:
            proc = subprocess.run(
                [_cli(), "predictions:structure-and-binding", "run",
                 "--input", f"@json://{payload}",
                 "--model", model,
                 "--run-dir", str(run_dir),
                 "--format", "json"],
                capture_output=True, text=True, timeout=TIMEOUT,
            )
        finally:
            payload.unlink(missing_ok=True)
        if proc.returncode != 0:
            raise BoltzUnavailable(
                f"boltz run failed for {name}: "
                f"{(proc.stderr or proc.stdout).strip()[-400:]}"
            )

    if not cached.exists():
        raise BoltzUnavailable(f"run completed but {cached} is missing")
    run_json = json.loads(cached.read_text())

    mean, sd, best = _confidence(run_json)
    files = sorted((run_dir / "outputs" / "files").glob("*.cif"))
    per_residue = _mean_plddt_per_residue(files[0], cropped) if files else {}
    binding = (best.get("binding_metrics") or run_json.get("binding_metrics")
               or {})

    conf = best.get("confidence", best)
    return Structure(
        path=files[0] if files else run_dir,
        sequence=cropped,
        confidence=mean,
        pocket_residues=list(pocket_hint or []),
        per_residue_confidence=per_residue,
        binding_confidence=binding.get("binding_confidence"),
        ptm=conf.get("ptm"),
        iptm=conf.get("iptm"),
        confidence_sd=sd,
        num_samples=num_samples,
        crop=(CROP_START, CROP_END),
        cost_usd=cost,
        prediction_id=run_json.get("id"),
    )
