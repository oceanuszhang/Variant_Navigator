"""Describe predicted mutant/WT geometry without treating it as response accuracy."""
from __future__ import annotations
import json
from pathlib import Path
import shlex
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

ROOT = Path("output/structures/boltz_benchmark")


def atoms(path: Path) -> list[dict]:
    """Read the downloaded Boltz atom-site loop (single-line atom records)."""
    columns, result = [], []
    for raw in path.read_text().splitlines():
        if raw.startswith("_atom_site."):
            columns.append(raw.split(".", 1)[1])
            continue
        if not columns:
            continue
        if raw.startswith(("#", "loop_", "_", "data_")):
            break
        if not raw.strip():
            continue
        tokens = shlex.split(raw)
        if len(tokens) != len(columns):
            raise ValueError("Unsupported CIF atom-site layout")
        row = dict(zip(columns, tokens))
        result.append({"chain": row["label_asym_id"], "name": row["label_atom_id"],
                       "element": row["type_symbol"], "residue": row["label_comp_id"], "position": int(row["label_seq_id"])+241 if row["label_seq_id"] not in (".", "?") else None,
                       "xyz": [float(row[k]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")]})
    return result


def main() -> None:
    import numpy as np
    from pipeline.mutant import THREE_TO_ONE, apply_variant
    study_reference = Path("data/cancer-exp/ABL1_NM_005157.6_DMS_reference_aa.txt").read_text().strip()
    results = json.loads((ROOT / "results.json").read_text())
    reference = atoms(Path(results["WILDTYPE"]["path"]))
    ref_ca = {a["position"]: a["xyz"] for a in reference if a["chain"] == "A" and a["name"] == "CA"}
    rows = []
    for variant, run in results.items():
        if "error" in run:
            continue
        current = atoms(Path(run["path"]))
        expected = study_reference if variant == "WILDTYPE" else apply_variant(study_reference, variant).sequence
        if any(THREE_TO_ONE.get(a["residue"]) != expected[a["position"]-1] for a in current if a["chain"] == "A" and a["name"] == "CA"):
            raise ValueError(f"Downloaded structure residue identities disagree for {variant}")
        ca = {a["position"]: a["xyz"] for a in current if a["chain"] == "A" and a["name"] == "CA"}
        positions = sorted(set(ca) & set(ref_ca))
        if len(positions) != 271:
            raise ValueError("Expected 271 aligned crop residues")
        x = np.array([ca[p] for p in positions]); y = np.array([ref_ca[p] for p in positions])
        x -= x.mean(0); y -= y.mean(0)
        u, _, vt = np.linalg.svd(x.T @ y)
        correction = np.eye(3); correction[-1,-1] = np.linalg.det(u @ vt)
        deviations = np.linalg.norm(x @ (u @ correction @ vt)-y, axis=1)
        rmsd = float(np.sqrt(np.mean(deviations**2)))
        if variant == "WILDTYPE" and rmsd > 1e-8:
            raise AssertionError("WT self-alignment failed")
        distances = {}
        lig = np.array([a["xyz"] for a in current if a["chain"] == "B" and a["element"] != "H"])
        for pos in (315, 396):
            site = np.array([a["xyz"] for a in current if a["chain"] == "A" and a["position"] == pos and a["element"] != "H"])
            if not len(lig) or not len(site):
                raise ValueError("Missing ligand or site atoms")
            distances[str(pos)] = float(np.linalg.norm(site[:,None,:]-lig[None,:,:], axis=2).min())
        rows.append({"variant": variant, "prediction_id": run["prediction_id"], "matched_ca_atoms": len(positions), "ca_rmsd_to_predicted_wt_angstrom": rmsd,
                     "site_min_heavy_atom_distance_to_ligand_angstrom": distances})
    report = {"rows": rows, "method": "best-sample C-alpha Kabsch alignment; fixed crop numbering; all 271 residues; no trimming",
              "caveat": "Predicted-to-predicted geometry, not experimental accuracy. Minimum contact distance and backbone movement do not measure binding affinity or efficacy."}
    (ROOT / "structure_comparison.json").write_text(json.dumps(report, indent=2),encoding="utf-8")
    path = ROOT / "report.md"
    text = path.read_text(encoding="utf-8").split("\n## Predicted geometry relative to predicted wild type")[0]
    text += "\n## Predicted geometry relative to predicted wild type\n\n| Variant | C-alpha RMSD (A) | Site 315 to ligand (A) | Site 396 to ligand (A) |\n| --- | --- | --- | --- |\n"
    for row in rows:
        d = row["site_min_heavy_atom_distance_to_ligand_angstrom"]
        text += f"| {row['variant']} | {row['ca_rmsd_to_predicted_wt_angstrom']:.3f} | {d['315']:.3f} | {d['396']:.3f} |\n"
    text += "\n" + report["caveat"] + "\n"
    path.write_text(text, encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
