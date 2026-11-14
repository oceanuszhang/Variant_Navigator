"""Compare a Boltz ABL1 crop against experimental 1IEP using residue-matched C-alpha atoms.

This is a structural agreement check, not drug-response accuracy or a held-out
benchmark. The experimental structure may have been in the model training data.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shlex
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def main() -> None:
    import numpy as np
    from pipeline.pocket import load_structure
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path)
    args = ap.parse_args()
    cif = next((args.run_dir / "outputs/files/prediction").glob("*_predicted.cif"))
    cols, predicted = [], {}
    for raw in cif.read_text().splitlines():
        if raw.startswith("_atom_site."):
            cols.append(raw.split(".", 1)[1])
            continue
        if not cols:
            continue
        if raw.startswith(("#", "loop_", "_", "data_")):
            break
        if not raw.strip():
            continue
        tokens = shlex.split(raw)
        if len(tokens) != len(cols):
            raise ValueError("Unsupported atom-site row layout")
        row = dict(zip(cols, tokens))
        if row["label_asym_id"] == "A" and row["label_atom_id"] == "CA":
            pos = int(row["label_seq_id"]) + 241
            predicted[pos] = (row["label_comp_id"], [float(row[k]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")])
    crystal = {a.residue_number: a for a in load_structure("1IEP").protein if a.chain == "A" and a.name == "CA"}
    positions = sorted(set(predicted) & set(crystal))
    mismatches = [p for p in positions if predicted[p][0] != crystal[p].residue_name]
    excluded = [{"position": p, "predicted_residue": predicted[p][0], "reference_residue": crystal[p].residue_name} for p in mismatches]
    positions = [p for p in positions if p not in mismatches]
    if len(positions) < 3:
        raise ValueError("Insufficient matched C-alpha atoms")
    mobile = np.array([predicted[p][1] for p in positions])
    fixed = np.array([[crystal[p].x, crystal[p].y, crystal[p].z] for p in positions])
    x, y = mobile-mobile.mean(0), fixed-fixed.mean(0)
    u, _, vt = np.linalg.svd(x.T @ y)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    distances = np.linalg.norm(x @ rotation-y, axis=1)
    report = {"reference_pdb": "1IEP", "reference_chain": "A", "predicted_chain": "A",
              "matched_ca_atoms": len(positions), "sequence_mismatches_excluded": excluded, "predicted_ca_atoms": len(predicted),
              "matched_position_range": [positions[0], positions[-1]],
              "ca_rmsd_angstrom": float(np.sqrt(np.mean(distances**2))),
              "median_ca_deviation_angstrom": float(np.median(distances)),
              "method": "residue identity verified; Kabsch rigid alignment; no outlier trimming",
              "caveat": "One known wild-type complex; possible training overlap. Not ligand-pose or drug-response accuracy."}
    (args.run_dir / "structural_validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
