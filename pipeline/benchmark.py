"""Read completed Boltz benchmark controls without submitting jobs or training models."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path("output/structures/boltz_benchmark")


def lookup(variant: str) -> dict[str, Any] | None:
    """Return one fixed-panel comparison with source hashes and wild-type delta."""
    if not (ROOT / "summary.json").exists() or not (ROOT / "manifest.json").exists():
        return None
    manifest = json.loads((ROOT / "manifest.json").read_text(encoding="utf-8"))
    if manifest["reference_sha256"] != hashlib.sha256(Path("data/cancer-exp/ABL1_NM_005157.6_DMS_reference_aa.txt").read_bytes()).hexdigest():
        raise ValueError("Cached Boltz benchmark reference changed; rerun benchmark")
    if manifest["dataset_sha256"] != hashlib.sha256(Path("data/datasets/abl_dms_netgrowth_1200nM.csv").read_bytes()).hexdigest():
        raise ValueError("Cached Boltz benchmark assay dataset changed; rerun benchmark")
    summary = json.loads((ROOT / "summary.json").read_text(encoding="utf-8"))
    if summary.get("manifest_sha256") != hashlib.sha256((ROOT / "manifest.json").read_bytes()).hexdigest():
        raise ValueError("Cached Boltz summary does not match its execution manifest")
    row = next((r for r in summary["rows"] if r["variant"] == variant), None)
    if row is None or "error" in row:
        return None
    geometry = None
    if (ROOT / "structure_comparison.json").exists():
        structural = json.loads((ROOT / "structure_comparison.json").read_text(encoding="utf-8"))
        geometry = next((r for r in structural["rows"] if r["variant"] == variant and r.get("prediction_id") == row["prediction_id"]), None)
    return {"geometry": geometry, "variant": variant, "model": manifest["model"], "crop": manifest["crop"],
            "num_samples": manifest["num_samples"], "compound": manifest["compound"],
            "constraints": manifest["constraints"], "binding_confidence": row["binding_confidence"],
            "delta_binding_confidence_vs_wt": row["delta_binding_confidence_vs_wt"],
            "confidence": row["confidence"], "confidence_sd": row["confidence_sd"],
            "prediction_id": row["prediction_id"], "caveats": summary["caveats"],
            "source": "cached fixed-panel contrast, not drug-response inference"}
