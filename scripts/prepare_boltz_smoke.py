"""Prepare a small ABL1 WT/imatinib input and stable key without submitting."""
from pathlib import Path
import json
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline.boltz_structure import build_input, crop_sequence, idempotency_key


def main() -> None:
    reference = Path("data/cancer-exp/ABL1_NM_005157.6_DMS_reference_aa.txt").read_text().strip()
    body = build_input(crop_sequence(reference), compound="imatinib", num_samples=1)
    out = Path("output/structures/boltz_smoke")
    out.mkdir(parents=True, exist_ok=True)
    (out / "input.json").write_text(json.dumps(body, sort_keys=True, indent=2), encoding="utf-8")
    key = idempotency_key(body)
    (out / "idempotency-key.txt").write_text(key, encoding="utf-8")
    print(f"Prepared {out / 'input.json'}; idempotency key {key}; no job submitted")


if __name__ == "__main__":
    main()
