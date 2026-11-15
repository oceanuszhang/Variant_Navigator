"""Local VariantNavigator demo server. Run: python ui/server.py --port 8000."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import sys
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from threading import Lock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
WORK_LOCK = Lock()


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def ptarmigan_control() -> dict:
    base = ROOT / "output/ptarmigan"
    return {
        "scope": (
            "Cached Ptarmigan-1 wild-type target control. It ranks compounds "
            "for wild-type ABL1 and does not score the submitted mutant."
        ),
        "hits": _read_json(base / "abl1_wt_hits.json"),
        "screen": _read_json(base / "abl1_wt_screen.json"),
    }


def submit_prediction(body: dict) -> dict:
    from pipeline.model_store import predict_saved
    return predict_saved(
        str(body.get("variant", "")),
        str(body.get("sequence", "")),
        str(body.get("model", "xgboost")),
    )


def _mutant_sequence(body: dict) -> tuple[str, str]:
    from pipeline.model_store import reference_sequence, validate_input
    variant, position, _ref, alt = validate_input(
        str(body.get("variant", "")), str(body.get("sequence", ""))
    )
    sequence = reference_sequence()
    return variant, sequence[: position - 1] + alt + sequence[position:]


def boltz_estimate(body: dict) -> dict:
    from pipeline.boltz_structure import predict_structure
    variant, sequence = _mutant_sequence(body)
    with WORK_LOCK:
        result = predict_structure(
            sequence,
            name=variant.lower(),
            compound="imatinib",
            pocket_hint=[286, 315, 381],
            max_usd=0.25,
            dry_run=True,
        )
    return {
        "variant": variant,
        "estimated_cost_usd": result.cost_usd,
        "model": "boltz-2.1",
        "crop": result.crop,
        "num_samples": result.num_samples,
        "submitted": False,
        "message": "Free estimate only; no paid prediction was submitted.",
    }


def boltz_run(body: dict) -> dict:
    from pipeline.boltz_structure import predict_structure
    if body.get("confirm_paid_run") is not True:
        raise ValueError("Paid Boltz run requires explicit confirmation.")
    max_usd = float(body.get("max_usd", 0.25))
    if not 0 < max_usd <= 0.25:
        raise ValueError("max_usd must be greater than $0 and no more than $0.25.")
    variant, sequence = _mutant_sequence(body)
    with WORK_LOCK:
        result = predict_structure(
            sequence,
            name=variant.lower(),
            compound="imatinib",
            pocket_hint=[286, 315, 381],
            max_usd=max_usd,
        )
    return {
        "variant": variant,
        "model": "boltz-2.1",
        "structure_path": str(result.path),
        "confidence": result.confidence,
        "binding_confidence": result.binding_confidence,
        "cost_usd": result.cost_usd,
        "crop": result.crop,
        "num_samples": result.num_samples,
        "prediction_id": result.prediction_id,
    }


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(ROOT / "ui"), **kwargs)

    def respond(self, status: int, data: dict) -> None:
        payload = json.dumps(data, allow_nan=False).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def do_GET(self) -> None:
        from pipeline import model_store
        try:
            if self.path == "/api/reference":
                self.respond(200, {"sequence": model_store.reference_sequence()})
            elif self.path == "/api/validation":
                self.respond(200, model_store.validation_summary())
            elif self.path == "/api/ptarmigan/abl1":
                self.respond(200, ptarmigan_control())
            elif self.path.startswith("/api/"):
                self.respond(404, {"error": "Unknown API route."})
            else:
                super().do_GET()
        except FileNotFoundError as exc:
            self.respond(404, {"error": str(exc)})
        except Exception as exc:
            self.log_error("GET failed: %s", exc)
            self.respond(503, {"error": str(exc)})

    def do_POST(self) -> None:
        origin = self.headers.get("Origin")
        if origin and origin != f"http://{self.headers.get('Host')}":
            self.respond(403, {"error": "Cross-origin submissions are not allowed."})
            return
        try:
            size = int(self.headers.get("Content-Length", "0"))
            if not 0 < size <= 65536:
                raise ValueError("Invalid request size.")
            body = json.loads(self.rfile.read(size))
            if not isinstance(body, dict):
                raise ValueError("JSON object required.")
            routes = {
                "/api/predict": submit_prediction,
                "/api/boltz/estimate": boltz_estimate,
                "/api/boltz/run": boltz_run,
            }
            action = routes.get(self.path)
            if action is None:
                self.respond(404, {"error": "Unknown API route."})
                return
            self.respond(200, action(body))
        except (ValueError, TypeError, json.JSONDecodeError) as exc:
            self.respond(400, {"error": str(exc)})
        except Exception as exc:
            self.log_error("POST failed: %s", exc)
            self.respond(503, {"error": str(exc)})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args()
    os.chdir(ROOT)
    print(f"VariantNavigator: http://127.0.0.1:{args.port}", flush=True)
    ThreadingHTTPServer(("127.0.0.1", args.port), Handler).serve_forever()
