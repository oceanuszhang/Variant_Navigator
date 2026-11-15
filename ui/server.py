"""Local UI and pipeline API. Run: python ui/server.py --port 8000."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import sys
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from threading import Lock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
INFERENCE_LOCK = Lock()


def validation() -> dict:
    path = ROOT / "output/model/validation_predictions.json"
    raw = path.read_bytes()
    data = json.loads(raw)
    return {**data, "artifact_sha256": hashlib.sha256(raw).hexdigest(),
            "provenance_limit": "Saved artifact; training membership and sequence hash may be absent."}


def submit(variant: str, supplied_sequence: str | None = None) -> dict:
    from pipeline import predict
    match = re.fullmatch(r"([ACDEFGHIKLMNPQRSTVWY])(\d{3})([ACDEFGHIKLMNPQRSTVWY])", variant)
    if not match:
        raise ValueError("Enter one substitution, e.g. T315I, at ABL1 position 242–512.")
    ref, position, alt = match.groups()
    sequence = predict.REFERENCE.read_text().strip()
    if supplied_sequence is not None:
        lines = supplied_sequence.strip().splitlines()
        if lines and lines[0].startswith('>'):
            lines = lines[1:]
        pasted = ''.join(''.join(lines).split()).upper()
        if not pasted or not re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+', pasted):
            raise ValueError('Paste one amino-acid sequence using the 20 standard residue letters; a single FASTA header is allowed.')
        if pasted != sequence:
            raise ValueError(f'This model supports only the exact {len(sequence)}-aa ABL1 study reference. Use Load ABL1 reference; other proteins and sequence backgrounds are not supported.')
    index = int(position) - 1
    if not 242 <= int(position) <= 512:
        raise ValueError("Response model supports ABL1 positions 242–512 only.")
    if ref == alt:
        raise ValueError("Choose a substitution that changes the amino acid.")
    if sequence[index] != ref:
        raise ValueError(f"ABL1 position {position} is {sequence[index]}, not {ref}.")
    with INFERENCE_LOCK:
        result = predict.predict_variant(int(position), ref, alt)
    held = predict.held_out(variant)
    return {"result": result, "held_out": held,
            "sequence": {"before": sequence[max(0, index-3):index],
                         "ref": ref, "alt": alt, "after": sequence[index+1:index+4]}}


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(ROOT / "ui"), **kwargs)

    def respond(self, status: int, data: dict):
        body = json.dumps(data, allow_nan=False).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path == "/api/reference":
            from pipeline import predict
            try:
                self.respond(200, {"sequence": predict.REFERENCE.read_text().strip()})
            except FileNotFoundError:
                self.respond(503, {"error": "Study reference file is missing."})
        elif self.path == "/api/validation":
            try:
                self.respond(200, validation())
            except FileNotFoundError:
                self.respond(404, {"error": "No validation artifact. Run scripts/train_embeddings.py first."})
            except (ValueError, KeyError):
                self.respond(503, {"error": "Validation artifact is invalid."})
        elif self.path.startswith("/api/"):
            self.respond(404, {"error": "Unknown API route"})
        else:
            super().do_GET()

    def do_POST(self):
        # Loopback service: reject cross-origin requests that could trigger work.
        origin = self.headers.get("Origin")
        if origin and origin != f"http://{self.headers.get('Host')}":
            self.respond(403, {"error": "Cross-origin submissions are not allowed."})
            return
        if self.path != "/api/predict":
            self.respond(404, {"error": "Unknown API route"})
            return
        try:
            size = int(self.headers.get("Content-Length", "0"))
            if not 0 < size <= 16384:
                raise ValueError("Invalid request size")
            body = json.loads(self.rfile.read(size))
            if not isinstance(body, dict) or not isinstance(body.get("variant"), str):
                raise ValueError("A variant string is required")
            if not isinstance(body.get("sequence"), str):
                raise ValueError("An amino-acid reference sequence is required")
            self.respond(200, submit(body["variant"].strip().upper(), body["sequence"]))
        except (ValueError, TypeError) as exc:
            self.respond(400, {"error": str(exc)})
        except Exception as exc:
            self.log_error("Prediction failed: %s", exc)
            self.respond(503, {"error": "Pipeline unavailable: check the server terminal, reference/data files and model dependencies."})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args()
    os.chdir(ROOT)
    print(f"VariantNavigator: http://127.0.0.1:{args.port}", flush=True)
    ThreadingHTTPServer(("127.0.0.1", args.port), Handler).serve_forever()
