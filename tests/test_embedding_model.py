"""Persistence contract tests using synthetic features, never biological metrics."""
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


@unittest.skipUnless(importlib.util.find_spec("sklearn"), "optional scikit-learn missing")
class PersistenceTests(unittest.TestCase):
    def test_roundtrip_and_contract_guards(self):
        import numpy as np
        from pipeline.embedding_model import predict_embedding
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "features.npz"
            reference = "ACDEFGHIK"
            sequence_hash = hashlib.sha256(reference.encode()).hexdigest()
            rng = np.random.default_rng(2)
            X = rng.normal(size=(60, 8)).astype(np.float32)
            y = X[:, 0] * 0.02
            np.savez(source, X=X, y=y, groups=np.repeat(np.arange(242, 262), 3),
                     rep1=y-0.002, rep2=y+0.002, names=np.array([f"variant{i}" for i in range(60)]),
                     model=np.array(["synthetic-test"]), crop=np.array([1, 9]),
                     reference=np.array([reference]), sequence_sha256=np.array([sequence_hash]),
                     dataset_sha256=np.array(["synthetic"]), feature_recipe=np.array(["synthetic"]),
                     partial=np.array([True]))
            run = root / "run"
            command = [sys.executable, "scripts/train_embeddings.py", "--npz", str(source),
                       "--out-dir", str(run), "--allow-partial", "--seed", "1"]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            manifest = json.loads((run / "manifest.json").read_text())
            self.assertFalse(set(manifest["train_positions"]) & set(manifest["test_positions"]))
            kwargs = dict(model="synthetic-test", sequence_sha256=sequence_hash, feature_recipe="synthetic")
            with self.assertRaisesRegex(ValueError, "Partial"):
                predict_embedding(run, X[0].tolist(), **kwargs)
            value = predict_embedding(run, X[0].tolist(), allow_partial=True, **kwargs)
            self.assertTrue(np.isfinite(value["predicted_netgrowth"]))
            with self.assertRaisesRegex(ValueError, "dimensions"):
                predict_embedding(run, [0.0], allow_partial=True, **kwargs)
            from pipeline.embedding_model import predict_variant
            with self.assertRaisesRegex(ValueError, "feature recipe"):
                predict_variant(run, "T315I", allow_partial=True)
            kwargs["model"] = "wrong-model"
            with self.assertRaisesRegex(ValueError, "contract"):
                predict_embedding(run, X[0].tolist(), allow_partial=True, **kwargs)

            kwargs["model"] = "synthetic-test"
            artifact = run / manifest["artifacts"][manifest["head"]]["file"]
            artifact.write_bytes(artifact.read_bytes() + b"corrupted")
            with self.assertRaisesRegex(ValueError, "checksum"):
                predict_embedding(run, X[0].tolist(), allow_partial=True, **kwargs)
