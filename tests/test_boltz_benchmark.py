"""Tests for fixed benchmark setup and confidence comparisons; no service calls."""
import unittest
from scripts.boltz_benchmark import plan, summarize, PANEL


class BenchmarkTests(unittest.TestCase):
    def test_fixed_panel_settings_and_measurement_provenance(self):
        jobs = plan()
        self.assertEqual(tuple(j["variant"] for j in jobs), PANEL)
        self.assertEqual(len({j["idempotency_key"] for j in jobs}), 4)
        for job in jobs:
            self.assertEqual(job["body"]["num_samples"], 5)
            self.assertEqual(len(job["body"]["entities"][0]["value"]), 271)
            self.assertNotIn("constraints", job["body"])
            self.assertEqual(job["body"]["entities"][1], jobs[0]["body"]["entities"][1])
        self.assertIsNone(jobs[0]["measured_netgrowth"])
        self.assertGreater(jobs[1]["measured_netgrowth"], jobs[2]["measured_netgrowth"])
        self.assertGreater(jobs[3]["measured_netgrowth"], jobs[2]["measured_netgrowth"])

    def test_contrasts_do_not_invent_missing_scores(self):
        results = {"WILDTYPE": {"binding_confidence": 0.8},
                   "T315I": {"binding_confidence": 0.2},
                   "T315A": {"binding_confidence": 0.9}}
        summary = summarize(results)
        self.assertTrue(summary["prespecified_contrasts"][0]["binding_confidence_lower_than_control"])
        self.assertIsNone(summary["prespecified_contrasts"][1]["binding_confidence_lower_than_control"])
        self.assertAlmostEqual(summary["rows"][1]["delta_binding_confidence_vs_wt"], -0.6)
        self.assertNotIn("accuracy", summary)

    def test_cached_lookup_checks_source_hashes(self):
        import hashlib
        import json
        from pathlib import Path
        import tempfile
        from unittest.mock import patch
        from pipeline import benchmark
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = {"reference_sha256": hashlib.sha256(Path("data/cancer-exp/ABL1_NM_005157.6_DMS_reference_aa.txt").read_bytes()).hexdigest(),
                        "dataset_sha256": hashlib.sha256(Path("data/datasets/abl_dms_netgrowth_1200nM.csv").read_bytes()).hexdigest(),
                        "model": "boltz-2.1", "crop": [242,512], "num_samples": 5,
                        "compound": "imatinib", "constraints": "none"}
            (root/"manifest.json").write_text(json.dumps(manifest))
            (root/"summary.json").write_text(json.dumps({"manifest_sha256": hashlib.sha256((root/"manifest.json").read_bytes()).hexdigest(), "rows": [{"variant": "T315I", "binding_confidence": 0.3,
                "delta_binding_confidence_vs_wt": -0.1, "confidence": 90, "confidence_sd": 1,
                "prediction_id": "mock"}], "caveats": ["test only"]}))
            with patch.object(benchmark, "ROOT", root):
                self.assertEqual(benchmark.lookup("T315I")["prediction_id"], "mock")
                self.assertIsNone(benchmark.lookup("Y253H"))
                manifest["reference_sha256"] = "wrong"
                (root/"manifest.json").write_text(json.dumps(manifest))
                with self.assertRaisesRegex(ValueError, "reference changed"):
                    benchmark.lookup("T315I")

    def test_estimate_does_not_relabel_completed_results(self):
        import contextlib
        import io
        from pathlib import Path
        import tempfile
        from unittest.mock import patch
        from scripts import boltz_benchmark as runner
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "manifest.json").write_text("original execution manifest")
            (root / "summary.json").write_text("original completed summary")
            with patch.object(runner, "OUT", root), patch.object(runner.bs, "estimate_cost", return_value={"estimated_cost_usd": "0.1"}), patch("sys.argv", ["benchmark"]), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(runner.main(), 0)
            self.assertEqual((root / "manifest.json").read_text(), "original execution manifest")
            self.assertEqual((root / "summary.json").read_text(), "original completed summary")
            self.assertTrue((root / "estimate.json").exists())
