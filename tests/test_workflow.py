"""Deterministic workflow checks; external services never run."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import main
from pipeline import predict, boltz_structure as bs, agent_report
from pipeline.fetch import Protein


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.seq = predict.REFERENCE.read_text().strip()
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.report_patch = patch.object(agent_report, "REPORT_DIR", Path(self.tmp.name))
        self.report_patch.start()
        self.addCleanup(self.report_patch.stop)

    def test_offline_measurement_and_separate_validation(self):
        with patch.object(main, "fetch_protein", side_effect=AssertionError("network")), patch.object(predict, "_fit", side_effect=AssertionError("fit")):
            result = main.run("ABL1", "T315I", "oncology", offline=True, use_agent=False)
        self.assertFalse(result["prediction"]["predicted"])
        self.assertEqual(result["validation"]["kind"], "held_out")
        report = next(Path(self.tmp.name).glob("*.md")).read_text(encoding="utf-8")
        self.assertIn("Skipped or unsupported stages", report)
        self.assertIn("Saved held-out model estimate", report)

    def test_optional_boltz_failure_still_reports(self):
        protein = Protein("P00519", "ABL1", "ABL1", "Homo sapiens", self.seq)
        for error in (bs.BoltzUnavailable("no auth"), bs.BoltzTooExpensive("cap")):
            with self.subTest(error=type(error).__name__), patch.object(main, "fetch_protein", return_value=protein), patch.object(main.esm2_score, "available", return_value=False), patch.object(bs, "predict_structure", side_effect=error), patch.object(main.baseline, "predict_panel", return_value=[]):
                result = main.run("ABL1", "T315I", "oncology", run_boltz=True, use_agent=False)
            self.assertTrue(any(x["stage"] == "boltz" for x in result["skipped_stages"]))
            self.assertIn("prediction", result)

    def test_egfr_never_uses_abl_crop(self):
        seq = "A" * 789 + "T" + "A" * 420
        protein = Protein("P00533", "EGFR", "EGFR", "Homo sapiens", seq)
        with patch.object(main, "fetch_protein", return_value=protein), patch.object(main.esm2_score, "available", return_value=False), patch.object(bs, "predict_structure") as call, patch.object(main.baseline, "predict_panel", return_value=[]):
            result = main.run("EGFR", "T790M", "oncology", run_boltz=True, use_agent=False)
        call.assert_not_called()
        self.assertNotIn("prediction", result)

    def test_corrupt_validation_does_not_lose_measurement(self):
        bad = Path(self.tmp.name) / "bad.json"
        bad.write_text("not json")
        with patch.object(predict, "VALIDATION", bad):
            result = main.run("ABL1", "T315I", "oncology", offline=True, use_agent=False)
        self.assertEqual(result["prediction"]["kind"], "measurement")
        self.assertTrue(any(x["stage"] == "validation" for x in result["skipped_stages"]))

    def test_boltz_metadata_survives_report(self):
        protein = Protein("P00519", "ABL1", "ABL1", "Homo sapiens", self.seq)
        structure = bs.Structure(Path("mock.cif"), self.seq, 85.0, binding_confidence=0.7,
                                 num_samples=5, crop=(242, 512), cost_usd=0.1)
        with patch.object(main, "fetch_protein", return_value=protein), patch.object(main.esm2_score, "available", return_value=False), patch.object(bs, "predict_structure", return_value=structure), patch.object(main.baseline, "predict_panel", return_value=[]):
            result = main.run("ABL1", "T315I", "oncology", run_boltz=True, use_agent=False)
        self.assertEqual(result["structure"]["binding_confidence"], 0.7)
        self.assertEqual(result["structure"]["num_samples"], 5)
        self.assertNotIn("sequence", result["structure"])

    def test_esm_service_failure_preserves_measurement(self):
        protein = Protein("P00519", "ABL1", "ABL1", "Homo sapiens", self.seq)
        with patch.object(main, "fetch_protein", return_value=protein), patch.object(main.esm2_score, "available", return_value=False), patch.object(main.baseline, "predict_panel", return_value=[]), patch("pipeline.embedding_model.predict_variant", side_effect=RuntimeError("service unavailable")):
            result = main.run("ABL1", "T315I", "oncology", esm_run_dir=Path("mock"), use_agent=False)
        self.assertEqual(result["prediction"]["kind"], "measurement")
        self.assertTrue(any(x["stage"] == "esm_inference" for x in result["skipped_stages"]))

    def test_scope_and_residue_guards(self):
        for args in (("RPOB", "S450L", "oncology"), ("ABL1", "T315I", "agriculture")):
            with self.assertRaises(ValueError):
                main.run(*args)
        for args in ((999, "A", "V"), (315, "A", "I"), (315, "T", "T")):
            with self.assertRaises(predict.PredictError):
                predict.predict_variant(*args)
        with self.assertRaises(predict.PredictError):
            predict.predict_variant(315, "T", "I", compound="erlotinib")
        self.assertTrue(all(json.loads(p.read_text())["mode"] == "oncology" for p in main.oncology_cases()))

    def test_held_out_metrics_are_bound_to_artifact_head(self):
        row = predict.held_out("T315I")
        self.assertEqual(row["head"], "hist_gbm")
        self.assertEqual(len(row["artifact_sha256"]), 64)
        self.assertNotIn("cv_pearson_r", predict.predict_variant(315, "T", "I"))


if __name__ == "__main__":
    unittest.main()
