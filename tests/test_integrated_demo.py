from __future__ import annotations

import unittest

from pipeline import boltz_structure
from pipeline.model_store import (
    predict_saved,
    reference_sequence,
    validation_summary,
)
from ui.server import ptarmigan_control, submit_prediction


class IntegratedDemoTests(unittest.TestCase):
    def test_saved_xgboost_demo_is_tied_to_held_out_artifact(self):
        result = predict_saved("E450R", reference_sequence(), "xgboost")
        self.assertEqual(result["measured_class"], "increase")
        self.assertEqual(result["held_out"]["prediction"], "increase")
        self.assertTrue(result["held_out"]["correct"])
        self.assertEqual(
            set(result["probabilities"]),
            {"decrease", "unchanged", "increase"},
        )

    def test_server_prediction_rejects_unknown_experimental_variant(self):
        with self.assertRaisesRegex(ValueError, "not one of"):
            submit_prediction({
                "variant": "E001A",
                "sequence": reference_sequence(),
                "model": "ridge",
            })

    def test_validation_artifact_is_site_grouped(self):
        data = validation_summary()
        split = data["metrics"]["split"]
        self.assertEqual(split["validation_rows"], 944)
        self.assertEqual(split["training_positions"], 216)
        self.assertEqual(split["validation_positions"], 54)
        self.assertFalse(
            set(split["training_position_values"])
            & set(split["validation_position_values"])
        )

    def test_ptarmigan_is_labeled_wild_type_control(self):
        data = ptarmigan_control()
        self.assertIn("wild-type", data["scope"].lower())
        self.assertEqual(data["hits"]["data"][0]["name"].lower(), "imatinib")

    def test_boltz_request_is_cropped_without_submission(self):
        crop = boltz_structure.crop_sequence(reference_sequence())
        body = boltz_structure.build_input(
            crop, compound="imatinib", pocket_hint=[286, 315, 381]
        )
        self.assertEqual(len(crop), 271)
        self.assertEqual(body["num_samples"], 5)
        self.assertEqual(body["binding"]["type"], "ligand_protein_binding")
        self.assertFalse(body["constraints"][0]["force"])


if __name__ == "__main__":
    unittest.main()
