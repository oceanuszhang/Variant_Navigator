"""Submission retry identity checks; no paid service calls."""
import unittest
from pipeline.boltz_structure import idempotency_key


class IdempotencyTests(unittest.TestCase):
    def test_payload_and_model_define_identity(self):
        first = {"num_samples": 1, "entities": [{"type": "protein", "value": "ACD"}]}
        reordered = {"entities": first["entities"], "num_samples": 1}
        self.assertEqual(idempotency_key(first), idempotency_key(reordered))
        self.assertNotEqual(idempotency_key(first), idempotency_key({**first, "num_samples": 5}))
        self.assertNotEqual(idempotency_key(first), idempotency_key(first, "another-model"))
