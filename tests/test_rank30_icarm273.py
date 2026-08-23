from __future__ import annotations

import importlib.util
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "baseline" / "verify_rank30_icarm273.py"
CERTIFICATE = ROOT / "baseline" / "rank30_icarm273_mod2_certificate.json"


def load_module():
    spec = importlib.util.spec_from_file_location("verify_rank30_icarm273", SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ICARMRank30CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_module()
        cls.fresh = cls.module.compute_certificate()
        cls.committed = json.loads(CERTIFICATE.read_text(encoding="utf-8"))

    def test_committed_certificate_matches_recomputation(self):
        self.assertEqual(self.committed, self.fresh)

    def test_exact_lower_bound_is_thirty(self):
        self.assertEqual(self.fresh["point_count"], 30)
        self.assertEqual(self.fresh["binary_matrix"]["rank"], 30)
        self.assertEqual(
            self.fresh["binary_matrix"]["pivot_columns_zero_based"],
            list(range(30)),
        )

    def test_torsion_certificate_is_unconditional(self):
        self.assertEqual(
            self.fresh["torsion_certificate"]["reductions"],
            [
                {"prime": 11, "group_order": 18},
                {"prime": 59, "group_order": 73},
            ],
        )
        self.assertEqual(self.fresh["conditional_assumptions"], [])

    def test_conditional_upper_bound_is_separated(self):
        note = self.fresh["conditional_exact_rank_note"]
        self.assertIn("GRH+BSD", note)
        self.assertNotIn("rank E(Q)=30", self.fresh["claim"])


if __name__ == "__main__":
    unittest.main()
