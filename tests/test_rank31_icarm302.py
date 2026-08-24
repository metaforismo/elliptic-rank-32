from __future__ import annotations

import importlib.util
import hashlib
import json
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "baseline" / "verify_rank31_icarm302.py"
INPUT = ROOT / "baseline" / "icarm_curve_302.json"
CERTIFICATE = ROOT / "baseline" / "rank31_icarm302_mod2_certificate.json"


def load_module():
    spec = importlib.util.spec_from_file_location("verify_rank31_icarm302", SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ICARMRank31CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_module()
        cls.fresh = cls.module.compute_certificate()
        cls.committed = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        cls.input_data = json.loads(INPUT.read_text(encoding="utf-8"))

    def test_committed_certificate_matches_recomputation(self):
        self.assertEqual(self.committed, self.fresh)

    def test_certificate_binds_both_verification_implementations(self):
        implementation = self.fresh["implementation"]
        self.assertEqual(
            implementation["script_sha256"], hashlib.sha256(SCRIPT.read_bytes()).hexdigest()
        )
        shared = ROOT / implementation["shared_exact_group_arithmetic"]
        self.assertEqual(
            implementation["shared_exact_group_arithmetic_sha256"],
            hashlib.sha256(shared.read_bytes()).hexdigest(),
        )

    def test_exact_unconditional_lower_bound_is_thirty_one(self):
        matrix = self.fresh["binary_matrix"]
        self.assertEqual(self.fresh["point_count"], 31)
        self.assertEqual(matrix["row_count"], 31)
        self.assertEqual(matrix["column_count"], 34)
        self.assertEqual(matrix["rank"], 31)
        self.assertEqual(len(matrix["pivot_columns_zero_based"]), 31)
        self.assertEqual(self.fresh["conditional_assumptions"], [])

    def test_generalized_model_and_short_model_formulas_are_pinned(self):
        self.assertEqual(self.input_data["a_invariants"][:3], ["1", "1", "1"])
        curve = self.fresh["curve"]
        invariants = {key: int(value) for key, value in curve["invariants"].items()}
        short = curve["short_model"]
        self.assertEqual(int(short["A"]), -27 * invariants["c4"])
        self.assertEqual(int(short["B"]), -54 * invariants["c6"])
        self.assertEqual(
            short["formula"],
            {
                "X": "36*x+3*b2",
                "Y": "108*(2*y+a1*x+a3)",
                "A": "-27*c4",
                "B": "-54*c6",
            },
        )
        self.assertEqual(curve["invariants"]["discriminant"], self.input_data["reported_discriminant"])

    def test_deterministic_local_prime_prefix_reaches_rank_thirty_one(self):
        records = self.fresh["local_quotients"]
        self.assertEqual(
            [record["prime"] for record in records],
            [
                113,
                241,
                263,
                281,
                283,
                337,
                347,
                373,
                389,
                467,
                491,
                541,
                601,
                607,
                617,
                653,
                701,
            ],
        )
        self.assertTrue(all(record["quotient_dimension"] == 2 for record in records))
        ranks = [record["matrix_rank_after_prime"] for record in records]
        self.assertEqual(ranks[-1], 31)
        self.assertLess(ranks[-2], 31)
        self.assertEqual(ranks, sorted(ranks))
        search = self.fresh["local_prime_search"]
        self.assertEqual(search["last_prime_considered"], 701)
        self.assertEqual(
            search["termination"],
            "first selected-prime prefix with binary row rank 31",
        )

    def test_torsion_certificate_is_exact_and_unconditional(self):
        torsion = self.fresh["torsion_certificate"]
        self.assertEqual(
            torsion["reductions"],
            [
                {"prime": 17, "group_order": 26},
                {"prime": 31, "group_order": 43},
            ],
        )
        first, second = torsion["reductions"]
        self.assertEqual(math.gcd(first["group_order"], second["group_order"]), 1)
        self.assertNotEqual(second["group_order"] % first["prime"], 0)
        self.assertNotEqual(first["group_order"] % second["prime"], 0)
        self.assertTrue(torsion["cross_characteristic_exclusions_verified"])

    def test_conditional_exact_rank_claim_is_separated(self):
        note = self.fresh["conditional_exact_rank_note"]
        self.assertIn("GRH+BSD", note)
        self.assertNotIn("rank E(Q)=31", self.fresh["claim"])
        self.assertEqual(
            self.input_data["reported_exact_rank_boundary"],
            {
                "claim": "rank equals 31 under GRH+BSD",
                "conditional_assumptions": ["GRH", "BSD"],
            },
        )


if __name__ == "__main__":
    unittest.main()
