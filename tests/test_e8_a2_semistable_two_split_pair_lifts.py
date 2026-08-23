from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    ROOT / "research" / "analyze_e8_a2_semistable_two_split_pair_lifts.py"
)
INPUT_PATH = (
    ROOT / "certificates" / "e8_a2_semistable_two_split_target_small_primes.json"
)
CERTIFICATE_PATH = (
    ROOT / "certificates"
    / "e8_a2_semistable_two_split_pair_incidence_lifts.json"
)

SPEC = importlib.util.spec_from_file_location(
    "analyze_e8_a2_semistable_two_split_pair_lifts", MODULE_PATH
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class E8A2SemistableTwoSplitPairLiftTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.source = json.loads(INPUT_PATH.read_text(encoding="utf-8"))

    def test_all_representative_seeds_are_smooth_dimension_two(self) -> None:
        expected_determinants = [4, 7, 4, 6, 2, 4, 9, 3]
        observed_determinants: list[int] = []
        for seed in self.source["representative_pair_seeds"]["seeds"]:
            prime = seed["prime"]
            vector, witness = MODULE.seed_vector(seed)
            residual, derivative = MODULE.jacobian(vector, prime)
            self.assertFalse(any(residual))
            self.assertEqual(MODULE.matrix_rank(derivative, prime), 28)
            determinant = MODULE.determinant_mod(
                [row[:28] for row in derivative[:28]], prime
            )
            self.assertNotEqual(determinant, 0)
            observed_determinants.append(determinant)
            self.assertEqual(witness["gcd_h_tV1_plus_V2"], [1])
            self.assertEqual(MODULE.poly_gcd_mod(
                witness["Cx"], witness["Cy"], prime
            ), [1])
        self.assertEqual(observed_determinants, expected_determinants)

    def test_first_seed_has_deterministic_formal_lift(self) -> None:
        seed = self.source["representative_pair_seeds"]["seeds"][0]
        vector, _ = MODULE.seed_vector(seed)
        lifted, trace, determinant = MODULE.hensel_lift_smooth_slice(
            vector, seed["prime"], 5
        )
        modulus = seed["prime"] ** 5
        self.assertEqual(determinant, 4)
        self.assertEqual(len(trace), 5)
        self.assertFalse(any(MODULE.evaluate_system(lifted, modulus)))
        self.assertEqual(lifted[-2:], vector[-2:])

    def test_certificate_keeps_claim_boundary(self) -> None:
        certificate = MODULE.compute_certificate(INPUT_PATH, digits=8)
        committed = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
        self.assertEqual(certificate, committed)
        self.assertTrue(
            certificate["claim_boundary"]["formal_p_adic_pair_lift_found"]
        )
        self.assertFalse(certificate["claim_boundary"]["rational_pair_lift_found"])
        self.assertFalse(certificate["claim_boundary"]["P3_lift_found"])
        self.assertFalse(certificate["claim_boundary"]["rank31_curve_found"])


if __name__ == "__main__":
    unittest.main()
