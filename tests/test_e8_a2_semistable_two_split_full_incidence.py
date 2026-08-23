from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    ROOT / "research" / "derive_e8_a2_semistable_two_split_full_incidence.py"
)
CERTIFICATE_PATH = (
    ROOT / "certificates" / "e8_a2_semistable_two_split_full_incidence.json"
)
SYMPY_AVAILABLE = importlib.util.find_spec("sympy") is not None

if SYMPY_AVAILABLE:
    SPEC = importlib.util.spec_from_file_location(
        "derive_e8_a2_semistable_two_split_full_incidence", MODULE_PATH
    )
    assert SPEC is not None and SPEC.loader is not None
    MODULE = importlib.util.module_from_spec(SPEC)
    sys.modules[SPEC.name] = MODULE
    SPEC.loader.exec_module(MODULE)
else:
    MODULE = None


@unittest.skipUnless(SYMPY_AVAILABLE, "SymPy exact algebra unavailable")
class E8A2SemistableTwoSplitFullIncidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.computed = MODULE.compute_certificate()

    def test_certificate_recomputes(self) -> None:
        committed = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
        self.assertEqual(self.computed, committed)

    def test_full_incidence_dimension_and_jacobian_boundary(self) -> None:
        incidence = self.computed["exact_full_incidence"]
        self.assertEqual(incidence["raw_variable_count"], 62)
        self.assertEqual(incidence["raw_equation_count"], 67)
        self.assertEqual(incidence["syzygetic_equation_count"], 6)
        self.assertEqual(incidence["reduced_local_equation_count"], 61)
        self.assertEqual(incidence["reduced_local_jacobian_shape"], [61, 62])
        self.assertEqual(incidence["expected_dimension"], 1)
        self.assertIsNone(incidence["observed_full_triple_jacobian_rank"])

    def test_triangular_P3_reduction_and_modular_replay(self) -> None:
        triangular = self.computed["P3_triangular_elimination"]
        self.assertEqual(
            [entry["section_variable_total_degree"]
             for entry in triangular["remainder_equations"]],
            [21, 19, 17, 15, 14, 12],
        )
        replay = self.computed["modular_evidence"][
            "independent_first_seed_P3_remainder_replay"
        ]
        self.assertEqual(replay["raw_r_Xlower_tuples"], 146410)
        self.assertEqual(replay["node_admissible_tuples"], 110000)
        self.assertEqual(replay["P3_hit_count"], 0)

    def test_claim_boundary(self) -> None:
        boundary = self.computed["claim_boundary"]
        self.assertFalse(boundary["full_modular_triple_found"])
        self.assertFalse(boundary["Noether_Lefschetz_polynomial_expanded"])
        self.assertFalse(boundary["characteristic_zero_triple_found"])
        self.assertFalse(boundary["rank31_curve_found"])


if __name__ == "__main__":
    unittest.main()
