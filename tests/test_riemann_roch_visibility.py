from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "riemann_roch_visibility_certificate.py"
CERTIFICATE_PATH = ROOT / "certificates" / "riemann_roch_visibility_barrier.json"

SPEC = importlib.util.spec_from_file_location(
    "riemann_roch_visibility_certificate", MODULE_PATH
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class RiemannRochVisibilityTests(unittest.TestCase):
    def test_certificate_matches_committed_json(self) -> None:
        computed = MODULE.compute_certificate()
        committed = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
        self.assertEqual(computed, committed)

    def test_near_square_barrier(self) -> None:
        result = MODULE.compute_certificate()
        frontier = result["near_square_frontier"][
            "table_degrees_2_through_16"
        ][-2]
        self.assertEqual(frontier["degree_Q"], 15)
        self.assertEqual(frontier["pole_order"], 30)
        self.assertEqual(
            frontier["rank_upper_bound_from_visible_points"], 29
        )

    def test_minimal_order_31_shape(self) -> None:
        result = MODULE.compute_certificate()
        target = result["minimal_single_function_rank30_target"]
        self.assertEqual(target["minimal_pole_order"], 31)
        self.assertEqual(
            target["canonical_choice"],
            {
                "deg_A_max": 15,
                "deg_B": 14,
                "deg_R": 3,
                "norm_form": "A(x)^2-R(x)*B(x)^2",
                "generic_norm_degree": 31,
            },
        )
        self.assertEqual(
            target["maximum_possible_span_after_forced_relation"], 30
        )
        self.assertEqual(
            result["theorem"]["conditional_assumptions"], []
        )

    def test_minimal_order_32_shape_and_universality(self) -> None:
        result = MODULE.compute_certificate()
        target = result["minimal_single_function_rank31_target"]
        self.assertEqual(target["minimal_pole_order"], 32)
        self.assertEqual(
            target["canonical_choice"],
            {
                "deg_A": 16,
                "deg_B_max": 14,
                "deg_R": 3,
                "norm_form": "A(x)^2-R(x)*B(x)^2",
                "generic_norm_degree": 32,
            },
        )
        self.assertEqual(
            target["maximum_possible_span_after_forced_relation"], 31
        )
        self.assertIn("P32=-(P1+...+P31)", target["universal_reverse_construction"])


if __name__ == "__main__":
    unittest.main()
