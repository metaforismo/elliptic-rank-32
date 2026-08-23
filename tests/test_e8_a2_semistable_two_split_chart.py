from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "derive_e8_a2_semistable_two_split_chart.py"
CERTIFICATE_PATH = ROOT / "certificates" / "e8_a2_semistable_two_split_chart.json"
SYMPY_AVAILABLE = importlib.util.find_spec("sympy") is not None

SPEC = importlib.util.spec_from_file_location(
    "derive_e8_a2_semistable_two_split_chart", MODULE_PATH
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class E8A2SemistableTwoSplitChartTests(unittest.TestCase):
    @unittest.skipUnless(SYMPY_AVAILABLE, "SymPy exact algebra unavailable")
    def test_certificate_recomputes(self) -> None:
        computed = MODULE.compute_certificate()
        committed = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
        self.assertEqual(computed, committed)
        self.assertEqual(computed["parameters"], ["lambda", "y", "z", "W"])
        self.assertEqual(computed["polynomials"]["H_degree"], 5)
        self.assertFalse(computed["claim_boundary"]["rank31_curve_found"])

    @unittest.skipUnless(SYMPY_AVAILABLE, "SymPy exact algebra unavailable")
    def test_third_split_fibre_is_not_assumed(self) -> None:
        computed = MODULE.compute_certificate()
        self.assertEqual(
            computed["gauge_and_node_values"]["third_fibre_split_iff"],
            "g is a nonzero square",
        )
        self.assertIn(
            "unnecessary for the target",
            computed["comparison_with_all_split_chart"],
        )


if __name__ == "__main__":
    unittest.main()
