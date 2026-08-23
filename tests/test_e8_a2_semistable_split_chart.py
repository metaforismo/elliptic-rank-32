from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "derive_e8_a2_semistable_split_chart.py"
CERTIFICATE_PATH = ROOT / "certificates" / "e8_a2_semistable_split_chart.json"
SYMPY_AVAILABLE = importlib.util.find_spec("sympy") is not None

SPEC = importlib.util.spec_from_file_location(
    "derive_e8_a2_semistable_split_chart", MODULE_PATH
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class E8A2SemistableSplitChartTests(unittest.TestCase):
    @unittest.skipUnless(SYMPY_AVAILABLE, "SymPy exact algebra unavailable")
    def test_certificate_recomputes_exactly(self) -> None:
        computed = MODULE.compute_certificate()
        committed = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
        self.assertEqual(computed, committed)
        self.assertEqual(computed["parameters"], ["lambda", "m", "R", "W"])
        self.assertEqual(computed["polynomials"]["residual_H_degree"], 5)
        self.assertFalse(computed["claim_boundary"]["rank31_curve_found"])

    @unittest.skipUnless(SYMPY_AVAILABLE, "SymPy exact algebra unavailable")
    def test_exact_relation_and_split_conic_are_recorded(self) -> None:
        computed = MODULE.compute_certificate()
        self.assertEqual(
            computed["polynomials"]["relation"],
            "a*gamma-3*beta^2=d*D",
        )
        self.assertEqual(
            computed["auxiliary_ratios"]["conic"],
            "y^2=(1-lambda)+lambda*x^2",
        )
        self.assertIn("dense birational chart", computed["coverage_boundary"])


if __name__ == "__main__":
    unittest.main()
