from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "derive_e8_a2_semistable_two_split_boundaries.py"
CERTIFICATE_PATH = ROOT / "certificates" / "e8_a2_semistable_two_split_boundaries.json"
SYMPY_AVAILABLE = importlib.util.find_spec("sympy") is not None

SPEC = importlib.util.spec_from_file_location(
    "derive_e8_a2_semistable_two_split_boundaries", MODULE_PATH
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class E8A2SemistableTwoSplitBoundaryTests(unittest.TestCase):
    @unittest.skipUnless(SYMPY_AVAILABLE, "SymPy exact algebra unavailable")
    def test_certificate_recomputes(self) -> None:
        computed = MODULE.compute_certificate()
        committed = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
        self.assertEqual(computed, committed)
        self.assertEqual(
            [chart["name"] for chart in computed["charts"]],
            ["s0_zero", "s_lambda_zero", "node1_beta_gamma_zero"],
        )
        self.assertTrue(all(chart["H_degree"] == 5 for chart in computed["charts"]))

    @unittest.skipUnless(SYMPY_AVAILABLE, "SymPy exact algebra unavailable")
    def test_no_rank31_claim(self) -> None:
        computed = MODULE.compute_certificate()
        self.assertFalse(computed["claim_boundary"]["rank31_curve_found"])
        self.assertFalse(
            computed["claim_boundary"]["finite_field_boundaries_searched"]
        )


if __name__ == "__main__":
    unittest.main()
