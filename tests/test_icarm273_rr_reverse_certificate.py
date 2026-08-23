from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "icarm273_rr_reverse_certificate.py"
CERTIFICATE_PATH = ROOT / "certificates" / "icarm273_rr_reverse_certificate.json"
SYMPY_AVAILABLE = importlib.util.find_spec("sympy") is not None

SPEC = importlib.util.spec_from_file_location(
    "icarm273_rr_reverse_certificate", MODULE_PATH
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class Icarm273ReverseRiemannRochTests(unittest.TestCase):
    @unittest.skipUnless(SYMPY_AVAILABLE, "SymPy exact QQ linear algebra unavailable")
    def test_reverse_certificate_matches(self) -> None:
        computed = MODULE.compute_certificate()
        committed = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
        self.assertEqual(computed, committed)
        self.assertEqual(computed["exact_linear_algebra"]["rank"], 30)
        self.assertEqual(computed["exact_linear_algebra"]["nullity"], 1)
        self.assertEqual(computed["norm_identity"]["degree"], 31)
        self.assertEqual(computed["norm_identity"]["distinct_rational_roots"], 31)
        self.assertTrue(
            computed["norm_identity"]["verified_by_exact_coefficient_comparison"]
        )


if __name__ == "__main__":
    unittest.main()
