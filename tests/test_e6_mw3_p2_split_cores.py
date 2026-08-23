from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "replay_e6_mw3_p2_split_cores.py"
CERTIFICATE_PATH = ROOT / "certificates" / "e6_mw3_p2_split_cores_gf31.json"
SYMPY_AVAILABLE = importlib.util.find_spec("sympy") is not None

SPEC = importlib.util.spec_from_file_location("e6_mw3_p2_replay", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class E6Mw3P2SplitCoreTests(unittest.TestCase):
    @unittest.skipUnless(SYMPY_AVAILABLE, "SymPy exact elimination unavailable")
    def test_certificate_matches_and_hit_is_dependent(self) -> None:
        computed = MODULE.compute_certificate()
        committed = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
        self.assertEqual(computed, committed)
        self.assertEqual(
            computed["summary"],
            {
                "declared_cores": 7,
                "reconstructed_surfaces": 2,
                "tested_pole_q0_pairs": 1674,
                "signed_canonical_P2_hits": 2,
                "geometric_canonical_P2_hits": 1,
                "dependent_geometric_P2_hits": 1,
                "independent_rank3_seed_hits": 0,
            },
        )
        core7 = computed["core_records"][6]
        hit = core7["surfaces"][0]["geometric_P2_hits"][0]
        self.assertFalse(hit["independent_rank3_seed"])
        self.assertEqual(
            hit["exact_P1_P2_P3_relations"],
            [
                {
                    "sign_P2": 1,
                    "sign_P3": 1,
                    "relation": "P1+sign_P2*P2+sign_P3*P3=O",
                }
            ],
        )

    def test_polynomial_square_root_round_trip(self) -> None:
        root = [3, 0, 8, 5]
        square = MODULE.poly_pow(root, 2)
        recovered = MODULE.polynomial_sqrt_roots(square)
        self.assertIn(MODULE.trim(root), recovered)
        self.assertIn(MODULE.poly_neg(root), recovered)


if __name__ == "__main__":
    unittest.main()
