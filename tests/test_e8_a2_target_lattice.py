from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "certify_e8_a2_target_lattice.py"
CERTIFICATE_PATH = ROOT / "certificates" / "e8_a2_target_lattice.json"

SPEC = importlib.util.spec_from_file_location(
    "certify_e8_a2_target_lattice", MODULE_PATH
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class E8A2TargetLatticeTests(unittest.TestCase):
    def test_certificate_recomputes(self) -> None:
        computed = MODULE.compute_certificate()
        committed = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
        self.assertEqual(computed, committed)
        self.assertTrue(
            computed["saturation"]["generated_rank_three_lattice_is_primitive"]
        )
        self.assertEqual(computed["target_height_gram"]["scaled_determinant"], 948)

    def test_unique_profiles_and_intersections(self) -> None:
        computed = MODULE.compute_certificate()
        self.assertEqual(
            computed["unique_component_orbit"]["canonical_labels_at_three_I3"],
            {"P1": [1, 1, 0], "P2": [2, 0, 0], "P3": [0, 0, 0]},
        )
        self.assertEqual(
            set(
                computed["unique_component_orbit"][
                    "pairwise_section_intersections"
                ].values()
            ),
            {"2"},
        )
        self.assertEqual(computed["torsion"]["status"], "trivial")


if __name__ == "__main__":
    unittest.main()
