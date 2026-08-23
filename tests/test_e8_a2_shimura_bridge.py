from __future__ import annotations

import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "certify_e8_a2_shimura_bridge.py"
CERTIFICATE_PATH = ROOT / "certificates" / "e8_a2_shimura_bridge.json"
NOTE_PATH = ROOT / "research" / "e8_a2_shimura_bridge.md"

SPEC = importlib.util.spec_from_file_location(
    "certify_e8_a2_shimura_bridge", MODULE_PATH
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class E8A2ShimuraBridgeTests(unittest.TestCase):
    def test_certificate_recomputes_exactly(self) -> None:
        computed = MODULE.compute_certificate()
        committed = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
        self.assertEqual(computed, committed)
        self.assertTrue(
            computed["claim_boundary"]["finite_quadratic_form_bridge_proved"]
        )
        self.assertFalse(computed["claim_boundary"]["P3_found"])

    def test_target_root_system_and_discriminant_form(self) -> None:
        computed = MODULE.compute_certificate()
        target = computed["target_essential_lattice"]
        self.assertEqual(target["determinant"], 948)
        self.assertEqual(target["root_system"], "E8+A2^3")
        self.assertEqual(target["root_count"], 258)
        self.assertEqual(target["roots_using_adjusted_section_coordinates"], 0)
        self.assertEqual(target["discriminant_group"], "Z/948Z")
        self.assertEqual(
            computed["finite_quadratic_form_isometries"]["target_congruence"],
            "1709*67^2 == 485 (mod 1896)",
        )

    def test_period_clifford_ramification(self) -> None:
        computed = MODULE.compute_certificate()
        period = computed["x6_79_period_lattice"]
        self.assertEqual(period["even_clifford_quaternion"], "(6,2)")
        self.assertEqual(period["quaternion_finite_ramification"], [2, 3])
        self.assertEqual(period["remaining_squarefree_level_factor"], 79)
        self.assertEqual(period["hilbert_symbols"]["79"], 1)

    def test_public_replay_terminal_matrix(self) -> None:
        note = NOTE_PATH.read_text(encoding="utf-8")
        matrix_text = (
            note.split(
                "The terminal rootless Gram matrix recovered directly from the log is",
                1,
            )[1]
            .split("```text", 1)[1]
            .split("```", 1)[0]
            .strip()
        )
        gram = ast.literal_eval(matrix_text)
        self.assertEqual(len(gram), 17)
        self.assertTrue(all(len(row) == 17 for row in gram))
        self.assertEqual(gram, [list(row) for row in zip(*gram)])
        canonical = json.dumps(gram, separators=(",", ":")).encode()
        self.assertEqual(
            hashlib.sha256(canonical).hexdigest(),
            "620a5e06473684d3e8015c0172f63c09c901e742ec02e77ba0aa35a923aa0295",
        )
        self.assertEqual(MODULE.bareiss_determinant(gram), 948)
        self.assertEqual(MODULE.enumerate_norm_two(gram), [])


if __name__ == "__main__":
    unittest.main()
