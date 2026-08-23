from __future__ import annotations

import ast
import hashlib
import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "research/search_e8_a2_target_neighbor_bridge.py"
WORKFLOW = ROOT / ".github/workflows/probe-e8-a2-target-neighbor-bridge.yml"


def load_search_module():
    spec = importlib.util.spec_from_file_location("target_neighbor_bridge", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TargetNeighborBridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.module = load_search_module()

    def test_frozen_inputs_have_the_certified_hashes(self) -> None:
        self.assertEqual(
            self.module.payload_sha256(self.module.target_essential_lattice()),
            self.module.TARGET_FORMULA_HASH,
        )
        self.assertEqual(
            self.module.payload_sha256(self.module.transparent_seed_lattice()),
            self.module.TRANSPARENT_SEED_FORMULA_HASH,
        )
        self.assertEqual(
            self.module.payload_sha256(self.module.FROZEN_ROOTLESS_GRAM),
            self.module.FROZEN_ROOTLESS_HASH,
        )

    def test_prime_parser_is_fail_closed(self) -> None:
        self.assertEqual(self.module.parse_primes("2,5"), (2, 5))
        with self.assertRaises(Exception):
            self.module.parse_primes("")
        with self.assertRaises(Exception):
            self.module.parse_primes("2,4")
        with self.assertRaises(Exception):
            self.module.parse_primes("5,5")

    def test_script_preserves_exact_neighbor_and_isometry_matrices(self) -> None:
        source = SCRIPT.read_text()
        ast.parse(source)
        required_fragments = [
            "find_primitive_p_divisible_vector__next",
            "find_p_neighbor_from_vec(ZZ(prime), vector, return_matrix=True)",
            "raw_transform * lll_transform",
            "transition.transpose() * parent_gram * transition != child_gram",
            "rank_mod_p_of_p_times_transition",
            "is_globally_equivalent_to(right_form, return_matrix=True)",
            '"parent_to_child_rational_transform"',
            '"exact-bridge.json"',
            '"bounded_search_completed"',
        ]
        for fragment in required_fragments:
            self.assertIn(fragment, source)

    def test_workflow_is_branch_scoped_manual_pinned_and_uploads_even_on_failure(self) -> None:
        source = WORKFLOW.read_text()
        self.assertIn("workflow_dispatch:", source)
        self.assertIn("\n  push:", source)
        self.assertIn("- codex/rank31-research", source)
        self.assertIn('"research/search_e8_a2_target_neighbor_bridge.py"', source)
        self.assertNotIn("branches: [main]", source)
        self.assertNotIn("pull_request:", source)
        self.assertNotIn("schedule:", source)
        self.assertIn("image: sagemath/sagemath:10.9", source)
        self.assertIn("if: always()", source)
        self.assertIn(
            "actions/upload-artifact@ea165f8d65b6e75b540449e92b4886f43607fa02",
            source,
        )
        self.assertNotIn("actions/upload-artifact@v4", source)
        self.assertIn("SHA256SUMS", source)

    def test_script_digest_is_not_accidentally_empty(self) -> None:
        digest = hashlib.sha256(SCRIPT.read_bytes()).hexdigest()
        self.assertEqual(len(digest), 64)
        self.assertNotEqual(digest, hashlib.sha256(b"").hexdigest())


if __name__ == "__main__":
    unittest.main()
