from __future__ import annotations

import importlib.util
import sys
import tempfile
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "research" / "verify_rank17_neighbor_chain_artifact.py"
SPEC = importlib.util.spec_from_file_location("rank17_artifact_verifier", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
VERIFIER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VERIFIER)


class Rank17ArtifactVerifierTests(unittest.TestCase):
    def test_public_chain_pins_lll_reduced_seed_and_target(self) -> None:
        self.assertEqual(
            VERIFIER.EXPECTED_SEED_HASH,
            "d29e4bb24bcb1db376f142708ceb841f02c7eba4b27a7a04d141b0f06b8d40d1",
        )
        self.assertEqual(
            VERIFIER.EXPECTED_TARGET_HASH,
            "620a5e06473684d3e8015c0172f63c09c901e742ec02e77ba0aa35a923aa0295",
        )

    def test_bareiss_and_matrix_hash_are_exact(self) -> None:
        matrix = [[4, -2, 0], [-2, 4, -2], [0, -2, 4]]
        self.assertEqual(VERIFIER.determinant_bareiss(matrix), 32)
        self.assertEqual(
            VERIFIER.canonical_matrix_hash(matrix),
            "1e9fb8f7f48f3323c4f7f4a738de932275f3e5594db17fa7369e02f391a54b1c",
        )

    def test_duplicate_json_key_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "duplicate.json"
            path.write_text('{"a": 1, "a": 2}\n')
            with self.assertRaises(VERIFIER.VerificationError):
                VERIFIER.load_json(path)

    def test_manifest_detects_any_unmanifested_or_mutated_file(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            payload = root / "payload.txt"
            payload.write_text("exact evidence\n")
            for name in VERIFIER.REQUIRED_FILES:
                (root / name).write_text(f"fixture for {name}\n")
            VERIFIER.write_manifest(root, VERIFIER.EXPECTED_SEMANTIC_SUMMARY)
            VERIFIER.verify_manifest(root)
            payload.write_text("mutated evidence\n")
            with self.assertRaises(VERIFIER.VerificationError):
                VERIFIER.verify_manifest(root)

    def test_missing_artifact_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(VERIFIER.VerificationError):
                VERIFIER.verify_artifact(Path(directory))

    def test_basis_capture_patch_matches_current_discovery_program(self) -> None:
        fake_sage = types.ModuleType("sage")
        fake_all = types.ModuleType("sage.all")
        for name in ("Matrix", "QQ", "ZZ", "identity_matrix"):
            setattr(fake_all, name, object())
        saved_sage = sys.modules.get("sage")
        saved_sage_all = sys.modules.get("sage.all")
        sys.modules["sage"] = fake_sage
        sys.modules["sage.all"] = fake_all
        sys.path.insert(0, str(ROOT / "research"))
        try:
            map_script = ROOT / "research" / "replay_rank17_neighbor_chain_with_maps.py"
            spec = importlib.util.spec_from_file_location("rank17_map_replay", map_script)
            assert spec is not None and spec.loader is not None
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            source = ROOT / ".github" / "workflows" / "probe-rank29-height-lattice.yml"
            base = module.patch_script(
                module.extract_discovery_script(source), Path("/tmp/test-rank17-map-patch")
            )
            patched = module.patch_basis_capture(base)
            self.assertIn("return_matrix=True", patched)
            self.assertIn("'target_neighbor_maps': target_neighbor_maps", patched)
            self.assertIn("'parent_transform': parent_transform", patched)
        finally:
            sys.path.pop(0)
            if saved_sage is None:
                sys.modules.pop("sage", None)
            else:
                sys.modules["sage"] = saved_sage
            if saved_sage_all is None:
                sys.modules.pop("sage.all", None)
            else:
                sys.modules["sage.all"] = saved_sage_all

    def test_manual_workflow_has_dispatch_and_checkpoint_trigger(self) -> None:
        workflow = (
            ROOT / ".github" / "workflows" / "replay-rank17-neighbor-chain-manual.yml"
        ).read_text()
        self.assertIn("workflow_dispatch:", workflow)
        self.assertIn("\n  push:", workflow)
        self.assertIn("branches: [codex/rank31-research]", workflow)
        self.assertIn("github.actor != 'github-actions[bot]'", workflow)
        self.assertIn("sagemath/sagemath:10.9@sha256:", workflow)
        self.assertIn("retention-days: 90", workflow)
        self.assertIn("--write-manifest", workflow)
        self.assertIn("artifact-verifier.py", workflow)
        self.assertIn("exact-neighbor-maps.json", workflow)


if __name__ == "__main__":
    unittest.main()
