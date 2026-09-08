from __future__ import annotations

import ast
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "research/search_e8_a2_target_neighbor_bridge.py"
WORKFLOW = ROOT / ".github/workflows/probe-e8-a2-target-neighbor-bridge.yml"


class FakeGram:
    def __init__(self, label: str, rank: int = 2) -> None:
        self.label = label
        self.rank = rank

    def nrows(self) -> int:
        return self.rank

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, FakeGram)
            and self.label == other.label
            and self.rank == other.rank
        )


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
        self.assertIn("math.isqrt", SCRIPT.read_text())

    def test_only_exact_expected_sage_p2_value_error_is_classified(self) -> None:
        expected = ValueError(self.module.EXPECTED_P2_CONSTRUCTION_ERROR)
        record = self.module.classify_expected_neighbor_construction_error(2, expected)
        self.assertIsNotNone(record)
        assert record is not None
        self.assertEqual(
            record["error_category"],
            "expected_sage_p2_nonmaximal_even_boundary",
        )
        self.assertFalse(record["prime_neighbor_enumeration_complete"])
        self.assertIn("2 divides the determinant", record["enumeration_boundary"])
        self.assertIsNone(
            self.module.classify_expected_neighbor_construction_error(5, expected)
        )
        self.assertIsNone(
            self.module.classify_expected_neighbor_construction_error(
                2, ValueError(self.module.EXPECTED_P2_CONSTRUCTION_ERROR + ".")
            )
        )

        def raise_expected():
            raise ValueError(self.module.EXPECTED_P2_CONSTRUCTION_ERROR)

        built, classified = self.module.call_neighbor_builder_fail_closed(
            2, raise_expected
        )
        self.assertIsNone(built)
        self.assertIsNotNone(classified)

        with self.assertRaisesRegex(ValueError, "unexpected construction failure"):
            self.module.call_neighbor_builder_fail_closed(
                2,
                lambda: (_ for _ in ()).throw(
                    ValueError("unexpected construction failure")
                ),
            )
        with self.assertRaisesRegex(AssertionError, "Gram assertion"):
            self.module.call_neighbor_builder_fail_closed(
                2,
                lambda: (_ for _ in ()).throw(AssertionError("Gram assertion")),
            )

        source = SCRIPT.read_text()
        self.assertIn("except ValueError as exc:", source)
        self.assertIn("if classification is None:\n            raise", source)

    def test_output_directory_rejects_stale_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            new_output = base / "new"
            self.module.ensure_fresh_output_directory(new_output)
            self.assertTrue(new_output.is_dir())

            empty_output = base / "empty"
            empty_output.mkdir()
            self.module.ensure_fresh_output_directory(empty_output)

            stale_output = base / "stale"
            stale_output.mkdir()
            (stale_output / "exact-bridge.json").write_text("{}\n")
            with self.assertRaises(self.module.OutputDirectoryError):
                self.module.ensure_fresh_output_directory(stale_output)

    @staticmethod
    def _state(side: str, uid: str, gram_hash: str, gram: FakeGram):
        return {
            "side": side,
            "uid": uid,
            "gram_sha256": gram_hash,
            "gram_object": gram,
            "checks": {
                "root_count_exact_theta": 10,
                "norm4_count_exact_theta": 20,
            },
        }

    def test_identical_gram_is_promoted_outside_zero_qfisom_budget(self) -> None:
        gram_hash = "a" * 64
        origin = self._state("origin", "origin:a", gram_hash, FakeGram("same"))
        endpoint = self._state(
            "transparent", "transparent:a", gram_hash, FakeGram("same")
        )
        indexes = {
            "origin": {(10, 20): [origin]},
            "transparent": {(10, 20): [endpoint]},
            "rootless": {},
        }
        stats = self.module.new_isometry_stats()
        with tempfile.TemporaryDirectory() as temporary:
            log = Path(temporary) / "isometry.jsonl"
            with mock.patch.object(
                self.module,
                "exact_isometry",
                side_effect=AssertionError("qfisom must not run"),
            ):
                meeting = self.module.find_meeting(
                    endpoint,
                    indexes,
                    ("transparent", "rootless"),
                    log,
                    stats,
                    0,
                )
            self.assertIsNotNone(meeting)
            record = json.loads(log.read_text().strip())
        self.assertEqual(record["comparison_method"], "identical_serialized_gram")
        self.assertEqual(stats["pairs_seen"], 1)
        self.assertEqual(stats["qfisom_checked"], 0)
        self.assertEqual(stats["identical_hash_matches"], 1)
        self.assertFalse(stats["budget_exhausted"])

    def test_first_unchecked_qfisom_pair_terminates_inconclusively(self) -> None:
        origin = self._state("origin", "origin:a", "a" * 64, FakeGram("origin"))
        endpoint = self._state(
            "transparent", "transparent:b", "b" * 64, FakeGram("endpoint")
        )
        indexes = {
            "origin": {(10, 20): [origin]},
            "transparent": {(10, 20): [endpoint]},
            "rootless": {},
        }
        stats = self.module.new_isometry_stats()
        with tempfile.TemporaryDirectory() as temporary:
            log = Path(temporary) / "isometry.jsonl"
            with mock.patch.object(
                self.module,
                "exact_isometry",
                side_effect=AssertionError("qfisom must not run beyond budget"),
            ):
                meeting = self.module.find_meeting(
                    endpoint,
                    indexes,
                    ("transparent", "rootless"),
                    log,
                    stats,
                    0,
                )
            record = json.loads(log.read_text().strip())
        self.assertIsNone(meeting)
        self.assertEqual(record["comparison_method"], "qfisom_skipped_budget_exhausted")
        self.assertEqual(stats["pairs_seen"], 1)
        self.assertEqual(stats["qfisom_checked"], 0)
        self.assertEqual(stats["skipped_due_to_budget"], 1)
        self.assertTrue(stats["budget_exhausted"])

    def test_beam_targets_dynamically_close_opposing_profiles(self) -> None:
        def state(uid: str, roots: int, norm4: int):
            return {
                "uid": uid,
                "gram_sha256": uid.rjust(64, "0"),
                "checks": {
                    "root_count_exact_theta": roots,
                    "norm4_count_exact_theta": norm4,
                },
            }

        close = state("1", 54, 2726)
        same_roots_but_far = state("2", 58, 3100)
        far = state("3", 84, 3376)
        selected = self.module.select_next_beam(
            [far, same_roots_but_far, close], 2, {"transparent": [(58, 2734)]}
        )
        self.assertEqual(selected[0]["uid"], close["uid"])
        self.assertEqual(
            self.module.fingerprint_distance((54, 2726), (58, 2734)),
            (12, 4, 8),
        )

    def test_origin_beam_has_endpoint_aware_exploitation_quotas(self) -> None:
        def state(uid: str, roots: int, norm4: int):
            return {
                "uid": uid,
                "gram_sha256": uid.encode().hex().rjust(64, "0"),
                "checks": {
                    "root_count_exact_theta": roots,
                    "norm4_count_exact_theta": norm4,
                },
            }

        transparent_states = [state(f"t{index}", index, index) for index in range(3)]
        transparent_duplicate = state("t-duplicate", 0, 0)
        rootless_states = [
            state(f"r{index}", 100 - index, 100 - index) for index in range(3)
        ]
        rootless_duplicate = state("r-duplicate", 100, 100)
        diversity = [state("d0", 45, 45), state("d1", 55, 55)]
        selected = self.module.select_next_beam(
            transparent_states
            + [transparent_duplicate]
            + rootless_states
            + [rootless_duplicate]
            + diversity,
            8,
            {
                "transparent": [(0, 0)],
                "rootless": [(100, 100)],
            },
        )
        self.assertEqual(
            self.module.exploitation_quotas(8, ["transparent", "rootless"]),
            {"transparent": 3, "rootless": 3},
        )
        self.assertTrue(all(item["uid"].startswith("t") for item in selected[:3]))
        self.assertTrue(all(item["uid"].startswith("r") for item in selected[3:6]))
        self.assertEqual(
            {self.module.state_fingerprint(item) for item in selected[:3]},
            {(0, 0), (1, 1), (2, 2)},
        )
        self.assertEqual(
            {self.module.state_fingerprint(item) for item in selected[3:6]},
            {(100, 100), (99, 99), (98, 98)},
        )
        self.assertEqual(
            len({self.module.state_fingerprint(item) for item in selected[:6]}),
            6,
        )
        self.assertEqual(len({item["uid"] for item in selected}), 8)

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
            "for opposing_side in opposing_sides",
            "origin_composite * meeting_transform * endpoint_composite.inverse()",
            "end_to_end.transpose() * origin_initial_gram * end_to_end",
            '"endpoint_initial_to_origin_initial_transform_M_sha256"',
            '"termination_reason"',
        ]
        for fragment in required_fragments:
            self.assertIn(fragment, source)

    def test_workflow_is_branch_scoped_manual_pinned_and_uploads_even_on_failure(self) -> None:
        source = WORKFLOW.read_text()
        self.assertIn("workflow_dispatch:", source)
        self.assertIn("\n  push:", source)
        self.assertIn("- codex/rank32-research", source)
        self.assertIn('default: "20000"', source)
        self.assertIn('echo "PRIMES=${INPUT_PRIMES:-2}"', source)
        self.assertIn('echo "ROUNDS=${INPUT_ROUNDS:-8}"', source)
        self.assertIn("rank32-20260908-p2-v1", source)
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
        self.assertIn("cd /tmp/e8-a2-target-neighbor-bridge", source)
        self.assertIn("retention-days: 90", source)
        self.assertIn("timeout --foreground --signal=TERM --kill-after=60s 350m", source)
        self.assertIn("SEARCH_CONTAINER_DIGEST:", source)
        self.assertIn("record['claim_boundary']['rank32_curve_found'] is False", source)

    def test_runtime_provenance_builds_github_run_url_without_secrets(self) -> None:
        provenance = self.module.runtime_provenance(
            {
                "GITHUB_SERVER_URL": "https://github.example",
                "GITHUB_REPOSITORY": "owner/repo",
                "GITHUB_RUN_ID": "123",
                "GITHUB_SHA": "abc",
                "UNRELATED_SECRET": "must-not-appear",
            }
        )
        self.assertEqual(
            provenance["github_run_url"],
            "https://github.example/owner/repo/actions/runs/123",
        )
        self.assertNotIn("unrelated_secret", provenance)

    def test_script_digest_is_not_accidentally_empty(self) -> None:
        digest = hashlib.sha256(SCRIPT.read_bytes()).hexdigest()
        self.assertEqual(len(digest), 64)
        self.assertNotEqual(digest, hashlib.sha256(b"").hexdigest())


if __name__ == "__main__":
    unittest.main()
