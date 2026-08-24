from __future__ import annotations

import copy
import importlib.util
import json
import shutil
import subprocess
import tempfile
import types
import unittest
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
IMPORTER_PATH = ROOT / "research" / "import_e8_a2_target_neighbor_bridge_artifact.py"
SEARCH_PATH = ROOT / "research" / "search_e8_a2_target_neighbor_bridge.py"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


IMPORTER = load_module("target_neighbor_artifact_importer", IMPORTER_PATH)
SEARCH = load_module("target_neighbor_search_fixture_source", SEARCH_PATH)


def write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_jsonl(path: Path, records: list[dict[str, object]]) -> None:
    path.write_text(
        "".join(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n" for record in records),
        encoding="utf-8",
    )


def checks(roots: int, norm4: int) -> dict[str, object]:
    return {
        "rank": 17,
        "determinant": 948,
        "even": True,
        "positive_definite": True,
        "target_genus_verified": True,
        "root_count_exact_theta": roots,
        "norm4_count_exact_theta": norm4,
    }


def state(
    side: str,
    gram: list[list[int]],
    source: str | None,
    round_discovered: int,
    parent_uid: str | None,
    roots: int,
    norm4: int,
) -> dict[str, object]:
    digest = IMPORTER.payload_sha256(gram)
    return {
        "uid": f"{side}:{digest}",
        "side": side,
        "source": source,
        "round_discovered": round_discovered,
        "gram_sha256": digest,
        "parent_uid": parent_uid,
        "incoming_move_sha256": None,
        "checks": checks(roots, norm4),
        "gram": gram,
    }


def conjugate_first_basis_vector(gram: list[list[int]]) -> list[list[int]]:
    return [
        [value * (-1 if row == 0 else 1) * (-1 if column == 0 else 1) for column, value in enumerate(line)]
        for row, line in enumerate(gram)
    ]


def source_files() -> dict[str, str]:
    return {
        "search_script": IMPORTER.SEARCH_SCRIPT_PATH,
        "search_script_sha256": IMPORTER.file_sha256(ROOT / IMPORTER.SEARCH_SCRIPT_PATH),
        "target_formula_script": IMPORTER.TARGET_SCRIPT_PATH,
        "target_formula_script_sha256": IMPORTER.file_sha256(ROOT / IMPORTER.TARGET_SCRIPT_PATH),
    }


def claim_boundary(found: bool, modern: bool = True) -> dict[str, object]:
    result: dict[str, object] = {
        "exact_lattice_bridge_found": found,
        "explicit_stable_NS_transport_completed": False,
        "explicit_K3_fibration_switch_completed": False,
        "explicit_moduli_map_completed": False,
        "rational_P3_found": False,
        "rank31_curve_found": False,
        "negative_result_scope": (
            "Only successfully constructed neighbors in the serialized finite window "
            "were excluded. Classified p=2 boundary-failure lines were attempted but "
            "are not excluded; this artifact does not settle any elliptic-rank record problem."
        ),
    }
    if modern:
        result.update(
            {
                "rank32_curve_found": False,
                "p2_neighbor_enumeration_complete": False,
                "p2_neighbor_enumeration_boundary": IMPORTER.P2_ENUMERATION_BOUNDARY,
                "bounded_negative_result_importable": True,
            }
        )
    return result


def initial_states() -> list[dict[str, object]]:
    records = [
        state("origin", SEARCH.target_essential_lattice(), "target", 0, None, 258, 4000),
        state("transparent", SEARCH.transparent_seed_lattice(), "seed", 0, None, 150, 5424),
        state("rootless", copy.deepcopy(SEARCH.FROZEN_ROOTLESS_GRAM), "rootless", 0, None, 0, 2622),
    ]
    assert {record["side"]: record["gram_sha256"] for record in records} == IMPORTER.INITIAL_HASHES
    return records


def make_move(parent: dict[str, object], child: dict[str, object], attempt: dict[str, object]) -> dict[str, object]:
    prime = 2
    move: dict[str, object] = {
        "parent_uid": parent["uid"],
        "parent_gram_sha256": parent["gram_sha256"],
        "prime": prime,
        "projective_line_ordinal_one_based": 1,
        "projective_isotropic_vector": [1] + [0] * 16,
        "q_of_vector": 2,
        "q_of_vector_mod_p": 0,
        "raw_neighbor_gram": child["gram"],
        "raw_neighbor_transform": [],
        "lll_unimodular_transform": [],
        "parent_to_child_rational_transform": [],
        "child_gram_sha256": child["gram_sha256"],
        "child_gram": child["gram"],
        "p_neighbor_relation": {
            "gram_identity_verified": True,
            "transition_determinant": "-1",
            "transition_denominator": prime,
            "inverse_transition_denominator": prime,
            "rank_mod_p_of_p_times_transition": 1,
            "rank_mod_p_of_p_times_inverse": 1,
            "intersection_index_in_parent": prime,
            "intersection_index_in_child": prime,
        },
        "child_exact_checks": child["checks"],
    }
    move["move_sha256"] = IMPORTER.payload_sha256(move)
    return move


def exact_raw_two_neighbor_fixture() -> tuple[list[list[int]], dict[str, object]]:
    """A sparse exact p=2 move frozen from the target lattice's ordinal-2 line."""

    parent = SEARCH.target_essential_lattice()
    transform = [
        [Fraction(int(row == column)) for column in range(17)]
        for row in range(17)
    ]
    transform[0][0] = Fraction(1, 2)
    transform[2][0] = Fraction(1, 2)
    transform[16][0] = Fraction(1, 2)
    transform[7][3] = Fraction(1)
    transform[7][7] = Fraction(2)
    child_fraction = [
        [
            sum(
                transform[left][row]
                * parent[left][right]
                * transform[right][column]
                for left in range(17)
                for right in range(17)
            )
            for column in range(17)
        ]
        for row in range(17)
    ]
    assert all(value.denominator == 1 for row in child_fraction for value in row)
    child = [[int(value) for value in row] for row in child_fraction]
    identity = [[int(row == column) for column in range(17)] for row in range(17)]
    vector = [0] * 17
    vector[0] = vector[2] = vector[16] = 1
    move: dict[str, object] = {
        "parent_uid": f"origin:{IMPORTER.payload_sha256(parent)}",
        "parent_gram_sha256": IMPORTER.payload_sha256(parent),
        "prime": 2,
        "projective_line_ordinal_one_based": 2,
        "projective_isotropic_vector": vector,
        "q_of_vector": 4,
        "q_of_vector_mod_p": 0,
        "raw_neighbor_gram": child,
        "raw_neighbor_transform": IMPORTER.serial_fraction_matrix(transform),
        "lll_unimodular_transform": identity,
        "parent_to_child_rational_transform": IMPORTER.serial_fraction_matrix(transform),
        "child_gram_sha256": IMPORTER.payload_sha256(child),
        "child_gram": child,
        "p_neighbor_relation": {
            "gram_identity_verified": True,
            "transition_determinant": "1",
            "transition_denominator": 2,
            "inverse_transition_denominator": 2,
            "rank_mod_p_of_p_times_transition": 1,
            "rank_mod_p_of_p_times_inverse": 1,
            "intersection_index_in_parent": 2,
            "intersection_index_in_child": 2,
        },
        "child_exact_checks": checks(162, 4758),
    }
    move["move_sha256"] = IMPORTER.payload_sha256(move)
    return parent, move


def refresh_manifest(root: Path) -> None:
    names = sorted(path.name for path in root.iterdir() if path.name != IMPORTER.MANIFEST_NAME)
    (root / IMPORTER.MANIFEST_NAME).write_text(
        "".join(f"{IMPORTER.file_sha256(root / name)}  ./{name}\n" for name in names),
        encoding="utf-8",
    )


def refresh_result_after_evidence_change(root: Path, name: str) -> None:
    result = json.loads((root / "result.json").read_text())
    result["evidence_sha256"][name] = IMPORTER.file_sha256(root / name)
    result.pop("result_payload_sha256", None)
    result["result_payload_sha256"] = IMPORTER.payload_sha256(result)
    write_json(root / "result.json", result)
    refresh_manifest(root)


def write_rehashed_result(root: Path, result: dict[str, object]) -> None:
    result.pop("result_payload_sha256", None)
    result["result_payload_sha256"] = IMPORTER.payload_sha256(result)
    write_json(root / "result.json", result)
    refresh_manifest(root)


def write_rehashed_bridge(root: Path, bridge: dict[str, object]) -> None:
    transport = bridge.get("end_to_end_transport")
    if isinstance(transport, dict):
        transport.pop("transport_sha256", None)
        transport["transport_sha256"] = IMPORTER.payload_sha256(transport)
    bridge.pop("bridge_sha256", None)
    bridge["bridge_sha256"] = IMPORTER.payload_sha256(bridge)
    write_json(root / "exact-bridge.json", bridge)
    result = json.loads((root / "result.json").read_text())
    result["bridge"] = bridge
    result["evidence_sha256"]["exact-bridge.json"] = IMPORTER.file_sha256(
        root / "exact-bridge.json"
    )
    write_rehashed_result(root, result)


def provenance() -> dict[str, object]:
    return {
        "repository": "metaforismo/elliptic-rank-31",
        "workflow": IMPORTER.WORKFLOW_PATH,
        "event_name": "workflow_dispatch",
        "git_ref": "refs/heads/codex/rank31-research",
        "commit_sha": "a" * 40,
        "run_id": "31393796753",
        "run_attempt": 2,
        "run_url": "https://github.com/metaforismo/elliptic-rank-31/actions/runs/31393796753",
        "artifact_name": "e8-a2-target-exact-neighbor-bridge",
        "artifact_id": "9064805239",
        "artifact_digest": "sha256:" + "b" * 64,
        "downloaded_archive_sha256": "b" * 64,
        "extra": {"download_tool": "gh-2.80.0", "operator_note": "kept exactly"},
    }


def embedded_runtime(supplied: dict[str, object]) -> dict[str, str]:
    return {
        "github_actions": "true",
        "github_server_url": "https://github.com",
        "github_repository": str(supplied["repository"]),
        "github_run_id": str(supplied["run_id"]),
        "github_run_attempt": str(supplied["run_attempt"]),
        "github_sha": str(supplied["commit_sha"]),
        "github_ref": str(supplied["git_ref"]),
        "github_ref_name": "codex/rank31-research",
        "github_workflow": IMPORTER.EXPECTED_WORKFLOW_NAME,
        "github_workflow_ref": (
            f'{supplied["repository"]}/{IMPORTER.WORKFLOW_PATH}@{supplied["git_ref"]}'
        ),
        "github_job": IMPORTER.EXPECTED_JOB_NAME,
        "runner_os": "Linux",
        "runner_arch": "X64",
        "search_container_image": IMPORTER.EXPECTED_CONTAINER_IMAGE,
        "search_container_digest": IMPORTER.EXPECTED_CONTAINER_DIGEST,
        "github_run_url": str(supplied["run_url"]),
    }


def build_fixture(
    root: Path,
    found: bool,
    *,
    modern: bool = True,
    supplied: dict[str, object] | None = None,
) -> str | None:
    supplied = provenance() if supplied is None else supplied
    initial = initial_states()
    states = list(initial)
    moves: list[dict[str, object]] = []
    isometries: list[dict[str, object]] = []
    attempt_base = {
        "attempt_index": 1,
        "round": 1,
        "side": "origin",
        "parent_uid": initial[0]["uid"],
        "prime": 2,
        "projective_line_ordinal_one_based": 1,
        "projective_isotropic_vector": [1] + [0] * 16,
    }
    bridge = None
    initial_hash_override = None
    if found:
        # A zero-edge synthetic bridge exercises the independent meeting-
        # isometry and end-to-end composition gates without fabricating a
        # p-neighbor. The test temporarily pins this synthetic endpoint hash.
        endpoint = copy.deepcopy(initial[0])
        endpoint["side"] = "rootless"
        endpoint["source"] = "synthetic identical endpoint"
        endpoint["uid"] = f'rootless:{endpoint["gram_sha256"]}'
        initial[2] = endpoint
        states = list(initial)
        initial_hash_override = endpoint["gram_sha256"]  # type: ignore[assignment]
        attempts = []
        identity = [
            [int(row == column) for column in range(17)]
            for row in range(17)
        ]
        meeting_isometry = {
            "from_uid": initial[0]["uid"],
            "to_uid": endpoint["uid"],
            "integral_unimodular_transform": identity,
            "determinant": 1,
            "gram_identity_verified": True,
            "fingerprint": [258, 4000],
        }
        if modern:
            meeting_isometry["promotion_method"] = (
                "identical_serialized_gram_sha256_and_matrix"
            )
            isometries = [
                {
                    "pair_index": 1,
                    "qfisom_check_index": None,
                    "check_index": None,
                    "origin_uid": initial[0]["uid"],
                    "endpoint_uid": endpoint["uid"],
                    "endpoint_side": "rootless",
                    "fingerprint": [258, 4000],
                    "comparison_method": "identical_serialized_gram",
                    "isometric": True,
                    "isometry": meeting_isometry,
                }
            ]
        else:
            isometries = [
                {
                    "check_index": 1,
                    "origin_uid": initial[0]["uid"],
                    "endpoint_uid": endpoint["uid"],
                    "endpoint_side": "rootless",
                    "fingerprint": [258, 4000],
                    "isometric": True,
                    "isometry": meeting_isometry,
                }
            ]
        bridge = {
            "schema_version": 1,
            "status": "exact_neighbor_bridge_found",
            "endpoint": "rootless",
            "orientation": "synthetic exact fixture",
            "origin_chain": [initial[0]],
            "origin_forward_moves": [],
            "meeting_isometry": meeting_isometry,
            "endpoint_chain_endpoint_to_meeting": [endpoint],
            "endpoint_forward_moves_to_invert": [],
            "neighbor_step_count": 0,
            "all_parent_child_gram_identities_verified": True,
            "all_neighbor_intersection_indices_verified": True,
            "endpoint_identification_verified_by_initial_source": True,
            "claim_boundary": (
                "This is an exact integral-lattice p-neighbor bridge. It is not yet an "
                "explicit K3 map, rational P3 section, or elliptic rank-record curve."
            ),
        }
        if modern:
            identity_fraction = IMPORTER.serial_fraction_matrix(
                IMPORTER.identity_matrix()
            )
            transport_core = {
                "matrix_convention": (
                    "Every edge T satisfies T^t G_parent T = G_child.  "
                    "M = C_origin Q C_endpoint^-1 maps endpoint initial "
                    "coordinates into origin initial coordinates."
                ),
                "origin_forward_composite_C_origin": identity_fraction,
                "endpoint_forward_composite_C_endpoint": identity_fraction,
                "endpoint_initial_to_origin_initial_transform_M": identity_fraction,
                "endpoint_initial_to_origin_initial_transform_M_sha256": (
                    IMPORTER.payload_sha256(identity_fraction)
                ),
                "M_determinant": "1",
                "M_determinant_is_plus_or_minus_one": True,
                "M_transpose_G_origin_initial_M_equals_G_endpoint_initial": True,
            }
            bridge["end_to_end_transport"] = {
                **transport_core,
                "transport_sha256": IMPORTER.payload_sha256(transport_core),
            }
        bridge["bridge_sha256"] = IMPORTER.payload_sha256(bridge)
        write_json(root / "exact-bridge.json", bridge)
        side_stats = {
            "parents": 0,
            "attempted_lines": 0,
            "successful_moves": 0,
            "failed_moves": 0,
            "new_unique_states": 0,
            "duplicate_child_grams": 0,
            "next_beam": [],
        }
        counters = {
            "attempted_lines": 0,
            "successful_moves": 0,
            "failed_moves": 0,
            "duplicate_child_grams": 0,
            "unique_states_discovered": 3,
            "move_budget_exhausted": False,
        }
        if modern:
            counters["expected_p2_construction_boundary_failures"] = 0
    else:
        if modern:
            attempts = [
                {
                    **attempt_base,
                    "status": "construction_error",
                    "construction_error_expected_and_bounded": True,
                    "error_category": "expected_sage_p2_nonmaximal_even_boundary",
                    "error_type": "ValueError",
                    "error_message": IMPORTER.EXPECTED_P2_ERROR,
                    "prime_neighbor_enumeration_complete": False,
                    "enumeration_boundary": IMPORTER.P2_ENUMERATION_BOUNDARY,
                }
            ]
        else:
            attempts = [
                {
                    **attempt_base,
                    "status": "construction_error",
                    "error": IMPORTER.EXPECTED_P2_ERROR_REPR,
                }
            ]
        side_stats = {
            "parents": 1,
            "attempted_lines": 1,
            "successful_moves": 0,
            "failed_moves": 1,
            "new_unique_states": 0,
            "duplicate_child_grams": 0,
            "next_beam": [],
        }
        counters = {
            "attempted_lines": 1,
            "successful_moves": 0,
            "failed_moves": 1,
            "duplicate_child_grams": 0,
            "unique_states_discovered": 3,
            "move_budget_exhausted": False,
        }
        if modern:
            counters["expected_p2_construction_boundary_failures"] = 1

    qfisom_checks = 0 if modern else len(isometries)
    isometry_stats = {
        "pairs_seen": len(isometries),
        "qfisom_checked": qfisom_checks,
        "identical_hash_matches": 1 if modern and found else 0,
        "skipped_due_to_budget": 0,
        "budget_exhausted": False,
    }
    round_record: dict[str, object] = {
        "round": 1,
        "sides": {"origin": side_stats},
        "cumulative": counters,
        "isometry_checks": qfisom_checks,
        "bridge_found": found,
    }
    if modern:
        round_record["isometry_comparisons"] = copy.deepcopy(isometry_stats)
    rounds = [
        {
            **round_record,
        }
    ]
    write_jsonl(root / "states.jsonl", states)
    write_jsonl(root / "moves.jsonl", moves)
    write_jsonl(root / "attempts.jsonl", attempts)
    write_jsonl(root / "rounds.jsonl", rounds)
    write_jsonl(root / "isometry-checks.jsonl", isometries)
    checkpoint: dict[str, object] = {
        "status": "bridge_found" if found else "searching",
        "rounds_completed": 1,
        "counters": counters,
        "isometry_checks": qfisom_checks,
        "current_beams": {"origin": []} if modern else {},
    }
    if modern:
        checkpoint["isometry_comparisons"] = copy.deepcopy(isometry_stats)
    write_json(root / "checkpoint.json", checkpoint)
    (root / "run.log").write_text("synthetic Sage run\n", encoding="utf-8")
    evidence_names = set(IMPORTER.RESULT_EVIDENCE_FILES)
    if found:
        evidence_names.add("exact-bridge.json")
    configuration: dict[str, object] = {
        "mode": "forward",
        "rounds_requested": 1,
        "rounds_completed": 1,
        "beam_size": 1,
        "primes": [2],
        "projective_line_offset": 0,
        "projective_lines_per_state_prime": 1,
        "max_successful_moves": 10,
        "max_isometry_checks": 10,
        "enumeration_order": "synthetic",
        "beam_selection": "synthetic",
    }
    if modern:
        configuration.update(
            {
                "side_expansion_order": ["origin"],
                "sequential_order_semantics": "synthetic sequential fixture",
                "prime_enumeration_boundaries": [
                    {
                        "prime": 2,
                        "complete": False,
                        "reason": IMPORTER.P2_ENUMERATION_BOUNDARY,
                    }
                ],
            }
        )
    result: dict[str, object] = {
        "schema_version": 1,
        "status": "bridge_found" if found else "bounded_search_completed",
        "truth_status": "Synthetic exact bridge." if found else "Synthetic finite negative result only.",
        "sage_version": "SageMath version 10.9, Release Date: synthetic",
        "configuration": configuration,
        "initial_states": initial,
        "target_genus": "synthetic determinant-948 genus",
        "counters": counters,
        "isometry_checks": qfisom_checks,
        "bridge": bridge,
        "claim_boundary": claim_boundary(found, modern),
        "source_files": source_files(),
        "evidence_sha256": {
            name: IMPORTER.file_sha256(root / name) for name in evidence_names
        },
    }
    if modern:
        result.update(
            {
                "termination_reason": (
                    "exact_bridge_found" if found else "requested_rounds_completed"
                ),
                "artifact_importable": True,
                "runtime_provenance": embedded_runtime(supplied),
                "isometry_comparisons": copy.deepcopy(isometry_stats),
            }
        )
    result["result_payload_sha256"] = IMPORTER.payload_sha256(result)
    write_json(root / "result.json", result)
    refresh_manifest(root)
    return initial_hash_override


def verify_fixture(
    root: Path,
    supplied: dict[str, object] | None = None,
    source_root: Path | None = None,
):
    supplied = provenance() if supplied is None else supplied
    result = json.loads((root / "result.json").read_text())
    rootless_hash = next(
        state["gram_sha256"]
        for state in result["initial_states"]
        if state["side"] == "rootless"
    )
    original_hashes = dict(IMPORTER.INITIAL_HASHES)
    IMPORTER.INITIAL_HASHES["rootless"] = rootless_hash
    try:
        return IMPORTER.verify_artifact(
            root, supplied, source_root=source_root
        )
    finally:
        IMPORTER.INITIAL_HASHES.clear()
        IMPORTER.INITIAL_HASHES.update(original_hashes)


class TargetNeighborArtifactImporterTests(unittest.TestCase):
    def test_expected_p2_compatibility_flag_is_counted_and_must_be_exact(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "artifact"
            root.mkdir()
            build_fixture(root, False)
            certificate = IMPORTER.verify_artifact(root, provenance())
            self.assertEqual(
                certificate["search"]["structured_expected_p2_boundary_failures"], 1
            )
            self.assertEqual(
                certificate["search"]["total_classified_p2_boundary_failures"], 1
            )

            attempts = [
                json.loads(line)
                for line in (root / "attempts.jsonl").read_text().splitlines()
            ]
            attempts[0]["error_message"] = "unexpected drift"
            write_jsonl(root / "attempts.jsonl", attempts)
            refresh_result_after_evidence_change(root, "attempts.jsonl")
            with self.assertRaises(IMPORTER.VerificationError):
                IMPORTER.verify_artifact(root, provenance())

    def test_modern_p2_type_and_enumeration_boundary_are_exact(self) -> None:
        for field, replacement in (
            ("error_type", "RuntimeError"),
            ("enumeration_boundary", IMPORTER.P2_ENUMERATION_BOUNDARY + " drift"),
            ("error_category", "wrong-category"),
        ):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as directory:
                root = Path(directory) / "artifact"
                root.mkdir()
                build_fixture(root, False)
                attempts = [
                    json.loads(line)
                    for line in (root / "attempts.jsonl").read_text().splitlines()
                ]
                attempts[0][field] = replacement
                write_jsonl(root / "attempts.jsonl", attempts)
                refresh_result_after_evidence_change(root, "attempts.jsonl")
                with self.assertRaises(IMPORTER.VerificationError):
                    IMPORTER.verify_artifact(root, provenance())

    def test_modern_p2_counter_and_claim_boundary_are_mandatory(self) -> None:
        for mutation in (
            "missing_counter",
            "missing_claim_boundary",
            "claim_boundary_drift",
        ):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as directory:
                root = Path(directory) / "artifact"
                root.mkdir()
                build_fixture(root, False)
                result = json.loads((root / "result.json").read_text())
                if mutation == "missing_counter":
                    result["counters"].pop(
                        "expected_p2_construction_boundary_failures"
                    )
                elif mutation == "missing_claim_boundary":
                    result["claim_boundary"].pop(
                        "p2_neighbor_enumeration_boundary"
                    )
                else:
                    result["claim_boundary"][
                        "p2_neighbor_enumeration_boundary"
                    ] += " drift"
                write_rehashed_result(root, result)
                with self.assertRaises(IMPORTER.VerificationError):
                    IMPORTER.verify_artifact(root, provenance())

    def test_legacy_exact_p2_repr_and_unclassified_errors_are_separate(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "artifact"
            root.mkdir()
            build_fixture(root, False, modern=False)
            certificate = IMPORTER.verify_artifact(root, provenance())
            self.assertEqual(certificate["search"]["legacy_exact_p2_boundary_failures"], 1)
            self.assertEqual(certificate["search"]["total_classified_p2_boundary_failures"], 1)
            self.assertEqual(certificate["search"]["legacy_unclassified_construction_errors"], 0)
            self.assertIn(
                "not excluded by the bounded negative window",
                certificate["audit_conclusion"],
            )

            attempts = [
                json.loads(line)
                for line in (root / "attempts.jsonl").read_text().splitlines()
            ]
            attempts[0]["error"] = IMPORTER.EXPECTED_P2_ERROR_REPR + "."
            write_jsonl(root / "attempts.jsonl", attempts)
            refresh_result_after_evidence_change(root, "attempts.jsonl")
            certificate = IMPORTER.verify_artifact(root, provenance())
            self.assertEqual(certificate["search"]["legacy_exact_p2_boundary_failures"], 0)
            self.assertEqual(certificate["search"]["legacy_unclassified_construction_errors"], 1)
            self.assertEqual(
                certificate["search"]["artifact_schema_generation"], "legacy"
            )
            self.assertIn("unclassified construction error", certificate["audit_conclusion"])

    def test_exact_fraction_replay_accepts_real_sparse_two_neighbor_and_rejects_drift(self) -> None:
        parent, move = exact_raw_two_neighbor_fixture()
        exact = IMPORTER.verify_move_exact_arithmetic(move, parent, "sparse p=2 fixture")
        self.assertEqual(exact["move_sha256"], move["move_sha256"])

        tampered = copy.deepcopy(move)
        tampered["raw_neighbor_transform"][0][1] = "1"  # type: ignore[index]
        tampered.pop("move_sha256")
        tampered["move_sha256"] = IMPORTER.payload_sha256(tampered)
        with self.assertRaises(IMPORTER.VerificationError):
            IMPORTER.verify_move_exact_arithmetic(tampered, parent, "tampered p=2 fixture")

    def test_imports_bounded_and_found_results_without_rank31_promotion(self) -> None:
        for found in (False, True):
            with self.subTest(found=found), tempfile.TemporaryDirectory() as directory:
                root = Path(directory) / "artifact"
                root.mkdir()
                override = build_fixture(root, found)
                original_hashes = dict(IMPORTER.INITIAL_HASHES)
                if override is not None:
                    IMPORTER.INITIAL_HASHES["rootless"] = override
                supplied = provenance()
                try:
                    certificate = IMPORTER.verify_artifact(root, supplied)
                finally:
                    IMPORTER.INITIAL_HASHES.clear()
                    IMPORTER.INITIAL_HASHES.update(original_hashes)
                self.assertEqual(certificate["provenance"], supplied)
                self.assertEqual(certificate["search"]["status"], "bridge_found" if found else "bounded_search_completed")
                self.assertIs(certificate["claim_boundary"]["rank31_curve_found"], False)
                self.assertIs(certificate["claim_boundary"]["rank32_curve_found"], False)
                self.assertIn("no explicit k3", certificate["audit_conclusion"].lower())
                core = dict(certificate)
                digest = core.pop("record_sha256")
                self.assertEqual(digest, IMPORTER.payload_sha256(core))
                self.assertEqual(
                    certificate["search"]["jsonl"]["states.jsonl"]["record_count"],
                    3,
                )
                self.assertEqual(
                    certificate["verification_implementation"]["path"],
                    IMPORTER.IMPORTER_SCRIPT_PATH,
                )
                self.assertEqual(
                    certificate["verification_implementation"]["sha256"],
                    IMPORTER.file_sha256(IMPORTER_PATH),
                )
                self.assertEqual(
                    certificate["search"]["artifact_schema_generation"],
                    "modern",
                )

    def test_modern_identical_isometry_and_promoted_bridge_are_linked(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "artifact"
            root.mkdir()
            build_fixture(root, True)
            certificate = verify_fixture(root)
            self.assertEqual(certificate["search"]["isometry_checks"], 0)
            self.assertEqual(
                certificate["search"]["isometry_comparisons"]["identical_hash_matches"],
                1,
            )

            isometries = [
                json.loads(line)
                for line in (root / "isometry-checks.jsonl").read_text().splitlines()
            ]
            isometries[0]["isometry"]["untrusted_extra"] = True
            write_jsonl(root / "isometry-checks.jsonl", isometries)
            refresh_result_after_evidence_change(root, "isometry-checks.jsonl")
            with self.assertRaises(IMPORTER.VerificationError):
                verify_fixture(root)

    def test_obsolete_canonical_gram_method_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "artifact"
            root.mkdir()
            build_fixture(root, True)
            isometries = [
                json.loads(line)
                for line in (root / "isometry-checks.jsonl").read_text().splitlines()
            ]
            isometries[0]["comparison_method"] = "identical_canonical_gram"
            isometries[0]["isometry"]["promotion_method"] = (
                "identical_canonical_gram_sha256_and_matrix"
            )
            write_jsonl(root / "isometry-checks.jsonl", isometries)
            refresh_result_after_evidence_change(root, "isometry-checks.jsonl")
            with self.assertRaisesRegex(
                IMPORTER.VerificationError, "unknown comparison method"
            ):
                verify_fixture(root)

    def test_every_positive_qfisom_matrix_is_recomputed_exactly(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "artifact"
            root.mkdir()
            build_fixture(root, True)
            isometries = [
                json.loads(line)
                for line in (root / "isometry-checks.jsonl").read_text().splitlines()
            ]
            record = isometries[0]
            record["comparison_method"] = "pari_qfisom"
            record["qfisom_check_index"] = 1
            record["check_index"] = 1
            record["isometry"].pop("promotion_method")
            record["isometry"]["integral_unimodular_transform"][0][0] = -1
            record["isometry"]["determinant"] = -1
            write_jsonl(root / "isometry-checks.jsonl", isometries)
            refresh_result_after_evidence_change(root, "isometry-checks.jsonl")
            with self.assertRaisesRegex(
                IMPORTER.VerificationError, r"T\^t G_parent T != G_child"
            ):
                verify_fixture(root)

    def test_modern_bridge_transport_determinant_is_recomputed(self) -> None:
        for mutation in ("determinant", "missing_transport"):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as directory:
                root = Path(directory) / "artifact"
                root.mkdir()
                build_fixture(root, True)
                bridge = json.loads((root / "exact-bridge.json").read_text())
                if mutation == "determinant":
                    bridge["end_to_end_transport"]["M_determinant"] = "-1"
                else:
                    bridge.pop("end_to_end_transport")
                write_rehashed_bridge(root, bridge)
                with self.assertRaises(IMPORTER.VerificationError):
                    verify_fixture(root)

    def test_modern_round_stats_runtime_and_termination_tampering_fail(self) -> None:
        for mutation in (
            "round_stats",
            "missing_round_stats",
            "runtime",
            "termination",
        ):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as directory:
                root = Path(directory) / "artifact"
                root.mkdir()
                build_fixture(root, False)
                if mutation in {"round_stats", "missing_round_stats"}:
                    rounds = [
                        json.loads(line)
                        for line in (root / "rounds.jsonl").read_text().splitlines()
                    ]
                    if mutation == "round_stats":
                        rounds[0]["isometry_comparisons"]["pairs_seen"] = 1
                    else:
                        rounds[0].pop("isometry_comparisons")
                    write_jsonl(root / "rounds.jsonl", rounds)
                    refresh_result_after_evidence_change(root, "rounds.jsonl")
                else:
                    result = json.loads((root / "result.json").read_text())
                    if mutation == "runtime":
                        result["runtime_provenance"]["search_container_image"] = "untrusted"
                    else:
                        result["termination_reason"] = "successful_move_budget_exhausted"
                    write_rehashed_result(root, result)
                with self.assertRaises(IMPORTER.VerificationError):
                    IMPORTER.verify_artifact(root, provenance())

    def test_modern_result_fields_are_mandatory_without_schema_downgrade(self) -> None:
        for field in IMPORTER.MODERN_RESULT_FIELDS:
            with self.subTest(field=field), tempfile.TemporaryDirectory() as directory:
                root = Path(directory) / "artifact"
                root.mkdir()
                build_fixture(root, False)
                result = json.loads((root / "result.json").read_text())
                result.pop(field)
                write_rehashed_result(root, result)
                with self.assertRaisesRegex(
                    IMPORTER.VerificationError, "mixes legacy and modern schemas"
                ):
                    IMPORTER.verify_artifact(root, provenance())

    def test_inconclusive_artifact_is_never_imported(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "artifact"
            root.mkdir()
            build_fixture(root, False)
            result = json.loads((root / "result.json").read_text())
            result["status"] = "inconclusive_isometry_budget_exhausted"
            result["artifact_importable"] = False
            result["termination_reason"] = "first_required_qfisom_pair_exceeded_budget"
            write_rehashed_result(root, result)
            with self.assertRaises(IMPORTER.VerificationError):
                IMPORTER.verify_artifact(root, provenance())

    def test_duplicate_json_key_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "duplicate.json"
            path.write_text('{"status":"ok","status":"overwritten"}\n')
            with self.assertRaises(IMPORTER.VerificationError):
                IMPORTER.load_json(path)

    def test_manifest_rejects_absolute_and_unmanifested_paths(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "artifact"
            root.mkdir()
            build_fixture(root, False)
            lines = (root / IMPORTER.MANIFEST_NAME).read_text().splitlines()
            lines[0] = lines[0].split("  ", 1)[0] + "  /tmp/attempts.jsonl"
            (root / IMPORTER.MANIFEST_NAME).write_text("\n".join(lines) + "\n")
            with self.assertRaises(IMPORTER.VerificationError):
                IMPORTER.verify_artifact(root, provenance())

            refresh_manifest(root)
            (root / "unmanifested.txt").write_text("not evidence\n")
            with self.assertRaises(IMPORTER.VerificationError):
                IMPORTER.verify_artifact(root, provenance())

    def test_result_payload_hash_is_checked_after_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "artifact"
            root.mkdir()
            build_fixture(root, False)
            result = json.loads((root / "result.json").read_text())
            result["truth_status"] = "tampered after hashing"
            write_json(root / "result.json", result)
            refresh_manifest(root)
            with self.assertRaises(IMPORTER.VerificationError):
                IMPORTER.verify_artifact(root, provenance())

    def test_round_cumulative_mismatch_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "artifact"
            root.mkdir()
            build_fixture(root, False)
            rounds = [json.loads(line) for line in (root / "rounds.jsonl").read_text().splitlines()]
            rounds[0]["cumulative"]["attempted_lines"] = 2
            write_jsonl(root / "rounds.jsonl", rounds)
            refresh_result_after_evidence_change(root, "rounds.jsonl")
            with self.assertRaises(IMPORTER.VerificationError):
                IMPORTER.verify_artifact(root, provenance())

    def test_rank_record_claim_cannot_be_promoted_by_artifact(self) -> None:
        for gate in ("rank31_curve_found", "rank32_curve_found"):
            with self.subTest(gate=gate), tempfile.TemporaryDirectory() as directory:
                root = Path(directory) / "artifact"
                root.mkdir()
                build_fixture(root, False)
                result = json.loads((root / "result.json").read_text())
                result["claim_boundary"][gate] = True
                result.pop("result_payload_sha256")
                result["result_payload_sha256"] = IMPORTER.payload_sha256(result)
                write_json(root / "result.json", result)
                refresh_manifest(root)
                with self.assertRaises(IMPORTER.VerificationError):
                    IMPORTER.verify_artifact(root, provenance())

    def test_duplicate_jsonl_key_fails_even_with_updated_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "artifact"
            root.mkdir()
            build_fixture(root, False)
            (root / "attempts.jsonl").write_text(
                '{"attempt_index":1,"attempt_index":2,"status":"construction_error"}\n'
            )
            refresh_result_after_evidence_change(root, "attempts.jsonl")
            with self.assertRaises(IMPORTER.VerificationError):
                IMPORTER.verify_artifact(root, provenance())

    def test_cli_provenance_preserves_values_and_rejects_duplicate_extras(self) -> None:
        args = types.SimpleNamespace(
            repository="metaforismo/elliptic-rank-31",
            workflow=IMPORTER.WORKFLOW_PATH,
            event_name="push",
            git_ref="refs/heads/codex/rank31-research",
            commit_sha="d" * 40,
            run_id="42",
            run_attempt=3,
            run_url="https://github.com/metaforismo/elliptic-rank-31/actions/runs/42",
            artifact_name="artifact-name",
            artifact_id="99",
            artifact_digest="sha256:" + "e" * 64,
            archive_sha256="e" * 64,
            provenance=["download_tool=gh", "note=value=with=equals"],
        )
        record = IMPORTER.cli_provenance(args)
        self.assertEqual(record["run_id"], "42")
        self.assertEqual(record["run_attempt"], 3)
        self.assertEqual(record["extra"]["note"], "value=with=equals")
        args.archive_sha256 = "f" * 64
        with self.assertRaises(IMPORTER.VerificationError):
            IMPORTER.cli_provenance(args)
        args.archive_sha256 = "e" * 64
        args.provenance = ["note=first", "note=second"]
        with self.assertRaises(IMPORTER.VerificationError):
            IMPORTER.cli_provenance(args)

    def test_explicit_historical_source_checkout_is_commit_locked(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            checkout = base / "checkout"
            (checkout / "research").mkdir(parents=True)
            shutil.copy2(SEARCH_PATH, checkout / IMPORTER.SEARCH_SCRIPT_PATH)
            shutil.copy2(
                ROOT / IMPORTER.TARGET_SCRIPT_PATH,
                checkout / IMPORTER.TARGET_SCRIPT_PATH,
            )
            subprocess.run(["git", "init", "-q", str(checkout)], check=True)
            subprocess.run(
                ["git", "-C", str(checkout), "add", "research"], check=True
            )
            subprocess.run(
                [
                    "git",
                    "-C",
                    str(checkout),
                    "-c",
                    "user.name=Artifact Test",
                    "-c",
                    "user.email=artifact@example.invalid",
                    "commit",
                    "-qm",
                    "historical source fixture",
                ],
                check=True,
            )
            head = subprocess.run(
                ["git", "-C", str(checkout), "rev-parse", "HEAD"],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
            supplied = provenance()
            supplied["commit_sha"] = head
            artifact = base / "artifact"
            artifact.mkdir()
            build_fixture(artifact, False, supplied=supplied)
            certificate = IMPORTER.verify_artifact(
                artifact, supplied, source_root=checkout
            )
            self.assertTrue(
                certificate["search"]["source_checkout"]["git_head_verified"]
            )
            self.assertEqual(
                certificate["search"]["source_checkout"]["git_head"], head
            )
            self.assertNotIn(
                "source_root", certificate["search"]["source_checkout"]
            )
            self.assertNotIn(
                "rank-31 problem", certificate["claim_boundary"]["negative_result_scope"]
            )

            source = checkout / IMPORTER.SEARCH_SCRIPT_PATH
            source.write_text(source.read_text() + "\n# dirty historical checkout\n")
            with self.assertRaises(IMPORTER.VerificationError):
                IMPORTER.verify_artifact(artifact, supplied, source_root=checkout)

    def test_certificate_writer_refuses_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "certificate.json"
            IMPORTER.write_certificate(path, {"verified": True})
            with self.assertRaises(IMPORTER.VerificationError):
                IMPORTER.write_certificate(path, {"verified": False})


if __name__ == "__main__":
    unittest.main()
