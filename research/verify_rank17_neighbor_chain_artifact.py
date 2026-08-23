#!/usr/bin/env python3
"""Fail-closed, standard-library audit of the rank-17 chain artifact.

The Sage replay remains responsible for constructing p-neighbors, enumerating
short vectors, and classifying ADE root systems.  This audit pins the
previously reproduced chain byte-for-byte, checks every Gram matrix and every
serialized rational basis identity by exact arithmetic, cross-checks the Sage
records, and verifies a closed file manifest.  It does not promote the lattice
chain to a geometric Neron-Severi marking, an explicit K3 map, a rational
section, or a rank-31 curve.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import re
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any


MANIFEST_NAME = "artifact-manifest.json"
ARTIFACT_KIND = "rank17-exact-neighbor-chain-sage-10.9"
EXPECTED_CHAIN_SHA256 = (
    "f549651c06190e97f53947475ef9ee3149ddb175b2db2a289dd7db02d60b9d4b"
)
EXPECTED_TARGET_HASH = (
    "620a5e06473684d3e8015c0172f63c09c901e742ec02e77ba0aa35a923aa0295"
)
EXPECTED_SEED_HASH = (
    "d29e4bb24bcb1db376f142708ceb841f02c7eba4b27a7a04d141b0f06b8d40d1"
)
EXPECTED_SOURCE_WORKFLOW_SHA256 = (
    "118eab82cd0ce85bf53970a8b1cfe7df9de784cc7e03910429da71f5df8d0af1"
)
EXPECTED_REPLAYED_DISCOVERY_SHA256 = (
    "a621cb29563b207236a505d211daece5471c0ec808b8135e53f3511fd7d994c9"
)
EXPECTED_CLASSIFICATION_RECORD_SHA256 = (
    "b2ce2416e0268cfdefe5c8003660379b64d164943c70aac5769e00777f28b048"
)
EXPECTED_PRIMES = [5, 5, 2, 2, 2, 2, 2]
EXPECTED_ROOT_COUNTS = [150, 30, 10, 6, 4, 2, 2, 0]
EXPECTED_NORM4_COUNTS = [5424, 2654, 2548, 2576, 2590, 2602, 2602, 2622]
EXPECTED_ROOT_SYSTEMS = [
    "A11 + A3 + A2",
    "A3 + A3 + A1 + A1 + A1",
    "A2 + A1 + A1",
    "A1 + A1 + A1",
    "A1 + A1",
    "A1",
    "A1",
    "0",
]
EXPECTED_MW_RANKS = [1, 8, 13, 14, 15, 16, 16, 17]
EXPECTED_SEMANTIC_SUMMARY = {
    "basis_transition_count": 7,
    "chain_sha256": EXPECTED_CHAIN_SHA256,
    "classification_record_sha256": EXPECTED_CLASSIFICATION_RECORD_SHA256,
    "composed_basis_identity_verified": True,
    "neighbor_primes": EXPECTED_PRIMES,
    "norm4_counts": EXPECTED_NORM4_COUNTS,
    "root_counts": EXPECTED_ROOT_COUNTS,
    "root_systems": EXPECTED_ROOT_SYSTEMS,
    "mordell_weil_ranks": EXPECTED_MW_RANKS,
    "seed_hash": EXPECTED_SEED_HASH,
    "source_workflow_sha256": EXPECTED_SOURCE_WORKFLOW_SHA256,
    "target_hash": EXPECTED_TARGET_HASH,
}
REQUIRED_FILES = {
    "ade-classification.json",
    "artifact-verifier.py",
    "chain-summary.json",
    "classification.log",
    "exact-neighbor-chain.json",
    "exact-neighbor-maps.json",
    "map-replay-harness.py",
    "map-replay.log",
    "map-replayed-discovery.py",
    "map-summary.json",
    "manual-workflow.yml",
    "provenance.json",
    "rank17-root-elimination.json",
    "replay.log",
    "replayed-discovery.py",
    "stdout-summary.json",
    f"target-{EXPECTED_TARGET_HASH}.json",
}
HEX_40 = re.compile(r"[0-9a-f]{40}\Z")
HEX_64 = re.compile(r"[0-9a-f]{64}\Z")


class VerificationError(RuntimeError):
    """Raised whenever artifact evidence is absent, ambiguous, or inconsistent."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise VerificationError(f"duplicate JSON key: {key!r}")
        result[key] = value
    return result


def load_json(path: Path) -> Any:
    try:
        return json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=reject_duplicate_keys,
            parse_constant=lambda value: (_ for _ in ()).throw(
                VerificationError(f"non-finite JSON number: {value}")
            ),
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise VerificationError(f"cannot parse {path}: {exc}") from exc


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_matrix_hash(matrix: list[list[int]]) -> str:
    raw = json.dumps(matrix, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def exact_record_hash(record: dict[str, Any], field: str = "record_sha256") -> str:
    payload = copy.deepcopy(record)
    require(field in payload, f"record is missing {field}")
    payload.pop(field)
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def determinant_bareiss(matrix: list[list[int]]) -> int:
    n = len(matrix)
    require(n > 0 and all(len(row) == n for row in matrix), "matrix is not square")
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(n - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next(
                (row for row in range(pivot_index + 1, n) if work[row][pivot_index]),
                None,
            )
            if swap is None:
                return 0
            work[pivot_index], work[swap] = work[swap], work[pivot_index]
            sign *= -1
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, n):
            for column in range(pivot_index + 1, n):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                require(
                    numerator % previous == 0,
                    "non-exact division in Bareiss determinant",
                )
                work[row][column] = numerator // previous
        previous = pivot
    return sign * work[-1][-1]


def is_integer(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def validate_gram(matrix: Any, label: str) -> list[list[int]]:
    require(isinstance(matrix, list) and len(matrix) == 17, f"{label}: expected 17 rows")
    require(
        all(isinstance(row, list) and len(row) == 17 for row in matrix),
        f"{label}: expected a 17 by 17 matrix",
    )
    require(
        all(is_integer(value) for row in matrix for value in row),
        f"{label}: Gram entries must be JSON integers",
    )
    require(
        all(matrix[i][j] == matrix[j][i] for i in range(17) for j in range(17)),
        f"{label}: Gram matrix is not symmetric",
    )
    require(
        all(matrix[i][i] % 2 == 0 for i in range(17)),
        f"{label}: lattice is not even",
    )
    for rank in range(1, 18):
        minor = [row[:rank] for row in matrix[:rank]]
        require(
            determinant_bareiss(minor) > 0,
            f"{label}: Gram matrix is not positive definite",
        )
    require(determinant_bareiss(matrix) == 948, f"{label}: determinant is not 948")
    return matrix


def quadratic_numerator(vector: list[int], gram: list[list[int]]) -> int:
    return sum(
        vector[i] * gram[i][j] * vector[j]
        for i in range(17)
        for j in range(17)
    )


def parse_rational_matrix(data: Any, label: str) -> list[list[Fraction]]:
    require(isinstance(data, list) and len(data) == 17, f"{label}: expected 17 rows")
    matrix: list[list[Fraction]] = []
    for row_index, row in enumerate(data):
        require(isinstance(row, list) and len(row) == 17, f"{label}: row {row_index} is not length 17")
        parsed_row: list[Fraction] = []
        for column_index, entry in enumerate(row):
            require(
                isinstance(entry, list)
                and len(entry) == 2
                and all(is_integer(value) for value in entry),
                f"{label}: bad rational at ({row_index},{column_index})",
            )
            numerator, denominator = entry
            require(denominator > 0, f"{label}: non-positive denominator")
            require(math.gcd(numerator, denominator) == 1, f"{label}: non-reduced rational")
            parsed_row.append(Fraction(numerator, denominator))
        matrix.append(parsed_row)
    return matrix


def determinant_fraction(matrix: list[list[Fraction]]) -> Fraction:
    n = len(matrix)
    require(n > 0 and all(len(row) == n for row in matrix), "rational matrix is not square")
    work = [row[:] for row in matrix]
    determinant = Fraction(1)
    for column in range(n):
        pivot_row = next((row for row in range(column, n) if work[row][column]), None)
        if pivot_row is None:
            return Fraction(0)
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            determinant *= -1
        pivot = work[column][column]
        determinant *= pivot
        for row in range(column + 1, n):
            if not work[row][column]:
                continue
            factor = work[row][column] / pivot
            for index in range(column + 1, n):
                work[row][index] -= factor * work[column][index]
            work[row][column] = Fraction(0)
    return determinant


def multiply_rational_matrices(
    left: list[list[Fraction]], right: list[list[Fraction]]
) -> list[list[Fraction]]:
    require(len(left) == 17 and len(right) == 17, "basis composition has wrong rank")
    return [
        [sum(left[i][k] * right[k][j] for k in range(17)) for j in range(17)]
        for i in range(17)
    ]


def verify_gram_transport(
    transform: list[list[Fraction]],
    parent_gram: list[list[int]],
    child_gram: list[list[int]],
    label: str,
) -> None:
    for row in range(17):
        for column in range(17):
            value = sum(
                transform[i][row] * parent_gram[i][j] * transform[j][column]
                for i in range(17)
                for j in range(17)
            )
            require(value == child_gram[row][column], f"{label}: T^t G_parent T != G_child")


def verify_chain(root: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    path = root / "exact-neighbor-chain.json"
    require(path.is_file() and not path.is_symlink(), "exact chain is missing or unsafe")
    require(sha256_file(path) == EXPECTED_CHAIN_SHA256, "exact chain SHA-256 drifted")
    chain = load_json(path)
    require(isinstance(chain, list) and len(chain) == 8, "chain must contain seed plus seven moves")
    expected_keys = {
        "gram",
        "hash",
        "move",
        "norm4_count",
        "parent_archive_index",
        "parent_hash",
        "root_count",
        "score",
        "source",
    }
    hashes: set[str] = set()
    for index, entry in enumerate(chain):
        require(isinstance(entry, dict), f"chain[{index}] is not an object")
        require(set(entry) == expected_keys, f"chain[{index}] schema drifted")
        label = f"chain[{index}]"
        gram = validate_gram(entry["gram"], label)
        entry_hash = entry["hash"]
        require(isinstance(entry_hash, str) and HEX_64.fullmatch(entry_hash), f"{label}: bad hash")
        require(entry_hash not in hashes, f"{label}: duplicate lattice hash")
        hashes.add(entry_hash)
        require(canonical_matrix_hash(gram) == entry_hash, f"{label}: Gram hash mismatch")
        require(entry["root_count"] == EXPECTED_ROOT_COUNTS[index], f"{label}: root count drifted")
        require(entry["norm4_count"] == EXPECTED_NORM4_COUNTS[index], f"{label}: norm-4 count drifted")
        require(
            entry["score"]
            == [entry["root_count"], abs(entry["norm4_count"] - 2622), entry["norm4_count"]],
            f"{label}: score is inconsistent",
        )
        if index == 0:
            require(entry_hash == EXPECTED_SEED_HASH, "transparent seed hash drifted")
            require(entry["source"] == "A11-plus-K6-det79", "unexpected chain seed")
            require(entry["move"] is None and entry["parent_hash"] is None, "seed has a parent move")
            require(entry["parent_archive_index"] is None, "seed has a parent archive index")
            continue
        require(entry["source"] is None, f"{label}: non-seed source must be null")
        require(entry["parent_hash"] == chain[index - 1]["hash"], f"{label}: broken parent link")
        require(
            is_integer(entry["parent_archive_index"]) and entry["parent_archive_index"] >= 0,
            f"{label}: invalid parent archive index",
        )
        move = entry["move"]
        require(isinstance(move, dict) and set(move) == {"prime", "vector"}, f"{label}: move schema drifted")
        prime = move["prime"]
        vector = move["vector"]
        require(prime == EXPECTED_PRIMES[index - 1], f"{label}: neighbor prime drifted")
        require(
            isinstance(vector, list)
            and len(vector) == 17
            and all(is_integer(value) for value in vector),
            f"{label}: invalid neighbor vector",
        )
        require(any(value % prime for value in vector), f"{label}: vector vanishes modulo p")
        require(
            quadratic_numerator(vector, chain[index - 1]["gram"]) % (2 * prime) == 0,
            f"{label}: vector is not p-divisible for the parent quadratic form",
        )
    require(chain[-1]["hash"] == EXPECTED_TARGET_HASH, "terminal target hash drifted")
    return chain, {
        "chain_sha256": EXPECTED_CHAIN_SHA256,
        "neighbor_primes": EXPECTED_PRIMES,
        "norm4_counts": EXPECTED_NORM4_COUNTS,
        "root_counts": EXPECTED_ROOT_COUNTS,
        "seed_hash": EXPECTED_SEED_HASH,
        "target_hash": EXPECTED_TARGET_HASH,
    }


def verify_summary(root: Path) -> dict[str, Any]:
    summary = load_json(root / "chain-summary.json")
    require(isinstance(summary, dict), "chain summary is not an object")
    expected = {
        "status": "completed",
        "source_workflow": ".github/workflows/probe-rank29-height-lattice.yml",
        "source_workflow_sha256": EXPECTED_SOURCE_WORKFLOW_SHA256,
        "replay_script_sha256": EXPECTED_REPLAYED_DISCOVERY_SHA256,
        "chain_length_including_seed": 8,
        "neighbor_step_count": 7,
        "root_counts": EXPECTED_ROOT_COUNTS,
        "norm4_counts": EXPECTED_NORM4_COUNTS,
        "neighbor_primes": EXPECTED_PRIMES,
        "target_hash": EXPECTED_TARGET_HASH,
        "target_gram_sha256": EXPECTED_TARGET_HASH,
        "chain_sha256": EXPECTED_CHAIN_SHA256,
    }
    require(summary == expected, "chain summary schema or values drifted")
    require(
        sha256_file(root / "replayed-discovery.py") == EXPECTED_REPLAYED_DISCOVERY_SHA256,
        "replayed discovery program drifted",
    )
    return summary


def verify_classification(root: Path, chain: list[dict[str, Any]]) -> dict[str, Any]:
    record = load_json(root / "ade-classification.json")
    require(isinstance(record, dict), "ADE classification is not an object")
    require(
        record.get("record_sha256") == EXPECTED_CLASSIFICATION_RECORD_SHA256,
        "ADE classification record hash drifted",
    )
    require(exact_record_hash(record) == record["record_sha256"], "ADE record hash is internally invalid")
    require(record.get("status") == "completed", "ADE classification did not complete")
    require(record.get("chain_sha256") == EXPECTED_CHAIN_SHA256, "ADE record names another chain")
    records = record.get("records")
    require(isinstance(records, list) and len(records) == 8, "ADE record count drifted")
    for index, (classified, entry) in enumerate(zip(records, chain)):
        require(isinstance(classified, dict), f"ADE record {index} is not an object")
        for field in ("hash", "gram", "move", "source", "root_count", "norm4_count"):
            require(classified.get(field) == entry[field], f"ADE record {index} disagrees on {field}")
        require(classified.get("root_system") == EXPECTED_ROOT_SYSTEMS[index], f"ADE type {index} drifted")
        require(classified.get("mordell_weil_rank") == EXPECTED_MW_RANKS[index], f"MW rank {index} drifted")
    seed = record.get("seed_fibration_data")
    require(isinstance(seed, dict), "seed fibration data are missing")
    require(seed.get("essential_root_system") == "A11 + A3 + A2", "seed ADE type drifted")
    require(seed.get("generator_height") == "79/12", "seed height drifted")
    require(
        seed.get("kodaira_configuration_possibilities")
        == [["I12", "I4", "I3"], ["I12", "I4", "IV"]],
        "I3-or-IV claim boundary was lost",
    )
    terminal = record.get("terminal_data")
    require(
        isinstance(terminal, dict)
        and terminal.get("root_system") == "0"
        and terminal.get("mordell_weil_rank") == 17
        and terminal.get("mordell_weil_regulator") == "948",
        "terminal lattice classification drifted",
    )
    return record


def verify_basis_maps(root: Path, chain: list[dict[str, Any]]) -> dict[str, Any]:
    path = root / "exact-neighbor-maps.json"
    record = load_json(path)
    require(isinstance(record, dict), "basis-map record is not an object")
    require(exact_record_hash(record) == record.get("record_sha256"), "basis-map record hash is invalid")
    require(record.get("schema_version") == 1 and record.get("status") == "completed", "basis-map replay did not complete")
    require(record.get("chain_sha256") == EXPECTED_CHAIN_SHA256, "basis maps name another chain")
    require(record.get("target_hash") == EXPECTED_TARGET_HASH, "basis maps name another target")
    require(record.get("source_workflow_sha256") == EXPECTED_SOURCE_WORKFLOW_SHA256, "basis-map source drifted")
    require(
        record.get("map_replayed_discovery_sha256")
        == sha256_file(root / "map-replayed-discovery.py"),
        "generated basis-map discovery script drifted",
    )
    require(record.get("transition_count") == 7, "basis-map transition count drifted")
    transitions = record.get("transitions")
    require(isinstance(transitions, list) and len(transitions) == 7, "seven basis maps are required")
    composed = [
        [Fraction(int(row == column)) for column in range(17)]
        for row in range(17)
    ]
    expected_transition_keys = {
        "basis_identity_verified",
        "basis_transform_child_in_parent",
        "child_hash",
        "determinant",
        "lll_composition",
        "parent_hash",
        "prime",
        "sage_api",
        "vector",
    }
    for index, transition in enumerate(transitions):
        require(isinstance(transition, dict) and set(transition) == expected_transition_keys, f"basis map {index} schema drifted")
        parent = chain[index]
        child = chain[index + 1]
        require(transition["parent_hash"] == parent["hash"], f"basis map {index} parent drifted")
        require(transition["child_hash"] == child["hash"], f"basis map {index} child drifted")
        require(transition["prime"] == child["move"]["prime"], f"basis map {index} prime drifted")
        require(transition["vector"] == child["move"]["vector"], f"basis map {index} vector drifted")
        require(transition["basis_identity_verified"] is True, f"basis map {index} lacks Sage identity gate")
        require(
            transition["sage_api"]
            == "QuadraticForm.find_p_neighbor_from_vec(return_matrix=True)",
            f"basis map {index} API drifted",
        )
        require(transition["lll_composition"] == "T = B * U", f"basis map {index} lost its LLL composition")
        transform = parse_rational_matrix(
            transition["basis_transform_child_in_parent"], f"basis map {index}"
        )
        verify_gram_transport(transform, parent["gram"], child["gram"], f"basis map {index}")
        determinant = determinant_fraction(transform)
        require(abs(determinant) == 1, f"basis map {index} determinant is not a unit")
        require(
            transition["determinant"]
            == [determinant.numerator, determinant.denominator],
            f"basis map {index} determinant record drifted",
        )
        composed = multiply_rational_matrices(composed, transform)
    recorded_composed = parse_rational_matrix(
        record.get("target_basis_in_seed_coordinates"),
        "composed target basis in seed coordinates",
    )
    require(composed == recorded_composed, "recorded composed basis does not equal product of transitions")
    verify_gram_transport(composed, chain[0]["gram"], chain[-1]["gram"], "composed basis map")
    determinant = determinant_fraction(composed)
    require(abs(determinant) == 1, "composed basis determinant is not a unit")
    require(
        record.get("composed_determinant") == [determinant.numerator, determinant.denominator],
        "composed determinant record drifted",
    )
    require(record.get("composed_basis_identity_verified") is True, "composed Gram identity gate is absent")

    summary = load_json(root / "map-summary.json")
    require(isinstance(summary, dict), "map summary is not an object")
    require(summary.get("status") == "completed", "map summary did not complete")
    require(summary.get("chain_sha256") == EXPECTED_CHAIN_SHA256, "map summary names another chain")
    require(summary.get("map_record_sha256") == record["record_sha256"], "map summary record hash drifted")
    require(summary.get("map_file_sha256") == sha256_file(path), "map summary file hash drifted")
    require(
        summary.get("map_replayed_discovery_sha256")
        == sha256_file(root / "map-replayed-discovery.py"),
        "map summary generated-script hash drifted",
    )
    require(summary.get("transition_count") == 7, "map summary transition count drifted")
    require(summary.get("target_hash") == EXPECTED_TARGET_HASH, "map summary target drifted")
    require(summary.get("composed_basis_identity_verified") is True, "map summary identity gate is absent")
    return record


def verify_discovery_record(root: Path) -> None:
    record = load_json(root / "rank17-root-elimination.json")
    require(isinstance(record, dict), "discovery record is not an object")
    require(exact_record_hash(record) == record.get("record_sha256"), "discovery record hash is invalid")
    require(record.get("status") == "target_found" and record.get("target_found") is True, "replay found no target")
    require(record.get("random_seed") == 2026081001, "random seed drifted")
    target = record.get("target")
    require(
        target == {"rank": 17, "determinant": 948, "minimum": 4, "kissing_number": 2622},
        "discovery target drifted",
    )
    verified = record.get("verified_targets")
    require(isinstance(verified, list) and len(verified) == 1, "expected exactly one verified target")
    target_record = verified[0]
    require(target_record.get("hash") == EXPECTED_TARGET_HASH, "discovery verified another target")
    require(target_record.get("target_verified") is True, "target lacks independent short-vector check")
    require(target_record.get("independent_short_vector_root_count") == 0, "terminal roots reappeared")
    require(target_record.get("independent_short_vector_norm4_count") == 2622, "terminal theta count drifted")


def verify_provenance(root: Path) -> None:
    provenance = load_json(root / "provenance.json")
    require(isinstance(provenance, dict), "provenance is not an object")
    require(provenance.get("schema_version") == 1, "provenance schema drifted")
    require(provenance.get("artifact_kind") == ARTIFACT_KIND, "artifact kind drifted")
    require(
        provenance.get("event_name") in {"workflow_dispatch", "push"},
        "artifact was not created by an authorized dispatch or checkpoint push",
    )
    require(provenance.get("repository") == "metaforismo/elliptic-rank-31", "unexpected repository")
    require(isinstance(provenance.get("commit_sha"), str) and HEX_40.fullmatch(provenance["commit_sha"]), "bad commit SHA")
    require(
        provenance.get("container_image")
        == "sagemath/sagemath:10.9@sha256:e068670ae5863b54b2550e72437ec637b0283acb0dc712c8584c124dbf44e667",
        "Sage container pin drifted",
    )
    require(str(provenance.get("sage_version", "")).startswith("SageMath version 10.9"), "Sage 10.9 was not recorded")
    source_hashes = provenance.get("source_sha256")
    require(isinstance(source_hashes, dict), "source hashes are missing")
    require(
        source_hashes.get(".github/workflows/probe-rank29-height-lattice.yml")
        == EXPECTED_SOURCE_WORKFLOW_SHA256,
        "embedded discovery source drifted",
    )
    require(
        source_hashes.get("research/replay_rank17_neighbor_chain.py")
        == "9726288ab1ea647f4c7f8e556fa886f1bd2c5cb97155bb404b9b977b3e3514a5",
        "replay harness drifted",
    )
    require(
        source_hashes.get("research/classify_rank17_neighbor_chain.py")
        == "9a0f549855f5cb22d6d7a6088215c2086cbb06ab6874b6cb0ba4d3c6b6bc7f67",
        "classification harness drifted",
    )
    require(
        source_hashes.get("research/verify_rank17_neighbor_chain_artifact.py")
        == sha256_file(root / "artifact-verifier.py"),
        "artifact verifier copy disagrees with provenance",
    )
    require(
        source_hashes.get("research/replay_rank17_neighbor_chain_with_maps.py")
        == sha256_file(root / "map-replay-harness.py"),
        "basis-map replay harness disagrees with provenance",
    )
    require(
        source_hashes.get(".github/workflows/replay-rank17-neighbor-chain-manual.yml")
        == sha256_file(root / "manual-workflow.yml"),
        "manual workflow copy disagrees with provenance",
    )
    claims = provenance.get("claim_gates")
    require(
        claims
        == {
            "exact_ade_classification_replayed": True,
            "exact_lattice_neighbor_chain_replayed": True,
            "exact_rational_parent_child_basis_maps_replayed": True,
            "explicit_integral_chain_from_e8_a2_target": False,
            "explicit_k3_neighbor_map": False,
            "geometric_neron_severi_marking": False,
            "rank31_curve_over_q": False,
            "rational_p3_section": False,
        },
        "claim gates drifted or were overpromoted",
    )


def regular_flat_files(root: Path, include_manifest: bool) -> dict[str, Path]:
    require(root.is_dir() and not root.is_symlink(), "artifact root is missing or a symlink")
    files: dict[str, Path] = {}
    for path in root.iterdir():
        require(not path.is_symlink(), f"artifact contains symlink: {path.name}")
        require(path.is_file(), f"artifact contains non-file entry: {path.name}")
        if not include_manifest and path.name == MANIFEST_NAME:
            continue
        files[path.name] = path
    return files


def make_manifest(root: Path, semantic_summary: dict[str, Any]) -> dict[str, Any]:
    files = regular_flat_files(root, include_manifest=False)
    entries = [
        {"path": name, "size": path.stat().st_size, "sha256": sha256_file(path)}
        for name, path in sorted(files.items())
    ]
    payload: dict[str, Any] = {
        "schema_version": 1,
        "artifact_kind": ARTIFACT_KIND,
        "required_files": sorted(REQUIRED_FILES),
        "semantic_summary": semantic_summary,
        "files": entries,
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["record_sha256"] = hashlib.sha256(raw).hexdigest()
    return payload


def write_manifest(root: Path, semantic_summary: dict[str, Any]) -> None:
    path = root / MANIFEST_NAME
    require(not path.exists(), f"refusing to overwrite existing {MANIFEST_NAME}")
    path.write_text(json.dumps(make_manifest(root, semantic_summary), indent=2, sort_keys=True) + "\n")


def verify_manifest(root: Path) -> None:
    manifest = load_json(root / MANIFEST_NAME)
    require(isinstance(manifest, dict), "manifest is not an object")
    require(exact_record_hash(manifest) == manifest.get("record_sha256"), "manifest record hash is invalid")
    require(manifest.get("schema_version") == 1, "manifest schema drifted")
    require(manifest.get("artifact_kind") == ARTIFACT_KIND, "manifest artifact kind drifted")
    require(manifest.get("required_files") == sorted(REQUIRED_FILES), "manifest required-file gate drifted")
    require(manifest.get("semantic_summary") == EXPECTED_SEMANTIC_SUMMARY, "manifest semantic pins drifted")
    entries = manifest.get("files")
    require(isinstance(entries, list), "manifest file list is missing")
    names: set[str] = set()
    for entry in entries:
        require(isinstance(entry, dict) and set(entry) == {"path", "size", "sha256"}, "bad manifest entry")
        name = entry["path"]
        require(isinstance(name, str) and name and "/" not in name and "\\" not in name, "unsafe manifest path")
        require(name not in names and name != MANIFEST_NAME, "duplicate or self-referential manifest path")
        names.add(name)
        path = root / name
        require(path.is_file() and not path.is_symlink(), f"manifest file missing or unsafe: {name}")
        require(entry["size"] == path.stat().st_size, f"file size mismatch: {name}")
        require(isinstance(entry["sha256"], str) and HEX_64.fullmatch(entry["sha256"]), f"bad SHA-256: {name}")
        require(entry["sha256"] == sha256_file(path), f"file hash mismatch: {name}")
    actual = set(regular_flat_files(root, include_manifest=True))
    require(actual == names | {MANIFEST_NAME}, "artifact has missing or unmanifested files")
    require(REQUIRED_FILES <= names, "artifact omits required evidence files")


def verify_semantics(root: Path) -> dict[str, Any]:
    chain, _ = verify_chain(root)
    verify_summary(root)
    classification = verify_classification(root, chain)
    verify_basis_maps(root, chain)
    verify_discovery_record(root)
    verify_provenance(root)
    require((root / "replay.log").stat().st_size > 0, "replay log is empty")
    require((root / "classification.log").stat().st_size > 0, "classification log is empty")
    require((root / "map-replay.log").stat().st_size > 0, "basis-map replay log is empty")
    summary = copy.deepcopy(EXPECTED_SEMANTIC_SUMMARY)
    require(classification["record_sha256"] == summary["classification_record_sha256"], "classification pin mismatch")
    return summary


def verify_artifact(root: Path, require_manifest_file: bool = True) -> dict[str, Any]:
    if require_manifest_file:
        verify_manifest(root)
    return verify_semantics(root)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", type=Path, help="extracted artifact directory")
    parser.add_argument(
        "--write-manifest",
        action="store_true",
        help="validate semantics and create the closed manifest (CI finalization only)",
    )
    args = parser.parse_args()
    root = args.artifact.resolve()
    try:
        if args.write_manifest:
            summary = verify_artifact(root, require_manifest_file=False)
            write_manifest(root, summary)
            verify_artifact(root, require_manifest_file=True)
        else:
            summary = verify_artifact(root, require_manifest_file=True)
    except VerificationError as exc:
        print(f"FAILED: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    print(
        json.dumps(
            {
                "status": "VERIFIED",
                "artifact_kind": ARTIFACT_KIND,
                "chain_sha256": summary["chain_sha256"],
                "target_hash": summary["target_hash"],
                "claim_boundary": (
                    "exact lattice chain and rational basis maps only; no "
                    "E8+A2^3 integral bridge, geometric NS/K3 map, rational "
                    "P3, or rank-31 curve"
                ),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
