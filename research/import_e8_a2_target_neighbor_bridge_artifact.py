#!/usr/bin/env python3
"""Fail-closed importer for an E8+A2^3 target-neighbor search artifact.

The Sage search is the authority for lattice construction, theta series, and
global isometry tests.  This standard-library audit checks the downloaded
artifact's closed SHA-256 manifest, internally hashed JSON/JSONL records,
counter/round consistency, and conservative claim gates.  It then writes a
small, hash-locked certificate suitable for committing to the repository.

Neither a successful lattice bridge nor a bounded negative search is an
explicit K3 transport, a rational fourth section, or an elliptic rank-record
curve.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import re
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import urlparse

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from research import sample_projective_lines as projective_sampler  # noqa: E402


SCHEMA_VERSION = 1
CERTIFICATE_KIND = "e8-a2-target-neighbor-bridge-artifact-audit"
WORKFLOW_PATH = ".github/workflows/probe-e8-a2-target-neighbor-bridge.yml"
SEARCH_SCRIPT_PATH = "research/search_e8_a2_target_neighbor_bridge.py"
TARGET_SCRIPT_PATH = "research/certify_e8_a2_shimura_bridge.py"
SAMPLER_SCRIPT_PATH = "research/sample_projective_lines.py"
SAMPLING_LOG_NAME = "sampling-windows.jsonl"
IMPORTER_SCRIPT_PATH = "research/import_e8_a2_target_neighbor_bridge_artifact.py"
EXPECTED_P2_ERROR = "either y is not primitive or self is not even, maximal at 2"
EXPECTED_P2_ERROR_REPR = f"ValueError({EXPECTED_P2_ERROR!r})"
P2_ENUMERATION_BOUNDARY = (
    "The determinant is 948, so 2 divides the determinant.  Sage's projective-line "
    "construction can reject the expected non-maximal/even 2-neighbor case; therefore "
    "the recorded p=2 window is not a complete enumeration of all 2-neighbors.  Any "
    "projective line with this classified boundary failure is recorded as attempted but "
    "is not excluded by the bounded negative result."
)
EXPECTED_WORKFLOW_NAME = "probe-e8-a2-target-neighbor-bridge"
EXPECTED_JOB_NAME = "exact-neighbor-search"
EXPECTED_CONTAINER_IMAGE = "sagemath/sagemath:10.9"
EXPECTED_CONTAINER_DIGEST = (
    "sha256:e068670ae5863b54b2550e72437ec637b0283acb0dc712c8584c124dbf44e667"
)
MANIFEST_NAME = "SHA256SUMS"
JSONL_NAMES = (
    "states.jsonl",
    "moves.jsonl",
    "attempts.jsonl",
    "rounds.jsonl",
    "isometry-checks.jsonl",
)
BASE_ARTIFACT_FILES = {
    *JSONL_NAMES,
    "checkpoint.json",
    "result.json",
    "run.log",
}
RESULT_EVIDENCE_FILES = {
    *JSONL_NAMES,
    "checkpoint.json",
}
HEX40 = re.compile(r"[0-9a-f]{40}\Z")
HEX64 = re.compile(r"[0-9a-f]{64}\Z")
REPOSITORY = re.compile(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+\Z")
PORTABLE_SUM_LINE = re.compile(r"([0-9a-f]{64}) ([ *])(.+)\Z")
INITIAL_HASHES = {
    "origin": "dc4108ea4612195c2ba8350b0b1b35067f92ec2d49801cd2cb597cbb61129690",
    "transparent": "10270a49860bff642fe34bc3eb0fc213623c7dff8aca08b6df51e2832c600bff",
    "rootless": "620a5e06473684d3e8015c0172f63c09c901e742ec02e77ba0aa35a923aa0295",
}
FALSE_CLAIM_GATES = (
    "explicit_stable_NS_transport_completed",
    "explicit_K3_fibration_switch_completed",
    "explicit_moduli_map_completed",
    "rational_P3_found",
    "rank31_curve_found",
)
CURRENT_TARGET_CLAIM_GATE = "rank32_curve_found"
COUNTER_FIELDS = (
    "attempted_lines",
    "successful_moves",
    "failed_moves",
    "duplicate_child_grams",
    "unique_states_discovered",
)
SIDE_COUNTER_FIELDS = {
    "attempted_lines": "attempted_lines",
    "successful_moves": "successful_moves",
    "failed_moves": "failed_moves",
    "duplicate_child_grams": "duplicate_child_grams",
    "unique_states_discovered": "new_unique_states",
}
ISOMETRY_COUNT_FIELDS = (
    "pairs_seen",
    "qfisom_checked",
    "identical_hash_matches",
    "skipped_due_to_budget",
)
MODERN_RESULT_FIELDS = (
    "artifact_importable",
    "termination_reason",
    "runtime_provenance",
    "isometry_comparisons",
)


class VerificationError(RuntimeError):
    """Raised for missing, ambiguous, unsafe, or inconsistent evidence."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def is_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def is_prime_integer(value: Any) -> bool:
    return (
        is_int(value)
        and value >= 2
        and all(value % divisor for divisor in range(2, math.isqrt(value) + 1))
    )


def is_modern_result(result: dict[str, Any]) -> bool:
    """Treat any modern marker as modern, preventing a partial schema downgrade."""

    return any(field in result for field in MODERN_RESULT_FIELDS)


def verify_selection_schema(result: dict[str, Any]) -> bool:
    """Return whether v2 hash sampling is enabled; reject partial downgrades."""
    config = result.get("configuration")
    require(isinstance(config, dict), "result: configuration missing")
    if result["schema_version"] == 1:
        require(
            not ({"line_selection", "sampling"} & set(config))
            and "projective_sampler_script" not in result.get("source_files", {}),
            "sampling metadata cannot be downgraded to schema 1",
        )
        return False
    require(is_modern_result(result), "schema 2 requires modern evidence")
    selection = config.get("line_selection")
    require(selection in {"lexicographic", "hash-projective"}, "unknown line selection")
    require("sampling" in config, "schema 2 lacks sampling configuration")
    if selection == "lexicographic":
        require(config["sampling"] is None, "lexicographic selection contains sampling metadata")
        require(config.get("enumeration_order") == "Sage 10.9 find_primitive_p_divisible_vector__next", "lexicographic order drift")
        return False
    sampling = config["sampling"]
    require(isinstance(sampling, dict), "hash sampling configuration missing")
    require(set(sampling) == {"algorithm", "seed", "max_draws_per_window", "ordinal_semantics"}, "sampling configuration schema drift")
    require(sampling["algorithm"] == projective_sampler.ALGORITHM, "unknown sampling algorithm")
    require(config.get("enumeration_order") == projective_sampler.ALGORITHM, "sampling enumeration-order drift")
    require(sampling["ordinal_semantics"] == "accepted unique isotropic projective lines, before offset", "sampling ordinal semantics drift")
    try:
        require(isinstance(config.get("primes"), list) and config["primes"], "sampling primes missing")
        for prime in config["primes"]:
            projective_sampler.validate_parameters(
                prime, sampling["seed"], config.get("projective_line_offset"),
                config.get("projective_lines_per_state_prime"), sampling["max_draws_per_window"],
            )
    except ValueError as exc:
        raise VerificationError(f"invalid sampling parameters: {exc}") from exc
    return True


def verify_sampling_windows(
    result: dict[str, Any], windows: list[dict[str, Any]],
    attempts: list[dict[str, Any]], states: dict[str, dict[str, Any]],
    rounds: list[dict[str, Any]],
) -> dict[str, Any]:
    """Replay every raw draw and bind each attempted line to its ordered window."""
    config = result["configuration"]
    sampling = config["sampling"]
    terminal_prefix = result["termination_reason"] in {
        "exact_bridge_found", "successful_move_budget_exhausted",
    }
    beams = {
        side: [state["uid"] for state in result["initial_states"] if state["side"] == side]
        for side in config["side_expansion_order"]
    }
    expanded: set[str] = set()
    expected_schedule: list[tuple[Any, ...]] = []
    for round_record in rounds:
        for side in config["side_expansion_order"]:
            for uid in beams[side]:
                require(uid in states and states[uid]["side"] == side, "sampling schedule has unknown/wrong-side parent")
                require(states[uid]["round_discovered"] < round_record["round"], "sampling parent was not discovered before its expansion")
                if uid in expanded:
                    continue
                expanded.add(uid)
                expected_schedule.extend(
                    (round_record["round"], side, uid, prime) for prime in config["primes"]
                )
            side_record = round_record["sides"].get(side)
            require(side_record is not None or (terminal_prefix and round_record is rounds[-1]), "sampling round omits an active side")
            beams[side] = [state["uid"] for state in side_record["next_beam"]] if side_record else []
    actual_schedule: list[tuple[Any, ...]] = []
    attempt_cursor = 0
    totals: dict[str, int] = {}
    coordinate_counts = {str(prime): [0] * 17 for prime in config["primes"]}
    for index, record in enumerate(windows):
        label = f"sampling-windows.jsonl:{index + 1}"
        require(set(record) == {"round", "side", "parent_uid", "window"}, f"{label}: record schema drift")
        window = record["window"]
        require(isinstance(window, dict), f"{label}: window missing")
        parent = states.get(record["parent_uid"])
        require(parent is not None and parent["side"] == record["side"], f"{label}: invalid parent")
        prime = window.get("prime")
        require(prime in config["primes"], f"{label}: unconfigured prime")
        actual_schedule.append((record["round"], record["side"], record["parent_uid"], prime))
        try:
            replay = projective_sampler.sample_window(
                parent["gram"], prime, sampling["seed"], config["projective_line_offset"],
                config["projective_lines_per_state_prime"], sampling["max_draws_per_window"],
            )
        except (ValueError, projective_sampler.SamplingBudgetError) as exc:
            raise VerificationError(f"{label}: sampler replay failed: {exc}") from exc
        require(window == replay, f"{label}: deterministic sampler replay mismatch")
        for key, count in replay["counters"].items():
            totals[key] = totals.get(key, 0) + count
        consumed = 0
        for line in replay["selected_lines"]:
            if attempt_cursor >= len(attempts):
                break
            attempt = attempts[attempt_cursor]
            if attempt.get("sampling_window_sha256") != replay["window_sha256"]:
                break
            for key in ("round", "side", "parent_uid"):
                require(attempt.get(key) == record[key], f"{label}: attempt {key} mismatch")
            require(attempt.get("prime") == prime, f"{label}: attempt prime mismatch")
            require(attempt.get("projective_line_ordinal_one_based") == line["ordinal_one_based"], f"{label}: attempt ordinal mismatch")
            require(attempt.get("projective_isotropic_vector") == line["vector"], f"{label}: attempt vector mismatch")
            for column, value in enumerate(line["vector"]):
                coordinate_counts[str(prime)][column] += int(value != 0)
            attempt_cursor += 1
            consumed += 1
        require(consumed > 0, f"{label}: unused or mislinked sampling window")
        require(
            consumed == len(replay["selected_lines"])
            or (terminal_prefix and index == len(windows) - 1),
            f"{label}: truncated window without a terminal search budget/bridge",
        )
    require(actual_schedule == expected_schedule[:len(actual_schedule)], "sampling windows do not follow the configured parent/prime schedule")
    require(terminal_prefix or len(actual_schedule) == len(expected_schedule), "sampling schedule ended early without a terminal budget/bridge")
    require(attempt_cursor == len(attempts), "attempts exist outside the replayed sampling windows")
    selected_count = len(windows) * config["projective_lines_per_state_prime"]
    return {
        "algorithm": projective_sampler.ALGORITHM,
        "all_recorded_draw_streams_replayed": True,
        "all_attempted_vectors_match_replayed_windows": True,
        "windows_replayed": len(windows),
        "sampling_costs": totals,
        "selected_lines_prepared": selected_count,
        "lines_attempted": attempt_cursor,
        "selected_lines_not_attempted_after_terminal_stop": selected_count - attempt_cursor,
        "attempted_coordinate_nonzero_counts_by_prime": coordinate_counts,
        "coverage_claim": "Finite hash-based full-coordinate sample, not an exhaustive graph enumeration.",
    }


def duplicate_key_guard(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise VerificationError(f"duplicate JSON key: {key!r}")
        result[key] = value
    return result


def parse_json_text(text: str, label: str) -> Any:
    try:
        return json.loads(
            text,
            object_pairs_hook=duplicate_key_guard,
            parse_constant=lambda value: (_ for _ in ()).throw(
                VerificationError(f"{label}: non-finite JSON number {value}")
            ),
        )
    except json.JSONDecodeError as exc:
        raise VerificationError(f"{label}: invalid JSON: {exc}") from exc


def load_json(path: Path) -> Any:
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise VerificationError(f"cannot read {path.name}: {exc}") from exc
    return parse_json_text(text, path.name)


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    try:
        raw = path.read_bytes()
        text = raw.decode("utf-8")
    except (OSError, UnicodeError) as exc:
        raise VerificationError(f"cannot read {path.name}: {exc}") from exc
    require(not raw or raw.endswith(b"\n"), f"{path.name}: final line is not newline-terminated")
    records: list[dict[str, Any]] = []
    for number, line in enumerate(text.splitlines(), 1):
        require(line.strip() == line and line, f"{path.name}:{number}: blank or padded line")
        record = parse_json_text(line, f"{path.name}:{number}")
        require(isinstance(record, dict), f"{path.name}:{number}: record is not an object")
        records.append(record)
    return records


def canonical_bytes(payload: Any) -> bytes:
    return json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")


def payload_sha256(payload: Any) -> str:
    return hashlib.sha256(canonical_bytes(payload)).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise VerificationError(f"cannot hash {path.name}: {exc}") from exc
    return digest.hexdigest()


def determinant_bareiss(matrix: list[list[int]]) -> int:
    """Return an exact integer determinant without third-party packages."""

    size = len(matrix)
    require(size > 0 and all(len(row) == size for row in matrix), "matrix is not square")
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next(
                (row for row in range(pivot_index + 1, size) if work[row][pivot_index]),
                None,
            )
            if swap is None:
                return 0
            work[pivot_index], work[swap] = work[swap], work[pivot_index]
            sign *= -1
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                require(numerator % previous == 0, "non-exact Bareiss division")
                work[row][column] = numerator // previous
        previous = pivot
    return sign * work[-1][-1]


def parse_integer_matrix(data: Any, label: str) -> list[list[int]]:
    require(isinstance(data, list) and len(data) == 17, f"{label}: expected 17 rows")
    require(
        all(isinstance(row, list) and len(row) == 17 for row in data),
        f"{label}: expected a 17 by 17 matrix",
    )
    require(
        all(is_int(value) for row in data for value in row),
        f"{label}: entries are not JSON integers",
    )
    return data


def parse_rational_string(value: Any, label: str) -> Fraction:
    require(isinstance(value, str), f"{label}: rational entry is not a string")
    require(
        re.fullmatch(r"-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?", value) is not None,
        f"{label}: malformed rational string",
    )
    parsed = Fraction(value)
    require(str(parsed) == value, f"{label}: rational string is not reduced/canonical")
    return parsed


def parse_rational_matrix(data: Any, label: str) -> list[list[Fraction]]:
    require(isinstance(data, list) and len(data) == 17, f"{label}: expected 17 rows")
    matrix: list[list[Fraction]] = []
    for row_index, row in enumerate(data):
        require(isinstance(row, list) and len(row) == 17, f"{label}: row {row_index} is not length 17")
        matrix.append(
            [
                parse_rational_string(value, f"{label}[{row_index},{column_index}]")
                for column_index, value in enumerate(row)
            ]
        )
    return matrix


def fraction_matrix_from_integers(data: Any, label: str) -> list[list[Fraction]]:
    return [
        [Fraction(value) for value in row]
        for row in parse_integer_matrix(data, label)
    ]


def identity_matrix() -> list[list[Fraction]]:
    return [
        [Fraction(int(row == column)) for column in range(17)]
        for row in range(17)
    ]


def identity_integer_matrix() -> list[list[int]]:
    return [[int(row == column) for column in range(17)] for row in range(17)]


def multiply_matrices(
    left: list[list[Fraction]], right: list[list[Fraction]]
) -> list[list[Fraction]]:
    require(len(left) == 17 and len(right) == 17, "matrix product has wrong rank")
    require(all(len(row) == 17 for row in left + right), "matrix product has wrong shape")
    return [
        [sum(left[row][inner] * right[inner][column] for inner in range(17)) for column in range(17)]
        for row in range(17)
    ]


def invert_matrix(matrix: list[list[Fraction]], label: str) -> list[list[Fraction]]:
    require(len(matrix) == 17 and all(len(row) == 17 for row in matrix), f"{label}: wrong shape")
    identity = identity_matrix()
    work = [matrix[row][:] + identity[row] for row in range(17)]
    for column in range(17):
        pivot_row = next((row for row in range(column, 17) if work[row][column]), None)
        require(pivot_row is not None, f"{label}: singular matrix")
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
        pivot = work[column][column]
        work[column] = [value / pivot for value in work[column]]
        for row in range(17):
            if row == column or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                work[row][index] - factor * work[column][index]
                for index in range(34)
            ]
    return [row[17:] for row in work]


def determinant_fraction(matrix: list[list[Fraction]], label: str) -> Fraction:
    require(len(matrix) == 17 and all(len(row) == 17 for row in matrix), f"{label}: wrong shape")
    work = [row[:] for row in matrix]
    determinant = Fraction(1)
    for column in range(17):
        pivot_row = next((row for row in range(column, 17) if work[row][column]), None)
        if pivot_row is None:
            return Fraction(0)
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            determinant *= -1
        pivot = work[column][column]
        determinant *= pivot
        for row in range(column + 1, 17):
            if not work[row][column]:
                continue
            factor = work[row][column] / pivot
            for index in range(column + 1, 17):
                work[row][index] -= factor * work[column][index]
            work[row][column] = Fraction(0)
    return determinant


def verify_gram_transport(
    transform: list[list[Fraction]],
    parent_gram: list[list[int]],
    child_gram: list[list[int]],
    label: str,
) -> None:
    for row in range(17):
        for column in range(17):
            value = sum(
                transform[left][row]
                * parent_gram[left][right]
                * transform[right][column]
                for left in range(17)
                for right in range(17)
            )
            require(value == child_gram[row][column], f"{label}: T^t G_parent T != G_child")


def matrix_denominator(matrix: list[list[Fraction]]) -> int:
    result = 1
    for row in matrix:
        for value in row:
            result = math.lcm(result, value.denominator)
    return result


def rank_mod_prime(matrix: list[list[int]], prime: int) -> int:
    work = [[value % prime for value in row] for row in matrix]
    rank = 0
    for column in range(17):
        pivot = next((row for row in range(rank, 17) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, prime)
        work[rank] = [(value * inverse) % prime for value in work[rank]]
        for row in range(17):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                (work[row][index] - factor * work[rank][index]) % prime
                for index in range(17)
            ]
        rank += 1
    return rank


def serial_fraction_matrix(matrix: list[list[Fraction]]) -> list[list[str]]:
    return [[str(value) for value in row] for row in matrix]


def validate_exact_gram(data: Any, label: str) -> list[list[int]]:
    gram = parse_integer_matrix(data, label)
    require(
        all(gram[row][column] == gram[column][row] for row in range(17) for column in range(17)),
        f"{label}: Gram matrix is not symmetric",
    )
    require(all(gram[index][index] % 2 == 0 for index in range(17)), f"{label}: lattice is not even")
    for rank in range(1, 18):
        leading_minor = [row[:rank] for row in gram[:rank]]
        require(determinant_bareiss(leading_minor) > 0, f"{label}: Gram matrix is not positive definite")
    require(determinant_bareiss(gram) == 948, f"{label}: determinant is not 948")
    return gram


def verify_embedded_hash(record: dict[str, Any], field: str, label: str) -> str:
    claimed = record.get(field)
    require(isinstance(claimed, str) and HEX64.fullmatch(claimed), f"{label}: invalid {field}")
    core = copy.deepcopy(record)
    core.pop(field)
    require(payload_sha256(core) == claimed, f"{label}: {field} mismatch")
    return claimed


def safe_flat_files(root: Path) -> dict[str, Path]:
    require(root.is_dir() and not root.is_symlink(), "artifact directory is missing or is a symlink")
    files: dict[str, Path] = {}
    for path in root.iterdir():
        require(not path.is_symlink(), f"artifact contains symlink: {path.name}")
        require(path.is_file(), f"artifact contains non-file entry: {path.name}")
        files[path.name] = path
    return files


def portable_manifest_name(raw_name: str, line_number: int) -> str:
    require("\\" not in raw_name and "\x00" not in raw_name, f"SHA256SUMS:{line_number}: unsafe path")
    name = raw_name[2:] if raw_name.startswith("./") else raw_name
    require(name and "/" not in name and name not in {".", ".."}, f"SHA256SUMS:{line_number}: path is not a flat relative name")
    require(not raw_name.startswith("/"), f"SHA256SUMS:{line_number}: absolute path")
    return name


def verify_manifest(root: Path) -> dict[str, dict[str, Any]]:
    files = safe_flat_files(root)
    path = root / MANIFEST_NAME
    require(path.is_file() and not path.is_symlink(), f"missing {MANIFEST_NAME}")
    try:
        raw = path.read_bytes()
        text = raw.decode("utf-8")
    except (OSError, UnicodeError) as exc:
        raise VerificationError(f"cannot read {MANIFEST_NAME}: {exc}") from exc
    require(raw.endswith(b"\n"), f"{MANIFEST_NAME}: final line is not newline-terminated")
    entries: dict[str, dict[str, Any]] = {}
    for number, line in enumerate(text.splitlines(), 1):
        match = PORTABLE_SUM_LINE.fullmatch(line)
        require(match is not None, f"{MANIFEST_NAME}:{number}: non-portable sha256sum line")
        digest, _mode, raw_name = match.groups()
        name = portable_manifest_name(raw_name, number)
        require(name != MANIFEST_NAME, f"{MANIFEST_NAME}:{number}: self reference")
        require(name not in entries, f"{MANIFEST_NAME}:{number}: duplicate path {name}")
        target = root / name
        require(target.is_file() and not target.is_symlink(), f"manifested file missing or unsafe: {name}")
        actual = file_sha256(target)
        require(actual == digest, f"manifest hash mismatch: {name}")
        entries[name] = {
            "sha256": digest,
            "size_bytes": target.stat().st_size,
        }
    require(entries, f"{MANIFEST_NAME}: no entries")
    require(set(files) == set(entries) | {MANIFEST_NAME}, "artifact has missing or unmanifested files")
    return entries


def verify_state(record: dict[str, Any], label: str) -> None:
    side = record.get("side")
    gram_hash = record.get("gram_sha256")
    uid = record.get("uid")
    gram = record.get("gram")
    require(side in {"origin", "transparent", "rootless"}, f"{label}: invalid side")
    require(isinstance(gram_hash, str) and HEX64.fullmatch(gram_hash), f"{label}: bad Gram hash")
    require(uid == f"{side}:{gram_hash}", f"{label}: UID/hash mismatch")
    require(isinstance(gram, list) and len(gram) == 17, f"{label}: Gram is not 17 by 17")
    require(all(isinstance(row, list) and len(row) == 17 for row in gram), f"{label}: Gram is not 17 by 17")
    require(all(is_int(value) for row in gram for value in row), f"{label}: non-integral Gram entry")
    require(payload_sha256(gram) == gram_hash, f"{label}: Gram payload hash mismatch")
    checks = record.get("checks")
    require(isinstance(checks, dict), f"{label}: exact checks are absent")
    required_checks = {
        "rank": 17,
        "determinant": 948,
        "even": True,
        "positive_definite": True,
        "target_genus_verified": True,
    }
    for key, expected in required_checks.items():
        require(checks.get(key) == expected, f"{label}: exact check {key} failed")
    for key in ("root_count_exact_theta", "norm4_count_exact_theta"):
        require(is_int(checks.get(key)) and checks[key] >= 0, f"{label}: invalid {key}")
    require(is_int(record.get("round_discovered")) and record["round_discovered"] >= 0, f"{label}: invalid discovery round")


def verify_move(move: Any, label: str) -> dict[str, Any]:
    require(isinstance(move, dict), f"{label}: move is not an object")
    verify_embedded_hash(move, "move_sha256", label)
    prime = move.get("prime")
    require(is_int(prime) and prime >= 2, f"{label}: invalid neighbor prime")
    for key in ("parent_gram_sha256", "child_gram_sha256"):
        require(isinstance(move.get(key), str) and HEX64.fullmatch(move[key]), f"{label}: invalid {key}")
    require(isinstance(move.get("parent_uid"), str), f"{label}: missing parent UID")
    require(move["parent_uid"].endswith(move["parent_gram_sha256"]), f"{label}: parent UID/hash mismatch")
    gram = move.get("child_gram")
    require(isinstance(gram, list) and payload_sha256(gram) == move["child_gram_sha256"], f"{label}: child Gram hash mismatch")
    relation = move.get("p_neighbor_relation")
    require(isinstance(relation, dict), f"{label}: p-neighbor checks missing")
    expected = {
        "gram_identity_verified": True,
        "transition_denominator": prime,
        "inverse_transition_denominator": prime,
        "rank_mod_p_of_p_times_transition": 1,
        "rank_mod_p_of_p_times_inverse": 1,
        "intersection_index_in_parent": prime,
        "intersection_index_in_child": prime,
    }
    for key, value in expected.items():
        require(relation.get(key) == value, f"{label}: p-neighbor gate {key} failed")
    require(relation.get("transition_determinant") in {"1", "-1"}, f"{label}: non-unit transition")
    child_checks = move.get("child_exact_checks")
    require(isinstance(child_checks, dict), f"{label}: child exact checks missing")
    for key, value in {
        "rank": 17,
        "determinant": 948,
        "even": True,
        "positive_definite": True,
        "target_genus_verified": True,
    }.items():
        require(child_checks.get(key) == value, f"{label}: child check {key} failed")
    return move


def verify_move_exact_arithmetic(
    move: dict[str, Any], parent_gram_data: Any, label: str
) -> dict[str, Any]:
    """Independently replay one serialized neighbor identity with Fraction."""

    move = verify_move(move, label)
    prime = move["prime"]
    require(is_prime_integer(prime), f"{label}: neighbor modulus is not prime")
    parent_gram = validate_exact_gram(parent_gram_data, f"{label} parent Gram")
    raw_neighbor_gram = validate_exact_gram(move.get("raw_neighbor_gram"), f"{label} raw neighbor Gram")
    child_gram = validate_exact_gram(move.get("child_gram"), f"{label} child Gram")
    raw_transform = parse_rational_matrix(move.get("raw_neighbor_transform"), f"{label} raw transform")
    lll_transform = fraction_matrix_from_integers(
        move.get("lll_unimodular_transform"), f"{label} LLL transform"
    )
    total_transform = parse_rational_matrix(
        move.get("parent_to_child_rational_transform"), f"{label} total transform"
    )
    require(
        multiply_matrices(raw_transform, lll_transform) == total_transform,
        f"{label}: raw times LLL does not equal total transform",
    )
    verify_gram_transport(raw_transform, parent_gram, raw_neighbor_gram, f"{label} raw")
    verify_gram_transport(lll_transform, raw_neighbor_gram, child_gram, f"{label} LLL")
    verify_gram_transport(total_transform, parent_gram, child_gram, f"{label} total")
    require(abs(determinant_fraction(raw_transform, f"{label} raw transform")) == 1, f"{label}: raw determinant is not a unit")
    require(abs(determinant_fraction(lll_transform, f"{label} LLL transform")) == 1, f"{label}: LLL determinant is not a unit")
    require(abs(determinant_fraction(total_transform, f"{label} total transform")) == 1, f"{label}: total determinant is not a unit")

    inverse = invert_matrix(total_transform, f"{label} total transform")
    require(matrix_denominator(total_transform) == prime, f"{label}: total denominator is not p")
    require(matrix_denominator(inverse) == prime, f"{label}: inverse denominator is not p")
    scaled = [[value * prime for value in row] for row in total_transform]
    scaled_inverse = [[value * prime for value in row] for row in inverse]
    require(all(value.denominator == 1 for row in scaled for value in row), f"{label}: pT is not integral")
    require(
        all(value.denominator == 1 for row in scaled_inverse for value in row),
        f"{label}: pT^-1 is not integral",
    )
    require(
        rank_mod_prime([[int(value) for value in row] for row in scaled], prime) == 1,
        f"{label}: rank(pT mod p) is not one",
    )
    require(
        rank_mod_prime([[int(value) for value in row] for row in scaled_inverse], prime) == 1,
        f"{label}: rank(pT^-1 mod p) is not one",
    )

    vector = move.get("projective_isotropic_vector")
    require(
        isinstance(vector, list) and len(vector) == 17 and all(is_int(value) for value in vector),
        f"{label}: invalid projective vector",
    )
    numerator = sum(
        vector[row] * parent_gram[row][column] * vector[column]
        for row in range(17)
        for column in range(17)
    )
    require(numerator % 2 == 0, f"{label}: odd quadratic numerator")
    q_value = numerator // 2
    require(q_value == move.get("q_of_vector"), f"{label}: q(vector) mismatch")
    require(q_value % prime == move.get("q_of_vector_mod_p") == 0, f"{label}: vector is not p-isotropic")
    return {
        "total_transform": total_transform,
        "parent_gram": parent_gram,
        "child_gram": child_gram,
        "move_sha256": move["move_sha256"],
    }


def sum_side_field(sides: dict[str, Any], field: str, label: str) -> int:
    total = 0
    for side, stats in sides.items():
        require(side in {"origin", "transparent", "rootless"}, f"{label}: invalid side {side}")
        require(isinstance(stats, dict), f"{label}: side statistics are not an object")
        value = stats.get(field)
        require(is_int(value) and value >= 0, f"{label}: invalid side counter {field}")
        total += value
    return total


def validate_isometry_stats(value: Any, label: str) -> dict[str, Any]:
    require(isinstance(value, dict), f"{label}: isometry statistics are absent")
    for field in ISOMETRY_COUNT_FIELDS:
        require(
            is_int(value.get(field)) and value[field] >= 0,
            f"{label}: invalid {field}",
        )
    require(
        isinstance(value.get("budget_exhausted"), bool),
        f"{label}: invalid budget_exhausted",
    )
    return value


def verify_rounds(
    rounds: list[dict[str, Any]], result: dict[str, Any], initial_count: int
) -> None:
    counters = result.get("counters")
    require(isinstance(counters, dict), "result: counters are absent")
    for key in COUNTER_FIELDS:
        require(is_int(counters.get(key)) and counters[key] >= 0, f"result: invalid counter {key}")
    require(isinstance(counters.get("move_budget_exhausted"), bool), "result: invalid move-budget gate")
    previous = {key: 0 for key in COUNTER_FIELDS}
    previous["unique_states_discovered"] = initial_count
    previous_budget = False
    previous_isometries = 0
    result_isometry_stats = result.get("isometry_comparisons")
    modern_isometry_stats = result_isometry_stats is not None
    previous_stats = {
        **{field: 0 for field in ISOMETRY_COUNT_FIELDS},
        "budget_exhausted": False,
    }
    if modern_isometry_stats:
        validate_isometry_stats(result_isometry_stats, "result")
    for expected_round, record in enumerate(rounds, 1):
        label = f"rounds.jsonl:{expected_round}"
        require(record.get("round") == expected_round, f"{label}: non-sequential round")
        sides = record.get("sides")
        cumulative = record.get("cumulative")
        require(isinstance(sides, dict) and sides, f"{label}: missing side statistics")
        require(isinstance(cumulative, dict), f"{label}: missing cumulative counters")
        for key in COUNTER_FIELDS:
            require(is_int(cumulative.get(key)) and cumulative[key] >= previous[key], f"{label}: non-monotone {key}")
            side_key = SIDE_COUNTER_FIELDS[key]
            require(
                cumulative[key] - previous[key] == sum_side_field(sides, side_key, label),
                f"{label}: side totals disagree on {key}",
            )
        budget = cumulative.get("move_budget_exhausted")
        require(isinstance(budget, bool) and (not previous_budget or budget), f"{label}: move-budget gate regressed")
        iso = record.get("isometry_checks")
        require(is_int(iso) and iso >= previous_isometries, f"{label}: non-monotone isometry count")
        if modern_isometry_stats:
            current_stats = validate_isometry_stats(
                record.get("isometry_comparisons"), label
            )
            for field in ISOMETRY_COUNT_FIELDS:
                require(
                    current_stats[field] >= previous_stats[field],
                    f"{label}: non-monotone {field}",
                )
            require(
                not previous_stats["budget_exhausted"]
                or current_stats["budget_exhausted"],
                f"{label}: isometry-budget gate regressed",
            )
            require(
                current_stats["skipped_due_to_budget"] == 0
                and current_stats["budget_exhausted"] is False,
                f"{label}: importable round skipped an isometry pair",
            )
            previous_stats = dict(current_stats)
        require(isinstance(record.get("bridge_found"), bool), f"{label}: invalid bridge gate")
        previous = {key: cumulative[key] for key in COUNTER_FIELDS}
        previous_budget = budget
        previous_isometries = iso
    require(previous == {key: counters[key] for key in COUNTER_FIELDS}, "final round counters disagree with result")
    require(previous_budget == counters["move_budget_exhausted"], "final move-budget gate disagrees")
    require(previous_isometries == result.get("isometry_checks"), "final isometry count disagrees")
    if modern_isometry_stats:
        require(
            previous_stats == result_isometry_stats,
            "final round isometry statistics disagree with result",
        )


def verify_claim_boundary(result: dict[str, Any]) -> dict[str, Any]:
    status = result.get("status")
    modern = is_modern_result(result)
    claims = result.get("claim_boundary")
    require(isinstance(claims, dict), "result: claim boundary is absent")
    required = {"exact_lattice_bridge_found", *FALSE_CLAIM_GATES, "negative_result_scope"}
    if modern:
        required.add(CURRENT_TARGET_CLAIM_GATE)
    require(required <= set(claims), "result: required claim-boundary fields are absent")
    require(claims["exact_lattice_bridge_found"] is (status == "bridge_found"), "result: lattice-bridge claim disagrees with status")
    for key in FALSE_CLAIM_GATES:
        require(claims[key] is False, f"result: forbidden promotion {key}")
    if modern:
        require(
            claims[CURRENT_TARGET_CLAIM_GATE] is False,
            f"result: forbidden promotion {CURRENT_TARGET_CLAIM_GATE}",
        )
    require(isinstance(claims["negative_result_scope"], str) and claims["negative_result_scope"], "result: finite-search scope is absent")
    if modern:
        require(
            {
                "p2_neighbor_enumeration_complete",
                "p2_neighbor_enumeration_boundary",
                "bounded_negative_result_importable",
            }
            <= set(claims),
            "result: modern p=2/importability claim boundary is incomplete",
        )
        require(
            claims.get("bounded_negative_result_importable") is True,
            "result: modern bounded-result importability gate is absent or false",
        )
    elif "bounded_negative_result_importable" in claims:
        require(claims["bounded_negative_result_importable"] is True, "result: search marks bounded result non-importable")
    if "p2_neighbor_enumeration_complete" in claims:
        primes = result.get("configuration", {}).get("primes", [])
        expected = False if 2 in primes else None
        require(claims["p2_neighbor_enumeration_complete"] is expected, "result: p=2 completeness boundary drift")
        boundary = claims.get("p2_neighbor_enumeration_boundary")
        require(
            (2 not in primes and boundary is None)
            or (2 in primes and boundary == P2_ENUMERATION_BOUNDARY),
            "result: p=2 enumeration limitation is absent",
        )
        if 2 in primes:
            require(
                "not excluded" in boundary,
                "result: p=2 boundary-failure lines are not explicitly outside the exclusion",
            )
            require(
                "not excluded" in claims["negative_result_scope"],
                "result: finite-negative scope silently includes failed p=2 lines",
            )
    return claims


def certificate_claim_boundary(
    claims: dict[str, Any],
    *,
    successful_moves: int,
    classified_p2_boundary_failures: int,
    unclassified_construction_errors: int,
) -> dict[str, Any]:
    """Keep legacy evidence exact while avoiding a stale global frontier claim."""

    normalized = copy.deepcopy(claims)
    normalized[CURRENT_TARGET_CLAIM_GATE] = False
    normalized["negative_result_scope"] = normalized["negative_result_scope"].replace(
        "the p-neighbor graph and rank-31 problem remain open.",
        "the p-neighbor graph remains open; this artifact does not settle any "
        "elliptic-rank record problem.",
    )
    if (
        not normalized["exact_lattice_bridge_found"]
        and (classified_p2_boundary_failures or unclassified_construction_errors)
    ):
        scope = (
            f"Only the {successful_moves} successfully constructed moves in the "
            "serialized finite window are excluded. "
        )
        if classified_p2_boundary_failures:
            scope += (
                f"The {classified_p2_boundary_failures} classified p=2 boundary-failure "
                "attempt(s) are not excluded. "
            )
        if unclassified_construction_errors:
            scope += (
                f"The {unclassified_construction_errors} unclassified construction-error "
                "attempt(s) are not excluded. "
            )
        normalized["negative_result_scope"] = (
            scope
            + "The p-neighbor graph remains open; this artifact does not settle any "
            "elliptic-rank record problem."
        )
    return normalized


def verify_modern_configuration_and_termination(
    result: dict[str, Any], checkpoint: dict[str, Any]
) -> None:
    if not is_modern_result(result):
        return
    config = result.get("configuration")
    require(isinstance(config, dict), "modern result lacks configuration")
    mode = config.get("mode")
    expected_sides = (
        ["origin"]
        if mode == "forward"
        else ["origin", "transparent", "rootless"]
        if mode == "bidirectional"
        else None
    )
    require(expected_sides is not None, "modern result has invalid search mode")
    require(
        config.get("side_expansion_order") == expected_sides,
        "modern result side-expansion order drift",
    )
    primes = config.get("primes")
    require(
        isinstance(primes, list)
        and primes
        and all(is_prime_integer(prime) for prime in primes)
        and len(primes) == len(set(primes)),
        "modern result has invalid neighbor primes",
    )
    expected_boundaries = (
        [{"prime": 2, "complete": False, "reason": P2_ENUMERATION_BOUNDARY}]
        if 2 in primes
        else []
    )
    require(
        config.get("prime_enumeration_boundaries") == expected_boundaries,
        "modern result prime-enumeration boundary drift",
    )

    status = result["status"]
    reason = result.get("termination_reason")
    allowed = {
        "bridge_found": {"exact_bridge_found"},
        "bounded_search_completed": {
            "requested_rounds_completed",
            "successful_move_budget_exhausted",
            "all_active_beams_empty",
        },
    }
    require(reason in allowed[status], "modern result has an inconsistent termination reason")
    counters = result["counters"]
    max_successful_moves = config.get("max_successful_moves")
    max_isometry_checks = config.get("max_isometry_checks")
    require(
        is_int(max_successful_moves) and max_successful_moves >= 0,
        "modern result has an invalid successful-move budget",
    )
    require(
        is_int(max_isometry_checks) and max_isometry_checks >= 0,
        "modern result has an invalid qfisom budget",
    )
    require(
        result["isometry_comparisons"]["qfisom_checked"]
        <= max_isometry_checks,
        "modern result exceeded its qfisom budget",
    )
    if reason == "successful_move_budget_exhausted":
        require(
            counters.get("move_budget_exhausted") is True
            and counters.get("successful_moves") == max_successful_moves,
            "move-budget termination without an exactly exhausted budget",
        )
    elif status == "bounded_search_completed":
        require(
            counters.get("move_budget_exhausted") is False
            and counters.get("successful_moves") <= max_successful_moves,
            "bounded result hid or exceeded its move budget",
        )
    if reason == "requested_rounds_completed":
        require(
            config.get("rounds_completed") == config.get("rounds_requested"),
            "requested-round termination occurred early",
        )
    if reason == "all_active_beams_empty":
        beams = checkpoint.get("current_beams")
        require(
            isinstance(beams, dict)
            and set(beams) == set(expected_sides)
            and all(isinstance(beams[side], list) and not beams[side] for side in expected_sides),
            "empty-beam termination has a nonempty or missing active beam",
        )


def verify_bridge(
    root: Path,
    result: dict[str, Any],
    state_by_uid: dict[str, dict[str, Any]],
    logged_moves_by_hash: dict[str, dict[str, Any]],
) -> dict[str, Any] | None:
    status = result["status"]
    modern = is_modern_result(result)
    bridge = result.get("bridge")
    bridge_path = root / "exact-bridge.json"
    if status == "bounded_search_completed":
        require(bridge is None, "bounded result embeds a bridge")
        require(not bridge_path.exists(), "bounded result unexpectedly has exact-bridge.json")
        return None
    require(isinstance(bridge, dict), "bridge result is missing its exact bridge")
    require(bridge_path.is_file() and not bridge_path.is_symlink(), "exact-bridge.json is missing or unsafe")
    disk_bridge = load_json(bridge_path)
    require(disk_bridge == bridge, "embedded bridge and exact-bridge.json disagree")
    verify_embedded_hash(bridge, "bridge_sha256", "exact-bridge.json")
    require(bridge.get("status") == "exact_neighbor_bridge_found", "exact bridge status drift")
    require(bridge.get("all_parent_child_gram_identities_verified") is True, "bridge lacks Gram identities")
    require(bridge.get("all_neighbor_intersection_indices_verified") is True, "bridge lacks neighbor-index gates")
    require(bridge.get("endpoint_identification_verified_by_initial_source") is True, "bridge endpoint is not identified")
    origin_moves = bridge.get("origin_forward_moves")
    endpoint_moves = bridge.get("endpoint_forward_moves_to_invert")
    require(isinstance(origin_moves, list) and isinstance(endpoint_moves, list), "bridge move lists are missing")
    require(bridge.get("neighbor_step_count") == len(origin_moves) + len(endpoint_moves), "bridge step count mismatch")
    origin_chain = bridge.get("origin_chain")
    endpoint_chain = bridge.get("endpoint_chain_endpoint_to_meeting")
    require(isinstance(origin_chain, list) and origin_chain, "origin bridge chain is missing")
    require(isinstance(endpoint_chain, list) and endpoint_chain, "endpoint bridge chain is missing")
    require(len(origin_chain) == len(origin_moves) + 1, "origin bridge chain length mismatch")
    require(len(endpoint_chain) == len(endpoint_moves) + 1, "endpoint bridge chain length mismatch")
    require(origin_chain[0].get("side") == "origin", "origin chain does not start at origin")
    require(
        endpoint_chain[0].get("side") in {"transparent", "rootless"}
        and bridge.get("endpoint") == endpoint_chain[0].get("side"),
        "bridge endpoint label disagrees with its initial chain state",
    )

    def verify_path(
        chain: list[dict[str, Any]], moves: list[dict[str, Any]], label: str
    ) -> list[list[Fraction]]:
        composite = identity_matrix()
        for state_index, state in enumerate(chain):
            require(isinstance(state, dict) and isinstance(state.get("uid"), str), f"{label}: invalid state")
            logged_state = state_by_uid.get(state["uid"])
            require(logged_state == state, f"{label}: state is absent from or differs from states.jsonl")
            if state_index == 0:
                require(state.get("parent_uid") is None, f"{label}: initial state has a parent")
                continue
            move = moves[state_index - 1]
            move_hash = move.get("move_sha256") if isinstance(move, dict) else None
            require(
                isinstance(move_hash, str) and logged_moves_by_hash.get(move_hash) == move,
                f"{label}: bridge move is absent from or differs from moves.jsonl",
            )
            require(move.get("parent_uid") == chain[state_index - 1]["uid"], f"{label}: broken move parent")
            require(move.get("child_gram_sha256") == state.get("gram_sha256"), f"{label}: broken move child")
            exact = verify_move_exact_arithmetic(
                move,
                chain[state_index - 1].get("gram"),
                f"{label} move {state_index}",
            )
            composite = multiply_matrices(composite, exact["total_transform"])
        verify_gram_transport(
            composite,
            validate_exact_gram(chain[0].get("gram"), f"{label} initial Gram"),
            validate_exact_gram(chain[-1].get("gram"), f"{label} meeting Gram"),
            f"{label} composed path",
        )
        return composite

    origin_composite = verify_path(origin_chain, origin_moves, "origin bridge path")
    endpoint_composite = verify_path(endpoint_chain, endpoint_moves, "endpoint bridge path")
    isometry = bridge.get("meeting_isometry")
    require(isinstance(isometry, dict), "bridge meeting isometry is absent")
    require(isometry.get("gram_identity_verified") is True, "meeting isometry lacks Gram identity")
    require(isometry.get("determinant") in {-1, 1}, "meeting isometry is not unimodular")
    require(isometry.get("from_uid") == origin_chain[-1]["uid"], "meeting isometry origin UID mismatch")
    require(isometry.get("to_uid") == endpoint_chain[-1]["uid"], "meeting isometry endpoint UID mismatch")
    expected_meeting_fingerprint = [
        origin_chain[-1]["checks"]["root_count_exact_theta"],
        origin_chain[-1]["checks"]["norm4_count_exact_theta"],
    ]
    require(
        isometry.get("fingerprint") == expected_meeting_fingerprint
        == [
            endpoint_chain[-1]["checks"]["root_count_exact_theta"],
            endpoint_chain[-1]["checks"]["norm4_count_exact_theta"],
        ],
        "meeting isometry fingerprint mismatch",
    )
    isometry_transform = fraction_matrix_from_integers(
        isometry.get("integral_unimodular_transform"), "meeting isometry transform"
    )
    require(
        determinant_fraction(isometry_transform, "meeting isometry transform")
        == isometry["determinant"],
        "meeting isometry determinant mismatch",
    )
    verify_gram_transport(
        isometry_transform,
        validate_exact_gram(origin_chain[-1].get("gram"), "origin meeting Gram"),
        validate_exact_gram(endpoint_chain[-1].get("gram"), "endpoint meeting Gram"),
        "meeting isometry",
    )

    endpoint_inverse = invert_matrix(endpoint_composite, "endpoint path composite")
    initial_composite = multiply_matrices(
        multiply_matrices(origin_composite, isometry_transform), endpoint_inverse
    )
    origin_initial_gram = validate_exact_gram(origin_chain[0].get("gram"), "origin initial Gram")
    endpoint_initial_gram = validate_exact_gram(endpoint_chain[0].get("gram"), "endpoint initial Gram")
    verify_gram_transport(
        initial_composite,
        origin_initial_gram,
        endpoint_initial_gram,
        "complete bridge composite",
    )
    require(
        abs(determinant_fraction(initial_composite, "complete bridge composite")) == 1,
        "complete bridge determinant is not a unit",
    )
    serial_composite = serial_fraction_matrix(initial_composite)
    recorded_composite = bridge.get("composed_initial_basis_transform")
    if recorded_composite is not None:
        require(recorded_composite == serial_composite, "recorded complete bridge matrix mismatch")
        recorded_hash = bridge.get("composed_initial_basis_transform_sha256")
        require(recorded_hash == payload_sha256(serial_composite), "recorded complete bridge hash mismatch")
    transport = bridge.get("end_to_end_transport")
    if modern:
        require(isinstance(transport, dict), "modern bridge lacks end-to-end transport")
    if transport is not None:
        require(isinstance(transport, dict), "end-to-end transport is not an object")
        verify_embedded_hash(transport, "transport_sha256", "end-to-end transport")
        require(
            transport.get("origin_forward_composite_C_origin")
            == serial_fraction_matrix(origin_composite),
            "recorded origin composite mismatch",
        )
        require(
            transport.get("endpoint_forward_composite_C_endpoint")
            == serial_fraction_matrix(endpoint_composite),
            "recorded endpoint composite mismatch",
        )
        require(
            transport.get("endpoint_initial_to_origin_initial_transform_M")
            == serial_composite,
            "recorded end-to-end matrix mismatch",
        )
        require(
            transport.get("endpoint_initial_to_origin_initial_transform_M_sha256")
            == payload_sha256(serial_composite),
            "recorded end-to-end matrix hash mismatch",
        )
        require(
            transport.get("M_transpose_G_origin_initial_M_equals_G_endpoint_initial") is True
            and transport.get("M_determinant_is_plus_or_minus_one") is True,
            "recorded end-to-end gates failed",
        )
        complete_determinant = determinant_fraction(
            initial_composite, "complete bridge composite"
        )
        require(
            transport.get("M_determinant") == str(complete_determinant),
            "recorded end-to-end determinant mismatch",
        )
        convention = transport.get("matrix_convention")
        require(
            isinstance(convention, str)
            and "M = C_origin Q C_endpoint^-1" in convention,
            "end-to-end matrix convention is absent or ambiguous",
        )
    boundary = bridge.get("claim_boundary")
    require(
        isinstance(boundary, str)
        and "not yet" in boundary
        and (
            "elliptic rank-record curve" in boundary
            or "rank-31 curve" in boundary
        ),
        "bridge claim boundary was overpromoted",
    )
    return {
        "bridge_sha256": bridge["bridge_sha256"],
        "endpoint": bridge.get("endpoint"),
        "neighbor_step_count": bridge["neighbor_step_count"],
        "meeting_isometry_verified": True,
        "all_bridge_move_matrix_identities_recomputed_with_fraction": True,
        "complete_bridge_composite_verified": True,
        "complete_bridge_composite_sha256": payload_sha256(serial_composite),
    }


def run_git(source_root: Path, arguments: list[str]) -> str:
    try:
        completed = subprocess.run(
            ["git", "-C", str(source_root), *arguments],
            check=False,
            capture_output=True,
            text=True,
            timeout=15,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise VerificationError(f"cannot inspect source checkout with git: {exc}") from exc
    require(
        completed.returncode == 0,
        f"git source-checkout inspection failed: {completed.stderr.strip()}",
    )
    return completed.stdout.strip()


def resolve_source_checkout(
    source_root: Path | None, provenance: dict[str, Any]
) -> tuple[Path, dict[str, Any]]:
    importer_root = Path(__file__).resolve().parents[1]
    if source_root is None:
        return importer_root, {
            "mode": "executing_importer_repository",
            "git_head_verified": False,
        }
    require(not source_root.is_symlink(), "explicit source root must not be a symlink")
    resolved = source_root.resolve()
    require(resolved.is_dir(), "explicit source root is not a directory")
    top_level = Path(run_git(resolved, ["rev-parse", "--show-toplevel"])).resolve()
    require(top_level == resolved, "explicit source root is not the checkout top level")
    head = run_git(resolved, ["rev-parse", "--verify", "HEAD"])
    require(head == provenance.get("commit_sha"), "source checkout HEAD disagrees with provenance commit")
    tracked_status = run_git(
        resolved,
        [
            "status",
            "--porcelain=v1",
            "--untracked-files=all",
            "--",
            SEARCH_SCRIPT_PATH,
            TARGET_SCRIPT_PATH,
            SAMPLER_SCRIPT_PATH,
        ],
    )
    require(not tracked_status, "source checkout has modified or untracked source evidence")
    return resolved, {
        "mode": "explicit_historical_git_checkout",
        "git_head_verified": True,
        "git_head": head,
    }


def verify_source_files(
    result: dict[str, Any], source_root: Path
) -> dict[str, str]:
    sources = result.get("source_files")
    require(isinstance(sources, dict), "result: source-file provenance missing")
    require(sources.get("search_script") == SEARCH_SCRIPT_PATH, "result: unexpected search script")
    require(sources.get("target_formula_script") == TARGET_SCRIPT_PATH, "result: unexpected target script")
    verified: dict[str, str] = {}
    pairs = [
        ("search_script", "search_script_sha256"),
        ("target_formula_script", "target_formula_script_sha256"),
    ]
    if result["schema_version"] == 2:
        require(sources.get("projective_sampler_script") == SAMPLER_SCRIPT_PATH, "result: unexpected sampler source")
        pairs.append(("projective_sampler_script", "projective_sampler_script_sha256"))
    for path_key, hash_key in pairs:
        relative = sources[path_key]
        digest = sources.get(hash_key)
        require(isinstance(digest, str) and HEX64.fullmatch(digest), f"result: invalid {hash_key}")
        local = source_root / relative
        require(local.is_file() and not local.is_symlink(), f"local source missing or unsafe: {relative}")
        require(file_sha256(local) == digest, f"artifact was generated by a different {relative}")
        if relative == SAMPLER_SCRIPT_PATH:
            require(file_sha256(REPO_ROOT / SAMPLER_SCRIPT_PATH) == digest, "sampler replay implementation differs from the producer")
        verified[relative] = digest
    return verified


def verify_runtime_provenance(
    result: dict[str, Any], supplied: dict[str, Any], modern: bool
) -> dict[str, Any] | None:
    runtime = result.get("runtime_provenance")
    if runtime is None:
        require(not modern, "modern result lacks embedded runtime provenance")
        return None
    require(isinstance(runtime, dict), "result: runtime provenance is not an object")
    require(runtime.get("github_actions") == "true", "result: runtime was not identified as GitHub Actions")
    expected = {
        "github_repository": supplied.get("repository"),
        "github_run_id": supplied.get("run_id"),
        "github_run_attempt": str(supplied.get("run_attempt")),
        "github_sha": supplied.get("commit_sha"),
        "github_ref": supplied.get("git_ref"),
        "github_run_url": supplied.get("run_url"),
    }
    for key, value in expected.items():
        require(runtime.get(key) == value, f"result: runtime provenance mismatch for {key}")
    require(runtime.get("github_server_url") == "https://github.com", "result: unexpected GitHub server")
    require(runtime.get("github_workflow") == EXPECTED_WORKFLOW_NAME, "result: workflow name drift")
    require(runtime.get("github_job") == EXPECTED_JOB_NAME, "result: workflow job drift")
    require(runtime.get("runner_os") == "Linux", "result: unexpected runner OS")
    require(runtime.get("runner_arch") == "X64", "result: unexpected runner architecture")
    require(
        runtime.get("search_container_image") == EXPECTED_CONTAINER_IMAGE,
        "result: Sage container image drift",
    )
    workflow_ref = runtime.get("github_workflow_ref")
    expected_workflow_ref = (
        f'{supplied.get("repository")}/{WORKFLOW_PATH}@{supplied.get("git_ref")}'
    )
    require(
        workflow_ref == expected_workflow_ref,
        "result: runtime workflow ref does not identify the expected workflow",
    )
    require(
        runtime.get("search_container_digest") == EXPECTED_CONTAINER_DIGEST,
        "result: pinned Sage container digest drift",
    )
    return copy.deepcopy(runtime)


def verify_artifact(
    root: Path,
    provenance: dict[str, Any],
    source_root: Path | None = None,
) -> dict[str, Any]:
    root = root.resolve()
    manifest = verify_manifest(root)
    result = load_json(root / "result.json")
    require(isinstance(result, dict), "result.json is not an object")
    require(type(result.get("schema_version")) is int and result["schema_version"] in {1, 2}, "result schema version drift")
    hash_sampling = verify_selection_schema(result)
    require(result.get("status") in {"bridge_found", "bounded_search_completed"}, "search did not finish with an importable status")
    modern = is_modern_result(result)
    if modern:
        require(
            all(field in result for field in MODERN_RESULT_FIELDS),
            "result mixes legacy and modern schemas",
        )
        require(result["artifact_importable"] is True, "search result is explicitly non-importable")
        require(
            isinstance(result.get("termination_reason"), str)
            and result["termination_reason"],
            "modern result lacks a termination reason",
        )
    verify_embedded_hash(result, "result_payload_sha256", "result.json")
    require(str(result.get("sage_version", "")).startswith("SageMath version 10.9"), "result was not produced by pinned Sage 10.9")
    require(isinstance(result.get("truth_status"), str) and result["truth_status"], "result truth status missing")
    if modern:
        validate_isometry_stats(result.get("isometry_comparisons"), "result")

    required_files = set(BASE_ARTIFACT_FILES)
    evidence_names = set(RESULT_EVIDENCE_FILES)
    if hash_sampling:
        required_files.add(SAMPLING_LOG_NAME)
        evidence_names.add(SAMPLING_LOG_NAME)
    if result["status"] == "bridge_found":
        required_files.add("exact-bridge.json")
        evidence_names.add("exact-bridge.json")
    require(set(manifest) == required_files, "artifact file schema drift or incomplete manifest")
    evidence = result.get("evidence_sha256")
    require(isinstance(evidence, dict) and set(evidence) == evidence_names, "result evidence-file schema drift")
    for name, digest in evidence.items():
        require(isinstance(digest, str) and HEX64.fullmatch(digest), f"result: invalid evidence hash for {name}")
        require(manifest[name]["sha256"] == digest, f"result evidence hash mismatch: {name}")

    jsonl_names = (*JSONL_NAMES, SAMPLING_LOG_NAME) if hash_sampling else JSONL_NAMES
    jsonl = {name: load_jsonl(root / name) for name in jsonl_names}
    states = jsonl["states.jsonl"]
    moves = jsonl["moves.jsonl"]
    attempts = jsonl["attempts.jsonl"]
    rounds = jsonl["rounds.jsonl"]
    isometries = jsonl["isometry-checks.jsonl"]

    initial = result.get("initial_states")
    require(isinstance(initial, list) and len(initial) == 3, "result: exactly three initial states required")
    require(states[:3] == initial, "states.jsonl does not begin with the recorded initial states")
    for index, state in enumerate(states):
        verify_state(state, f"states.jsonl:{index + 1}")
    uids = [state["uid"] for state in states]
    require(len(uids) == len(set(uids)), "states.jsonl contains duplicate UIDs")
    state_by_uid = {state["uid"]: state for state in states}
    initial_by_side = {state["side"]: state for state in initial}
    require(set(initial_by_side) == set(INITIAL_HASHES), "initial sides drifted")
    for side, digest in INITIAL_HASHES.items():
        require(initial_by_side[side]["gram_sha256"] == digest, f"initial {side} lattice hash drifted")
        require(initial_by_side[side]["round_discovered"] == 0 and initial_by_side[side]["parent_uid"] is None, f"initial {side} has a parent")

    attempt_by_index: dict[int, dict[str, Any]] = {}
    expected_p2_boundary_attempts = 0
    legacy_exact_p2_boundary_attempts = 0
    legacy_unclassified_construction_errors = 0
    for index, attempt in enumerate(attempts, 1):
        require(attempt.get("attempt_index") == index, f"attempts.jsonl:{index}: non-sequential attempt")
        require(
            attempt.get("status")
            in {
                "construction_error",
                "expected_construction_boundary",
                "success_duplicate",
                "success_new_state",
            },
            f"attempts.jsonl:{index}: invalid status",
        )
        explicit_expected_p2_boundary = (
            attempt["status"] == "expected_construction_boundary"
            or attempt.get("construction_error_expected_and_bounded") is True
        )
        if "construction_error_expected_and_bounded" in attempt:
            require(
                isinstance(attempt["construction_error_expected_and_bounded"], bool),
                f"attempts.jsonl:{index}: invalid expected-error gate",
            )
            require(
                attempt["status"] == "construction_error",
                f"attempts.jsonl:{index}: expected-error compatibility flag has the wrong status",
            )
        if explicit_expected_p2_boundary:
            expected_p2_boundary_attempts += 1
            require(attempt.get("prime") == 2, f"attempts.jsonl:{index}: boundary failure is not p=2")
            require(
                attempt.get("error_category")
                == "expected_sage_p2_nonmaximal_even_boundary",
                f"attempts.jsonl:{index}: unexpected boundary category",
            )
            require(
                attempt.get("error_message")
                == EXPECTED_P2_ERROR,
                f"attempts.jsonl:{index}: unexpected p=2 error",
            )
            require(
                attempt.get("error_type") == "ValueError",
                f"attempts.jsonl:{index}: unexpected p=2 error type",
            )
            require(
                attempt.get("prime_neighbor_enumeration_complete") is False,
                f"attempts.jsonl:{index}: p=2 failure overclaims completeness",
            )
            require(
                attempt.get("enumeration_boundary") == P2_ENUMERATION_BOUNDARY,
                f"attempts.jsonl:{index}: p=2 enumeration boundary drift",
            )
        elif attempt["status"] == "construction_error":
            if (
                attempt.get("prime") == 2
                and attempt.get("error") == EXPECTED_P2_ERROR_REPR
            ):
                legacy_exact_p2_boundary_attempts += 1
            else:
                legacy_unclassified_construction_errors += 1
        attempt_by_index[index] = attempt
    successful_attempts = [attempt for attempt in attempts if attempt["status"].startswith("success_")]
    require(len(moves) == len(successful_attempts), "move count disagrees with successful attempts")
    new_state_move_hashes: set[str] = set()
    logged_moves_by_hash: dict[str, dict[str, Any]] = {}
    for index, record in enumerate(moves, 1):
        label = f"moves.jsonl:{index}"
        require(record.get("status") in {"success_duplicate", "success_new_state"}, f"{label}: invalid status")
        attempt = record.get("attempt")
        require(isinstance(attempt, dict) and is_int(attempt.get("attempt_index")), f"{label}: attempt link missing")
        logged = attempt_by_index.get(attempt["attempt_index"])
        require(logged is not None and logged["status"] == record["status"], f"{label}: attempt status mismatch")
        for key, value in attempt.items():
            require(logged.get(key) == value, f"{label}: attempt field {key} mismatch")
        move = verify_move(record.get("move"), label)
        require(move["move_sha256"] not in logged_moves_by_hash, f"{label}: duplicate move hash")
        logged_moves_by_hash[move["move_sha256"]] = move
        require(logged.get("move_sha256") == move["move_sha256"], f"{label}: move hash link mismatch")
        require(logged.get("child_gram_sha256") == move["child_gram_sha256"], f"{label}: child hash link mismatch")
        require(move["parent_uid"] == attempt.get("parent_uid"), f"{label}: parent link mismatch")
        require(move["parent_uid"] in state_by_uid, f"{label}: parent state missing")
        require(move["prime"] == attempt.get("prime"), f"{label}: neighbor-prime link mismatch")
        require(
            move.get("projective_line_ordinal_one_based")
            == attempt.get("projective_line_ordinal_one_based"),
            f"{label}: projective-line ordinal mismatch",
        )
        require(
            move.get("projective_isotropic_vector")
            == attempt.get("projective_isotropic_vector"),
            f"{label}: projective vector mismatch",
        )
        duplicate = record.get("duplicate_state_on_same_side")
        require(duplicate is (record["status"] == "success_duplicate"), f"{label}: duplicate gate mismatch")
        if not duplicate:
            child_uid = f'{attempt.get("side")}:{move["child_gram_sha256"]}'
            child = state_by_uid.get(child_uid)
            require(child is not None, f"{label}: new child state missing")
            require(child["parent_uid"] == move["parent_uid"], f"{label}: child parent link mismatch")
            require(child["incoming_move_sha256"] == move["move_sha256"], f"{label}: child move link mismatch")
            require(child["round_discovered"] == attempt.get("round"), f"{label}: child round mismatch")
            require(move["move_sha256"] not in new_state_move_hashes, f"{label}: duplicate new-state move hash")
            new_state_move_hashes.add(move["move_sha256"])

    for state_record in states[len(initial) :]:
        incoming = state_record.get("incoming_move_sha256")
        require(
            isinstance(incoming, str) and incoming in new_state_move_hashes,
            f'state {state_record["uid"]}: no unique successful incoming move',
        )

    modern_isometry_log = any("pair_index" in record for record in isometries)
    require(
        not modern_isometry_log or all("pair_index" in record for record in isometries),
        "isometry log mixes legacy and pair-index schemas",
    )
    if modern:
        require(
            not isometries or modern_isometry_log,
            "modern result contains a legacy isometry log",
        )
    qfisom_log_count = 0
    identical_log_count = 0
    for index, record in enumerate(isometries, 1):
        if modern_isometry_log:
            require(record.get("pair_index") == index, f"isometry-checks.jsonl:{index}: non-sequential pair")
            method = record.get("comparison_method")
            require(
                method
                in {
                    "pari_qfisom",
                    "identical_serialized_gram",
                },
                f"isometry-checks.jsonl:{index}: incomplete/unknown comparison method",
            )
            if method == "pari_qfisom":
                qfisom_log_count += 1
                require(
                    record.get("qfisom_check_index") == qfisom_log_count
                    and record.get("check_index") == qfisom_log_count,
                    f"isometry-checks.jsonl:{index}: non-sequential qfisom check",
                )
            else:
                identical_log_count += 1
                require(
                    record.get("qfisom_check_index") is None
                    and record.get("check_index") is None,
                    f"isometry-checks.jsonl:{index}: identical Gram consumed qfisom",
                )
        else:
            require(record.get("check_index") == index, f"isometry-checks.jsonl:{index}: non-sequential check")
            qfisom_log_count += 1
        require(isinstance(record.get("isometric"), bool), f"isometry-checks.jsonl:{index}: invalid gate")
        origin_uid = record.get("origin_uid")
        endpoint_uid = record.get("endpoint_uid")
        require(origin_uid in state_by_uid and endpoint_uid in state_by_uid, f"isometry-checks.jsonl:{index}: unknown state")
        require(state_by_uid[origin_uid]["side"] == "origin", f"isometry-checks.jsonl:{index}: left state is not origin")
        require(
            state_by_uid[endpoint_uid]["side"] in {"transparent", "rootless"},
            f"isometry-checks.jsonl:{index}: right state is not an endpoint",
        )
        require(
            record.get("endpoint_side") == state_by_uid[endpoint_uid]["side"],
            f"isometry-checks.jsonl:{index}: endpoint side mismatch",
        )
        expected_fingerprint = [
            state_by_uid[origin_uid]["checks"]["root_count_exact_theta"],
            state_by_uid[origin_uid]["checks"]["norm4_count_exact_theta"],
        ]
        require(record.get("fingerprint") == expected_fingerprint, f"isometry-checks.jsonl:{index}: fingerprint mismatch")
        require(
            expected_fingerprint
            == [
                state_by_uid[endpoint_uid]["checks"]["root_count_exact_theta"],
                state_by_uid[endpoint_uid]["checks"]["norm4_count_exact_theta"],
            ],
            f"isometry-checks.jsonl:{index}: states have different fingerprints",
        )
        if (
            modern_isometry_log
            and record.get("comparison_method") == "identical_serialized_gram"
        ):
            require(
                state_by_uid[origin_uid]["gram_sha256"]
                == state_by_uid[endpoint_uid]["gram_sha256"],
                f"isometry-checks.jsonl:{index}: identical-Gram method used for different Grams",
            )
            require(
                state_by_uid[origin_uid]["gram"]
                == state_by_uid[endpoint_uid]["gram"],
                f"isometry-checks.jsonl:{index}: identical hash has different serialized Gram",
            )
            require(record["isometric"] is True, f"isometry-checks.jsonl:{index}: identical Grams rejected")
        if record["isometric"]:
            iso = record.get("isometry")
            require(isinstance(iso, dict) and iso.get("gram_identity_verified") is True, f"isometry-checks.jsonl:{index}: missing exact isometry")
            require(iso.get("determinant") in {-1, 1}, f"isometry-checks.jsonl:{index}: non-unimodular isometry")
            require(
                iso.get("from_uid") == origin_uid
                and iso.get("to_uid") == endpoint_uid,
                f"isometry-checks.jsonl:{index}: isometry UID mismatch",
            )
            require(
                iso.get("fingerprint") == expected_fingerprint,
                f"isometry-checks.jsonl:{index}: isometry fingerprint mismatch",
            )
            transform = fraction_matrix_from_integers(
                iso.get("integral_unimodular_transform"),
                f"isometry-checks.jsonl:{index} transform",
            )
            require(
                determinant_fraction(
                    transform, f"isometry-checks.jsonl:{index} transform"
                )
                == iso["determinant"],
                f"isometry-checks.jsonl:{index}: isometry determinant mismatch",
            )
            verify_gram_transport(
                transform,
                validate_exact_gram(
                    state_by_uid[origin_uid]["gram"],
                    f"isometry-checks.jsonl:{index} origin Gram",
                ),
                validate_exact_gram(
                    state_by_uid[endpoint_uid]["gram"],
                    f"isometry-checks.jsonl:{index} endpoint Gram",
                ),
                f"isometry-checks.jsonl:{index}",
            )
            if (
                modern_isometry_log
                and record.get("comparison_method") == "identical_serialized_gram"
            ):
                require(
                    iso.get("promotion_method")
                    == "identical_serialized_gram_sha256_and_matrix",
                    f"isometry-checks.jsonl:{index}: identical-Gram promotion method drift",
                )
                require(
                    iso.get("integral_unimodular_transform")
                    == identity_integer_matrix(),
                    f"isometry-checks.jsonl:{index}: identical-Gram promotion is not identity",
                )
        else:
            require(record.get("isometry") is None, f"isometry-checks.jsonl:{index}: negative check embeds an isometry")

    comparison_stats = result.get("isometry_comparisons")
    if comparison_stats is not None:
        require(isinstance(comparison_stats, dict), "result: invalid isometry-comparison statistics")
        require(comparison_stats.get("pairs_seen") == len(isometries), "result: isometry pair count mismatch")
        require(comparison_stats.get("qfisom_checked") == qfisom_log_count, "result: qfisom count mismatch")
        require(
            comparison_stats.get("identical_hash_matches") == identical_log_count,
            "result: identical-Gram count mismatch",
        )
        require(comparison_stats.get("skipped_due_to_budget") == 0, "importable result skipped an isometry pair")
        require(comparison_stats.get("budget_exhausted") is False, "importable result exhausted isometry budget")
        require(result.get("isometry_checks") == qfisom_log_count, "result: compatibility qfisom count mismatch")
    else:
        require(result.get("isometry_checks") == len(isometries), "isometry JSONL count mismatch")

    origin_by_fingerprint: dict[tuple[int, int], list[str]] = {}
    endpoint_by_fingerprint: dict[tuple[int, int], list[str]] = {}
    for state in states:
        fingerprint = (
            state["checks"]["root_count_exact_theta"],
            state["checks"]["norm4_count_exact_theta"],
        )
        target = origin_by_fingerprint if state["side"] == "origin" else endpoint_by_fingerprint
        target.setdefault(fingerprint, []).append(state["uid"])
    expected_cross_fingerprint_pairs = sum(
        len(origin_uids) * len(endpoint_by_fingerprint.get(fingerprint, []))
        for fingerprint, origin_uids in origin_by_fingerprint.items()
    )
    logged_pairs = [(record["origin_uid"], record["endpoint_uid"]) for record in isometries]
    require(len(logged_pairs) == len(set(logged_pairs)), "isometry log repeats a cross-side pair")
    if result["status"] == "bounded_search_completed":
        require(
            len(isometries) == expected_cross_fingerprint_pairs,
            "bounded result left a cross-side fingerprint pair unchecked",
        )
        origin_hashes = {state["gram_sha256"] for state in states if state["side"] == "origin"}
        endpoint_hashes = {
            state["gram_sha256"] for state in states if state["side"] in {"transparent", "rootless"}
        }
        require(not (origin_hashes & endpoint_hashes), "bounded result contains an unpromoted identical Gram")

    config = result.get("configuration")
    require(isinstance(config, dict), "result: configuration is missing")
    require(config.get("rounds_completed") == len(rounds), "result: completed-round count mismatch")
    require(is_int(config.get("rounds_requested")) and config["rounds_requested"] >= len(rounds), "result: invalid requested rounds")
    verify_rounds(rounds, result, len(initial))
    counters = result["counters"]
    require(counters["attempted_lines"] == len(attempts), "attempt JSONL count mismatch")
    require(counters["successful_moves"] == len(moves), "move JSONL count mismatch")
    require(
        counters["failed_moves"]
        == sum(
            attempt["status"] in {"construction_error", "expected_construction_boundary"}
            for attempt in attempts
        ),
        "failed-move count mismatch",
    )
    if "expected_p2_construction_boundary_failures" in counters:
        require(
            counters["expected_p2_construction_boundary_failures"]
            == expected_p2_boundary_attempts,
            "expected p=2 boundary-failure count mismatch",
        )
        require(
            counters["failed_moves"] == expected_p2_boundary_attempts,
            "new fail-closed schema contains an unclassified construction failure",
        )
    if modern:
        require(
            "expected_p2_construction_boundary_failures" in counters,
            "modern result lacks the classified p=2 boundary counter",
        )
    require(counters["duplicate_child_grams"] == sum(record["status"] == "success_duplicate" for record in moves), "duplicate-move count mismatch")
    require(counters["unique_states_discovered"] == len(states), "state JSONL count mismatch")
    checkpoint = load_json(root / "checkpoint.json")
    require(isinstance(checkpoint, dict), "checkpoint.json is not an object")
    require(checkpoint.get("rounds_completed") == len(rounds), "checkpoint round count mismatch")
    require(checkpoint.get("counters") == counters, "checkpoint counters disagree with result")
    require(checkpoint.get("isometry_checks") == qfisom_log_count, "checkpoint qfisom count mismatch")
    if comparison_stats is not None:
        require(
            checkpoint.get("isometry_comparisons") == comparison_stats,
            "checkpoint isometry-comparison statistics mismatch",
        )
    verify_modern_configuration_and_termination(result, checkpoint)
    sampling_audit = (
        verify_sampling_windows(result, jsonl[SAMPLING_LOG_NAME], attempts, state_by_uid, rounds)
        if hash_sampling else None
    )
    if not hash_sampling:
        require(not any("sampling_window_sha256" in attempt for attempt in attempts), "sampling attempt metadata without hash sampling")

    claims = verify_claim_boundary(result)
    bridge_summary = verify_bridge(root, result, state_by_uid, logged_moves_by_hash)
    positive_isometries = [record for record in isometries if record["isometric"]]
    if result["status"] == "bridge_found":
        require(rounds and rounds[-1]["bridge_found"] is True, "found bridge is absent from final round")
        require(len(positive_isometries) == 1, "found bridge does not have exactly one positive isometry check")
        require(
            positive_isometries[0]["isometry"]
            == result["bridge"]["meeting_isometry"],
            "positive isometry log and promoted bridge disagree",
        )
        logged_move_hashes = set(logged_moves_by_hash)
        bridge_move_hashes = {
            move["move_sha256"]
            for move in [
                *result["bridge"]["origin_forward_moves"],
                *result["bridge"]["endpoint_forward_moves_to_invert"],
            ]
        }
        require(bridge_move_hashes <= logged_move_hashes, "exact bridge cites an unlogged neighbor move")
    else:
        require(not any(record["bridge_found"] for record in rounds), "bounded result has a bridge round")
        require(not positive_isometries, "bounded result has a positive isometry check")
    source_checkout_root, source_checkout = resolve_source_checkout(
        source_root, provenance
    )
    source_hashes = verify_source_files(result, source_checkout_root)
    embedded_runtime_provenance = verify_runtime_provenance(
        result, provenance, modern
    )

    importer_path = Path(__file__).resolve()
    importer_root = importer_path.parents[1]
    require(
        importer_path == importer_root / IMPORTER_SCRIPT_PATH,
        "verification implementation path drift",
    )
    legacy_scope_notes: list[str] = []
    if legacy_exact_p2_boundary_attempts:
        legacy_scope_notes.append(
            f"The legacy log contains {legacy_exact_p2_boundary_attempts} exact p=2 "
            "Sage boundary failure(s); those attempted lines are not excluded by the "
            "bounded negative window."
        )
    if legacy_unclassified_construction_errors:
        legacy_scope_notes.append(
            f"The legacy log contains {legacy_unclassified_construction_errors} "
            "unclassified construction error(s); those attempted lines are recorded "
            "but are not excluded by this audit."
        )
    legacy_error_note = (
        " " + " ".join(legacy_scope_notes) if legacy_scope_notes else ""
    )

    certificate: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "certificate_kind": CERTIFICATE_KIND,
        "verification_status": "verified",
        "provenance": copy.deepcopy(provenance),
        "verification_implementation": {
            "path": IMPORTER_SCRIPT_PATH,
            "sha256": file_sha256(importer_path),
        },
        "artifact": {
            "artifact_directory_name": root.name,
            "sha256sums_sha256": file_sha256(root / MANIFEST_NAME),
            "manifest_file_count": len(manifest),
            "manifest_total_bytes": sum(entry["size_bytes"] for entry in manifest.values()),
        },
        "search": {
            "artifact_schema_generation": "explicit-selection-v2" if result["schema_version"] == 2 else "modern" if modern else "legacy",
            "sampling_audit": sampling_audit,
            "status": result["status"],
            "truth_status": result["truth_status"],
            "sage_version": result["sage_version"],
            "configuration": copy.deepcopy(config),
            "counters": copy.deepcopy(counters),
            "isometry_checks": qfisom_log_count,
            "isometry_comparisons": copy.deepcopy(comparison_stats),
            "structured_expected_p2_boundary_failures": expected_p2_boundary_attempts,
            "legacy_exact_p2_boundary_failures": legacy_exact_p2_boundary_attempts,
            "total_classified_p2_boundary_failures": (
                expected_p2_boundary_attempts + legacy_exact_p2_boundary_attempts
            ),
            "legacy_unclassified_construction_errors": legacy_unclassified_construction_errors,
            "cross_fingerprint_pairs_expected": expected_cross_fingerprint_pairs,
            "cross_fingerprint_pairs_all_checked_for_bounded_result": (
                result["status"] != "bounded_search_completed"
                or len(isometries) == expected_cross_fingerprint_pairs
            ),
            "result_payload_sha256": result["result_payload_sha256"],
            "result_file_sha256": manifest["result.json"]["sha256"],
            "source_sha256": source_hashes,
            "source_checkout": source_checkout,
            "embedded_runtime_provenance": embedded_runtime_provenance,
            "jsonl": {
                name: {"record_count": len(jsonl[name]), "sha256": manifest[name]["sha256"]}
                for name in jsonl_names
            },
            "initial_state_gram_sha256": copy.deepcopy(INITIAL_HASHES),
            "bridge": bridge_summary,
            "independent_exact_arithmetic_scope": {
                "all_logged_search_moves_recomputed": False,
                "found_bridge_moves_recomputed": bridge_summary is not None,
                "found_bridge_end_to_end_composite_recomputed": bridge_summary is not None,
                "note": (
                    "The Sage run asserts every accepted local move. The standard-library importer "
                    "recomputes matrices only for the path promoted into an exact bridge."
                ),
            },
        },
        "claim_boundary": certificate_claim_boundary(
            claims,
            successful_moves=counters["successful_moves"],
            classified_p2_boundary_failures=(
                expected_p2_boundary_attempts + legacy_exact_p2_boundary_attempts
            ),
            unclassified_construction_errors=legacy_unclassified_construction_errors,
        ),
        "audit_conclusion": (
            "Artifact closure, hashes, internal links, finite-search counters, cross-fingerprint "
            "comparison coverage, and the conservative claim boundary were verified. Accepted "
            "local moves are Sage-asserted but are not all independently replayed by this importer. "
            "For a found bridge, every cited matrix identity and the end-to-end composite are also "
            "recomputed independently with Fraction. No explicit K3/moduli transport, rational P3, "
            "or elliptic rank-record curve is certified by this audit."
            + legacy_error_note
        ),
    }
    certificate["record_sha256"] = payload_sha256(certificate)
    return certificate


def parse_extra_provenance(values: Iterable[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for raw in values:
        key, separator, value = raw.partition("=")
        require(separator == "=" and key and value, f"invalid --provenance value: {raw!r}")
        require(re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*", key) is not None, f"invalid provenance key: {key!r}")
        require(key not in result, f"duplicate provenance key: {key}")
        result[key] = value
    return result


def cli_provenance(args: argparse.Namespace) -> dict[str, Any]:
    require(REPOSITORY.fullmatch(args.repository) is not None, "invalid repository owner/name")
    require(args.workflow == WORKFLOW_PATH, f"unexpected workflow path: {args.workflow}")
    require(HEX40.fullmatch(args.commit_sha) is not None, "commit SHA must be 40 lowercase hex characters")
    require(args.run_id.isdecimal() and int(args.run_id) > 0, "run ID must be a positive decimal string")
    require(args.run_attempt > 0, "run attempt must be positive")
    require(
        re.fullmatch(r"refs/(?:heads|tags)/[^\s]+", args.git_ref) is not None,
        "git ref must be the full GitHub refs/heads/... or refs/tags/... value",
    )
    parsed_url = urlparse(args.run_url)
    expected_path = f"/{args.repository}/actions/runs/{args.run_id}"
    require(
        parsed_url.scheme == "https"
        and parsed_url.netloc == "github.com"
        and parsed_url.path.rstrip("/") == expected_path
        and not parsed_url.query
        and not parsed_url.fragment,
        "run URL does not identify the declared GitHub repository/run",
    )
    require(args.artifact_name and "/" not in args.artifact_name, "artifact name is missing or unsafe")
    require(
        args.artifact_id is not None
        and args.artifact_id.isdecimal()
        and int(args.artifact_id) > 0,
        "artifact ID must be a positive decimal",
    )
    require(
        args.artifact_digest is not None
        and re.fullmatch(r"sha256:[0-9a-f]{64}", args.artifact_digest),
        "artifact digest must be sha256:<lowercase hex>",
    )
    require(args.archive_sha256 is None or HEX64.fullmatch(args.archive_sha256), "archive SHA-256 must be lowercase hex")
    if args.archive_sha256 is not None:
        require(
            args.archive_sha256 == args.artifact_digest.removeprefix("sha256:"),
            "downloaded archive hash disagrees with GitHub artifact digest",
        )
    return {
        "repository": args.repository,
        "workflow": args.workflow,
        "event_name": args.event_name,
        "git_ref": args.git_ref,
        "commit_sha": args.commit_sha,
        "run_id": args.run_id,
        "run_attempt": args.run_attempt,
        "run_url": args.run_url,
        "artifact_name": args.artifact_name,
        "artifact_id": args.artifact_id,
        "artifact_digest": args.artifact_digest,
        "downloaded_archive_sha256": args.archive_sha256,
        "extra": parse_extra_provenance(args.provenance),
    }


def write_certificate(path: Path, certificate: dict[str, Any]) -> None:
    require(not path.exists(), f"refusing to overwrite certificate: {path}")
    require(path.parent.is_dir() and not path.parent.is_symlink(), "certificate parent is missing or unsafe")
    path.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", type=Path, help="downloaded and extracted artifact directory")
    parser.add_argument("--certificate", type=Path, required=True, help="new compact certificate path")
    parser.add_argument("--repository", required=True)
    parser.add_argument("--workflow", default=WORKFLOW_PATH)
    parser.add_argument("--event-name", choices=("push", "workflow_dispatch"), required=True)
    parser.add_argument("--git-ref", required=True)
    parser.add_argument("--commit-sha", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--run-attempt", type=int, required=True)
    parser.add_argument("--run-url", required=True)
    parser.add_argument("--artifact-name", default="e8-a2-target-exact-neighbor-bridge")
    parser.add_argument("--artifact-id", required=True)
    parser.add_argument("--artifact-digest", required=True)
    parser.add_argument("--archive-sha256")
    parser.add_argument(
        "--source-root",
        type=Path,
        help=(
            "optional clean Git checkout whose HEAD equals --commit-sha; use it "
            "to verify historical artifacts against their producing source"
        ),
    )
    parser.add_argument(
        "--provenance",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="repeatable additional provenance string (duplicate keys rejected)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    try:
        provenance = cli_provenance(args)
        certificate = verify_artifact(
            args.artifact, provenance, source_root=args.source_root
        )
        write_certificate(args.certificate, certificate)
    except (VerificationError, OSError) as exc:
        print(f"FAILED: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    print(
        json.dumps(
            {
                "status": "VERIFIED",
                "certificate": str(args.certificate),
                "record_sha256": certificate["record_sha256"],
                "search_status": certificate["search"]["status"],
                "rank31_curve_found": False,
                "rank32_curve_found": False,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
