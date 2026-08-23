#!/usr/bin/env sage-python
"""Search exact p-neighbor paths out of the E8+A2^3 target lattice.

This is a deterministic, bounded Sage 10.9 computation.  It searches from
the positive essential lattice forced by the three-section E8+A2^3 chart and,
optionally, backwards from two already certified endpoint lattices:

* the transparent A11 + K6(det=79) seed; and
* the frozen rootless determinant-948 lattice.

Every successfully constructed neighbor move is written to ``moves.jsonl``.
In particular, the record contains the exact rational basis transformation T
with

    T.transpose() * G_parent * T == G_child.

Thus a successful run records more than a chain of Gram hashes: it preserves
the rational embeddings needed for a later Neron--Severi transport.  A
negative run proves only that no bridge occurred in the explicitly recorded
finite search window.

Run with Sage, not CPython::

    sage -python research/search_e8_a2_target_neighbor_bridge.py \
      --output /tmp/e8-a2-neighbor-bridge
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import traceback
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Iterator


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from research.certify_e8_a2_shimura_bridge import (  # noqa: E402
    PERIOD_LATTICE,
    target_essential_lattice,
    transparent_seed_lattice,
)


try:  # Keep the module importable by the standard-library unit tests.
    from sage.all import (  # type: ignore[import-not-found]
        GF,
        Genus,
        IntegralLattice,
        Matrix,
        QQ,
        QuadraticForm,
        ZZ,
        pari,
        version,
    )

    SAGE_AVAILABLE = True
except ModuleNotFoundError:  # pragma: no cover - exercised outside Sage.
    SAGE_AVAILABLE = False


SCHEMA_VERSION = 1
EXPECTED_SAGE_PREFIX = "SageMath version 10.9"
TARGET_FORMULA_HASH = "dc4108ea4612195c2ba8350b0b1b35067f92ec2d49801cd2cb597cbb61129690"
TRANSPARENT_SEED_FORMULA_HASH = "10270a49860bff642fe34bc3eb0fc213623c7dff8aca08b6df51e2832c600bff"
FROZEN_ROOTLESS_HASH = "620a5e06473684d3e8015c0172f63c09c901e742ec02e77ba0aa35a923aa0295"


FROZEN_ROOTLESS_GRAM: list[list[int]] = [
    [4, -2, 2, -2, -1, 2, -2, -2, 1, 1, 1, -2, 1, -2, -1, 0, 1],
    [-2, 4, 0, 0, -1, 0, 0, 1, -2, -2, -2, 0, 1, 0, -1, 1, 1],
    [2, 0, 4, 0, 0, 2, -1, 0, -1, -1, -1, -1, 0, -1, -2, -1, 0],
    [-2, 0, 0, 4, 1, -1, 1, 2, 0, 0, 0, 1, -1, 1, 0, -2, -1],
    [-1, -1, 0, 1, 4, 0, 2, 0, 1, 1, 1, 2, -1, 0, 2, -2, -1],
    [2, 0, 2, -1, 0, 4, 0, 0, -1, 1, -1, -2, 2, -1, -1, 0, 0],
    [-2, 0, -1, 1, 2, 0, 4, 2, 0, 1, 0, 1, -1, 1, 2, -1, -2],
    [-2, 1, 0, 2, 0, 0, 2, 4, -2, 0, -1, 0, -1, 2, 0, -1, -2],
    [1, -2, -1, 0, 1, -1, 0, -2, 4, 2, 2, 0, 0, -1, 2, -1, 1],
    [1, -2, -1, 0, 1, 1, 1, 0, 2, 4, 2, -1, 0, 0, 2, -1, 0],
    [1, -2, -1, 0, 1, -1, 0, -1, 2, 2, 4, 0, -1, -1, 1, -1, 0],
    [-2, 0, -1, 1, 2, -2, 1, 0, 0, -1, 0, 4, -2, 1, 1, 0, 0],
    [1, 1, 0, -1, -1, 2, -1, -1, 0, 0, -1, -2, 4, -2, -1, 1, 1],
    [-2, 0, -1, 1, 0, -1, 1, 2, -1, 0, -1, 1, -2, 4, 1, 0, -1],
    [-1, -1, -2, 0, 2, -1, 2, 0, 2, 2, 1, 1, -1, 1, 4, -1, 0],
    [0, 1, -1, -2, -2, 0, -1, -1, -1, -1, -1, 0, 1, 0, -1, 4, 1],
    [1, 1, 0, -1, -1, 0, -2, -2, 1, 0, 0, 0, 1, -1, 0, 1, 4],
]


def canonical_json_bytes(payload: Any) -> bytes:
    return json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()


def payload_sha256(payload: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(payload)).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_primes(raw: str) -> tuple[int, ...]:
    values = tuple(int(part.strip()) for part in raw.split(",") if part.strip())
    if not values:
        raise argparse.ArgumentTypeError("at least one prime is required")
    for value in values:
        if value < 2 or any(value % divisor == 0 for divisor in range(2, int(value**0.5) + 1)):
            raise argparse.ArgumentTypeError(f"not a prime: {value}")
    if len(set(values)) != len(values):
        raise argparse.ArgumentTypeError("neighbor primes must be distinct")
    return values


def serial_matrix(matrix: Any) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.ncols())]
        for row in range(matrix.nrows())
    ]


def rational_string(value: Any) -> str:
    numerator = int(value.numerator())
    denominator = int(value.denominator())
    return str(numerator) if denominator == 1 else f"{numerator}/{denominator}"


def serial_rational_matrix(matrix: Any) -> list[list[str]]:
    return [
        [rational_string(matrix[row, column]) for column in range(matrix.ncols())]
        for row in range(matrix.nrows())
    ]


def matrix_hash(matrix: Any) -> str:
    return payload_sha256(serial_matrix(matrix))


def append_jsonl(path: Path, payload: dict[str, Any]) -> None:
    with path.open("a", encoding="utf-8") as handle:
        handle.write(canonical_json_bytes(payload).decode())
        handle.write("\n")


def exact_theta_counts(gram: Any) -> tuple[int, int]:
    form = QuadraticForm(ZZ, gram)
    theta = form.theta_series(3)
    return int(theta[1]), int(theta[2])


def independent_short_counts(gram: Any) -> tuple[int, int]:
    short = IntegralLattice(gram).short_vectors(5)
    return len(short[2]), len(short[4])


def exact_state_checks(gram: Any, target_genus: Any) -> dict[str, Any]:
    if gram.nrows() != 17 or gram.ncols() != 17:
        raise AssertionError(("wrong rank", gram.dimensions()))
    if gram != gram.transpose():
        raise AssertionError("Gram matrix is not symmetric")
    determinant = int(gram.det())
    if determinant != 948:
        raise AssertionError(("wrong determinant", determinant))
    if any(int(gram[index, index]) % 2 for index in range(17)):
        raise AssertionError("Gram matrix is not even")
    if not gram.is_positive_definite():
        raise AssertionError("Gram matrix is not positive definite")
    genus = Genus(gram)
    if genus != target_genus:
        raise AssertionError(("left target genus", genus, target_genus))
    root_count, norm4_count = exact_theta_counts(gram)
    return {
        "rank": 17,
        "determinant": determinant,
        "even": True,
        "positive_definite": True,
        "target_genus_verified": True,
        "root_count_exact_theta": root_count,
        "norm4_count_exact_theta": norm4_count,
    }


def reduce_gram_with_transform(gram: Any) -> tuple[Any, Any]:
    transform = Matrix(ZZ, gram.LLL_gram())
    if abs(int(transform.det())) != 1:
        raise AssertionError(("LLL transform is not unimodular", transform.det()))
    reduced = Matrix(ZZ, transform.transpose() * gram * transform)
    return reduced, transform


def matrix_denominator(matrix: Any) -> int:
    denominator = ZZ(1)
    for entry in matrix.list():
        denominator = denominator.lcm(ZZ(entry.denominator()))
    return int(denominator)


def verify_p_neighbor_transition(parent_gram: Any, child_gram: Any, transition: Any, prime: int) -> dict[str, Any]:
    p = ZZ(prime)
    transition = Matrix(QQ, transition)
    if transition.transpose() * parent_gram * transition != child_gram:
        raise AssertionError("parent-to-child rational Gram identity failed")
    determinant = transition.det()
    if abs(determinant) != 1:
        raise AssertionError(("neighbor basis covolumes differ", determinant))
    inverse = transition.inverse()
    forward_denominator = matrix_denominator(transition)
    inverse_denominator = matrix_denominator(inverse)
    if forward_denominator != prime or inverse_denominator != prime:
        raise AssertionError(
            ("unexpected p-neighbor denominators", prime, forward_denominator, inverse_denominator)
        )
    scaled_forward = Matrix(ZZ, p * transition)
    scaled_inverse = Matrix(ZZ, p * inverse)
    forward_rank = Matrix(GF(p), scaled_forward).rank()
    inverse_rank = Matrix(GF(p), scaled_inverse).rank()
    if forward_rank != 1 or inverse_rank != 1:
        raise AssertionError(("not an index-p/index-p neighbor", forward_rank, inverse_rank))
    return {
        "gram_identity_verified": True,
        "transition_determinant": rational_string(determinant),
        "transition_denominator": forward_denominator,
        "inverse_transition_denominator": inverse_denominator,
        "rank_mod_p_of_p_times_transition": int(forward_rank),
        "rank_mod_p_of_p_times_inverse": int(inverse_rank),
        "intersection_index_in_parent": prime,
        "intersection_index_in_child": prime,
    }


def projective_isotropic_lines(form: Any, prime: int, offset: int, limit: int) -> Iterator[tuple[int, Any]]:
    """Yield an exact deterministic window in Sage's projective ordering."""

    previous = None
    ordinal = 0
    yielded = 0
    while yielded < limit:
        current = form.find_primitive_p_divisible_vector__next(ZZ(prime), previous)
        if current is None:
            return
        previous = current
        ordinal += 1
        if ordinal <= offset:
            continue
        q_value = ZZ(form(current))
        if q_value % prime:
            raise AssertionError(("Sage returned a non-isotropic line", prime, current, q_value))
        yield ordinal, current
        yielded += 1


def build_neighbor(parent: dict[str, Any], prime: int, line_ordinal: int, vector: Any, target_genus: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    parent_gram = parent["gram_object"]
    parent_form = QuadraticForm(ZZ, parent_gram)
    raw_transform = Matrix(QQ, parent_form.find_p_neighbor_from_vec(ZZ(prime), vector, return_matrix=True))
    raw_neighbor_form = parent_form(raw_transform)
    raw_gram_qq = Matrix(QQ, raw_neighbor_form.Hessian_matrix())
    if any(entry.denominator() != 1 for entry in raw_gram_qq.list()):
        raise AssertionError("Sage neighbor Hessian is not integral")
    raw_gram = Matrix(ZZ, raw_gram_qq)
    child_gram, lll_transform = reduce_gram_with_transform(raw_gram)
    total_transition = Matrix(QQ, raw_transform * lll_transform)
    relation = verify_p_neighbor_transition(parent_gram, child_gram, total_transition, prime)
    checks = exact_state_checks(child_gram, target_genus)
    child_hash = matrix_hash(child_gram)
    q_value = ZZ(parent_form(vector))

    move_core = {
        "parent_uid": parent["uid"],
        "parent_gram_sha256": parent["gram_sha256"],
        "prime": prime,
        "projective_line_ordinal_one_based": line_ordinal,
        "projective_isotropic_vector": [int(entry) for entry in vector],
        "q_of_vector": int(q_value),
        "q_of_vector_mod_p": int(q_value % prime),
        "raw_neighbor_gram": serial_matrix(raw_gram),
        "raw_neighbor_transform": serial_rational_matrix(raw_transform),
        "lll_unimodular_transform": serial_matrix(lll_transform),
        "parent_to_child_rational_transform": serial_rational_matrix(total_transition),
        "child_gram_sha256": child_hash,
        "child_gram": serial_matrix(child_gram),
        "p_neighbor_relation": relation,
        "child_exact_checks": checks,
    }
    move = {**move_core, "move_sha256": payload_sha256(move_core)}
    child = {
        "side": parent["side"],
        "gram_object": child_gram,
        "gram_sha256": child_hash,
        "parent_uid": parent["uid"],
        "move": move,
        "checks": checks,
        "source": None,
        "round_discovered": parent["round_discovered"] + 1,
    }
    child["uid"] = f'{child["side"]}:{child_hash}'
    return child, move


def state_public_record(state: dict[str, Any]) -> dict[str, Any]:
    return {
        "uid": state["uid"],
        "side": state["side"],
        "source": state["source"],
        "round_discovered": state["round_discovered"],
        "gram_sha256": state["gram_sha256"],
        "parent_uid": state["parent_uid"],
        "incoming_move_sha256": state["move"]["move_sha256"] if state["move"] else None,
        "checks": state["checks"],
        "gram": serial_matrix(state["gram_object"]),
    }


def make_initial_state(side: str, source: str, gram_data: list[list[int]], target_genus: Any) -> dict[str, Any]:
    gram = Matrix(ZZ, gram_data)
    checks = exact_state_checks(gram, target_genus)
    state = {
        "side": side,
        "source": source,
        "gram_object": gram,
        "gram_sha256": matrix_hash(gram),
        "parent_uid": None,
        "move": None,
        "checks": checks,
        "round_discovered": 0,
    }
    state["uid"] = f'{side}:{state["gram_sha256"]}'
    return state


def state_fingerprint(state: dict[str, Any]) -> tuple[int, int]:
    checks = state["checks"]
    return checks["root_count_exact_theta"], checks["norm4_count_exact_theta"]


def distance_key(state: dict[str, Any], profiles: Iterable[tuple[int, int]]) -> tuple[Any, ...]:
    roots, norm4 = state_fingerprint(state)
    distances = sorted((abs(roots - goal_roots), abs(norm4 - goal_norm4)) for goal_roots, goal_norm4 in profiles)
    return (distances[0], roots, norm4, state["gram_sha256"])


def select_next_beam(candidates: list[dict[str, Any]], size: int, profiles: list[tuple[int, int]]) -> list[dict[str, Any]]:
    if not candidates or size <= 0:
        return []
    selected: list[dict[str, Any]] = []
    used: set[str] = set()

    # Preserve one exact queue per opposing endpoint profile.
    for profile in profiles:
        ordered = sorted(
            candidates,
            key=lambda state: (
                abs(state_fingerprint(state)[0] - profile[0]),
                abs(state_fingerprint(state)[1] - profile[1]),
                state["gram_sha256"],
            ),
        )
        for state in ordered:
            if state["uid"] not in used:
                selected.append(state)
                used.add(state["uid"])
                break
            
    # Then exploit the combined score while retaining deterministic diversity.
    ordered = sorted(candidates, key=lambda state: distance_key(state, profiles))
    if ordered:
        for numerator, denominator in ((1, 4), (1, 2), (3, 4)):
            state = ordered[min(len(ordered) - 1, len(ordered) * numerator // denominator)]
            if state["uid"] not in used and len(selected) < size:
                selected.append(state)
                used.add(state["uid"])
    for state in ordered:
        if len(selected) >= size:
            break
        if state["uid"] not in used:
            selected.append(state)
            used.add(state["uid"])
    return selected[:size]


def exact_isometry(left: dict[str, Any], right: dict[str, Any]) -> dict[str, Any] | None:
    if state_fingerprint(left) != state_fingerprint(right):
        return None
    left_form = QuadraticForm(ZZ, left["gram_object"])
    right_form = QuadraticForm(ZZ, right["gram_object"])
    transform = left_form.is_globally_equivalent_to(right_form, return_matrix=True)
    if transform is False:
        return None
    transform = Matrix(ZZ, transform)
    if abs(int(transform.det())) != 1:
        raise AssertionError("global equivalence transform is not unimodular")
    if transform.transpose() * left["gram_object"] * transform != right["gram_object"]:
        raise AssertionError("global equivalence transform fails Gram identity")
    return {
        "from_uid": left["uid"],
        "to_uid": right["uid"],
        "integral_unimodular_transform": serial_matrix(transform),
        "determinant": int(transform.det()),
        "gram_identity_verified": True,
        "fingerprint": list(state_fingerprint(left)),
    }


def path_to(state: dict[str, Any], states: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    path = []
    cursor = state
    while cursor is not None:
        path.append(cursor)
        parent_uid = cursor["parent_uid"]
        cursor = states.get(parent_uid) if parent_uid is not None else None
    path.reverse()
    return path


def bridge_record(origin_state: dict[str, Any], endpoint_state: dict[str, Any], isometry: dict[str, Any], states: dict[str, dict[str, Any]], endpoint_name: str) -> dict[str, Any]:
    origin_path = path_to(origin_state, states)
    endpoint_path = path_to(endpoint_state, states)
    origin_moves = [state["move"] for state in origin_path[1:]]
    endpoint_moves = [state["move"] for state in endpoint_path[1:]]
    payload = {
        "schema_version": SCHEMA_VERSION,
        "status": "exact_neighbor_bridge_found",
        "endpoint": endpoint_name,
        "orientation": (
            "Follow origin_chain forward, apply meeting_isometry, then follow "
            "endpoint_chain in reverse using the inverse of each recorded rational transform."
        ),
        "origin_chain": [state_public_record(state) for state in origin_path],
        "origin_forward_moves": origin_moves,
        "meeting_isometry": isometry,
        "endpoint_chain_endpoint_to_meeting": [state_public_record(state) for state in endpoint_path],
        "endpoint_forward_moves_to_invert": endpoint_moves,
        "neighbor_step_count": len(origin_moves) + len(endpoint_moves),
        "all_parent_child_gram_identities_verified": True,
        "all_neighbor_intersection_indices_verified": True,
        "endpoint_identification_verified_by_initial_source": True,
        "claim_boundary": (
            "This is an exact integral-lattice p-neighbor bridge.  It is not yet an "
            "explicit K3 fibration switch, moduli map, rational P3 section, or rank-31 curve."
        ),
    }
    return {**payload, "bridge_sha256": payload_sha256(payload)}


def find_meeting(
    state: dict[str, Any],
    indexes: dict[str, dict[tuple[int, int], list[dict[str, Any]]]],
    endpoint_sides: tuple[str, ...],
    isometry_log: Path,
    counter: dict[str, int],
    limit: int,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], str] | None:
    if state["side"] == "origin":
        comparisons = endpoint_sides
        left_is_state = True
    elif state["side"] in endpoint_sides:
        comparisons = ("origin",)
        left_is_state = False
    else:
        return None

    fingerprint = state_fingerprint(state)
    for other_side in comparisons:
        for other in indexes[other_side].get(fingerprint, []):
            if counter["count"] >= limit:
                return None
            left, right = (state, other) if left_is_state else (other, state)
            counter["count"] += 1
            isometry = exact_isometry(left, right)
            check_record = {
                "check_index": counter["count"],
                "origin_uid": left["uid"],
                "endpoint_uid": right["uid"],
                "endpoint_side": right["side"],
                "fingerprint": list(fingerprint),
                "isometric": isometry is not None,
                "isometry": isometry,
            }
            append_jsonl(isometry_log, check_record)
            if isometry is not None:
                return left, right, isometry, right["side"]
    return None


def validate_configuration(args: argparse.Namespace) -> None:
    if args.rounds < 0:
        raise ValueError("rounds must be nonnegative")
    if args.beam_size < 1:
        raise ValueError("beam size must be positive")
    if args.lines_per_prime < 1:
        raise ValueError("lines per prime must be positive")
    if args.line_offset < 0:
        raise ValueError("line offset must be nonnegative")
    if args.max_successful_moves < 1:
        raise ValueError("max successful moves must be positive")
    if args.max_isometry_checks < 1:
        raise ValueError("max isometry checks must be positive")


def run_search(args: argparse.Namespace) -> dict[str, Any]:
    validate_configuration(args)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    for filename in ("states.jsonl", "moves.jsonl", "attempts.jsonl", "rounds.jsonl", "isometry-checks.jsonl"):
        (output / filename).write_text("", encoding="utf-8")

    sage_version = str(version())
    if not sage_version.startswith(EXPECTED_SAGE_PREFIX):
        raise RuntimeError(f"requires pinned Sage 10.9, found {sage_version}")
    pari.allocatemem(2 * 1024**3)

    target_data = target_essential_lattice()
    transparent_data = transparent_seed_lattice()
    if payload_sha256(target_data) != TARGET_FORMULA_HASH:
        raise AssertionError("target formula drift")
    if payload_sha256(transparent_data) != TRANSPARENT_SEED_FORMULA_HASH:
        raise AssertionError("transparent seed formula drift")
    if payload_sha256(FROZEN_ROOTLESS_GRAM) != FROZEN_ROOTLESS_HASH:
        raise AssertionError("frozen rootless matrix drift")

    period = Matrix(ZZ, PERIOD_LATTICE)
    target_genus = IntegralLattice(period).discriminant_group().genus((17, 0))
    origin = make_initial_state("origin", "E8+A2^3 three-section target", target_data, target_genus)
    transparent = make_initial_state("transparent", "A11+K6(det=79) certified seed", transparent_data, target_genus)
    rootless = make_initial_state("rootless", "frozen rootless target 620a5e...", FROZEN_ROOTLESS_GRAM, target_genus)
    if origin["gram_sha256"] != TARGET_FORMULA_HASH:
        raise AssertionError("unexpected origin hash")
    if transparent["gram_sha256"] != TRANSPARENT_SEED_FORMULA_HASH:
        raise AssertionError("unexpected transparent hash")
    if rootless["gram_sha256"] != FROZEN_ROOTLESS_HASH:
        raise AssertionError("unexpected frozen rootless hash")

    initial_states = [origin, transparent, rootless]
    for state in initial_states:
        exact_short = independent_short_counts(state["gram_object"])
        if exact_short != state_fingerprint(state):
            raise AssertionError(("initial theta/short-vector mismatch", state["uid"], exact_short, state_fingerprint(state)))
        state["checks"]["independent_short_vector_counts_verified"] = list(exact_short)
        append_jsonl(output / "states.jsonl", state_public_record(state))

    all_states = {state["uid"]: state for state in initial_states}
    seen_by_side: dict[str, set[str]] = {state["side"]: {state["gram_sha256"]} for state in initial_states}
    indexes: dict[str, dict[tuple[int, int], list[dict[str, Any]]]] = {
        state["side"]: defaultdict(list) for state in initial_states
    }
    for state in initial_states:
        indexes[state["side"]][state_fingerprint(state)].append(state)

    sides = ("origin",) if args.mode == "forward" else ("origin", "transparent", "rootless")
    endpoint_sides = ("transparent", "rootless")
    beams: dict[str, list[dict[str, Any]]] = {side: [all_states[next(uid for uid in all_states if uid.startswith(side + ":"))]] for side in sides}
    expanded: dict[str, set[str]] = {side: set() for side in sides}
    profile_by_side = {state["side"]: state_fingerprint(state) for state in initial_states}

    counters = {
        "attempted_lines": 0,
        "successful_moves": 0,
        "failed_moves": 0,
        "duplicate_child_grams": 0,
        "unique_states_discovered": len(initial_states),
        "move_budget_exhausted": False,
    }
    isometry_counter = {"count": 0}
    bridge = None
    rounds_completed = 0

    initial_checkpoint = {
        "status": "initialized",
        "rounds_completed": 0,
        "counters": counters,
        "isometry_checks": 0,
        "current_beams": {
            side: [state_public_record(state) for state in beam]
            for side, beam in beams.items()
        },
    }
    (output / "checkpoint.json").write_text(
        json.dumps(initial_checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    # The initial endpoint fingerprints differ, so any initial match would be
    # an input inconsistency.  Search only newly constructed states.
    for round_index in range(1, args.rounds + 1):
        round_record: dict[str, Any] = {"round": round_index, "sides": {}}
        for side in sides:
            candidates: list[dict[str, Any]] = []
            side_stats = {
                "parents": 0,
                "attempted_lines": 0,
                "successful_moves": 0,
                "failed_moves": 0,
                "new_unique_states": 0,
                "duplicate_child_grams": 0,
            }
            for parent in beams[side]:
                if parent["uid"] in expanded[side]:
                    continue
                expanded[side].add(parent["uid"])
                side_stats["parents"] += 1
                parent_form = QuadraticForm(ZZ, parent["gram_object"])
                for prime in args.primes:
                    for ordinal, vector in projective_isotropic_lines(
                        parent_form, prime, args.line_offset, args.lines_per_prime
                    ):
                        if counters["successful_moves"] >= args.max_successful_moves:
                            counters["move_budget_exhausted"] = True
                            break
                        counters["attempted_lines"] += 1
                        side_stats["attempted_lines"] += 1
                        attempt_base = {
                            "attempt_index": counters["attempted_lines"],
                            "round": round_index,
                            "side": side,
                            "parent_uid": parent["uid"],
                            "prime": prime,
                            "projective_line_ordinal_one_based": ordinal,
                            "projective_isotropic_vector": [int(entry) for entry in vector],
                        }
                        try:
                            child, move = build_neighbor(parent, prime, ordinal, vector, target_genus)
                        except Exception as exc:
                            counters["failed_moves"] += 1
                            side_stats["failed_moves"] += 1
                            append_jsonl(
                                output / "attempts.jsonl",
                                {**attempt_base, "status": "construction_error", "error": repr(exc)},
                            )
                            continue
                        counters["successful_moves"] += 1
                        side_stats["successful_moves"] += 1
                        duplicate = child["gram_sha256"] in seen_by_side[side]
                        move_record = {
                            "attempt": attempt_base,
                            "status": "success_duplicate" if duplicate else "success_new_state",
                            "duplicate_state_on_same_side": duplicate,
                            "move": move,
                        }
                        append_jsonl(output / "moves.jsonl", move_record)
                        append_jsonl(
                            output / "attempts.jsonl",
                            {
                                **attempt_base,
                                "status": move_record["status"],
                                "move_sha256": move["move_sha256"],
                                "child_gram_sha256": child["gram_sha256"],
                            },
                        )
                        if duplicate:
                            counters["duplicate_child_grams"] += 1
                            side_stats["duplicate_child_grams"] += 1
                            continue

                        seen_by_side[side].add(child["gram_sha256"])
                        all_states[child["uid"]] = child
                        indexes[side][state_fingerprint(child)].append(child)
                        candidates.append(child)
                        counters["unique_states_discovered"] += 1
                        side_stats["new_unique_states"] += 1
                        append_jsonl(output / "states.jsonl", state_public_record(child))

                        meeting = find_meeting(
                            child,
                            indexes,
                            endpoint_sides,
                            output / "isometry-checks.jsonl",
                            isometry_counter,
                            args.max_isometry_checks,
                        )
                        if meeting is not None:
                            left, right, isometry, endpoint_name = meeting
                            # Independently enumerate the two meeting lattices
                            # before promoting the exact global equivalence.
                            for meeting_state in (left, right):
                                counts = independent_short_counts(meeting_state["gram_object"])
                                if counts != state_fingerprint(meeting_state):
                                    raise AssertionError(("meeting theta/short mismatch", meeting_state["uid"], counts))
                            bridge = bridge_record(left, right, isometry, all_states, endpoint_name)
                            (output / "exact-bridge.json").write_text(
                                json.dumps(bridge, indent=2, sort_keys=True) + "\n", encoding="utf-8"
                            )
                            break
                    if bridge is not None or counters["move_budget_exhausted"]:
                        break
                if bridge is not None or counters["move_budget_exhausted"]:
                    break

            opposing_profiles = (
                [profile_by_side["transparent"], profile_by_side["rootless"]]
                if side == "origin"
                else [profile_by_side["origin"]]
            )
            beams[side] = select_next_beam(candidates, args.beam_size, opposing_profiles)
            side_stats["next_beam"] = [
                {
                    "uid": state["uid"],
                    "gram_sha256": state["gram_sha256"],
                    "fingerprint": list(state_fingerprint(state)),
                }
                for state in beams[side]
            ]
            round_record["sides"][side] = side_stats
            if bridge is not None or counters["move_budget_exhausted"]:
                break

        rounds_completed = round_index
        round_record["cumulative"] = dict(counters)
        round_record["isometry_checks"] = isometry_counter["count"]
        round_record["bridge_found"] = bridge is not None
        append_jsonl(output / "rounds.jsonl", round_record)
        checkpoint = {
            "status": "bridge_found" if bridge is not None else "searching",
            "rounds_completed": rounds_completed,
            "counters": counters,
            "isometry_checks": isometry_counter["count"],
            "current_beams": {
                side: [state_public_record(state) for state in beam]
                for side, beam in beams.items()
            },
        }
        (output / "checkpoint.json").write_text(
            json.dumps(checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(json.dumps(round_record, sort_keys=True), flush=True)
        if bridge is not None or counters["move_budget_exhausted"] or all(not beams[side] for side in sides):
            break

    result: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "status": "bridge_found" if bridge is not None else "bounded_search_completed",
        "truth_status": (
            "An exact p-neighbor bridge with every rational parent-to-child basis transform and an "
            "integral meeting isometry was independently verified."
            if bridge is not None
            else "Finite negative result only: no exact bridge was found among the recorded deterministic "
            "projective-line windows, beams, rounds, successful moves, and isometry checks."
        ),
        "sage_version": sage_version,
        "configuration": {
            "mode": args.mode,
            "rounds_requested": args.rounds,
            "rounds_completed": rounds_completed,
            "beam_size": args.beam_size,
            "primes": list(args.primes),
            "projective_line_offset": args.line_offset,
            "projective_lines_per_state_prime": args.lines_per_prime,
            "max_successful_moves": args.max_successful_moves,
            "max_isometry_checks": args.max_isometry_checks,
            "enumeration_order": "Sage 10.9 find_primitive_p_divisible_vector__next",
        },
        "initial_states": [state_public_record(state) for state in initial_states],
        "target_genus": str(target_genus),
        "counters": counters,
        "isometry_checks": isometry_counter["count"],
        "bridge": bridge,
        "claim_boundary": {
            "exact_lattice_bridge_found": bridge is not None,
            "explicit_stable_NS_transport_completed": False,
            "explicit_K3_fibration_switch_completed": False,
            "explicit_moduli_map_completed": False,
            "rational_P3_found": False,
            "rank31_curve_found": False,
            "negative_result_scope": (
                "Only the serialized finite window was excluded; the p-neighbor graph and rank-31 problem remain open."
            ),
        },
        "source_files": {
            "search_script": str(Path(__file__).relative_to(REPO_ROOT)),
            "search_script_sha256": file_sha256(Path(__file__)),
            "target_formula_script": "research/certify_e8_a2_shimura_bridge.py",
            "target_formula_script_sha256": file_sha256(
                REPO_ROOT / "research/certify_e8_a2_shimura_bridge.py"
            ),
        },
    }
    evidence_files = [
        "states.jsonl",
        "moves.jsonl",
        "attempts.jsonl",
        "rounds.jsonl",
        "isometry-checks.jsonl",
        "checkpoint.json",
    ]
    if bridge is not None:
        evidence_files.append("exact-bridge.json")
    result["evidence_sha256"] = {
        filename: file_sha256(output / filename) for filename in evidence_files
    }
    result_core = dict(result)
    result["result_payload_sha256"] = payload_sha256(result_core)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode", choices=("forward", "bidirectional"), default="bidirectional")
    parser.add_argument("--rounds", type=int, default=6)
    parser.add_argument("--beam-size", type=int, default=5)
    parser.add_argument("--primes", type=parse_primes, default=parse_primes("2,5"))
    parser.add_argument("--line-offset", type=int, default=0)
    parser.add_argument("--lines-per-prime", dest="lines_per_prime", type=int, default=32)
    parser.add_argument("--max-successful-moves", type=int, default=10000)
    parser.add_argument("--max-isometry-checks", type=int, default=250)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if not SAGE_AVAILABLE:
        parser.error("this exact computation requires `sage -python` (SageMath 10.9)")
    args.output.mkdir(parents=True, exist_ok=True)
    try:
        result = run_search(args)
    except Exception as exc:
        result = {
            "schema_version": SCHEMA_VERSION,
            "status": "error",
            "error": repr(exc),
            "traceback": traceback.format_exc(),
            "claim_boundary": "No mathematical search conclusion may be drawn from an error run.",
        }
    result_path = args.output / "result.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)
    if result["status"] == "error":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
