#!/usr/bin/env python3
"""Bounded exact Eichler reduction of an abstract U-stabilized marking.

Every accepted transformation is an integral isometry. Search termination
does not prove global optimality, geometric effectivity, or an elliptic rank.
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from research import stable_hyperbolic_neighbor_transport as stable  # noqa: E402

audit = stable.audit
require = audit.require
INPUT_PATH = "certificates/short_hyperbolic_rootless_marking_20260908.json"
INPUT_RECORD_SHA256 = "ca67417287c26fb4c1d6665559a0486609aa7cd72286875f1161aabc22923a71"
MAX_ITERATIONS = 2048
MAX_NODES_PER_SHIFT = 100000


def gram_schmidt(gram):
    size = len(gram)
    require(0 < size <= 20 and all(len(row) == size for row in gram), "invalid Gram dimensions")
    require(all(type(x) is int for row in gram for x in row), "Gram is not integral")
    require(all(gram[i][j] == gram[j][i] for i in range(size) for j in range(size)), "Gram is not symmetric")
    require(all(gram[i][i] % 2 == 0 for i in range(size)), "Gram is not even")
    mu = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    diagonal = []
    for i in range(size):
        for j in range(i):
            mu[i][j] = (gram[i][j] - sum(mu[i][k] * mu[j][k] * diagonal[k] for k in range(j))) / diagonal[j]
        diagonal.append(Fraction(gram[i][i]) - sum(mu[i][k] ** 2 * diagonal[k] for k in range(i)))
        require(diagonal[-1] > 0, "Gram is not positive definite")
    return mu, diagonal


def nearest_integer(value):
    value = Fraction(value)
    return (2 * value.numerator + value.denominator) // (2 * value.denominator)


def bounded_shift(gram, coordinates, divisor, gs=None, max_nodes=MAX_NODES_PER_SHIFT):
    require(type(divisor) is int and divisor > 0, "divisor must be a positive integer")
    require(type(max_nodes) is int and max_nodes >= 0, "node budget must be nonnegative")
    require(len(coordinates) == len(gram), "coordinate dimension mismatch")
    mu, diagonal = gs if gs is not None else gram_schmidt(gram)
    size = len(gram)
    target = [Fraction(x, divisor) for x in coordinates]
    shifted = target[:]
    babai = [0] * size
    for j in reversed(range(size)):
        babai[j] = -nearest_integer(shifted[j] + sum(mu[i][j] * shifted[i] for i in range(j + 1, size)))
        shifted[j] += babai[j]
    best = stable.pairing(shifted, gram, shifted)
    initial_bound = best
    solution = babai[:]
    reached_target = best < 2
    exhausted = False
    nodes = 0
    path = [0] * size
    values = target[:]

    def enumerate_suffix(j, spent):
        nonlocal best, solution, reached_target, exhausted, nodes
        if reached_target or exhausted:
            return
        if j < 0:
            if spent < best:
                best, solution = spent, path[:]
                reached_target = best < 2
            return
        center = target[j] + sum(mu[i][j] * values[i] for i in range(j + 1, size))
        left = (-center).numerator // (-center).denominator
        right = left + 1
        while True:
            if abs(left + center) <= abs(right + center):
                integer = left
                left -= 1
            else:
                integer = right
                right += 1
            contribution = diagonal[j] * (integer + center) ** 2
            # Integer choices are ordered by distance, so all later choices
            # also miss this exact strict-improvement radius.
            if spent + contribution >= best:
                return
            if nodes >= max_nodes:
                exhausted = True
                return
            nodes += 1
            path[j], values[j] = integer, target[j] + integer
            enumerate_suffix(j - 1, spent + contribution)
            if reached_target or exhausted:
                return

    enumerate_suffix(size - 1, Fraction(0))
    witness = [target[i] + solution[i] for i in range(size)]
    require(stable.pairing(witness, gram, witness) == best, "shift distance mismatch")
    return solution, {
        "babai_squared_norm": str(Fraction(initial_bound)),
        "selected_squared_norm": str(Fraction(best)),
        "nodes_visited": nodes,
        "node_budget_exhausted": exhausted,
        "strict_norm_below_two_found": reached_target,
        "enumeration_finished_without_early_target_or_budget": not reached_target and not exhausted,
    }


def eichler_matrix(gram, shift):
    require(len(shift) == len(gram) and all(type(x) is int for x in shift), "nonintegral Eichler shift")
    size = len(gram)
    matrix = stable.identity(size + 2)
    dual = [sum(a * b for a, b in zip(row, shift)) for row in gram]
    norm = sum(a * b for a, b in zip(shift, dual))
    require(norm % 2 == 0, "odd Eichler shift norm")
    matrix[0][1] = norm // 2
    for i in range(size):
        matrix[0][i + 2] = dual[i]
        matrix[i + 2][1] = shift[i]
    return matrix


def swap_hyperbolic_coordinates(matrix, side):
    result = [row[:] for row in matrix]
    if side == "left":
        result[0], result[1] = result[1], result[0]
    else:
        require(side == "right", "unknown swap side")
        for row in result:
            row[0], row[1] = row[1], row[0]
    return result


def reduce_matrix(matrix, origin, endpoint, max_iterations=MAX_ITERATIONS, max_nodes=MAX_NODES_PER_SHIFT):
    require(type(max_iterations) is int and max_iterations >= 0, "iteration budget must be nonnegative")
    stable.verify_stable_matrix(matrix, origin, endpoint)
    matrix = [row[:] for row in matrix]
    gs = {"left": gram_schmidt(origin), "right": gram_schmidt(endpoint)}
    endpoint_inverse = stable.inverse(endpoint)
    trace = []
    termination = "iteration_budget_exhausted"
    for iteration in range(max_iterations):
        improved = False
        for side, gram in (("left", origin), ("right", endpoint)):
            before = matrix[1][0]
            require(before > 0, "nonpositive fiber intersection requires a separate boundary analysis")
            if side == "left":
                coordinates = [row[0] for row in matrix[2:]]
            else:
                row = matrix[1][2:]
                coordinates = [sum(a * b for a, b in zip(inverse_row, row)) for inverse_row in endpoint_inverse]
            shift, search = bounded_shift(gram, coordinates, before, gs[side], max_nodes)
            isometry = eichler_matrix(gram, shift)
            matrix = stable.multiply(isometry, matrix) if side == "left" else stable.multiply(matrix, isometry)
            other = matrix[0][0] if side == "left" else matrix[1][1]
            require(other > 0, "isotropic boundary with zero opposite coordinate is not promoted")
            improved = other < before
            if improved:
                matrix = swap_hyperbolic_coordinates(matrix, side)
            trace.append({
                "iteration": iteration,
                "side": side,
                "shift": shift,
                "fiber_intersection_before": before,
                "opposite_hyperbolic_coordinate_after_translation": other,
                "swap_performed": improved,
                "fiber_intersection_after": matrix[1][0],
                "search": search,
            })
            if improved:
                break
        if not improved:
            termination = "no_improving_swap_from_the_recorded_bounded_shifts"
            break
    determinant = stable.verify_stable_matrix(matrix, origin, endpoint)
    return matrix, trace, determinant, termination


def replay_trace(initial, trace, origin, endpoint):
    matrix = [row[:] for row in initial]
    for index, event in enumerate(trace):
        side = event["side"]
        require(side in {"left", "right"}, "unknown trace side")
        require(matrix[1][0] == event["fiber_intersection_before"], f"trace counter drift at {index}")
        gram = origin if side == "left" else endpoint
        isometry = eichler_matrix(gram, event["shift"])
        require(stable.verify_stable_matrix(isometry, gram, gram) == 1, "Eichler determinant is not one")
        matrix = stable.multiply(isometry, matrix) if side == "left" else stable.multiply(matrix, isometry)
        other = matrix[0][0] if side == "left" else matrix[1][1]
        require(other == event["opposite_hyperbolic_coordinate_after_translation"], "trace opposite coordinate drift")
        require(event["swap_performed"] == (0 < other < event["fiber_intersection_before"]), "invalid descent claim")
        if event["swap_performed"]:
            matrix = swap_hyperbolic_coordinates(matrix, side)
        require(matrix[1][0] == event["fiber_intersection_after"], "trace result drift")
    stable.verify_stable_matrix(matrix, origin, endpoint)
    return matrix


def build_certificate():
    source = audit.load_json(ROOT / INPUT_PATH)
    source_core = {key: value for key, value in source.items() if key != "record_sha256"}
    require(source["record_sha256"] == INPUT_RECORD_SHA256 == audit.payload_sha256(source_core), "short-marking input drift")
    require(audit.file_sha256(ROOT / stable.BRIDGE_PATH) == stable.BRIDGE_FILE_SHA256, "frozen bridge drift")
    bridge = audit.load_json(ROOT / stable.BRIDGE_PATH)
    origin, endpoint = bridge["origin_chain"][0]["gram"], bridge["endpoint_chain_endpoint_to_meeting"][0]["gram"]
    initial = source["rootless_to_origin_stable_matrix"]
    matrix, trace, determinant, termination = reduce_matrix(initial, origin, endpoint)
    require(replay_trace(initial, trace, origin, endpoint) == matrix, "independent action replay mismatch")
    paths = [
        "research/reduce_hyperbolic_fiber_degree.py", "research/stable_hyperbolic_neighbor_transport.py",
        "research/import_e8_a2_target_neighbor_bridge_artifact.py", INPUT_PATH, stable.BRIDGE_PATH,
    ]
    core = {
        "schema_version": 1,
        "certificate_kind": "bounded-exact-Eichler-reduction-of-hyperbolic-fiber-degree",
        "source_certificate_record_sha256": source["record_sha256"],
        "configuration": {"maximum_iterations": MAX_ITERATIONS, "maximum_nodes_per_shift": MAX_NODES_PER_SHIFT},
        "termination": termination,
        "trace": trace,
        "trace_length": len(trace),
        "cvp_nodes_visited": sum(event["search"]["nodes_visited"] for event in trace),
        "cvp_searches_exhausting_budget": sum(event["search"]["node_budget_exhausted"] for event in trace),
        "strictly_decreasing_swaps": sum(event["swap_performed"] for event in trace),
        "reduced_rootless_to_origin_stable_matrix": matrix,
        "stable_matrix_sha256": audit.payload_sha256(matrix),
        "determinant": determinant,
        "comparison": {
            "input_fiber_intersection": initial[1][0], "reduced_fiber_intersection": matrix[1][0],
            "input_maximum_absolute_entry": max(abs(x) for row in initial for x in row),
            "reduced_maximum_absolute_entry": max(abs(x) for row in matrix for x in row),
        },
        "source_sha256": {path: audit.file_sha256(ROOT / path) for path in paths},
        "claim_boundary": {
            "integral_unimodular_gram_identity_and_action_replay_verified": True,
            "global_minimum_fiber_intersection_proved": False,
            "effective_or_nef_divisor_classes_certified": False,
            "marked_geometric_K3_transport_completed": False,
            "rational_P3_found": False, "rank32_curve_found": False,
        },
    }
    return {**core, "record_sha256": audit.payload_sha256(core)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-certificate", type=Path)
    parser.add_argument("--compare", type=Path)
    args = parser.parse_args()
    record = build_certificate()
    if args.compare:
        require(audit.load_json(args.compare) == record, "Eichler certificate mismatch")
    if args.write_certificate:
        with args.write_certificate.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: record[key] for key in (
        "record_sha256", "stable_matrix_sha256", "determinant", "comparison", "termination",
        "trace_length", "cvp_nodes_visited", "cvp_searches_exhausting_budget", "claim_boundary",
    )}, indent=2))


if __name__ == "__main__":
    main()
