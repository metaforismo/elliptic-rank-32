#!/usr/bin/env python3
"""Bounded exact height reduction by square-minus-two reflections.

The result is an abstract integral marking. Finite root searches do not
certify nefness, a geometric K3 map, or a new elliptic-rank record.
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
from research import reduce_hyperbolic_fiber_degree as eichler  # noqa: E402

stable = eichler.stable
audit = stable.audit
require = audit.require
INPUT_PATH = "certificates/reduced_hyperbolic_fiber_degree_20260908.json"
INPUT_RECORD_SHA256 = "a9460950a15e447885206e6cff6cef71a714c5e51885e2b3ebc22507e6e09055"
MAX_ROOT_HEIGHT = 64
MAX_REFLECTIONS = 128
MAX_SEARCH_NODES = 500000


def decreasing_root(gram, fiber, height, max_nodes, gs=None):
    require(type(height) is int and height > 0, "root height must be positive")
    require(type(max_nodes) is int and max_nodes >= 0, "invalid root-search node budget")
    mu, diagonal = gs if gs is not None else eichler.gram_schmidt(gram)
    require(len(fiber) == len(gram) + 2 and all(type(x) is int for x in fiber), "invalid integral fiber")
    sg = stable.stabilized_gram(gram)
    require(stable.pairing(fiber, sg, fiber) == 0 and fiber[1] > 0, "fiber must be isotropic with positive height")
    size, degree = len(gram), fiber[1]
    center = [Fraction(height * z, degree) for z in fiber[2:]]
    coordinates = [0] * size
    errors = [-x for x in center]
    nodes = 0
    leaves = 0
    exhausted = False
    solution = None

    def visit(j, cost):
        nonlocal nodes, leaves, exhausted, solution
        if exhausted or solution is not None:
            return
        if j < 0:
            leaves += 1
            norm = stable.pairing(coordinates, gram, coordinates)
            require(norm % 2 == 0, "odd norm in root search")
            q = norm // 2
            if q % height != 1 % height:
                return
            root = [(q - 1) // height, height, *coordinates]
            require(stable.pairing(root, sg, root) == -2, "root norm identity failed")
            pairing = stable.pairing(fiber, sg, root)
            require(pairing < 0, "short root is not decreasing")
            require(degree + pairing * height > 0, "zero-height cusp needs a separate analysis")
            solution = root
            return
        offset = -center[j] + sum(mu[i][j] * errors[i] for i in range(j + 1, size))
        left = (-offset).numerator // (-offset).denominator
        right = left + 1
        while True:
            if abs(left + offset) <= abs(right + offset):
                integer = left
                left -= 1
            else:
                integer = right
                right += 1
            term = diagonal[j] * (integer + offset) ** 2
            if cost + term >= 2:
                return
            if nodes >= max_nodes:
                exhausted = True
                return
            nodes += 1
            coordinates[j], errors[j] = integer, integer - center[j]
            visit(j - 1, cost + term)
            if exhausted or solution is not None:
                return

    visit(size - 1, Fraction(0))
    return solution, {
        "root_height": height, "nodes_visited": nodes, "leaves_checked": leaves,
        "node_budget_exhausted": exhausted, "decreasing_root_found": solution is not None,
    }


def reflection_matrix(gram, root):
    sg = stable.stabilized_gram(gram)
    require(len(root) == len(sg) and all(type(x) is int for x in root), "invalid integral root")
    require(stable.pairing(root, sg, root) == -2, "reflection vector is not a square-minus-two root")
    dual = [sum(a * b for a, b in zip(row, root)) for row in sg]
    return [[int(i == j) + root[i] * dual[j] for j in range(len(sg))] for i in range(len(sg))]


def reduce_roots(initial, origin, endpoint, maximum_height=MAX_ROOT_HEIGHT,
                 maximum_reflections=MAX_REFLECTIONS, maximum_nodes=MAX_SEARCH_NODES):
    require(all(type(x) is int and x >= 0 for x in (maximum_height, maximum_reflections, maximum_nodes)), "invalid root-search bounds")
    stable.verify_stable_matrix(initial, origin, endpoint)
    matrix = [row[:] for row in initial]
    gs = {"left": eichler.gram_schmidt(origin), "right": eichler.gram_schmidt(endpoint)}
    endpoint_inverse = stable.inverse(endpoint)
    trace, searches = [], []
    nodes = 0
    termination = "reflection_budget_exhausted"
    for iteration in range(maximum_reflections):
        selected = None
        for side, gram in (("left", origin), ("right", endpoint)):
            if side == "left":
                fiber = [row[0] for row in matrix]
            else:
                tail = matrix[1][2:]
                inverse_tail = [-sum(a * b for a, b in zip(row, tail)) for row in endpoint_inverse]
                fiber = stable.integer_vector([matrix[1][1], matrix[1][0], *inverse_tail], "inverse image of old fiber")
            for height in range(1, maximum_height + 1):
                root, search = decreasing_root(gram, fiber, height, maximum_nodes - nodes, gs[side])
                nodes += search["nodes_visited"]
                searches.append({"iteration": iteration, "side": side, **search})
                if root is not None:
                    selected = side, gram, fiber, root
                    break
                if search["node_budget_exhausted"]:
                    break
            if selected is not None or nodes >= maximum_nodes:
                break
        if selected is None:
            termination = "search_node_budget_exhausted" if nodes >= maximum_nodes else "no_decreasing_root_in_recorded_height_window"
            break
        side, gram, fiber, root = selected
        reflection = reflection_matrix(gram, root)
        require(stable.verify_stable_matrix(reflection, gram, gram) == -1, "root reflection determinant failed")
        before = matrix[1][0]
        pairing = stable.pairing(fiber, stable.stabilized_gram(gram), root)
        matrix = stable.multiply(reflection, matrix) if side == "left" else stable.multiply(matrix, reflection)
        after = matrix[1][0]
        require(after == before + root[1] * pairing and 0 < after < before, "root reflection did not decrease the degree")
        trace.append({
            "iteration": iteration, "side": side, "root": root,
            "fiber_root_pairing": pairing, "degree_before": before, "degree_after": after,
        })
    stable.verify_stable_matrix(matrix, origin, endpoint)
    return matrix, trace, searches, termination


def replay_reflections(initial, trace, origin, endpoint):
    matrix = [row[:] for row in initial]
    for event in trace:
        side = event["side"]
        require(side in {"left", "right"}, "unknown reflection side")
        gram = origin if side == "left" else endpoint
        require(matrix[1][0] == event["degree_before"], "reflection trace input drift")
        reflection = reflection_matrix(gram, event["root"])
        matrix = stable.multiply(reflection, matrix) if side == "left" else stable.multiply(matrix, reflection)
        require(matrix[1][0] == event["degree_after"], "reflection trace output drift")
        require(event["degree_after"] == event["degree_before"] + event["root"][1] * event["fiber_root_pairing"], "reflection pairing drift")
        require(0 < event["degree_after"] < event["degree_before"], "reflection trace is not decreasing")
    stable.verify_stable_matrix(matrix, origin, endpoint)
    return matrix


def build_certificate():
    source = audit.load_json(ROOT / INPUT_PATH)
    core = {key: value for key, value in source.items() if key != "record_sha256"}
    require(source["record_sha256"] == INPUT_RECORD_SHA256 == audit.payload_sha256(core), "Eichler input drift")
    require(audit.file_sha256(ROOT / stable.BRIDGE_PATH) == stable.BRIDGE_FILE_SHA256, "bridge input drift")
    bridge = audit.load_json(ROOT / stable.BRIDGE_PATH)
    origin, endpoint = bridge["origin_chain"][0]["gram"], bridge["endpoint_chain_endpoint_to_meeting"][0]["gram"]
    initial = source["reduced_rootless_to_origin_stable_matrix"]
    reflected, trace, searches, termination = reduce_roots(initial, origin, endpoint)
    require(replay_reflections(initial, trace, origin, endpoint) == reflected, "reflection replay mismatch")
    matrix, cleanup, determinant, cleanup_termination = eichler.reduce_matrix(reflected, origin, endpoint)
    require(eichler.replay_trace(reflected, cleanup, origin, endpoint) == matrix, "cleanup replay mismatch")
    paths = [
        "research/reduce_hyperbolic_root_reflections.py", "research/reduce_hyperbolic_fiber_degree.py",
        "research/stable_hyperbolic_neighbor_transport.py", "research/import_e8_a2_target_neighbor_bridge_artifact.py",
        INPUT_PATH, stable.BRIDGE_PATH,
    ]
    core = {
        "schema_version": 1,
        "certificate_kind": "bounded-square-minus-two-root-reflection-reduction",
        "source_certificate_record_sha256": source["record_sha256"],
        "configuration": {"maximum_root_height": MAX_ROOT_HEIGHT, "maximum_reflections": MAX_REFLECTIONS, "maximum_search_nodes": MAX_SEARCH_NODES},
        "termination": termination,
        "root_searches": searches,
        "search_nodes_visited": sum(item["nodes_visited"] for item in searches),
        "reflection_trace": trace,
        "reflection_count": len(trace),
        "postprocessing_eichler_trace": cleanup,
        "postprocessing_termination": cleanup_termination,
        "rootless_to_origin_stable_matrix": matrix,
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
            "all_negative_pairing_roots_excluded": False,
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
        require(audit.load_json(args.compare) == record, "root-reflection certificate mismatch")
    if args.write_certificate:
        with args.write_certificate.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: record[key] for key in (
        "record_sha256", "stable_matrix_sha256", "determinant", "comparison", "termination",
        "reflection_count", "search_nodes_visited", "claim_boundary",
    )}, indent=2))


if __name__ == "__main__":
    main()
