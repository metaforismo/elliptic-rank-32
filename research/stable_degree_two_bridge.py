#!/usr/bin/env python3
"""Integral U-stabilization of the frozen fifteen-edge, all-p=2 bridge.

Degree two here means adjacent abstract fiber intersections. It does not
mean a geometric map of degree two, or an end-to-end intersection of two.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from research import reduce_stable_hyperbolic_marking as short  # noqa: E402

stable = short.stable
audit = stable.audit
require = audit.require
BRIDGE_PATH = "certificates/e8_a2_transparent_p2_bridge_34174072896/exact-bridge.json"
BRIDGE_FILE_SHA256 = "54e2d92b6163f9339f26cec977b47ff09c2ad586f0b3df029ebd38c05d878e33"


def verified_committed_bridge():
    path = ROOT / BRIDGE_PATH
    require(audit.file_sha256(path) == BRIDGE_FILE_SHA256, "frozen p=2 bridge file hash drift")
    bridge = audit.load_json(path)
    require(bridge["endpoint"] == "transparent", "wrong endpoint for the p=2 bridge")
    states = {state["uid"]: state for state in bridge["origin_chain"] + bridge["endpoint_chain_endpoint_to_meeting"]}
    all_moves = bridge["origin_forward_moves"] + bridge["endpoint_forward_moves_to_invert"]
    require(len(all_moves) == bridge["neighbor_step_count"] == 15, "wrong p=2 bridge length")
    require(all(move["prime"] == 2 for move in all_moves), "bridge contains a non-two neighbor")
    moves = {move["move_sha256"]: move for move in all_moves}
    for uid, state in states.items():
        audit.verify_state(state, uid)
    require(bridge["origin_chain"][0]["gram_sha256"] == audit.INITIAL_HASHES["origin"], "origin input drift")
    require(bridge["endpoint_chain_endpoint_to_meeting"][0]["gram_sha256"] == audit.INITIAL_HASHES["transparent"], "transparent input drift")
    audit.verify_bridge(path.parent, {"status": "bridge_found", "termination_reason": "exact_bridge_found", "bridge": bridge}, states, moves)
    return bridge


def build_certificate():
    bridge = verified_committed_bridge()
    transports, composites = {}, {}
    for side, chain, moves in (
        ("origin", bridge["origin_chain"], bridge["origin_forward_moves"]),
        ("transparent", bridge["endpoint_chain_endpoint_to_meeting"], bridge["endpoint_forward_moves_to_invert"]),
    ):
        records, composite = [], stable.identity(19)
        for parent, move in zip(chain, moves):
            record = short.reduce_move(move, parent["gram"])
            require(record["abstract_fiber"][1] == 2, "adjacent fiber intersection is not two")
            records.append(record)
            composite = stable.multiply(composite, record["matrix"])
        stable.verify_stable_matrix(composite, chain[0]["gram"], chain[-1]["gram"])
        transports[side], composites[side] = records, composite
    meeting = stable.identity(19)
    for i, row in enumerate(bridge["meeting_isometry"]["integral_unimodular_transform"]):
        for j, value in enumerate(row):
            meeting[i + 2][j + 2] = value
    matrix = stable.integer_matrix(
        stable.multiply(stable.multiply(composites["origin"], meeting), stable.inverse(composites["transparent"])),
        "all-p=2 stabilized composite",
    )
    origin, endpoint = bridge["origin_chain"][0]["gram"], bridge["endpoint_chain_endpoint_to_meeting"][0]["gram"]
    determinant = stable.verify_stable_matrix(matrix, origin, endpoint)
    paths = [
        "research/stable_degree_two_bridge.py", "research/reduce_stable_hyperbolic_marking.py",
        "research/stable_hyperbolic_neighbor_transport.py", "research/import_e8_a2_target_neighbor_bridge_artifact.py",
        BRIDGE_PATH,
    ]
    core = {
        "schema_version": 1,
        "certificate_kind": "explicit-integral-hyperbolic-stabilization-of-fifteen-degree-two-neighbors",
        "source_run_id": 34174072896,
        "source_bridge_sha256": bridge["bridge_sha256"],
        "endpoint": "transparent_A11_direct_sum_K6",
        "matrix_convention": "Columns express child coordinates in parent coordinates; S=U direct-sum (-G). The endpoint chain is inverted at composition.",
        "per_edge_transports": transports,
        "composites": composites,
        "lifted_meeting_isometry": meeting,
        "transparent_to_origin_stable_matrix": matrix,
        "stable_matrix_sha256": audit.payload_sha256(matrix),
        "determinant": determinant,
        "neighbor_step_count": 15,
        "every_adjacent_fiber_intersection": 2,
        "end_to_end_fiber_intersection": matrix[1][0],
        "maximum_absolute_matrix_entry": max(abs(x) for row in matrix for x in row),
        "source_sha256": {path: audit.file_sha256(ROOT / path) for path in paths},
        "claim_boundary": {
            "all_fifteen_edges_extended_integrally": True,
            "integral_unimodular_stable_gram_identity_verified": True,
            "endpoint_is_the_rootless_lattice": False,
            "end_to_end_fiber_intersection_is_two": matrix[1][0] == 2,
            "effective_or_nef_divisor_classes_certified": False,
            "period_or_ample_cone_compatibility_certified": False,
            "marked_geometric_K3_transport_completed": False,
            "riemann_roch_pencils_or_birational_maps_constructed": False,
            "rational_P3_found": False,
            "rank32_curve_found": False,
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
        require(audit.load_json(args.compare) == record, "all-p=2 stable certificate mismatch")
    if args.write_certificate:
        with args.write_certificate.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: record[key] for key in (
        "record_sha256", "stable_matrix_sha256", "determinant", "neighbor_step_count",
        "end_to_end_fiber_intersection", "maximum_absolute_matrix_entry", "claim_boundary",
    )}, indent=2))


if __name__ == "__main__":
    main()
