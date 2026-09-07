#!/usr/bin/env python3
"""Construct exact integral U-stabilizations of the certified neighbor bridge.

This is an abstract lattice calculation. No nef/effective cone, period,
Weierstrass surface, birational map, or new rational point is certified.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from research import import_e8_a2_target_neighbor_bridge_artifact as audit  # noqa: E402

BRIDGE_PATH = "certificates/e8_a2_rootless_bridge_34164507132/exact-bridge.json"
BRIDGE_FILE_SHA256 = "0714dd21746e9be94e53e91b675d445487411e6c9b140af17523543daaf60139"
require = audit.require


def identity(size):
    return [[int(i == j) for j in range(size)] for i in range(size)]


def transpose(matrix):
    return [list(column) for column in zip(*matrix)]


def multiply(left, right):
    require(left and right and len(left[0]) == len(right), "matrix product shape mismatch")
    return [[sum(x * y for x, y in zip(row, column)) for column in zip(*right)] for row in left]


def inverse(matrix):
    size = len(matrix)
    require(all(len(row) == size for row in matrix), "inverse requires a square matrix")
    work = [[Fraction(value) for value in row] + identity(size)[i] for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((i for i in range(column, size) if work[i][column]), None)
        require(pivot is not None, "singular stable matrix")
        work[column], work[pivot] = work[pivot], work[column]
        factor = work[column][column]
        work[column] = [value / factor for value in work[column]]
        for i in range(size):
            if i != column:
                factor = work[i][column]
                work[i] = [a - factor * b for a, b in zip(work[i], work[column])]
    return [row[size:] for row in work]


def integer_vector(vector, label):
    require(all(Fraction(value).denominator == 1 for value in vector), f"{label} is not integral")
    return [int(value) for value in vector]


def integer_matrix(matrix, label):
    return [integer_vector(row, label) for row in matrix]


def stabilized_gram(gram):
    size = len(gram) + 2
    return [[
        int(i + j == 1) if i < 2 and j < 2
        else -gram[i - 2][j - 2] if i >= 2 and j >= 2 else 0
        for j in range(size)
    ] for i in range(size)]


def pairing(left, gram, right):
    return sum(left[i] * sum(value * x for value, x in zip(row, right)) for i, row in enumerate(gram))


def extended_gcd(a, b):
    if b == 0:
        return abs(a), 1 if a > 0 else -1, 0
    gcd, x, y = extended_gcd(b, a % b)
    return gcd, y, x - (a // b) * y


def bezout_one(values):
    coefficients = []
    gcd = 0
    for value in values:
        gcd, x, y = extended_gcd(gcd, value)
        coefficients = [x * old for old in coefficients] + [y]
    require(gcd == 1, "isotropic fiber has divisibility greater than one")
    require(sum(x * y for x, y in zip(coefficients, values)) == 1, "Bezout identity failed")
    return coefficients


def lift_isotropic_vector(gram, vector, prime):
    """Reproduce the even lift used by the pinned Sage neighbor construction."""
    require(audit.is_prime_integer(prime), "neighbor modulus is not prime")
    size = len(gram)
    require(len(vector) == size and all(type(value) is int for value in vector), "invalid isotropic vector")
    y = vector[:]
    dual = [sum(value * x for value, x in zip(row, y)) for row in gram]
    norm = sum(x * z for x, z in zip(y, dual))
    require(norm % (2 * prime) == 0, "vector is not isotropic for q modulo p")
    if prime != 2 and norm % prime**2:
        pivot = next((i for i, value in enumerate(dual) if value % prime), None)
        require(pivot is not None, "radical line: odd-prime lift is unavailable")
        y[pivot] -= norm * pow(2 * dual[pivot], -1, prime)
    elif prime == 2 and norm % 8 == 4:
        pivot = next((i for i, value in enumerate(dual) if value % 2), None)
        require(pivot is not None, "radical line: non-maximal even p=2 boundary")
        y[pivot] += 2
    q = pairing(y, gram, y) // 2
    require(q % prime**2 == 0, "lift does not have q divisible by p squared")
    require(all((a - b) % prime == 0 for a, b in zip(y, vector)), "lift changed the projective line")
    return y, q


def verify_stable_matrix(matrix, parent_gram, child_gram):
    size = len(parent_gram) + 2
    require(len(matrix) == size and all(len(row) == size for row in matrix), "wrong stable-matrix dimension")
    matrix = integer_matrix(matrix, "stable basis matrix")
    transported = multiply(transpose(matrix), multiply(stabilized_gram(parent_gram), matrix))
    require(transported == stabilized_gram(child_gram), "stable Gram identity failed")
    determinant = audit.determinant_bareiss(matrix)
    require(abs(determinant) == 1, "stable basis matrix is not unimodular")
    return determinant


def extend_move(move, parent_gram):
    exact = audit.verify_move_exact_arithmetic(move, parent_gram, "stable-extension input")
    gram = exact["parent_gram"]
    child = exact["child_gram"]
    transform = exact["total_transform"]
    prime = move["prime"]
    y, q = lift_isotropic_vector(gram, move["projective_isotropic_vector"], prime)
    size = len(gram)
    dual = [sum(value * x for value, x in zip(row, y)) for row in gram]
    membership = integer_vector(
        [row[0] for row in multiply(audit.invert_matrix(transform, "neighbor basis"), [[Fraction(x, prime)] for x in y])],
        "lift divided by p in the child basis",
    )
    stable = stabilized_gram(gram)
    fiber = [q // prime, prime, *y]
    dual_fiber = [prime, q // prime, *[-value for value in dual]]
    bezout = bezout_one(dual_fiber)
    bezout_norm = pairing(bezout, stable, bezout)
    require(bezout_norm % 2 == 0, "stabilized lattice is not even")
    companion = [z - (bezout_norm // 2) * f for z, f in zip(bezout, fiber)]
    pivot = next((i for i in reversed(range(size)) if y[i] % prime), None)
    require(pivot is not None, "lift is not primitive modulo p")
    unit = pow(y[pivot], -1, prime)
    frame_columns = []
    for column in range(size):
        w = [row[column] for row in transform]
        t = prime * w[pivot] * unit
        require(t.denominator == 1, "neighbor residue coefficient is not integral")
        lifted = [
            sum(x * z for x, z in zip(w, dual)) / prime - t * Fraction(q, prime**2),
            -t,
            *[x - t * Fraction(z, prime) for x, z in zip(w, y)],
        ]
        lifted = integer_vector(lifted, "lifted frame vector")
        projection = pairing(companion, stable, lifted)
        frame_columns.append([x - projection * f for x, f in zip(lifted, fiber)])
    matrix = transpose([fiber, companion, *frame_columns])
    determinant = verify_stable_matrix(matrix, gram, child)
    section = [z - f for z, f in zip(companion, fiber)]
    require(pairing(fiber, stable, fiber) == 0, "new fiber is not isotropic")
    require(pairing(fiber, stable, companion) == 1 and pairing(companion, stable, companion) == 0, "new U basis failed")
    require(pairing(section, stable, section) == -2 and pairing(section, stable, fiber) == 1, "abstract section-class identities failed")
    core = {
        "prime": prime,
        "move_sha256": move["move_sha256"],
        "parent_gram_sha256": audit.payload_sha256(gram),
        "child_gram_sha256": audit.payload_sha256(child),
        "lifted_isotropic_vector": y,
        "q_of_lifted_vector": q,
        "lift_divided_by_p_in_child_basis": membership,
        "abstract_new_fiber_class": fiber,
        "bezout_pairing_one_witness": bezout,
        "abstract_new_hyperbolic_companion": companion,
        "abstract_new_zero_section_class": section,
        "old_new_fiber_intersection": prime,
        "child_stabilized_basis_in_parent_coordinates": matrix,
        "stable_matrix_sha256": audit.payload_sha256(matrix),
        "determinant": determinant,
        "integral_unimodular_gram_identity_verified": True,
        "geometric_effectivity_or_nefness_verified": False,
    }
    return {**core, "transport_sha256": audit.payload_sha256(core)}


def verified_committed_bridge():
    path = ROOT / BRIDGE_PATH
    require(audit.file_sha256(path) == BRIDGE_FILE_SHA256, "frozen bridge file hash drift")
    bridge = audit.load_json(path)
    states = {state["uid"]: state for state in bridge["origin_chain"] + bridge["endpoint_chain_endpoint_to_meeting"]}
    moves = {move["move_sha256"]: move for move in bridge["origin_forward_moves"] + bridge["endpoint_forward_moves_to_invert"]}
    for uid, state in states.items():
        audit.verify_state(state, uid)
    require(bridge["origin_chain"][0]["gram_sha256"] == audit.INITIAL_HASHES["origin"], "origin input drift")
    require(bridge["endpoint_chain_endpoint_to_meeting"][0]["gram_sha256"] == audit.INITIAL_HASHES["rootless"], "rootless input drift")
    audit.verify_bridge(path.parent, {"status": "bridge_found", "termination_reason": "exact_bridge_found", "bridge": bridge}, states, moves)
    return bridge


def build_certificate():
    bridge = verified_committed_bridge()
    transports = {}
    composites = {}
    for name, chain, moves in (
        ("origin", bridge["origin_chain"], bridge["origin_forward_moves"]),
        ("rootless", bridge["endpoint_chain_endpoint_to_meeting"], bridge["endpoint_forward_moves_to_invert"]),
    ):
        composite = identity(19)
        records = []
        for parent, move in zip(chain, moves):
            record = extend_move(move, parent["gram"])
            records.append(record)
            composite = multiply(composite, record["child_stabilized_basis_in_parent_coordinates"])
        verify_stable_matrix(composite, chain[0]["gram"], chain[-1]["gram"])
        transports[name] = records
        composites[name] = composite
    meeting = identity(19)
    for i, row in enumerate(bridge["meeting_isometry"]["integral_unimodular_transform"]):
        for j, value in enumerate(row):
            meeting[i + 2][j + 2] = value
    matrix = integer_matrix(multiply(multiply(composites["origin"], meeting), inverse(composites["rootless"])), "end-to-end stable matrix")
    parent = bridge["origin_chain"][0]["gram"]
    child = bridge["endpoint_chain_endpoint_to_meeting"][0]["gram"]
    determinant = verify_stable_matrix(matrix, parent, child)
    paths = [
        "research/stable_hyperbolic_neighbor_transport.py",
        "research/import_e8_a2_target_neighbor_bridge_artifact.py",
        "research/sample_projective_lines.py",
        BRIDGE_PATH,
    ]
    core = {
        "schema_version": 1,
        "certificate_kind": "explicit-integral-hyperbolic-stabilization-of-seven-edge-rootless-bridge",
        "source_run_id": 34164507132,
        "source_bridge_sha256": bridge["bridge_sha256"],
        "matrix_convention": "Columns express child (e,f,N_child) coordinates in parent (e,f,N_parent); S=U direct-sum (-G), U=[[0,1],[1,0]].",
        "per_edge_transports": transports,
        "origin_forward_composite": composites["origin"],
        "rootless_forward_composite": composites["rootless"],
        "lifted_meeting_isometry": meeting,
        "rootless_initial_to_origin_initial_stable_matrix": matrix,
        "stable_matrix_sha256": audit.payload_sha256(matrix),
        "stable_matrix_determinant": determinant,
        "old_fiber_intersection_with_transported_rootless_fiber": matrix[1][0],
        "integral_unimodular_stable_gram_identity_verified": True,
        "source_sha256": {path: audit.file_sha256(ROOT / path) for path in paths},
        "claim_boundary": {
            "abstract_U_stabilized_lattice_isometry_constructed": True,
            "all_seven_neighbor_edges_extended_integrally": True,
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
        require(audit.load_json(args.compare) == record, "stable transport certificate mismatch")
    if args.write_certificate:
        with args.write_certificate.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: record[key] for key in (
        "record_sha256", "stable_matrix_sha256", "stable_matrix_determinant",
        "old_fiber_intersection_with_transported_rootless_fiber", "claim_boundary",
    )}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
