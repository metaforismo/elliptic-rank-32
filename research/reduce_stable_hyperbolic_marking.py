#!/usr/bin/env python3
"""Bounded smaller U-stabilizations of a frozen exact neighbor bridge.

This constructs abstract integral isometries, not effective/nef divisors,
K3 surface maps, rational points, or elliptic-rank records.
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
ORIGINAL_PATH = "certificates/stable_hyperbolic_rootless_bridge_34164507132.json"
ORIGINAL_RECORD_SHA256 = "40f65721819c5e4d4725814a7732f1bea77543b52afe778e0e6bdc84debad5a0"
MAX_LIFTS_PER_MOVE = 4


def centered_residue(value, prime):
    residue = value % prime
    return residue - prime if 2 * residue > prime else residue


def q_value(gram, vector):
    norm = stable.pairing(vector, gram, vector)
    require(norm % 2 == 0, "odd norm in an even lattice")
    return norm // 2


def short_lifts(gram, vector, prime):
    """Finite candidate family, not a shortest-vector or CVP oracle."""
    require(audit.is_prime_integer(prime), "neighbor modulus is not prime")
    require(2 <= prime <= 251, "short-lift search prime exceeds its finite bound")
    require(len(vector) == len(gram) and all(type(x) is int for x in vector), "invalid vector")
    require(any(x % prime for x in vector), "zero projective vector")
    require(q_value(gram, vector) % prime == 0, "vector is not isotropic")
    candidates = set()
    considered = 0
    for unit in range(1, prime):
        base = [centered_residue(unit * x, prime) for x in vector]
        dual = [sum(a * b for a, b in zip(row, base)) for row in gram]
        q = q_value(gram, base)
        require(q % prime == 0, "scalar multiple changed isotropy")
        if q % prime**2 == 0:
            candidates.add(tuple(base))
        for pivot, value in enumerate(dual):
            if value % prime == 0:
                continue
            correction = (-(q // prime) * pow(value, -1, prime)) % prime
            for coefficient in sorted({correction, correction - prime}):
                considered += 1
                lifted = base[:]
                lifted[pivot] += prime * coefficient
                require(q_value(gram, lifted) % prime**2 == 0, "short lift failed modulo p squared")
                candidates.add(tuple(lifted))
    require(candidates, "no admissible even lift; possible bad-prime radical boundary")
    ordered = sorted(candidates, key=lambda y: (q_value(gram, y), max(map(abs, y)), y))
    return [list(y) for y in ordered], considered


def companion_witnesses(gram, lifted, prime):
    dual = [sum(a * b for a, b in zip(row, lifted)) for row in gram]
    q = q_value(gram, lifted)
    fiber = [q // prime, prime, *lifted]
    stable_gram = stable.stabilized_gram(gram)
    witnesses = set()
    witnesses.add(tuple(stable.bezout_one([prime, q // prime, *[-x for x in dual]])))
    for pivot, value in enumerate(dual):
        if value % prime == 0:
            continue
        residue = pow(-value, -1, prime)
        for coefficient in sorted({residue, residue - prime}):
            horizontal = (1 + value * coefficient) // prime
            witness = [horizontal, 0, *[0] * len(gram)]
            witness[pivot + 2] = coefficient
            require(stable.pairing(fiber, stable_gram, witness) == 1, "pairing-one witness failed")
            witnesses.add(tuple(witness))
    return [list(z) for z in sorted(witnesses)]


def lifted_frame(gram, transform, lifted, prime):
    q = q_value(gram, lifted)
    require(q % prime**2 == 0, "lift not even modulo p squared")
    dual = [sum(a * b for a, b in zip(row, lifted)) for row in gram]
    membership = stable.integer_vector(
        [row[0] for row in stable.multiply(
            audit.invert_matrix(transform, "neighbor basis"),
            [[Fraction(x, prime)] for x in lifted],
        )], "short lift divided by p in the child basis",
    )
    pivot = next(i for i in reversed(range(len(gram))) if lifted[i] % prime)
    inverse = pow(lifted[pivot], -1, prime)
    columns = []
    for column in range(len(gram)):
        w = [row[column] for row in transform]
        t = prime * w[pivot] * inverse
        require(t.denominator == 1, "neighbor residue coefficient is not integral")
        columns.append(stable.integer_vector([
            sum(x * z for x, z in zip(w, dual)) / prime - t * Fraction(q, prime**2),
            -t,
            *[x - t * Fraction(z, prime) for x, z in zip(w, lifted)],
        ], "short lifted frame vector"))
    return columns, membership


def matrix_for_witness(gram, lifted, prime, frame, witness):
    q = q_value(gram, lifted)
    fiber = [q // prime, prime, *lifted]
    stable_gram = stable.stabilized_gram(gram)
    require(stable.pairing(fiber, stable_gram, witness) == 1, "witness does not pair to one")
    norm = stable.pairing(witness, stable_gram, witness)
    require(norm % 2 == 0, "odd witness norm")
    companion = [z - (norm // 2) * f for z, f in zip(witness, fiber)]
    columns = []
    for vector in frame:
        coefficient = stable.pairing(companion, stable_gram, vector)
        columns.append([x - coefficient * f for x, f in zip(vector, fiber)])
    return stable.transpose([fiber, companion, *columns]), fiber, companion


def reduce_move(move, parent_gram):
    exact = audit.verify_move_exact_arithmetic(move, parent_gram, "short-marking input")
    gram, child, transform = exact["parent_gram"], exact["child_gram"], exact["total_transform"]
    prime = move["prime"]
    lifts, corrections = short_lifts(gram, move["projective_isotropic_vector"], prime)
    selected_lifts = lifts[:MAX_LIFTS_PER_MOVE]
    best = None
    tested = 0
    for lifted in selected_lifts:
        frame, membership = lifted_frame(gram, transform, lifted, prime)
        for witness in companion_witnesses(gram, lifted, prime):
            tested += 1
            matrix, fiber, companion = matrix_for_witness(gram, lifted, prime, frame, witness)
            score = (max(abs(x) for row in matrix for x in row),
                     sum(abs(x) for row in matrix for x in row),
                     q_value(gram, lifted), tuple(lifted), tuple(witness))
            if best is None or score < best[0]:
                best = (score, matrix, lifted, witness, fiber, companion, membership)
    require(best is not None, "finite marking search produced no candidate")
    score, matrix, lifted, witness, fiber, companion, membership = best
    determinant = stable.verify_stable_matrix(matrix, gram, child)
    section = [z - f for z, f in zip(companion, fiber)]
    sg = stable.stabilized_gram(gram)
    require(stable.pairing(fiber, sg, fiber) == 0, "fiber is not isotropic")
    require(stable.pairing(section, sg, section) == -2 and stable.pairing(section, sg, fiber) == 1,
            "section-class identities failed")
    record = {
        "prime": prime,
        "move_sha256": move["move_sha256"],
        "parent_gram_sha256": audit.payload_sha256(gram),
        "child_gram_sha256": audit.payload_sha256(child),
        "coordinate_corrections_considered": corrections,
        "unique_admissible_lifts": len(lifts),
        "lifts_tested": len(selected_lifts),
        "pairing_one_witnesses_tested": tested,
        "selected_lift": lifted,
        "selected_lift_q": q_value(gram, lifted),
        "lift_divided_by_p_in_child_basis": membership,
        "selected_pairing_one_witness": witness,
        "abstract_fiber": fiber,
        "abstract_companion": companion,
        "abstract_zero_section": section,
        "matrix": matrix,
        "maximum_absolute_entry": score[0],
        "determinant": determinant,
        "integral_unimodular_gram_identity_verified": True,
        "geometric_effectivity_or_nefness_verified": False,
    }
    return {**record, "record_sha256": audit.payload_sha256(record)}


def build_certificate():
    original = audit.load_json(ROOT / ORIGINAL_PATH)
    original_core = {key: value for key, value in original.items() if key != "record_sha256"}
    require(original["record_sha256"] == ORIGINAL_RECORD_SHA256 == audit.payload_sha256(original_core),
            "original stable-certificate hash drift")
    bridge = stable.verified_committed_bridge()
    transports, composites = {}, {}
    for side, chain, moves in (
        ("origin", bridge["origin_chain"], bridge["origin_forward_moves"]),
        ("rootless", bridge["endpoint_chain_endpoint_to_meeting"], bridge["endpoint_forward_moves_to_invert"]),
    ):
        records = []
        composite = stable.identity(19)
        for parent, move in zip(chain, moves):
            record = reduce_move(move, parent["gram"])
            records.append(record)
            composite = stable.multiply(composite, record["matrix"])
        stable.verify_stable_matrix(composite, chain[0]["gram"], chain[-1]["gram"])
        transports[side], composites[side] = records, composite
    meeting = stable.identity(19)
    for i, row in enumerate(bridge["meeting_isometry"]["integral_unimodular_transform"]):
        for j, value in enumerate(row):
            meeting[i + 2][j + 2] = value
    matrix = stable.integer_matrix(stable.multiply(stable.multiply(composites["origin"], meeting),
                                                   stable.inverse(composites["rootless"])), "short stable composite")
    parent, child = bridge["origin_chain"][0]["gram"], bridge["endpoint_chain_endpoint_to_meeting"][0]["gram"]
    determinant = stable.verify_stable_matrix(matrix, parent, child)
    old_matrix = original["rootless_initial_to_origin_initial_stable_matrix"]
    stable.verify_stable_matrix(old_matrix, parent, child)
    paths = [
        "research/reduce_stable_hyperbolic_marking.py",
        "research/stable_hyperbolic_neighbor_transport.py",
        "research/import_e8_a2_target_neighbor_bridge_artifact.py",
        stable.BRIDGE_PATH, ORIGINAL_PATH,
    ]
    core = {
        "schema_version": 1,
        "certificate_kind": "bounded-short-lift-integral-hyperbolic-marking",
        "original_certificate_record_sha256": original["record_sha256"],
        "source_bridge_sha256": bridge["bridge_sha256"],
        "maximum_lifts_per_move": MAX_LIFTS_PER_MOVE,
        "selection_rule": "For the first four lifts ordered by (q,maxabs,vector), minimize (matrix maxabs,matrix entry L1,q,vector,witness).",
        "per_edge_transports": transports,
        "composites": composites,
        "lifted_meeting_isometry": meeting,
        "rootless_to_origin_stable_matrix": matrix,
        "stable_matrix_sha256": audit.payload_sha256(matrix),
        "determinant": determinant,
        "comparison": {
            "original_maximum_absolute_entry": max(abs(x) for row in old_matrix for x in row),
            "new_maximum_absolute_entry": max(abs(x) for row in matrix for x in row),
            "original_fiber_intersection": old_matrix[1][0],
            "new_fiber_intersection": matrix[1][0],
        },
        "source_sha256": {path: audit.file_sha256(ROOT / path) for path in paths},
        "claim_boundary": {
            "integral_unimodular_stable_gram_identity_verified": True,
            "global_minimum_marking_or_intersection_proved": False,
            "effective_or_nef_divisor_classes_certified": False,
            "marked_geometric_K3_transport_completed": False,
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
        require(audit.load_json(args.compare) == record, "short marking certificate mismatch")
    if args.write_certificate:
        with args.write_certificate.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: record[key] for key in ("record_sha256", "stable_matrix_sha256", "determinant", "comparison", "claim_boundary")}, indent=2))


if __name__ == "__main__":
    main()
