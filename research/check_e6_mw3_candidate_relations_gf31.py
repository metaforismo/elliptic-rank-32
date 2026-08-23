#!/usr/bin/env python3
"""Exact bounded relation search for a modular E6/MW3 section triple.

This is a discovery/collision filter, not a substitute for a Shioda height
pairing.  It loads an exact core-slice certificate, reconstructs the canonical
P2 and polynomial P3 searches, and evaluates the group law in GF(31)(t).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import replay_e6_mw3_core_exhaustion_gf31 as core_module
import replay_e6_mw3_p2_split_cores as p2_module

from replay_e6_mw3_p2_split_cores import (
    PRIME,
    canonical_sha256,
    poly_add,
    poly_mul,
    poly_neg,
    poly_pow,
    poly_scale,
    poly_sub,
    search_p2,
    search_p3_polynomial,
    trim,
)
from replay_e6_mw3_core_exhaustion_gf31 import poly_exact_divide, poly_gcd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = Path("/private/tmp/e6_slice_3_1_1.json")
DEFAULT_OUTPUT = ROOT / "certificates" / "e6_mw3_candidate_r3_s1_x1_core_1_0_14_23_relations.json"


def configure_prime(prime):
    global PRIME, ZERO, ONE
    PRIME = prime
    p2_module.PRIME = prime
    core_module.PRIME = prime
    ZERO = rational([0])
    ONE = rational([1])


def rational(numerator, denominator=(1,)):
    numerator = trim(list(numerator))
    denominator = trim(list(denominator))
    if denominator == [0]:
        raise ZeroDivisionError("zero rational-function denominator")
    if numerator == [0]:
        return ((0,), (1,))
    common = poly_gcd(numerator, denominator)
    numerator = poly_exact_divide(numerator, common)
    denominator = poly_exact_divide(denominator, common)
    scalar = pow(denominator[-1], -1, PRIME)
    return (
        tuple(poly_scale(numerator, scalar)),
        tuple(poly_scale(denominator, scalar)),
    )


ZERO = rational([0])
ONE = rational([1])


def radd(left, right):
    return rational(
        poly_add(
            poly_mul(list(left[0]), list(right[1])),
            poly_mul(list(right[0]), list(left[1])),
        ),
        poly_mul(list(left[1]), list(right[1])),
    )


def rneg(value):
    return rational(poly_neg(list(value[0])), value[1])


def rsub(left, right):
    return radd(left, rneg(right))


def rmul(left, right):
    return rational(
        poly_mul(list(left[0]), list(right[0])),
        poly_mul(list(left[1]), list(right[1])),
    )


def rdiv(left, right):
    if right[0] == (0,):
        raise ZeroDivisionError("rational-function division by zero")
    return rational(
        poly_mul(list(left[0]), list(right[1])),
        poly_mul(list(left[1]), list(right[0])),
    )


def rsquare(value):
    return rational(poly_pow(list(value[0]), 2), poly_pow(list(value[1]), 2))


def polynomial(coefficients):
    return rational(coefficients)


def equal(left, right):
    return left == right


def negate(point):
    if point is None:
        return None
    return point[0], rneg(point[1])


def add(left, right, A):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if equal(x1, x2):
        if equal(y1, rneg(y2)):
            return None
        if equal(y1, ZERO):
            return None
        slope = rdiv(
            radd(rmul(polynomial([3]), rsquare(x1)), A),
            rmul(polynomial([2]), y1),
        )
    else:
        slope = rdiv(rsub(y2, y1), rsub(x2, x1))
    x3 = rsub(rsub(rsquare(slope), x1), x2)
    y3 = rsub(rmul(slope, rsub(x1, x3)), y1)
    return x3, y3


def multiply(coefficient, point, A):
    if coefficient < 0:
        return multiply(-coefficient, negate(point), A)
    result = None
    addend = point
    value = coefficient
    while value:
        if value & 1:
            result = add(result, addend, A)
        addend = add(addend, addend, A)
        value >>= 1
    return result


def point_on_curve(point, A, B):
    if point is None:
        return True
    x, y = point
    return equal(
        rsub(
            rsub(rsquare(y), rmul(rsquare(x), x)),
            radd(rmul(A, x), B),
        ),
        ZERO,
    )


def point_fingerprint(point):
    if point is None:
        return "O"
    return [
        {"numerator": list(point[0][0]), "denominator": list(point[0][1])},
        {"numerator": list(point[1][0]), "denominator": list(point[1][1])},
    ]


def search_relation_boxes(
    prime,
    A_coefficients,
    B_coefficients,
    P1_data,
    P2_data,
    P3_data,
    bound,
):
    """Search one coefficient box for several P3 candidates on one surface."""
    configure_prime(prime)
    A = polynomial(A_coefficients)
    B = polynomial(B_coefficients)
    P1 = polynomial(P1_data["X"]), polynomial(P1_data["Y"])
    z = [-P2_data["pole"], 1]
    P2 = rational(P2_data["X2"], poly_pow(z, 2)), rational(
        P2_data["Y2"], poly_pow(z, 3)
    )
    candidates = [
        (polynomial(item["X3"]), polynomial(item["Y3_up_to_sign"]))
        for item in P3_data
    ]
    if not all(point_on_curve(point, A, B) for point in [P1, P2, *candidates]):
        raise AssertionError("one reconstructed point failed the curve equation")

    coefficients = list(range(-bound, bound + 1))
    p1_multiples = {
        coefficient: multiply(coefficient, P1, A) for coefficient in coefficients
    }
    p2_multiples = {
        coefficient: multiply(coefficient, P2, A) for coefficient in coefficients
    }
    partials = {
        (first, second): add(p1_multiples[first], p2_multiples[second], A)
        for first in coefficients
        for second in coefficients
    }
    results = []
    checked = (2 * bound + 1) ** 3 - 1
    for data, candidate in zip(P3_data, candidates, strict=True):
        p3_multiples = {
            coefficient: multiply(coefficient, candidate, A)
            for coefficient in coefficients
        }
        relations = []
        for first in coefficients:
            for second in coefficients:
                partial = partials[(first, second)]
                for third in coefficients:
                    if first == second == third == 0:
                        continue
                    if partial == negate(p3_multiples[third]):
                        relations.append([first, second, third])
        results.append(
            {
                "P3": data,
                "nonzero_triples_checked": checked,
                "relations": relations,
                "no_relation_in_box": not relations,
                "point_fingerprints": [
                    point_fingerprint(point) for point in (P1, P2, candidate)
                ],
            }
        )
    return results


def compute(input_path: Path, core, bound: int):
    configure_prime(31)
    source = json.loads(input_path.read_text())
    record = next(
        item for item in source["valid_two_I2_surfaces"] if item["core"] == list(core)
    )
    coordinate_slice = source["coordinate_slice"]
    s0 = coordinate_slice["s0"]
    A_coefficients = record["A"]
    B_coefficients = record["B"]
    roots = record["valid_I2_roots"]
    if len(roots) != 2:
        raise AssertionError(f"expected exactly two I2 roots, got {roots}")
    nodes = {int(key): value for key, value in record["singular_nodes"].items()}
    lam, mu = roots
    p2_tested, p2_hits = search_p2(
        A_coefficients,
        B_coefficients,
        lam,
        mu,
        nodes[lam][0],
        nodes[mu][0],
        s0=s0,
    )
    p3_tested, p3_hits = search_p3_polynomial(
        A_coefficients,
        B_coefficients,
        core[3],
        lam,
        mu,
        nodes[lam][0],
        nodes[mu][0],
        s0=s0,
    )
    if len(p2_hits) != 2 or len(p3_hits) != 1:
        raise AssertionError(
            f"expected one geometric P2 and one P3, got {len(p2_hits)}, {len(p3_hits)}"
        )
    p2_data = p2_hits[0]
    p3_data = p3_hits[0]
    relation_result = search_relation_boxes(
        31,
        A_coefficients,
        B_coefficients,
        {"X": record["X1"], "Y": record["Y1"]},
        p2_data,
        [p3_data],
        bound,
    )[0]
    relations = relation_result["relations"]
    checked = relation_result["nonzero_triples_checked"]

    payload = {
        "schema": "elliptic-rank30/e6-mw3-bounded-relations-gf31/v1",
        "field": "GF(31)(t)",
        "source_certificate": str(input_path),
        "source_certificate_sha256": source["certificate_sha256"],
        "coordinate_slice": coordinate_slice,
        "core_a1_a2_a4_s1": list(core),
        "A": A_coefficients,
        "B": B_coefficients,
        "P1": {"X": record["X1"], "Y": record["Y1"]},
        "P2": p2_data,
        "P3": p3_data,
        "P2_cases_tested": p2_tested,
        "P3_quadratic_X_cases_tested": p3_tested,
        "coefficient_box": [-bound, bound],
        "nonzero_triples_checked": checked,
        "relations": relations,
        "no_relation_in_box": not relations,
        "point_fingerprints": relation_result["point_fingerprints"],
        "claim_boundary": (
            "An exact bounded collision filter only. Full MW-rank-3 independence "
            "still requires local component verification and a positive Shioda regulator."
        ),
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--core", default="1,0,14,23")
    parser.add_argument("--bound", type=int, default=8)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    core = tuple(int(value) for value in args.core.split(","))
    if len(core) != 4:
        raise SystemExit("--core must contain a1,a2,a4,s1")
    payload = compute(args.input, core, args.bound)
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "relations": payload["relations"],
                "nonzero_triples_checked": payload["nonzero_triples_checked"],
                "certificate_sha256": payload["certificate_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
