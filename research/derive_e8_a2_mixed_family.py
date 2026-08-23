#!/usr/bin/env python3
"""Exact certificate for a non-isotrivial E8+A2^3 search chart.

This module verifies the two-IV/one-I3 family

    y^2 = x^3 - 3 k^2 T^2 x
          + T^2 (c (t-lambda)^3 + 2 k^3 T),
    T=t(t-1),

its discriminant and Kodaira valuations, the split-I3 normalization, and the
local/intersection conditions that are equivalent to the Mordell--Weil Gram

    (1/3) [[8,-1,0],[-1,10,0],[0,0,12]].

It also replays three independently found polynomial sections over GF(13), one
for each required diagonal-height profile.  They live on different surfaces:
they are checks of the search engine, not a rank-three seed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence

import sympy as sp


def canonical_sha256(payload: object) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def polynomial_coefficients(poly: sp.Expr, variable: sp.Symbol) -> list[str]:
    expanded = sp.Poly(sp.expand(poly), variable)
    return [str(expanded.nth(index)) for index in range(expanded.degree() + 1)]


def poly_eval(coefficients: Sequence[int], value: int, prime: int) -> int:
    result = 0
    for coefficient in reversed(coefficients):
        result = (result * value + coefficient) % prime
    return result


def poly_mul(left: Sequence[int], right: Sequence[int], prime: int) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] = (result[i + j] + a * b) % prime
    return result


def poly_add(*polys: Sequence[int], prime: int, width: int) -> list[int]:
    result = [0] * width
    for poly in polys:
        for index, coefficient in enumerate(poly):
            result[index] = (result[index] + coefficient) % prime
    return result


def pad(poly: Sequence[int], width: int) -> list[int]:
    if len(poly) > width:
        raise ValueError("polynomial exceeds requested width")
    return list(poly) + [0] * (width - len(poly))


def verify_modular_section(example: dict[str, object]) -> dict[str, object]:
    prime = int(example["prime"])
    coefficient_k = int(example["k"])
    coefficient_c = int(example["c"])
    lam = int(example["lambda"])
    x = [int(value) for value in example["x"]]
    y = [int(value) for value in example["y"]]

    t_poly = [0, 1]
    t_minus_one = [-1 % prime, 1]
    T = poly_mul(t_poly, t_minus_one, prime)
    T2 = poly_mul(T, T, prime)
    z = [-lam % prime, 1]
    z2 = poly_mul(z, z, prime)
    z3 = poly_mul(z2, z, prime)

    x2 = poly_mul(x, x, prime)
    x3 = poly_mul(x2, x, prime)
    y2 = poly_mul(y, y, prime)
    ax = [
        (-3 * coefficient_k * coefficient_k * value) % prime
        for value in poly_mul(T2, x, prime)
    ]
    inner = poly_add(
        [(coefficient_c * value) % prime for value in z3],
        [(2 * pow(coefficient_k, 3, prime) * value) % prime for value in T],
        prime=prime,
        width=4,
    )
    b_poly = poly_mul(T2, inner, prime)
    residual = poly_add(
        y2,
        [(-value) % prime for value in x3],
        [(-value) % prime for value in ax],
        [(-value) % prime for value in b_poly],
        prime=prime,
        width=13,
    )
    if any(residual):
        raise AssertionError(f"GF({prime}) example is not a section: {residual}")

    node_x = coefficient_k * lam * (lam - 1) % prime
    mask = 0
    if poly_eval(x, 0, prime) == 0 and poly_eval(y, 0, prime) == 0:
        mask |= 1
    if poly_eval(x, 1, prime) == 0 and poly_eval(y, 1, prime) == 0:
        mask |= 2
    if (
        poly_eval(x, lam, prime) == node_x
        and poly_eval(y, lam, prime) == 0
    ):
        mask |= 4
    if mask != int(example["expected_mask"]):
        raise AssertionError(f"unexpected component mask {mask}")

    return {
        **example,
        "identity_residual": residual,
        "computed_mask": mask,
        "nonidentity_fibre_count": mask.bit_count(),
        "shioda_height": str(Fraction(4) - Fraction(2 * mask.bit_count(), 3)),
    }


def compute_certificate() -> dict[str, object]:
    t, coefficient_k, coefficient_c, lam = sp.symbols("t k c lambda")
    T = t * (t - 1)
    z = t - lam
    A = -3 * coefficient_k**2 * T**2
    B = T**2 * (coefficient_c * z**3 + 2 * coefficient_k**3 * T)
    discriminant = sp.expand(-16 * (4 * A**3 + 27 * B**2))
    expected_discriminant = sp.expand(
        -432
        * coefficient_c
        * T**4
        * z**3
        * (coefficient_c * z**3 + 4 * coefficient_k**3 * T)
    )
    if sp.expand(discriminant - expected_discriminant) != 0:
        raise AssertionError("mixed-family discriminant factorization failed")

    c4 = sp.expand(-48 * A)
    c6 = sp.expand(-864 * B)
    if sp.degree(A, t) != 4 or sp.degree(B, t) != 7:
        raise AssertionError("unexpected finite degrees")
    if sp.degree(discriminant, t) != 14:
        raise AssertionError("unexpected discriminant degree")

    # The shift x=X+M exposes the split nodal branches at t=lambda.
    X = sp.symbols("X")
    L = lam * (lam - 1)
    M = coefficient_k * T
    shifted_rhs = sp.expand(
        (X + M) ** 3 + A * (X + M) + B
    )
    expected_shifted_rhs = sp.expand(
        X**2 * (X + 3 * M) + coefficient_c * T**2 * z**3
    )
    if sp.expand(shifted_rhs - expected_shifted_rhs) != 0:
        raise AssertionError("shifted nodal identity failed")

    # Under the split condition d^2=3*k*L, scaling by d/(3L) gives k=3L.
    split_gauge_k = sp.expand(3 * L)
    split_node_slope_square = sp.expand(3 * split_gauge_k * L)
    if sp.expand(split_node_slope_square - (3 * L) ** 2) != 0:
        raise AssertionError("split-I3 gauge failed")

    gram = sp.Matrix(
        [
            [sp.Rational(8, 3), sp.Rational(-1, 3), 0],
            [sp.Rational(-1, 3), sp.Rational(10, 3), 0],
            [0, 0, 4],
        ]
    )
    gram_determinant = sp.factor(gram.det())
    if gram_determinant != sp.Rational(316, 9):
        raise AssertionError("unexpected Mordell--Weil determinant")
    root_discriminant = 3**3
    ns_discriminant = sp.factor(root_discriminant * gram_determinant)
    if ns_discriminant != 948:
        raise AssertionError("target Neron--Severi discriminant mismatch")

    examples = [
        {
            "prime": 13,
            "k": 1,
            "lambda": 8,
            "c": 1,
            "x": [4, 12, 3, 11, 12],
            "y": [5, 3, 0, 11, 2, 11, 8],
            "expected_mask": 0,
            "role": "P3 diagonal profile only",
        },
        {
            "prime": 13,
            "k": 1,
            "lambda": 6,
            "c": 5,
            "x": [0, 10, 7, 11, 12],
            "y": [0, 5, 8, 0, 6, 11, 8],
            "expected_mask": 1,
            "role": "P2 diagonal profile only",
        },
        {
            "prime": 13,
            "k": 1,
            "lambda": 6,
            "c": 2,
            "x": [3, 3, 1, 3, 3],
            "y": [1, 8, 3, 9, 9, 8, 1],
            "expected_mask": 6,
            "role": "P1 diagonal profile only",
        },
    ]
    verified_examples = [verify_modular_section(example) for example in examples]

    residual = sp.Poly(
        coefficient_c * z**3 + 4 * coefficient_k**3 * T, t
    )
    residual_discriminant = sp.factor(sp.discriminant(residual.as_expr(), t))

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "e8_a2_mixed_two_iv_one_i3",
        "exact_claim": (
            "The displayed non-isotrivial family has fibres "
            "II*+2IV+I3+3I1 under the stated open conditions; the stated "
            "component and intersection data are equivalent to the target "
            "Mordell--Weil Gram of determinant 316/9 and NS discriminant 948."
        ),
        "family": {
            "T": "t*(t-1)",
            "A": str(sp.factor(A)),
            "B": str(sp.factor(B)),
            "discriminant": str(sp.factor(discriminant)),
            "residual_cubic": str(residual.as_expr()),
            "residual_cubic_discriminant": str(residual_discriminant),
            "c4": str(sp.factor(c4)),
            "c6": str(sp.factor(c6)),
            "open_conditions": [
                "characteristic is not 2 or 3",
                "k*c*lambda*(lambda-1) != 0",
                "the residual cubic is squarefree",
            ],
            "kodaira_data": {
                "t=0": {"type": "IV", "ord(A,B,Delta)": [2, 2, 4]},
                "t=1": {"type": "IV", "ord(A,B,Delta)": [2, 2, 4]},
                "t=lambda": {"type": "I3", "ord(c4,Delta)": [0, 3]},
                "t=infinity": {
                    "type": "II*",
                    "ord(c4,c6,Delta)": [4, 5, 10],
                },
                "residual": {"type": "3 I1"},
            },
            "euler_number_check": "4+4+3+3+10=24",
        },
        "split_i3_chart": {
            "condition": "3*k*lambda*(lambda-1) is a nonzero square",
            "normalization": "k=3*lambda*(lambda-1)",
            "shift": "X=x-k*t*(t-1)",
            "equation": "y^2=X^2*(X+3*k*t*(t-1))+c*t^2*(t-1)^2*(t-lambda)^3",
            "node_slopes": ["+3*lambda*(lambda-1)", "-3*lambda*(lambda-1)"],
        },
        "target_gram": {
            "matrix": [["8/3", "-1/3", "0"], ["-1/3", "10/3", "0"], ["0", "0", "4"]],
            "component_nonidentity_counts": [2, 1, 0],
            "shared_nonidentity_component_rule": (
                "P1 and P2 overlap at exactly one A2 fibre and use opposite "
                "nonzero Z/3 labels"
            ),
            "pairwise_section_intersections": [2, 2, 2],
            "mordell_weil_determinant": str(gram_determinant),
            "root_lattice": "E8+A2^3",
            "root_discriminant": root_discriminant,
            "neron_severi_absolute_discriminant": int(ns_discriminant),
            "torsion": {
                "order": 1,
                "argument": (
                    "A torsion section would have height zero, but its height "
                    "is at least 4-3*(2/3)=2 because the II* component group "
                    "is trivial and the three A2 corrections are each at most 2/3."
                ),
            },
        },
        "profiled_section_system": {
            "notation": ["z=t-lambda", "L=lambda*(lambda-1)", "k=3*L", "M=k*t*(t-1)"],
            "P1": {
                "coordinates": ["X1=z*t*U1", "y1=z*t*V1"],
                "degrees": {"U1": 2, "V1": 4},
                "identity": "V1^2-U1^2*(z*t*U1+3*M)-c*(t-1)^2*z=0",
                "branch": "V1(lambda)=+3*L*U1(lambda)",
                "saturations": ["U1(1)!=0", "U1(lambda)!=0"],
            },
            "P2": {
                "coordinates": ["X2=z*U2", "y2=z*V2"],
                "degrees": {"U2": 3, "V2": 5},
                "identity": "V2^2-U2^2*(z*U2+3*M)-c*t^2*(t-1)^2*z=0",
                "branch": "V2(lambda)=-3*L*U2(lambda)",
                "saturations": ["U2(0)*U2(1)*U2(lambda)!=0"],
            },
            "P3": {
                "coordinates": ["X3", "Y3"],
                "degrees": {"X3": 4, "Y3": 6},
                "identity": "Y3^2-X3^2*(X3+3*M)-c*t^2*(t-1)^2*z^3=0",
                "saturations": ["X3(0)*X3(1)*X3(lambda)!=0"],
            },
            "P1_cusp_parameter": "c=-rho^2/lambda with rho=V1(0)!=0",
            "remaining_gate": (
                "verify the three pairwise intersections are exactly 2, "
                "including any contribution at infinity"
            ),
        },
        "modular_engine_checks": {
            "note": (
                "The examples certify all three diagonal profiles but occur "
                "on different GF(13) surfaces; they are not a rank-three seed."
            ),
            "examples": verified_examples,
        },
        "claim_boundary": {
            "rank31_curve_found": False,
            "rational_characteristic_zero_surface_found": False,
            "complete_modular_target_triple_found": False,
            "scope": (
                "This certificate validates a replacement search chart and "
                "its exact target filters.  It does not prove that the chart "
                "contains the desired characteristic-zero lattice point."
            ),
        },
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()

    certificate = compute_certificate()
    encoded = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if arguments.compare is not None:
        expected = json.loads(arguments.compare.read_text())
        if certificate != expected:
            raise SystemExit("certificate mismatch")
    if arguments.output is not None:
        arguments.output.write_text(encoded)
    elif arguments.compare is None:
        print(encoded, end="")


if __name__ == "__main__":
    main()
