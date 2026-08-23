#!/usr/bin/env python3
"""Derive the full marked split-II* + 3 I3 search chart.

The earlier two-IV/one-I3 chart is a useful small experiment but is a proper
subfamily.  This module records a four-dimensional chart for the full
E8+A2^3 lattice-polarized moduli problem with all three A2 fibres semistable.
It is the appropriate ambient family before imposing the three target
Mordell--Weil sections.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


def canonical_sha256(payload: object) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def compute_certificate() -> dict[str, object]:
    t, lam = sp.symbols("t lambda")
    a0, a1, a2 = sp.symbols("a0 a1 a2")
    beta0, beta1 = sp.symbols("beta0 beta1")
    gamma0, gamma1 = sp.symbols("gamma0 gamma1")
    quotient_d = sp.symbols("d")

    D = t * (t - 1) * (t - lam)
    a = a0 + a1 * t + a2 * t**2
    beta = beta0 + beta1 * t
    gamma = gamma0 + gamma1 * t
    b = D * beta
    c_error = D**2 * gamma

    A = sp.expand(-3 * (a**2 + 2 * b))
    B = sp.expand(2 * (a**3 + 3 * a * b) + c_error)
    discriminant = sp.expand(-16 * (4 * A**3 + 27 * B**2))

    relation = sp.expand(a * gamma - 3 * beta**2 - quotient_d * D)
    quotient_Q = sp.expand(
        4 * a**2 * quotient_d
        + 12 * a * beta * gamma
        - 32 * beta**3
        + D * gamma**2
    )
    expected_on_relation = sp.expand(-432 * D**3 * quotient_Q)
    error_identity = sp.factor(discriminant - expected_on_relation)
    expected_error = sp.factor(-1728 * a**2 * D**2 * relation)
    if sp.expand(error_identity - expected_error) != 0:
        raise AssertionError("semistable discriminant identity failed")

    relation_poly = sp.Poly(relation, t)
    if relation_poly.degree() > 3:
        raise AssertionError("relation unexpectedly exceeds degree three")
    relation_equations = [
        str(sp.factor(relation_poly.nth(index))) for index in range(4)
    ]

    X = sp.symbols("X")
    shifted_rhs = sp.expand((X + a) ** 3 + A * (X + a) + B)
    expected_shifted_rhs = sp.expand(
        X**2 * (X + 3 * a) - 6 * D * beta * X + D**2 * gamma
    )
    if sp.expand(shifted_rhs - expected_shifted_rhs) != 0:
        raise AssertionError("shifted full-family identity failed")

    # Prove that the mixed family is the beta=0 specialization.
    coefficient_k, coefficient_c = sp.symbols("k c")
    T = t * (t - 1)
    z = t - lam
    mixed_substitution = {
        a0: 0,
        a1: -coefficient_k,
        a2: coefficient_k,
        beta0: 0,
        beta1: 0,
        gamma0: -coefficient_c * lam,
        gamma1: coefficient_c,
        quotient_d: coefficient_k * coefficient_c,
    }
    if sp.expand(relation.subs(mixed_substitution)) != 0:
        raise AssertionError("mixed subfamily does not satisfy the I3 relation")
    mixed_A = sp.factor(A.subs(mixed_substitution))
    mixed_B = sp.factor(B.subs(mixed_substitution))
    expected_mixed_A = -3 * coefficient_k**2 * T**2
    expected_mixed_B = T**2 * (
        coefficient_c * z**3 + 2 * coefficient_k**3 * T
    )
    if sp.expand(mixed_A - expected_mixed_A) != 0:
        raise AssertionError("mixed A specialization failed")
    if sp.expand(mixed_B - expected_mixed_B) != 0:
        raise AssertionError("mixed B specialization failed")

    # Section identities after extracting the forced nodal factors.
    U, V = sp.symbols("U V")
    S1 = t * z
    S2 = z
    p1_reduced_rhs = sp.expand(
        S1 * U**3
        + 3 * a * U**2
        - 6 * (t - 1) * beta * U
        + (t - 1) ** 2 * gamma
    )
    p2_reduced_rhs = sp.expand(
        S2 * U**3
        + 3 * a * U**2
        - 6 * T * beta * U
        + T**2 * gamma
    )

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "e8_a2_full_semistable_three_i3_chart",
        "exact_claim": (
            "On the marked-node cover, the displayed four-dimensional chart "
            "parametrizes short Weierstrass K3 models with II* at infinity "
            "and three prescribed I3 fibres, subject to explicit open "
            "conditions.  The previous mixed chart is its beta=0 subfamily."
        ),
        "notation": {
            "D": "t*(t-1)*(t-lambda)",
            "a": "a0+a1*t+a2*t^2",
            "beta": "beta0+beta1*t",
            "gamma": "gamma0+gamma1*t",
            "b": "D*beta",
            "c_error": "D^2*gamma",
        },
        "weierstrass_model": {
            "A": "-3*(a^2+2*D*beta)",
            "B": "2*(a^3+3*a*D*beta)+D^2*gamma",
            "equation": "y^2=x^3+A*x+B",
            "shifted_coordinate": "X=x-a",
            "shifted_equation": (
                "y^2=X^2*(X+3*a)-6*D*beta*X+D^2*gamma"
            ),
        },
        "three_i3_condition": {
            "polynomial_identity": "a*gamma-3*beta^2=d*D",
            "coefficient_equations_low_to_high": relation_equations,
            "discriminant_on_identity": "-432*D^3*Q",
            "Q": "4*a^2*d+12*a*beta*gamma-32*beta^3+D*gamma^2",
            "unreduced_error_identity": str(error_identity),
        },
        "kodaira_open_conditions": {
            "characteristic": "not 2 or 3",
            "marked_fibres": [
                "lambda not in {0,1}",
                "a(0)*a(1)*a(lambda) != 0",
                "Q(0)*Q(1)*Q(lambda) != 0",
            ],
            "residual_fibres": "Q has degree 5, is squarefree, and is coprime to D",
            "infinity": "gamma1 != 0 gives ord_infinity(c4,c6,Delta)=(4,5,10), hence II*",
            "configuration": "II* + I3 + I3 + I3 + 5 I1",
            "euler_number_check": "10+3+3+3+5=24",
        },
        "split_cover": {
            "condition": (
                "3*a(0), 3*a(1), and 3*a(lambda) are nonzero squares"
            ),
            "reason": (
                "At a root r of D the nodal cubic is "
                "y^2=(x-a(r))^2*(x+2*a(r)); its tangent slopes are "
                "+/-sqrt(3*a(r))."
            ),
        },
        "parameter_count": {
            "raw_variables": 9,
            "variables": [
                "a0",
                "a1",
                "a2",
                "beta0",
                "beta1",
                "gamma0",
                "gamma1",
                "d",
                "lambda",
            ],
            "relation_coefficient_equations": 4,
            "weierstrass_homothety_dimension": 1,
            "moduli_dimension": 4,
            "interpretation": (
                "This matches 20-rank(U+E8+A2^3)=4; unlike the mixed "
                "two-dimensional chart, it is not dimensionally forced to "
                "miss a general discriminant-948 Noether--Lefschetz curve."
            ),
        },
        "target_section_chart": {
            "profile": (
                "P1 nonidentity at 0 and lambda; P2 nonidentity only at "
                "lambda on the opposite branch; P3 identity at all three I3 fibres"
            ),
            "P1": {
                "coordinates": ["X1=t*(t-lambda)*U1", "y1=t*(t-lambda)*V1"],
                "degrees": {"U1": 2, "V1": 4},
                "reduced_identity": f"V1^2={str(p1_reduced_rhs)}",
            },
            "P2": {
                "coordinates": ["X2=(t-lambda)*U2", "y2=(t-lambda)*V2"],
                "degrees": {"U2": 3, "V2": 5},
                "reduced_identity": f"V2^2={str(p2_reduced_rhs)}",
            },
            "P3": {
                "coordinates": ["X3", "Y3"],
                "degrees": {"X3": 4, "Y3": 6},
                "identity": (
                    "Y3^2=X3^2*(X3+3*a)-6*D*beta*X3+D^2*gamma"
                ),
            },
            "remaining_exact_filters": [
                "chosen nonzero component labels and opposite overlap branch",
                "all three pairwise section intersections equal 2",
                "open-condition saturations",
            ],
        },
        "mixed_subfamily": {
            "specialization": [
                "beta=0",
                "a=k*t*(t-1)",
                "gamma=c*(t-lambda)",
                "d=k*c",
            ],
            "resulting_A": str(mixed_A),
            "resulting_B": str(mixed_B),
            "dimension": 2,
            "ambient_position": (
                "boundary beta=0 with gcd(a,D)!=1 at t=0,1; there the "
                "corresponding I3+I1 fibres specialize to IV"
            ),
            "consequence": (
                "Failure to find a target triple in the mixed chart is not a "
                "global obstruction; the full four-dimensional chart remains open."
            ),
        },
        "claim_boundary": {
            "full_target_triple_found": False,
            "characteristic_zero_rank17_fibration_reconstructed": False,
            "rank31_curve_found": False,
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
        expected = json.loads(arguments.compare.read_text(encoding="utf-8"))
        if certificate != expected:
            raise SystemExit("certificate mismatch")
    if arguments.output is not None:
        arguments.output.write_text(encoded, encoding="utf-8")
    elif arguments.compare is None:
        print(encoded, end="")


if __name__ == "__main__":
    main()
