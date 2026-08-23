#!/usr/bin/env python3
"""Certify the four-parameter chart forced by the target section profiles.

For the target Mordell--Weil Gram, P1 meets nonidentity components at two I3
fibres and P2 meets one of those two.  Those two fibres must be split over the
ground field.  The third I3 fibre need not be split because every target
section meets its identity component.  This module therefore parametrizes
the necessary two-split chart, which is arithmetically broader than the
all-three-split chart.
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


def is_zero(expression: sp.Expr) -> bool:
    return sp.factor(sp.cancel(expression)) == 0


def compute_certificate() -> dict[str, object]:
    t, lam, y_ratio, z_ratio, W = sp.symbols("t lambda y z W")
    g_ratio = sp.cancel((y_ratio**2 + lam - 1) / lam)
    q_ratio = sp.cancel((1 - lam + lam * z_ratio) / y_ratio)

    a2 = sp.factor(
        (-lam * z_ratio + lam + y_ratio**2 - 1) ** 2
        / (3 * lam * y_ratio**2 * (lam + y_ratio**2 - 1))
    )
    a_at_1 = sp.cancel(z_ratio**2 / (3 * g_ratio))
    a = sp.cancel(
        sp.Rational(1, 3)
        + (a_at_1 - sp.Rational(1, 3) - a2) * t
        + a2 * t**2
    )
    beta = sp.cancel(W * (1 + (z_ratio - 1) * t) / 3)
    gamma = sp.cancel(W**2 * (1 + (g_ratio - 1) * t))
    D = t * (t - 1) * (t - lam)
    gamma1 = sp.factor(sp.diff(gamma, t))
    quotient_d = sp.factor(a2 * gamma1)

    interpolation_checks = {
        "3a(0)-1": 3 * a.subs(t, 0) - 1,
        "3a(1)-z^2/g": 3 * a.subs(t, 1) - z_ratio**2 / g_ratio,
        "3a(lambda)-q^2": 3 * a.subs(t, lam) - q_ratio**2,
        "3beta(0)-W": 3 * beta.subs(t, 0) - W,
        "3beta(1)-zW": 3 * beta.subs(t, 1) - z_ratio * W,
        "3beta(lambda)-qyW": (
            3 * beta.subs(t, lam) - q_ratio * y_ratio * W
        ),
        "gamma(0)-W^2": gamma.subs(t, 0) - W**2,
        "gamma(1)-gW^2": gamma.subs(t, 1) - g_ratio * W**2,
        "gamma(lambda)-y^2W^2": (
            gamma.subs(t, lam) - y_ratio**2 * W**2
        ),
    }
    failures = [
        name for name, residual in interpolation_checks.items()
        if not is_zero(residual)
    ]
    if failures:
        raise AssertionError(f"interpolation failed: {failures}")

    relation = sp.factor(a * gamma - 3 * beta**2)
    if not is_zero(relation - quotient_d * D):
        raise AssertionError("a*gamma-3*beta^2=dD failed")
    if t in quotient_d.free_symbols:
        raise AssertionError("d unexpectedly depends on t")

    expected_gamma1 = sp.factor(W**2 * (y_ratio**2 - 1) / lam)
    if not is_zero(gamma1 - expected_gamma1):
        raise AssertionError("gamma1 formula failed")

    residual_H = (
        4 * a**2 * quotient_d
        + 12 * a * beta * gamma
        - 32 * beta**3
        + D * gamma**2
    )
    residual_poly = sp.Poly(residual_H, t)
    if residual_poly.degree() != 5:
        raise AssertionError("H does not have degree five")
    if not is_zero(residual_poly.LC() - gamma1**2):
        raise AssertionError("leading coefficient of H is not gamma1^2")

    H_at_nodes = {
        "0": sp.factor(residual_H.subs(t, 0)),
        "1": sp.factor(residual_H.subs(t, 1)),
        "lambda": sp.factor(residual_H.subs(t, lam)),
    }
    expected_H_at_nodes = {
        "0": sp.factor(sp.Rational(4, 27) * (3 * quotient_d + W**3)),
        "1": sp.factor(
            4
            * z_ratio**3
            * (3 * z_ratio * quotient_d + g_ratio**2 * W**3)
            / (27 * g_ratio**2)
        ),
        "lambda": sp.factor(
            sp.Rational(4, 27)
            * q_ratio**3
            * (3 * q_ratio * quotient_d + y_ratio**3 * W**3)
        ),
    }
    for node in H_at_nodes:
        if not is_zero(H_at_nodes[node] - expected_H_at_nodes[node]):
            raise AssertionError(f"H({node}) formula failed")

    # Universal identities, checked without expanding the parameter chart.
    X, av, bv, gv, dv, Dv = sp.symbols("X av bv gv dv Dv")
    universal_A = -3 * (av**2 + 2 * Dv * bv)
    universal_B = 2 * (av**3 + 3 * av * Dv * bv) + Dv**2 * gv
    shifted_error = sp.expand(
        (X + av) ** 3
        + universal_A * (X + av)
        + universal_B
        - (X**2 * (X + 3 * av) - 6 * Dv * bv * X + Dv**2 * gv)
    )
    if shifted_error != 0:
        raise AssertionError("universal shifted equation failed")
    universal_H = (
        4 * av**2 * dv
        + 12 * av * bv * gv
        - 32 * bv**3
        + Dv * gv**2
    )
    universal_relation = av * gv - 3 * bv**2 - dv * Dv
    universal_discriminant = -16 * (
        4 * universal_A**3 + 27 * universal_B**2
    )
    if sp.expand(
        universal_discriminant
        + 432 * Dv**3 * universal_H
        + 1728 * av**2 * Dv**2 * universal_relation
    ) != 0:
        raise AssertionError("universal discriminant identity failed")

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "e8_a2_semistable_two_split_rational_chart",
        "exact_claim": (
            "The formulas give a rational four-parameter dense chart of the "
            "marked semistable family in which precisely the two I3 fibres "
            "forced to split by the target component profiles are split."
        ),
        "parameters": ["lambda", "y", "z", "W"],
        "auxiliary_ratios": {
            "g=gamma(1)/gamma(0)": str(g_ratio),
            "y=s_lambda/s0": "y",
            "z=3*beta(1)/W": "z",
            "q=r_lambda/r0": str(q_ratio),
            "interpolation_equations": [
                "y^2=(1-lambda)+lambda*g",
                "q*y=(1-lambda)+lambda*z",
            ],
        },
        "gauge_and_node_values": {
            "homothety_gauge": "r0=sqrt(3*a(0))=1",
            "at_0": ["3a=1", "3beta=W", "gamma=W^2"],
            "at_1": ["3a=z^2/g", "3beta=zW", "gamma=gW^2"],
            "at_lambda": [
                "3a=q^2",
                "3beta=qyW",
                "gamma=y^2W^2",
            ],
            "third_fibre_split_iff": "g is a nonzero square",
        },
        "polynomials": {
            "D": "t*(t-1)*(t-lambda)",
            "a2": str(a2),
            "a": "1/3+(z^2/(3g)-1/3-a2)*t+a2*t^2",
            "beta": "W*(1+(z-1)*t)/3",
            "gamma": "W^2*(1+(g-1)*t)",
            "gamma1": str(expected_gamma1),
            "d": str(quotient_d),
            "relation": "a*gamma-3*beta^2=d*D",
            "H": "4*a^2*d+12*a*beta*gamma-32*beta^3+D*gamma^2",
            "H_degree": residual_poly.degree(),
        },
        "weierstrass_model": {
            "A": "-3*(a^2+2*D*beta)",
            "B": "2*(a^3+3*a*D*beta)+D^2*gamma",
            "shifted_equation": (
                "Y^2=X^2*(X+3*a)-6*D*beta*X+D^2*gamma"
            ),
            "discriminant": "-432*D^3*H",
        },
        "node_coprimality_factors": {
            "H(0)_nonzero_iff": "3*d+W^3 != 0",
            "H(1)_nonzero_iff": "3*z*d+g^2*W^3 != 0",
            "H(lambda)_nonzero_iff": "3*q*d+y^3*W^3 != 0",
        },
        "birational_coverage": {
            "target_forces_split_at": ["0", "lambda"],
            "target_does_not_force_split_at": ["1"],
            "inverse_on_dense_open": [
                "use the homothety to set r0=1",
                "W=3*beta(0)",
                "g=gamma(1)/W^2",
                "z=3*beta(1)/W",
                "y=s_lambda/W",
                "q=r_lambda",
            ],
            "relation_at_1": "a(1)=3*beta(1)^2/gamma(1)",
            "dimension": 4,
        },
        "open_conditions": [
            "characteristic is not 2 or 3",
            "lambda*(lambda-1) != 0",
            "W*y*z*q*g != 0",
            "y^2 != 1 (equivalently gamma1 != 0)",
            "3*d+W^3 != 0",
            "3*z*d+g^2*W^3 != 0",
            "3*q*d+y^3*W^3 != 0",
            "H is squarefree",
        ],
        "coverage_boundary": (
            "The chart is dense but excludes divisors where the chosen node "
            "coordinates W, gamma(1), or s_lambda vanish, plus the usual "
            "Kodaira collision loci. Separate boundary charts are required "
            "before any global nonexistence conclusion."
        ),
        "comparison_with_all_split_chart": (
            "The earlier all-split chart is recovered when g=x^2 and z=R*x. "
            "Requiring that square condition is unnecessary for the target."
        ),
        "claim_boundary": {
            "target_sections_found": False,
            "characteristic_zero_lift_found": False,
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
