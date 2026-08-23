#!/usr/bin/env python3
"""Certify a rational four-parameter split chart for II* + 3 I3 K3s.

The full marked semistable family is described by

    D = t(t-1)(t-lambda),
    a*gamma - 3*beta^2 = d*D,

where deg(a)<=2 and deg(beta),deg(gamma)<=1.  On the split cover write

    3*a(i) = r_i^2,   3*beta(i) = r_i*s_i,

at i=0,1,lambda.  Then gamma(i)=s_i^2.  The two interpolation
conditions are a conic plus a linear equation, and yield a birational
four-parameter chart after the Weierstrass homothety r_0=1.

This script checks every identity with exact rational-function arithmetic and
emits a deterministic JSON certificate.  It does not assert the existence of
the three extra Mordell--Weil sections or of a rank-31 curve over Q.
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
    t, lam, m, R, W = sp.symbols("t lambda m R W")

    denominator = m**2 - lam
    x_ratio = sp.cancel((m**2 - 2 * m + lam) / denominator)
    y_ratio = sp.cancel((-m**2 + 2 * m * lam - lam) / denominator)
    q_ratio = sp.cancel(
        ((1 - lam) + lam * R * x_ratio) / y_ratio
    )

    # The compact quadratic interpolation formula is equivalent to the
    # three-node Lagrange formula but keeps the exact replay inexpensive.
    a2 = sp.cancel((R - x_ratio) ** 2 / (3 * y_ratio**2))
    a = sp.cancel(
        sp.Rational(1, 3)
        + ((R**2 - 1) / 3 - a2) * t
        + a2 * t**2
    )
    beta = sp.cancel(W * (1 - t + R * x_ratio * t) / 3)
    gamma = sp.cancel(W**2 * (1 - t + x_ratio**2 * t))
    D = t * (t - 1) * (t - lam)

    conic_residual = sp.factor(
        y_ratio**2 - ((1 - lam) + lam * x_ratio**2)
    )
    if conic_residual != 0:
        raise AssertionError("conic parametrization failed")

    interpolation_checks = {
        "3a(0)-1": 3 * a.subs(t, 0) - 1,
        "3a(1)-R^2": 3 * a.subs(t, 1) - R**2,
        "3a(lambda)-Q^2": 3 * a.subs(t, lam) - q_ratio**2,
        "3beta(0)-W": 3 * beta.subs(t, 0) - W,
        "3beta(1)-R*x*W": 3 * beta.subs(t, 1) - R * x_ratio * W,
        "3beta(lambda)-Q*y*W": (
            3 * beta.subs(t, lam) - q_ratio * y_ratio * W
        ),
        "gamma(0)-W^2": gamma.subs(t, 0) - W**2,
        "gamma(1)-x^2*W^2": gamma.subs(t, 1) - x_ratio**2 * W**2,
        "gamma(lambda)-y^2*W^2": (
            gamma.subs(t, lam) - y_ratio**2 * W**2
        ),
    }
    failed = [name for name, value in interpolation_checks.items() if not is_zero(value)]
    if failed:
        raise AssertionError(f"interpolation checks failed: {failed}")

    relation = sp.factor(sp.cancel(a * gamma - 3 * beta**2))
    relation_poly = sp.Poly(relation, t)
    quotient_d = sp.factor(relation_poly.LC())
    if not is_zero(relation - quotient_d * D):
        raise AssertionError("a*gamma-3*beta^2 is not d*D")
    if t in quotient_d.free_symbols:
        raise AssertionError("d unexpectedly depends on t")

    expected_d = sp.factor(a2 * W**2 * (x_ratio**2 - 1))
    if not is_zero(quotient_d - expected_d):
        raise AssertionError("closed formula for d failed")

    gamma1 = sp.factor(sp.diff(gamma, t))
    expected_gamma1 = sp.factor(
        -4 * W**2 * m * (m - lam) * (m - 1) / denominator**2
    )
    if not is_zero(gamma1 - expected_gamma1):
        raise AssertionError("gamma leading coefficient formula failed")

    # Recheck the Weierstrass and discriminant identities universally, before
    # substitution.  Expanding the four-parameter rational functions would be
    # algebraically redundant and dramatically slower.
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
        raise AssertionError("universal shifted Weierstrass identity failed")

    residual_H = (
        4 * a**2 * quotient_d
        + 12 * a * beta * gamma
        - 32 * beta**3
        + D * gamma**2
    )
    residual_poly = sp.Poly(residual_H, t)
    residual_degree = residual_poly.degree()
    if residual_degree != 5 or not is_zero(residual_poly.LC() - gamma1**2):
        raise AssertionError("residual H is not a degree-five polynomial in t")

    universal_H = (
        4 * av**2 * dv
        + 12 * av * bv * gv
        - 32 * bv**3
        + Dv * gv**2
    )
    universal_relation = av * gv - 3 * bv**2 - dv * Dv
    universal_discriminant = -16 * (4 * universal_A**3 + 27 * universal_B**2)
    discriminant_error = sp.expand(
        universal_discriminant
        + 432 * Dv**3 * universal_H
        + 1728 * av**2 * Dv**2 * universal_relation
    )
    if discriminant_error != 0:
        raise AssertionError("universal discriminant identity failed")

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "e8_a2_semistable_split_rational_chart",
        "exact_claim": (
            "The displayed formulas give a rational four-parameter dense "
            "chart of the split marked II*+3I3 semistable family on the "
            "stated open set. They satisfy a*gamma-3*beta^2=dD exactly."
        ),
        "parameters": ["lambda", "m", "R", "W"],
        "auxiliary_ratios": {
            "denominator": "m^2-lambda",
            "x=s1/s0": str(x_ratio),
            "y=s_lambda/s0": str(y_ratio),
            "Q=r_lambda/r0": str(q_ratio),
            "conic": "y^2=(1-lambda)+lambda*x^2",
            "linear_interpolation": (
                "Q*y=(1-lambda)+lambda*R*x"
            ),
        },
        "gauge_and_values": {
            "homothety_gauge": "r0=1",
            "R": "r1/r0",
            "W": "s0/r0^3",
            "at_0": ["3a=1", "3beta=W", "gamma=W^2"],
            "at_1": ["3a=R^2", "3beta=R*x*W", "gamma=x^2*W^2"],
            "at_lambda": [
                "3a=Q^2",
                "3beta=Q*y*W",
                "gamma=y^2*W^2",
            ],
        },
        "polynomials": {
            "D": "t*(t-1)*(t-lambda)",
            "a": (
                "1/3+((R^2-1)/3-a2)*t+a2*t^2"
            ),
            "a2": "(R-x)^2/(3*y^2)",
            "beta": "W*(1-t+R*x*t)/3",
            "gamma": "W^2*(1-t+x^2*t)",
            "d": str(expected_d),
            "gamma1": str(expected_gamma1),
            "relation": "a*gamma-3*beta^2=d*D",
            "residual_H": (
                "4*a^2*d+12*a*beta*gamma-32*beta^3+D*gamma^2"
            ),
            "residual_H_degree": residual_degree,
        },
        "weierstrass_model": {
            "A": "-3*(a^2+2*D*beta)",
            "B": "2*(a^3+3*a*D*beta)+D^2*gamma",
            "shifted_equation": (
                "y^2=X^2*(X+3*a)-6*D*beta*X+D^2*gamma"
            ),
            "discriminant": "-432*D^3*H",
        },
        "birational_derivation": {
            "split_variables": (
                "3a(i)=r_i^2 and 3beta(i)=r_i*s_i at i=0,1,lambda"
            ),
            "forced_values": "gamma(i)=s_i^2",
            "interpolation_equations": [
                "s_lambda^2=(1-lambda)*s0^2+lambda*s1^2",
                "r_lambda*s_lambda=(1-lambda)*r0*s0+lambda*r1*s1",
            ],
            "conic_line": "y=1+m*(x-1) through (x,y)=(1,1)",
            "inverse_on_dense_open": [
                "R=r1/r0",
                "W=s0/r0^3 after the homothety r0=1",
                "x=s1/s0",
                "y=s_lambda/s0",
                "m=(y-1)/(x-1)",
            ],
            "dimension": "7 split variables - 2 equations - 1 homothety = 4",
        },
        "open_conditions": [
            "characteristic is not 2 or 3",
            "lambda*(lambda-1) != 0",
            "m^2-lambda != 0",
            "-m^2+2*m*lambda-lambda != 0",
            "R*Q*W != 0",
            "m*(m-lambda)*(m-1) != 0 (equivalently gamma1 != 0 here)",
            "3*d+W^3 != 0",
            "3*R*d+x^3*W^3 != 0",
            "3*Q*d+y^3*W^3 != 0",
            "H is degree 5, squarefree, and coprime to D",
        ],
        "coverage_boundary": (
            "This is a dense birational chart, not every boundary point of "
            "the split cover: it divides by r0, s0, s_lambda, x-1, the conic "
            "denominator, and the listed open factors."
        ),
        "claim_boundary": {
            "three_target_sections_found": False,
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
