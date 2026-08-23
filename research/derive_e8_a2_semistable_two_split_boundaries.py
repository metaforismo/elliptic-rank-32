#!/usr/bin/env python3
"""Certify the three boundary charts missing from the dense two-split chart.

Together with ``derive_e8_a2_semistable_two_split_chart.py``, these charts
cover the marked locus relevant to the target component profiles before the
Kodaira open conditions are imposed.  The three cases are s_0=0,
s_lambda=0, and beta(1)=gamma(1)=0.
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


def interpolate_quadratic(
    t: sp.Symbol,
    lam: sp.Expr,
    value_0: sp.Expr,
    value_1: sp.Expr,
    value_lam: sp.Expr,
) -> tuple[sp.Expr, sp.Expr]:
    coefficient_2 = sp.factor(
        (value_lam - value_0 - lam * (value_1 - value_0))
        / (lam * (lam - 1))
    )
    polynomial = sp.cancel(
        value_0
        + (value_1 - value_0 - coefficient_2) * t
        + coefficient_2 * t**2
    )
    return polynomial, coefficient_2


def certify_chart(
    *,
    name: str,
    parameters: list[str],
    t: sp.Symbol,
    lam: sp.Expr,
    a: sp.Expr,
    beta: sp.Expr,
    gamma: sp.Expr,
    expected_a2: sp.Expr,
    expected_d: sp.Expr,
    node_values: dict[str, list[str]],
    open_conditions: list[str],
) -> dict[str, object]:
    D = sp.expand(t * (t - 1) * (t - lam))
    a_poly = sp.Poly(a, t)
    beta_poly = sp.Poly(beta, t)
    gamma_poly = sp.Poly(gamma, t)
    if a_poly.degree() > 2 or beta_poly.degree() > 1 or gamma_poly.degree() > 1:
        raise AssertionError(f"degree bounds failed in {name}")
    if not is_zero(a_poly.nth(2) - expected_a2):
        raise AssertionError(f"a2 formula failed in {name}")
    relation = sp.factor(a * gamma - 3 * beta**2)
    if not is_zero(relation - expected_d * D):
        raise AssertionError(f"semistable relation failed in {name}")
    gamma1 = sp.factor(gamma_poly.nth(1))
    if gamma1 == 0:
        raise AssertionError(f"gamma1 vanishes identically in {name}")
    H = (
        4 * a**2 * expected_d
        + 12 * a * beta * gamma
        - 32 * beta**3
        + D * gamma**2
    )
    H_poly = sp.Poly(H, t)
    if H_poly.degree() != 5 or not is_zero(H_poly.LC() - gamma1**2):
        raise AssertionError(f"residual quintic failed in {name}")
    return {
        "name": name,
        "parameters": parameters,
        "lambda": str(lam),
        "a": str(sp.factor(a)),
        "beta": str(sp.factor(beta)),
        "gamma": str(sp.factor(gamma)),
        "a2": str(sp.factor(expected_a2)),
        "d": str(sp.factor(expected_d)),
        "relation": "a*gamma-3*beta^2=d*t*(t-1)*(t-lambda)",
        "H_degree": H_poly.degree(),
        "H_leading_coefficient": str(sp.factor(gamma1**2)),
        "node_values": node_values,
        "open_conditions": open_conditions,
        "remaining_kodaira_filter": "H is squarefree and coprime to D",
    }


def compute_certificate() -> dict[str, object]:
    t, lam, q, S, W, h = sp.symbols("t lambda q S W h")

    # Boundary A: s_0=0.  The nonzero value s_lambda=S determines the
    # linear beta and gamma, and the relation fixes a(1).
    a_a, a2_a = interpolate_quadratic(
        t, lam, sp.Rational(1, 3), q**2 / (3 * lam), q**2 / 3
    )
    beta_a = t * q * S / (3 * lam)
    gamma_a = t * S**2 / lam
    d_a = S**2 / (3 * lam**2)

    # Boundary B: s_lambda=0.  Linear interpolation forces beta(lambda)
    # and gamma(lambda) to vanish while a(lambda) remains free and split.
    a_b, a2_b = interpolate_quadratic(
        t,
        lam,
        sp.Rational(1, 3),
        (lam - 1) / (3 * lam),
        q**2 / 3,
    )
    beta_b = W * (1 - t / lam) / 3
    gamma_b = W**2 * (1 - t / lam)
    d_b = -W**2 * q**2 / (3 * lam**2 * (lam - 1))

    # Boundary C: beta(1)=gamma(1)=0 with s_0,s_lambda nonzero.  Splitting
    # at lambda forces 1-lambda=q^2, so lambda=1-q^2.
    lam_c = 1 - q**2
    a_c, a2_c = interpolate_quadratic(
        t, lam_c, sp.Rational(1, 3), h / 3, q**2 / 3
    )
    beta_c = W * (1 - t) / 3
    gamma_c = W**2 * (1 - t)
    d_c = -W**2 * h / (3 * q**2)

    charts = [
        certify_chart(
            name="s0_zero",
            parameters=["lambda", "q", "S"],
            t=t,
            lam=lam,
            a=a_a,
            beta=beta_a,
            gamma=gamma_a,
            expected_a2=sp.Rational(1, 3) / lam,
            expected_d=d_a,
            node_values={
                "0": ["3a=1", "3beta=0", "gamma=0", "r=1", "s=0"],
                "1": ["3a=q^2/lambda", "3beta=qS/lambda", "gamma=S^2/lambda"],
                "lambda": ["3a=q^2", "3beta=qS", "gamma=S^2", "r=q", "s=S"],
            },
            open_conditions=[
                "lambda*(lambda-1)*q*S != 0",
                "H is squarefree and H(0)*H(1)*H(lambda) != 0",
            ],
        ),
        certify_chart(
            name="s_lambda_zero",
            parameters=["lambda", "q", "W"],
            t=t,
            lam=lam,
            a=a_b,
            beta=beta_b,
            gamma=gamma_b,
            expected_a2=q**2 / (3 * lam * (lam - 1)),
            expected_d=d_b,
            node_values={
                "0": ["3a=1", "3beta=W", "gamma=W^2", "r=1", "s=W"],
                "1": ["3a=(lambda-1)/lambda", "3beta=W*(lambda-1)/lambda", "gamma=W^2*(lambda-1)/lambda"],
                "lambda": ["3a=q^2", "3beta=0", "gamma=0", "r=q", "s=0"],
            },
            open_conditions=[
                "lambda*(lambda-1)*q*W != 0",
                "H is squarefree and H(0)*H(1)*H(lambda) != 0",
            ],
        ),
        certify_chart(
            name="node1_beta_gamma_zero",
            parameters=["q", "W", "h"],
            t=t,
            lam=lam_c,
            a=a_c,
            beta=beta_c,
            gamma=gamma_c,
            expected_a2=h / (3 * q**2),
            expected_d=d_c,
            node_values={
                "0": ["3a=1", "3beta=W", "gamma=W^2", "r=1", "s=W"],
                "1": ["3a=h", "3beta=0", "gamma=0"],
                "lambda=1-q^2": ["3a=q^2", "3beta=q^2*W", "gamma=q^2*W^2", "r=q", "s=qW"],
            },
            open_conditions=[
                "q*W*h*(q^2-1) != 0",
                "H is squarefree and H(0)*H(1)*H(lambda) != 0",
            ],
        ),
    ]

    if [sp.factor(a2_a), sp.factor(a2_b), sp.factor(a2_c)] != [
        sp.Rational(1, 3) / lam,
        q**2 / (3 * lam * (lam - 1)),
        h / (3 * q**2),
    ]:
        raise AssertionError("compact boundary a2 formulas failed")

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "e8_a2_semistable_two_split_boundary_charts",
        "exact_claim": (
            "The three displayed three-dimensional charts exactly cover the "
            "intrinsic divisors omitted by the dense two-split chart, before "
            "imposing the residual-quintic Kodaira open conditions."
        ),
        "ambient_model": {
            "D": "t*(t-1)*(t-lambda)",
            "relation": "a*gamma-3*beta^2=d*D",
            "discriminant": "-432*D^3*H",
            "H": "4*a^2*d+12*a*beta*gamma-32*beta^3+D*gamma^2",
        },
        "charts": charts,
        "coverage_argument": [
            "If s0=0, gamma is nonzero linear and s_lambda cannot also vanish; use s0_zero.",
            "If s0!=0 but s_lambda=0, use s_lambda_zero.",
            "If both are nonzero but gamma(1)=0, the node relation gives beta(1)=0 and use node1_beta_gamma_zero.",
            "Otherwise W, y, g and z are defined and the dense two-split chart applies.",
        ],
        "arithmetic_scope": (
            "After relabelling, the target component orbit forces splitting "
            "only at 0 and lambda. The dense chart plus these boundaries "
            "therefore covers the marked two-split target locus."
        ),
        "claim_boundary": {
            "target_sections_found": False,
            "finite_field_boundaries_searched": False,
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
