#!/usr/bin/env python3
"""Audit and correct the proposed E8+A2^3 rank-jump section system.

The public ansatz x=q^2+r, y=q^3+s with deg(r),deg(s)<=1 omits two
nonzero high coefficients: t^9=-3*q2^4*r1 and, after r1=0,
t^8=-3*q2^4*r0.  On q2!=0 it therefore collapses to r=0.

The complete degree-balanced ansatz is

    x = q^2+r,
    y = q^3+(3/2)q*r+h,

with deg(q)=2 and deg(r),deg(h)<=1.  Its residual has degree at most seven
identically.  The t^7 and t^6 equations eliminate h1 and mu, leaving six
exact equations in seven variables on q2!=0.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

from replay_e6_mw3_p2_split_cores import canonical_sha256


VARIABLE_NAMES = ("lam", "q0", "q1", "q2", "r0", "r1", "h0")


def is_prime(value):
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, int(value**0.5) + 1))


def symbolic_system():
    t = sp.symbols("t")
    lam, mu, q0, q1, q2, r0, r1, h0, h1 = sp.symbols(
        "lam mu q0 q1 q2 r0 r1 h0 h1"
    )
    q = q0 + q1 * t + q2 * t**2
    r = r0 + r1 * t
    h = h0 + h1 * t
    x = q**2 + r
    y = q**3 + sp.Rational(3, 2) * q * r + h
    fiber = (t * (t - 1) * (t - lam)) ** 2 * (t - mu)
    residual = sp.Poly(sp.expand(y**2 - x**3 - fiber), t)
    if residual.degree() != 7:
        raise AssertionError("corrected residual did not have degree seven")
    h1_formula = 1 / (2 * q2**3)
    mu_formula = (
        -2 * h0 * q2**3
        - 6 * h1_formula * q1 * q2**2
        - 2 * lam
        + sp.Rational(3, 4) * q2**2 * r1**2
        - 2
    )
    substitutions = {h1: h1_formula, mu: mu_formula}
    residual_substituted = sp.Poly(
        sp.together(residual.as_expr().subs(substitutions)), t
    )
    if any(
        sp.cancel(residual_substituted.coeff_monomial(t**degree)) != 0
        for degree in (6, 7)
    ):
        raise AssertionError("triangular h1/mu elimination failed")
    variables = (lam, q0, q1, q2, r0, r1, h0)
    equations = []
    for degree in range(6):
        coefficient = residual_substituted.coeff_monomial(t**degree)
        equations.append(sp.together(coefficient).as_numer_denom()[0])
    legacy_s0, legacy_s1 = sp.symbols("legacy_s0 legacy_s1")
    legacy_y = q**3 + legacy_s0 + legacy_s1 * t
    legacy = sp.Poly(sp.expand(legacy_y**2 - x**3 - fiber), t)
    legacy_high = {
        "t9": sp.factor(legacy.coeff_monomial(t**9)),
        "t8_after_r1_zero": sp.factor(
            legacy.coeff_monomial(t**8).subs({r1: 0})
        ),
    }
    return variables, h1_formula, mu_formula, equations, legacy_high


def derive(prime):
    if not is_prime(prime) or prime in (2, 3, 79):
        raise ValueError("choose a prime different from 2, 3, and 79")
    variables, h1, mu, equations, legacy = symbolic_system()
    return (
        variables,
        h1,
        mu,
        [sp.Poly(equation, variables, modulus=prime) for equation in equations],
        legacy,
    )


def poly_record(poly):
    return {
        "total_degree": poly.total_degree(),
        "terms": [
            {"monomial": list(monomial), "coefficient": int(coefficient)}
            for monomial, coefficient in poly.terms()
        ],
    }


def singular_expression(poly, prime):
    terms = []
    for monomial, coefficient in poly.terms():
        coefficient = int(coefficient) % prime
        factors = []
        if coefficient != 1 or not any(monomial):
            factors.append(str(coefficient))
        for name, exponent in zip(VARIABLE_NAMES, monomial, strict=True):
            if exponent == 1:
                factors.append(name)
            elif exponent:
                factors.append(f"{name}^{exponent}")
        terms.append("*".join(factors) if factors else "0")
    return "+".join(terms) if terms else "0"


def singular_probe(prime, equations, slice_seed):
    coefficients = [
        (slice_seed * (index + 2) ** 2 + 3 * index + 1) % prime
        for index in range(len(VARIABLE_NAMES))
    ]
    constant = (7 * slice_seed + 1) % prime
    linear = "+".join(
        [
            f"{coefficient}*{name}"
            for name, coefficient in zip(VARIABLE_NAMES, coefficients, strict=True)
            if coefficient
        ]
        + [str(-constant)]
    )
    generators = [singular_expression(poly, prime) for poly in equations]
    generators += ["z*q2-1", linear]
    script = "\n".join(
        [
            f"ring R={prime},(lam,q0,q1,q2,r0,r1,h0,z),dp;",
            "option(redSB);",
            "ideal I=" + ",\n".join(generators) + ";",
            "ideal G=std(I);",
            'print("R3CORRECTED|unit="+string(reduce(1,G)==0));',
            'print("R3CORRECTED|dim="+string(dim(G)));',
            'print("R3CORRECTED|size="+string(size(G)));',
            'if (dim(G)==0) { print("R3CORRECTED|vdim="+string(vdim(G))); }',
            "quit;",
            "",
        ]
    )
    return script, coefficients, constant


def compute(prime, slice_seed):
    _variables, h1, mu, equations, legacy = derive(prime)
    script, coefficients, constant = singular_probe(prime, equations, slice_seed)
    payload = {
        "schema": "elliptic-rank30/e8-a2-rank3-jump-corrected/v1",
        "field_prime": prime,
        "variables": list(VARIABLE_NAMES),
        "chart": "q2_nonzero",
        "corrected_ansatz": "x=q^2+r; y=q^3+(3/2)q*r+h; deg(r),deg(h)<=1",
        "legacy_obstruction_over_Q": {
            name: str(expression) for name, expression in legacy.items()
        },
        "eliminated_formulas_over_Q": {"h1": str(h1), "mu": str(mu)},
        "equations": [poly_record(poly) for poly in equations],
        "slice": {"seed": slice_seed, "coefficients": coefficients, "constant": constant},
        "proof": (
            "Direct expansion gives degree at most seven for the corrected "
            "residual. Its t7 and t6 coefficients give the recorded h1 and mu. "
            "The six stored numerators are exactly coefficients t0 through t5."
        ),
        "claim_boundary": (
            "This repairs the one-section equations. The j=0 family still "
            "cannot realize an odd geometric Mordell-Weil rank, so it is not "
            "by itself the target rank-3 K3 family."
        ),
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload, script


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=31)
    parser.add_argument("--slice-seed", type=int, default=1)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--singular-output", type=Path)
    args = parser.parse_args()
    payload, script = compute(args.prime, args.slice_seed)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if args.singular_output:
        args.singular_output.parent.mkdir(parents=True, exist_ok=True)
        args.singular_output.write_text(script)
    print(
        json.dumps(
            {
                "legacy_obstruction": payload["legacy_obstruction_over_Q"],
                "equation_degrees": [e["total_degree"] for e in payload["equations"]],
                "equation_terms": [len(e["terms"]) for e in payload["equations"]],
                "certificate_sha256": payload["certificate_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
