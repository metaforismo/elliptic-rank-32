#!/usr/bin/env python3
"""Exact replay of the 893,730-point E6/MW3 core exhaustion over GF(31).

The public reduction fixes the smooth P1 chart

    (r0, s0, x1) = (4, 18, 27)

and leaves four coordinates (a1,a2,a4,s1), with s1 nonzero.  The two exact
I4 conditions Delta''(1)=Delta'''(1)=0 cut out the finite core set.  This
script derives those equations from the Weierstrass and section identities,
enumerates the whole declared branch, reconstructs every solution, and then
computes the repeated-root polynomial of the residual discriminant.

It is deliberately independent of Sage and msolve.  SymPy is used only once
to derive the two equations; the exhaustive finite-field scan and all
classification arithmetic use the small exact routines in the companion
P2 replay.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from replay_e6_mw3_p2_split_cores import (
    CORES,
    PRIME,
    build_p1_candidate,
    canonical_sha256,
    derive_triangular_expressions,
    evaluate_expression,
    poly_derivative,
    poly_eval,
    poly_pow,
    poly_add,
    poly_scale,
    section_residual,
    trim,
    valuation_at,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERTIFICATE = ROOT / "certificates" / "e6_mw3_core_exhaustion_gf31.json"
UPSTREAM_COMMIT = "e821b0c258bf635f38a1612ff8f2f769bf4ed04b"
UPSTREAM_NOTE = (
    "https://github.com/royvanrijn/jacobian-research/blob/"
    + UPSTREAM_COMMIT
    + "/elkies-k3/E6_P2_REDUCTION_2026-08-20.md"
)


def poly_divmod(left: list[int], right: list[int]) -> tuple[list[int], list[int]]:
    left = trim(left)
    right = trim(right)
    if right == [0]:
        raise ZeroDivisionError("polynomial division by zero")
    if len(left) < len(right):
        return [0], left
    quotient = [0] * (len(left) - len(right) + 1)
    remainder = left[:]
    inverse_leading = pow(right[-1], -1, PRIME)
    while remainder != [0] and len(remainder) >= len(right):
        shift = len(remainder) - len(right)
        coefficient = remainder[-1] * inverse_leading % PRIME
        quotient[shift] = coefficient
        for index, value in enumerate(right):
            remainder[index + shift] = (
                remainder[index + shift] - coefficient * value
            ) % PRIME
        remainder = trim(remainder)
    return trim(quotient), remainder


def poly_exact_divide(left: list[int], right: list[int]) -> list[int]:
    quotient, remainder = poly_divmod(left, right)
    if remainder != [0]:
        raise AssertionError(f"nonexact polynomial division, remainder={remainder}")
    return quotient


def poly_monic(poly: list[int]) -> list[int]:
    poly = trim(poly)
    if poly == [0]:
        return poly
    return poly_scale(poly, pow(poly[-1], -1, PRIME))


def poly_gcd(left: list[int], right: list[int]) -> list[int]:
    left, right = trim(left), trim(right)
    while right != [0]:
        _, remainder = poly_divmod(left, right)
        left, right = right, remainder
    return poly_monic(left)


def compile_polynomial(expression, variables) -> list[tuple[tuple[int, ...], int]]:
    import sympy as sp

    numerator, denominator = sp.cancel(expression, modulus=PRIME).as_numer_denom()
    if int(denominator) % PRIME == 0:
        raise AssertionError("derived equation has zero scalar denominator")
    polynomial = sp.Poly(numerator, *variables, modulus=PRIME)
    return [
        (tuple(int(exponent) for exponent in monomial), int(coefficient) % PRIME)
        for monomial, coefficient in polynomial.terms()
    ]


def compile_rational(expression, variables):
    import sympy as sp

    numerator, denominator = sp.cancel(expression, modulus=PRIME).as_numer_denom()
    return compile_polynomial(numerator, variables), compile_polynomial(
        denominator, variables
    )


def evaluate_rational(compiled, values):
    numerator = evaluate_terms(compiled[0], values)
    denominator = evaluate_terms(compiled[1], values)
    if denominator == 0:
        raise ZeroDivisionError("triangular chart denominator vanished")
    return numerator * pow(denominator, -1, PRIME) % PRIME


def compile_specialized_triangular_expressions(
    expressions, symbol_map, r0_value=4, s0_value=18, x1_value=27
):
    import sympy as sp

    variables = tuple(symbol_map[name] for name in ("a1", "a2", "a4", "s1"))
    fixed_values = {
        symbol_map["s0"]: s0_value,
        symbol_map["x0"]: r0_value**2 - 2 * s0_value,
        symbol_map["x1"]: x1_value,
        symbol_map["y0"]: r0_value * (r0_value**2 - 3 * s0_value),
        symbol_map["y1"]: (
            symbol_map["a1"]
            + 3 * (r0_value**2 - s0_value) * x1_value
        ) / sp.Integer(2 * r0_value),
    }
    composed = {}
    for key in ("a3", "y2", "b5", "b4"):
        expression = expressions[key]
        substitutions = {
            symbol_map[previous]: value for previous, value in composed.items()
        }
        expression = sp.cancel(expression.subs(substitutions).subs(fixed_values), modulus=PRIME)
        unexpected = expression.free_symbols - set(variables)
        if unexpected:
            raise AssertionError(
                f"{key} retained eliminated symbols: {sorted(map(str, unexpected))}"
            )
        composed[key] = expression
    return {key: compile_rational(expression, variables) for key, expression in composed.items()}


POWER_TABLE = [
    [pow(value, exponent, PRIME) for exponent in range(16)]
    for value in range(PRIME)
]


def evaluate_terms(
    terms: list[tuple[tuple[int, ...], int]], values: tuple[int, ...]
) -> int:
    total = 0
    for monomial, coefficient in terms:
        term = coefficient
        for index, exponent in enumerate(monomial):
            term = term * POWER_TABLE[values[index]][exponent] % PRIME
        total += term
    return total % PRIME


def enumerate_cores(core_equations, symbol_map):
    variables = tuple(symbol_map[name] for name in ("a1", "a2", "a4", "s1"))
    d2_terms = compile_polynomial(core_equations["D2_at_1"], variables)
    d3_terms = compile_polynomial(core_equations["D3_at_1"], variables)

    # D2 is only quadratic in a4.  Group its terms by that exponent so the
    # complete scan tests D3 only at D2 roots.
    d2_by_a4: dict[int, list[tuple[tuple[int, ...], int]]] = {0: [], 1: [], 2: []}
    for monomial, coefficient in d2_terms:
        a4_exponent = monomial[2]
        reduced = (monomial[0], monomial[1], monomial[3])
        d2_by_a4[a4_exponent].append((reduced, coefficient))

    solutions = []
    d2_roots_tested_by_d3 = 0
    examined = 0
    for a1 in range(PRIME):
        for a2 in range(PRIME):
            for s1 in range(1, PRIME):
                coefficients = [
                    evaluate_terms(d2_by_a4[degree], (a1, a2, s1))
                    for degree in range(3)
                ]
                for a4 in range(PRIME):
                    examined += 1
                    if sum(coefficients[k] * pow(a4, k, PRIME) for k in range(3)) % PRIME:
                        continue
                    d2_roots_tested_by_d3 += 1
                    core = (a1, a2, a4, s1)
                    if evaluate_terms(d3_terms, core) == 0:
                        solutions.append(core)
    return examined, d2_roots_tested_by_d3, solutions, d2_terms, d3_terms


def reconstruct_surface(core, compiled_expressions, coordinates=(4, 18, 27)):
    r0, s0, x1 = coordinates
    a1, a2, a4, s1 = core
    x0 = (r0**2 - 2 * s0) % PRIME
    y0 = r0 * (r0**2 - 3 * s0) % PRIME
    y1 = (a1 + 3 * (r0**2 - s0) * x1) * pow(2 * r0, -1, PRIME) % PRIME
    values = {}
    for key in ("a3", "y2", "b5", "b4"):
        values[key] = evaluate_rational(compiled_expressions[key], core)
    A, B, X1, Y1 = build_p1_candidate(
        core,
        values["a3"],
        values["y2"],
        values["b5"],
        values["b4"],
        coordinates,
    )
    if section_residual(A, B, X1, Y1) != [0]:
        raise AssertionError(f"section reconstruction failed for core {core}")
    discriminant = poly_scale(
        poly_add(poly_scale(poly_pow(A, 3), 4), poly_scale(poly_pow(B, 2), 27)),
        -16,
    )
    return A, B, X1, Y1, discriminant


def classify_solution(core, compiled_expressions, coordinates=(4, 18, 27)):
    try:
        A, B, X1, Y1, discriminant = reconstruct_surface(
            core, compiled_expressions, coordinates
        )
    except (AssertionError, ZeroDivisionError):
        return {
            "core": list(core),
            "status": "outside_direct_triangular_chart",
        }
    multiplicity_zero = valuation_at(discriminant, 0)
    multiplicity_one = valuation_at(discriminant, 1)
    record = {
        "core": list(core),
        "multiplicity_at_0": multiplicity_zero,
        "multiplicity_at_1": multiplicity_one,
    }
    if multiplicity_zero != 4 or multiplicity_one != 4:
        record["status"] = "boundary_or_wrong_I4_multiplicity"
        return record

    known_i4 = poly_pow([0, 1], 4)
    known_i4 = poly_add([0], known_i4)  # make a fresh normalized list
    from replay_e6_mw3_p2_split_cores import poly_mul

    known_i4 = poly_mul(known_i4, poly_pow([-1, 1], 4))
    residual = poly_exact_divide(discriminant, known_i4)
    repeated = poly_gcd(residual, poly_derivative(residual))
    repeated_squarefree = poly_gcd(repeated, poly_derivative(repeated)) == [1]
    roots = [point for point in range(PRIME) if poly_eval(repeated, point) == 0]
    derivative_a = poly_derivative(A)
    derivative_b = poly_derivative(B)
    node_map = {}
    for point in roots:
        nodes = [
            node
            for node in range(PRIME)
            if poly_eval(A, point) == (-3 * node**2) % PRIME
            and poly_eval(B, point) == (2 * node**3) % PRIME
            and (
                poly_eval(derivative_b, point)
                + node * poly_eval(derivative_a, point)
            )
            % PRIME
            == 0
        ]
        if nodes:
            node_map[str(point)] = nodes
    valid_i2_roots = [
        point
        for point in roots
        if valuation_at(discriminant, point) == 2 and str(point) in node_map
    ]
    record.update(
        {
            "status": "saturated",
            "A": A,
            "B": B,
            "X1": X1,
            "Y1": Y1,
            "residual_discriminant": residual,
            "repeated_root_polynomial": repeated,
            "repeated_root_degree": len(repeated) - 1,
            "repeated_root_polynomial_squarefree": repeated_squarefree,
            "roots_in_GF31": roots,
            "singular_nodes": node_map,
            "valid_I2_roots": valid_i2_roots,
        }
    )
    return record


def compute_certificate(r0=4, s0=18, x1=27):
    coordinates = (r0 % PRIME, s0 % PRIME, x1 % PRIME)
    r0, s0, x1 = coordinates
    expressions, symbol_map, core_equations = derive_triangular_expressions(
        r0, s0, x1
    )
    compiled_expressions = compile_specialized_triangular_expressions(
        expressions, symbol_map, r0, s0, x1
    )
    examined, d2_roots, solutions, d2_terms, d3_terms = enumerate_cores(
        core_equations, symbol_map
    )
    classified = [
        classify_solution(core, compiled_expressions, coordinates)
        for core in solutions
    ]
    gcd_candidates = [
        record
        for record in classified
        if record.get("status") == "saturated"
        and record.get("repeated_root_degree", 0) >= 2
        and record.get("repeated_root_polynomial_squarefree")
    ]
    split_candidates = [
        record
        for record in gcd_candidates
        if len(record.get("roots_in_GF31", []))
        == record.get("repeated_root_degree")
    ]
    declared = {tuple(core) for core in CORES}
    recovered = {tuple(record["core"]) for record in split_candidates}
    y0 = r0 * (r0**2 - 3 * s0) % PRIME
    constant_in_y1 = 3 * (r0**2 - s0) * x1 % PRIME
    degenerate_a1 = (-constant_in_y1 - 2 * r0 * y0) % PRIME
    valid_two_i2_surfaces = [
        record
        for record in classified
        if len(record.get("valid_I2_roots", [])) >= 2
        and record["core"][0] != degenerate_a1
    ]
    payload = {
        "schema": "elliptic-rank30/e6-mw3-core-exhaustion-gf31/v1",
        "field_prime": PRIME,
        "coordinate_slice": {"r0": r0, "s0": s0, "x1": x1},
        "branch_condition": "s1 != 0",
        "upstream_commit": UPSTREAM_COMMIT,
        "upstream_note": UPSTREAM_NOTE,
        "equation_variable_order": ["a1", "a2", "a4", "s1"],
        "D2_terms": [[list(monomial), coefficient] for monomial, coefficient in d2_terms],
        "D3_terms": [[list(monomial), coefficient] for monomial, coefficient in d3_terms],
        "counts": {
            "quadruples_examined": examined,
            "D2_roots_tested_by_D3": d2_roots,
            "simultaneous_core_solutions": len(solutions),
            "squarefree_repeated_gcd_candidates": len(gcd_candidates),
            "fully_split_repeated_gcd_candidates": len(split_candidates),
            "valid_surfaces_with_at_least_two_I2_roots": len(valid_two_i2_surfaces),
            "degenerate_a1_pivot_value": degenerate_a1,
            "solutions_on_degenerate_a1_pivot_plane": sum(
                1 for core in solutions if core[0] == degenerate_a1
            ),
            "direct_chart_exact_I4_surfaces": sum(
                1
                for record in classified
                if record.get("status") == "saturated"
                and record["core"][0] != degenerate_a1
            ),
            "direct_chart_surfaces_with_repeated_factor_degree_at_least_1": sum(
                1
                for record in classified
                if record.get("repeated_root_degree", 0) >= 1
                and record["core"][0] != degenerate_a1
            ),
            "direct_chart_surfaces_with_at_least_one_valid_I2_root": sum(
                1
                for record in classified
                if len(record.get("valid_I2_roots", [])) >= 1
                and record["core"][0] != degenerate_a1
            ),
        },
        "all_core_solutions_sha256": canonical_sha256([list(core) for core in solutions]),
        "gcd_candidates": gcd_candidates,
        "valid_two_I2_surfaces": valid_two_i2_surfaces,
        "declared_split_cores_recovered_exactly": (
            recovered == declared if coordinates == (4, 18, 27) else None
        ),
        "claim_boundary": (
            "Complete only for the fixed GF(31) coordinate slice and s1!=0 branch; "
            "it is not a proof that another slice, prime, or K3 neighbor has no rank-3 seed."
        ),
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--r0", type=int, default=4)
    parser.add_argument("--s0", type=int, default=18)
    parser.add_argument("--x1", type=int, default=27)
    args = parser.parse_args()
    certificate = compute_certificate(args.r0, args.s0, args.x1)
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "counts": certificate["counts"],
        "declared_split_cores_recovered_exactly": certificate[
            "declared_split_cores_recovered_exactly"
        ],
        "certificate_sha256": certificate["certificate_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
