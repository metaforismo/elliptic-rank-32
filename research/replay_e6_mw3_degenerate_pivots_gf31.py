#!/usr/bin/env python3
"""Exhaust the exceptional P1 triangular-pivot planes over GF(31).

For each selected coordinate slice, the generic chain

    P1_7 -> b4, P1_6 -> b5, P1_5 -> y2, P1_4 -> a3

fails at the unique value of a1 for which y0+y1=0.  On that plane the exact
replacement chain is

    P1_7 -> b4, P1_6 -> b5, P1_5 -> a3, P1_4 -> y2.

After these eliminations only (a2,a4,s1) remain, so the declared s1!=0 branch
has exactly 31^2*30 = 28,830 triples per coordinate slice.  This script
exhausts them, reconstructs every surface, verifies the two split I2 fibers,
and runs the canonical P2/P3 searches before emitting a certificate.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import sympy as sp

from replay_e6_mw3_core_exhaustion_gf31 import (
    compile_polynomial,
    compile_rational,
    evaluate_rational,
    evaluate_terms,
    poly_exact_divide,
    poly_gcd,
)
from replay_e6_mw3_p2_split_cores import (
    PRIME,
    build_p1_candidate,
    canonical_sha256,
    canonical_sign,
    derive_degenerate_branch_equations,
    poly_add,
    poly_derivative,
    poly_eval,
    poly_mul,
    poly_pow,
    poly_scale,
    search_p2,
    search_p3_polynomial,
    section_residual,
    valuation_at,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "certificates" / "e6_mw3_degenerate_pivots_gf31.json"
SLICES = [
    (4, 18, 27),
    (1, 1, 1),
    (2, 1, 1),
    (3, 1, 1),
    (1, 2, 1),
    (1, 1, 2),
]


def solve_affine(equation, variable):
    coefficient = sp.cancel(sp.diff(equation, variable), modulus=PRIME)
    if coefficient == 0 or sp.diff(coefficient, variable) != 0:
        raise AssertionError(f"equation is not affine-linear in {variable}")
    return sp.cancel(
        -equation.subs(variable, 0) / coefficient,
        modulus=PRIME,
    ), coefficient


def compile_slice(coordinates):
    expressions, symbol_map, equations, metadata = derive_degenerate_branch_equations(
        *coordinates
    )
    a2, a3, a4, s1, y2 = (
        symbol_map[name] for name in ("a2", "a3", "a4", "s1", "y2")
    )
    variables = (a2, a4, s1)
    a3_expression, a3_pivot = solve_affine(equations["P1_5"], a3)
    p14 = sp.cancel(equations["P1_4"].subs(a3, a3_expression), modulus=PRIME)
    y2_expression, y2_pivot = solve_affine(p14, y2)
    substitutions = {a3: a3_expression, y2: y2_expression}

    core_equations = {
        name: sp.cancel(
            equation.subs(a3, a3_expression).subs(y2, y2_expression),
            modulus=PRIME,
        )
        for name, equation in equations.items()
        if name in ("D2_at_1", "D3_at_1")
    }
    compiled_core_equations = {
        name: compile_rational(equation, variables)
        for name, equation in core_equations.items()
    }

    b5_expression = sp.cancel(
        expressions["b5"].subs(substitutions), modulus=PRIME
    )
    b4_expression = sp.cancel(
        expressions["b4"].subs(substitutions), modulus=PRIME
    )
    coordinate_expressions = {
        "a3": a3_expression,
        "y2": y2_expression,
        "b5": b5_expression,
        "b4": b4_expression,
    }
    compiled_coordinates = {
        name: compile_rational(expression, variables)
        for name, expression in coordinate_expressions.items()
    }
    equation_terms = {}
    for name, expression in core_equations.items():
        numerator, denominator = sp.cancel(expression, modulus=PRIME).as_numer_denom()
        equation_terms[name] = {
            "numerator": [
                [list(monomial), coefficient]
                for monomial, coefficient in compile_polynomial(numerator, variables)
            ],
            "denominator": [
                [list(monomial), coefficient]
                for monomial, coefficient in compile_polynomial(denominator, variables)
            ],
        }
    return {
        "symbol_map": symbol_map,
        "metadata": metadata,
        "a3_pivot": str(a3_pivot),
        "y2_pivot": str(y2_pivot),
        "compiled_core_equations": compiled_core_equations,
        "compiled_coordinates": compiled_coordinates,
        "equation_terms": equation_terms,
    }


def enumerate_core_solutions(compiled):
    solutions = []
    denominator_failures = 0
    examined = 0
    equations = compiled["compiled_core_equations"]
    for a2 in range(PRIME):
        for a4 in range(PRIME):
            for s1 in range(1, PRIME):
                examined += 1
                values = (a2, a4, s1)
                try:
                    if all(
                        evaluate_rational(equation, values) == 0
                        for equation in equations.values()
                    ):
                        solutions.append(values)
                except ZeroDivisionError:
                    denominator_failures += 1
    return examined, denominator_failures, solutions


def reconstruct_surface(coordinates, compiled, triple):
    a2, a4, s1 = triple
    metadata = compiled["metadata"]
    a1 = metadata["a1"]
    values = {
        name: evaluate_rational(expression, triple)
        for name, expression in compiled["compiled_coordinates"].items()
    }
    core = (a1, a2, a4, s1)
    A, B, X1, Y1 = build_p1_candidate(
        core,
        values["a3"],
        values["y2"],
        values["b5"],
        values["b4"],
        coordinates,
    )
    if section_residual(A, B, X1, Y1) != [0]:
        raise AssertionError(("P1 identity failed", coordinates, triple, values))
    discriminant = poly_scale(
        poly_add(poly_scale(poly_pow(A, 3), 4), poly_scale(poly_pow(B, 2), 27)),
        -16,
    )
    multiplicity_zero = valuation_at(discriminant, 0)
    multiplicity_one = valuation_at(discriminant, 1)
    record = {
        "core_a1_a2_a4_s1": list(core),
        "triangular_coordinates_a3_y2_b5_b4": [
            values[name] for name in ("a3", "y2", "b5", "b4")
        ],
        "multiplicity_at_0": multiplicity_zero,
        "multiplicity_at_1": multiplicity_one,
    }
    if multiplicity_zero != 4 or multiplicity_one != 4:
        record["status"] = "boundary_or_wrong_I4_multiplicity"
        return record

    known_i4 = poly_mul(poly_pow([0, 1], 4), poly_pow([-1, 1], 4))
    residual_discriminant = poly_exact_divide(discriminant, known_i4)
    repeated = poly_gcd(
        residual_discriminant, poly_derivative(residual_discriminant)
    )
    repeated_roots = [
        point for point in range(PRIME) if poly_eval(repeated, point) == 0
    ]
    derivative_a = poly_derivative(A)
    derivative_b = poly_derivative(B)
    nodes = {}
    for point in repeated_roots:
        matches = [
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
        if matches:
            nodes[str(point)] = matches
    valid_i2_roots = [
        point
        for point in repeated_roots
        if valuation_at(discriminant, point) == 2 and str(point) in nodes
    ]
    record.update(
        {
            "status": "exact_I4",
            "A": A,
            "B": B,
            "X1": X1,
            "Y1": Y1,
            "residual_discriminant": residual_discriminant,
            "repeated_root_polynomial": repeated,
            "repeated_roots_in_GF31": repeated_roots,
            "singular_nodes": nodes,
            "valid_I2_roots": valid_i2_roots,
        }
    )
    return record


def p2_key(hit):
    return (
        hit["pole"],
        hit["q0"],
        tuple(hit["X2"]),
        tuple(canonical_sign(hit["Y2"])),
    )


def search_sections(surface, s0):
    roots = surface["valid_I2_roots"]
    nodes = {int(key): value for key, value in surface["singular_nodes"].items()}
    p2_cases = 0
    p2_hits = {}
    p3_cases = 0
    p3_hits = {}
    for left, right in itertools.combinations(roots, 2):
        for lam, mu in ((left, right), (right, left)):
            for sl in nodes[lam]:
                for sm in nodes[mu]:
                    tested, hits = search_p2(
                        surface["A"],
                        surface["B"],
                        lam,
                        mu,
                        sl,
                        sm,
                        s0=s0,
                    )
                    p2_cases += tested
                    for hit in hits:
                        p2_hits[p2_key(hit)] = hit
        tested, hits = search_p3_polynomial(
            surface["A"],
            surface["B"],
            surface["core_a1_a2_a4_s1"][3],
            left,
            right,
            nodes[left][0],
            nodes[right][0],
            s0=s0,
        )
        p3_cases += tested
        for hit in hits:
            key = (tuple(hit["X3"]), tuple(hit["Y3_up_to_sign"]))
            p3_hits[key] = hit
    return {
        "ordered_P2_cases": p2_cases,
        "geometric_P2_hits": list(p2_hits.values()),
        "P3_cases": p3_cases,
        "P3_candidates": list(p3_hits.values()),
    }


def compute_certificate():
    slice_records = []
    totals = {
        "coordinate_slices": len(SLICES),
        "degenerate_core_triples_examined": 0,
        "simultaneous_core_solutions": 0,
        "valid_two_I2_surfaces": 0,
        "ordered_P2_cases": 0,
        "geometric_P2_hits": 0,
        "P3_cases": 0,
        "P3_candidates": 0,
        "joint_P2_P3_candidates": 0,
    }
    joint_candidates = []
    for coordinates in SLICES:
        compiled = compile_slice(coordinates)
        examined, denominator_failures, solutions = enumerate_core_solutions(compiled)
        reconstructed = [
            reconstruct_surface(coordinates, compiled, triple) for triple in solutions
        ]
        valid_surfaces = [
            record for record in reconstructed if len(record.get("valid_I2_roots", [])) >= 2
        ]
        for surface in valid_surfaces:
            searches = search_sections(surface, coordinates[1])
            surface.update(searches)
            if searches["geometric_P2_hits"] and searches["P3_candidates"]:
                joint_candidates.append(
                    {
                        "coordinate_slice": list(coordinates),
                        "core_a1_a2_a4_s1": surface["core_a1_a2_a4_s1"],
                        "surface": surface,
                    }
                )
            totals["ordered_P2_cases"] += searches["ordered_P2_cases"]
            totals["geometric_P2_hits"] += len(searches["geometric_P2_hits"])
            totals["P3_cases"] += searches["P3_cases"]
            totals["P3_candidates"] += len(searches["P3_candidates"])
        totals["degenerate_core_triples_examined"] += examined
        totals["simultaneous_core_solutions"] += len(solutions)
        totals["valid_two_I2_surfaces"] += len(valid_surfaces)
        slice_records.append(
            {
                "coordinate_slice": list(coordinates),
                "degenerate_a1": compiled["metadata"]["a1"],
                "replacement_pivots": {
                    "P1_5_to_a3": compiled["a3_pivot"],
                    "P1_4_to_y2": compiled["y2_pivot"],
                },
                "core_equations": compiled["equation_terms"],
                "counts": {
                    "triples_examined": examined,
                    "denominator_failures": denominator_failures,
                    "simultaneous_core_solutions": len(solutions),
                    "exact_I4_surfaces": sum(
                        record.get("status") == "exact_I4" for record in reconstructed
                    ),
                    "valid_two_I2_surfaces": len(valid_surfaces),
                },
                "all_core_solutions_sha256": canonical_sha256(
                    [list(solution) for solution in solutions]
                ),
                "valid_two_I2_surfaces": valid_surfaces,
            }
        )
    totals["joint_P2_P3_candidates"] = len(joint_candidates)
    payload = {
        "schema": "elliptic-rank30/e6-mw3-degenerate-pivots-gf31/v1",
        "field": "GF(31)",
        "branch": "y0+y1=0 with s1!=0",
        "slice_records": slice_records,
        "totals": totals,
        "joint_candidates": joint_candidates,
        "claim_boundary": (
            "Complete for the exceptional y0+y1=0 planes of the six listed "
            "GF(31) coordinate slices. Other slices, primes, and K3 neighbors remain open."
        ),
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    payload = compute_certificate()
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "totals": payload["totals"],
                "certificate_sha256": payload["certificate_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
