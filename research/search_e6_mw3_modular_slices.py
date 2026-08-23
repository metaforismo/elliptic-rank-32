#!/usr/bin/env python3
"""Search exact regular and exceptional E6/MW3 coordinate slices at a good prime.

The original replay modules intentionally freeze p=31 for stable certificates.
This driver reconfigures their small finite-field kernels at runtime and emits
a prime-specific scouting certificate.  Every retained surface is reconstructed
from the original section identity and filtered by exact Kodaira multiplicities
and rational nodal I2 conditions before P2/P3 searches are attempted.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import sympy as sp

import replay_e6_mw3_core_exhaustion_gf31 as core_search
import replay_e6_mw3_degenerate_pivots_gf31 as degenerate_search
import replay_e6_mw3_p2_split_cores as p2_search
import check_e6_mw3_candidate_relations_gf31 as group
import classify_e6_mw3_component_profiles as profiles


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SLICES = [(1, 1, 1)]


def is_prime(value):
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, int(value**0.5) + 1))


def configure_prime(prime):
    if not is_prime(prime) or prime in (2, 3, 79):
        raise ValueError("choose a good prime different from 2, 3, and 79")
    p2_search.PRIME = prime
    core_search.PRIME = prime
    degenerate_search.PRIME = prime
    group.configure_prime(prime)
    core_search.POWER_TABLE = [
        [pow(value, exponent, prime) for exponent in range(24)]
        for value in range(prime)
    ]


def canonical_p2_key(hit):
    return (
        hit["pole"],
        hit["q0"],
        tuple(hit["X2"]),
        tuple(p2_search.canonical_sign(hit["Y2"])),
    )


def search_sections(surface, core_key, s0, target_profile_only=False):
    roots = surface["valid_I2_roots"]
    nodes = {int(key): value for key, value in surface["singular_nodes"].items()}
    p2_cases = 0
    p2_hits = {}
    p3_cases = 0
    p3_hits = {}
    profile_pairs = []
    A_function = group.polynomial(surface["A"])
    B_function = group.polynomial(surface["B"])
    P1 = profiles.polynomial_point(surface["X1"], surface["Y1"])
    s1 = surface[core_key][3]
    for left, right in itertools.combinations(roots, 2):
        finite_fibers = {
            "I4_0": {"point": 0, "node": s0},
            "I4_1": {"point": 1, "node": s1},
            "I2_lambda": {"point": left, "node": nodes[left][0]},
            "I2_mu": {"point": right, "node": nodes[right][0]},
        }
        p1_classification = profiles.classify_point(
            P1, A_function, finite_fibers, p2_search.PRIME
        )
        p1_vector = profiles.component_vector(p1_classification)
        p1_matches = p1_vector == profiles.TARGET_FINITE["P1"]
        if target_profile_only and not p1_matches:
            profile_pairs.append(
                {
                    "I2_roots": [left, right],
                    "P1": {
                        "finite_component_vector": p1_vector,
                        "matches_target": False,
                        "fibers": p1_classification,
                    },
                    "P2_records": [],
                    "P3_records": [],
                    "matching_target_P2": 0,
                    "matching_target_P3": 0,
                    "matching_target_triples": 0,
                    "search_skipped_after_P1_filter": True,
                }
            )
            continue

        pair_p2_hits = {}
        for lam, mu in ((left, right), (right, left)):
            for sl in nodes[lam]:
                for sm in nodes[mu]:
                    tested, hits = p2_search.search_p2(
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
                        p2_hits[canonical_p2_key(hit)] = hit
                        pair_p2_hits[canonical_p2_key(hit)] = hit
        p2_records = []
        for hit in pair_p2_hits.values():
            classification = profiles.classify_point(
                profiles.p2_point(hit),
                A_function,
                finite_fibers,
                p2_search.PRIME,
            )
            vector = profiles.component_vector(classification)
            p2_records.append(
                {
                    "section": hit,
                    "finite_component_vector": vector,
                    "matches_target": vector == profiles.TARGET_FINITE["P2"],
                    "fibers": classification,
                }
            )
        matching_p2 = sum(record["matches_target"] for record in p2_records)

        pair_p3_hits = {}
        if not target_profile_only or matching_p2:
            tested, hits = p2_search.search_p3_polynomial(
                surface["A"],
                surface["B"],
                surface[core_key][3],
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
                pair_p3_hits[key] = hit
        p3_records = []
        for hit in pair_p3_hits.values():
            point = profiles.polynomial_point(hit["X3"], hit["Y3_up_to_sign"])
            classification = profiles.classify_point(
                point, A_function, finite_fibers, p2_search.PRIME
            )
            vector = profiles.component_vector(classification)
            p3_records.append(
                {
                    "section": hit,
                    "finite_component_vector": vector,
                    "matches_target": vector == profiles.TARGET_FINITE["P3"],
                    "fibers": classification,
                }
            )
        matching_p3 = sum(record["matches_target"] for record in p3_records)
        profile_pairs.append(
            {
                "I2_roots": [left, right],
                "P1": {
                    "finite_component_vector": p1_vector,
                    "matches_target": p1_matches,
                    "fibers": p1_classification,
                },
                "P2_records": p2_records,
                "P3_records": p3_records,
                "matching_target_P2": matching_p2,
                "matching_target_P3": matching_p3,
                "matching_target_triples": (
                    matching_p2 * matching_p3
                    if p1_matches
                    else 0
                ),
                "search_skipped_after_P2_filter": bool(
                    target_profile_only and not matching_p2
                ),
            }
        )
    return {
        "ordered_P2_cases": p2_cases,
        "geometric_P2_hits": list(p2_hits.values()),
        "P3_cases": p3_cases,
        "P3_candidates": list(p3_hits.values()),
        "finite_component_profile_pairs": profile_pairs,
    }


def equation_terms(equations, symbol_map):
    variables = tuple(symbol_map[name] for name in ("a1", "a2", "a4", "s1"))
    result = {}
    for name, expression in equations.items():
        numerator, denominator = sp.cancel(
            expression, modulus=p2_search.PRIME
        ).as_numer_denom()
        result[name] = {
            "numerator": [
                [list(monomial), coefficient]
                for monomial, coefficient in core_search.compile_polynomial(
                    numerator, variables
                )
            ],
            "denominator": [
                [list(monomial), coefficient]
                for monomial, coefficient in core_search.compile_polynomial(
                    denominator, variables
                )
            ],
        }
    return result


def regular_slice(coordinates, target_profile_only=False):
    expressions, symbol_map, equations = p2_search.derive_triangular_expressions(
        *coordinates
    )
    compiled_expressions = core_search.compile_specialized_triangular_expressions(
        expressions, symbol_map, *coordinates
    )
    examined, d2_roots, solutions, _d2_terms, _d3_terms = core_search.enumerate_cores(
        equations, symbol_map
    )
    classified = [
        core_search.classify_solution(solution, compiled_expressions, coordinates)
        for solution in solutions
    ]
    for record in classified:
        if "roots_in_GF31" in record:
            record["roots_in_field"] = record.pop("roots_in_GF31")
    r0, s0, x1 = coordinates
    prime = p2_search.PRIME
    y0 = r0 * (r0**2 - 3 * s0) % prime
    constant = 3 * (r0**2 - s0) * x1 % prime
    degenerate_a1 = (-constant - 2 * r0 * y0) % prime
    surfaces = [
        record
        for record in classified
        if len(record.get("valid_I2_roots", [])) >= 2
        and record["core"][0] != degenerate_a1
    ]
    for surface in surfaces:
        surface.update(
            search_sections(surface, "core", s0, target_profile_only)
        )
    return {
        "coordinate_slice": list(coordinates),
        "degenerate_a1": degenerate_a1,
        "core_equations": equation_terms(equations, symbol_map),
        "counts": {
            "quadruples_examined": examined,
            "D2_roots_tested_by_D3": d2_roots,
            "simultaneous_core_solutions": len(solutions),
            "solutions_on_degenerate_plane": sum(
                solution[0] == degenerate_a1 for solution in solutions
            ),
            "valid_two_I2_surfaces": len(surfaces),
        },
        "all_core_solutions_sha256": p2_search.canonical_sha256(
            [list(solution) for solution in solutions]
        ),
        "valid_two_I2_surfaces": surfaces,
    }


def degenerate_slice(coordinates, target_profile_only=False):
    compiled = degenerate_search.compile_slice(coordinates)
    examined, denominator_failures, solutions = (
        degenerate_search.enumerate_core_solutions(compiled)
    )
    reconstructed = [
        degenerate_search.reconstruct_surface(coordinates, compiled, solution)
        for solution in solutions
    ]
    for record in reconstructed:
        if "repeated_roots_in_GF31" in record:
            record["repeated_roots_in_field"] = record.pop(
                "repeated_roots_in_GF31"
            )
    surfaces = [
        record
        for record in reconstructed
        if len(record.get("valid_I2_roots", [])) >= 2
    ]
    for surface in surfaces:
        surface.update(
            search_sections(
                surface,
                "core_a1_a2_a4_s1",
                coordinates[1],
                target_profile_only,
            )
        )
    return {
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
            "valid_two_I2_surfaces": len(surfaces),
        },
        "all_core_solutions_sha256": p2_search.canonical_sha256(
            [list(solution) for solution in solutions]
        ),
        "valid_two_I2_surfaces": surfaces,
    }


def compute_certificate(
    prime,
    slices,
    target_profile_only=False,
    regular_only=False,
    skip_chart_errors=False,
):
    configure_prime(prime)
    records = []
    totals = {
        "coordinate_slices": len(slices),
        "regular_quadruples_examined": 0,
        "degenerate_triples_examined": 0,
        "regular_core_solutions": 0,
        "degenerate_core_solutions": 0,
        "valid_two_I2_surfaces": 0,
        "ordered_P2_cases": 0,
        "geometric_P2_hits": 0,
        "P3_cases": 0,
        "P3_candidates": 0,
        "joint_P2_P3_candidates": 0,
        "target_finite_profile_triples": 0,
        "slice_chart_errors": 0,
    }
    joint_candidates = []
    for coordinates in slices:
        coordinates = tuple(value % prime for value in coordinates)
        try:
            regular = regular_slice(coordinates, target_profile_only)
        except (AssertionError, ZeroDivisionError) as error:
            if not skip_chart_errors:
                raise
            totals["slice_chart_errors"] += 1
            records.append(
                {
                    "coordinate_slice": list(coordinates),
                    "status": "chart_error_skipped",
                    "error_type": type(error).__name__,
                    "error": str(error),
                }
            )
            continue
        if regular_only:
            degenerate = {
                "coordinate_slice": list(coordinates),
                "status": "not_run_regular_only",
                "counts": {
                    "triples_examined": 0,
                    "simultaneous_core_solutions": 0,
                },
                "valid_two_I2_surfaces": [],
            }
        else:
            degenerate = degenerate_slice(coordinates, target_profile_only)
            if (
                regular["counts"]["solutions_on_degenerate_plane"]
                != degenerate["counts"]["simultaneous_core_solutions"]
            ):
                raise AssertionError("regular/degenerate core counts disagree")
        surfaces = (
            regular["valid_two_I2_surfaces"]
            + degenerate["valid_two_I2_surfaces"]
        )
        for branch, branch_record in (("regular", regular), ("degenerate", degenerate)):
            core_key = "core" if branch == "regular" else "core_a1_a2_a4_s1"
            for surface in branch_record["valid_two_I2_surfaces"]:
                totals["ordered_P2_cases"] += surface["ordered_P2_cases"]
                totals["geometric_P2_hits"] += len(surface["geometric_P2_hits"])
                totals["P3_cases"] += surface["P3_cases"]
                totals["P3_candidates"] += len(surface["P3_candidates"])
                if surface["geometric_P2_hits"] and surface["P3_candidates"]:
                    joint_candidates.append(
                        {
                            "branch": branch,
                            "coordinate_slice": list(coordinates),
                            "core": surface[core_key],
                            "surface": surface,
                        }
                    )
                totals["target_finite_profile_triples"] += sum(
                    pair["matching_target_triples"]
                    for pair in surface["finite_component_profile_pairs"]
                )
        totals["regular_quadruples_examined"] += regular["counts"][
            "quadruples_examined"
        ]
        totals["degenerate_triples_examined"] += degenerate["counts"][
            "triples_examined"
        ]
        totals["regular_core_solutions"] += regular["counts"][
            "simultaneous_core_solutions"
        ]
        totals["degenerate_core_solutions"] += degenerate["counts"][
            "simultaneous_core_solutions"
        ]
        totals["valid_two_I2_surfaces"] += len(surfaces)
        records.append(
            {
                "coordinate_slice": list(coordinates),
                "regular": regular,
                "degenerate": degenerate,
            }
        )
    totals["joint_P2_P3_candidates"] = len(joint_candidates)
    payload = {
        "schema": "elliptic-rank30/e6-mw3-modular-slices/v2",
        "field_prime": prime,
        "target_profile_only": target_profile_only,
        "regular_only": regular_only,
        "skip_chart_errors": skip_chart_errors,
        "slice_records": records,
        "totals": totals,
        "joint_candidates": joint_candidates,
        "claim_boundary": (
            ("Complete for the regular branch only of the listed " if regular_only
             else "Complete for the regular and exceptional branches of the listed ")
            + f"GF({prime}) coordinate slices only. A modular seed is not a "
            "characteristic-zero rank certificate. Finite component profiles "
            "are exact up to the sign ambiguity between odd I4 classes 1 and 3. "
            + ("Any listed chart errors were skipped and are not negative results."
               if skip_chart_errors else "")
        ),
    }
    payload["certificate_sha256"] = p2_search.canonical_sha256(payload)
    return payload


def parse_slices(text):
    slices = []
    for item in text.split(";"):
        values = tuple(int(value) for value in item.split(","))
        if len(values) != 3:
            raise ValueError("each slice must be r0,s0,x1")
        slices.append(values)
    return slices


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=29)
    parser.add_argument("--slices", default="1,1,1")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--target-profile-only", action="store_true")
    parser.add_argument("--regular-only", action="store_true")
    parser.add_argument("--skip-chart-errors", action="store_true")
    args = parser.parse_args()
    slices = parse_slices(args.slices)
    payload = compute_certificate(
        args.prime,
        slices,
        args.target_profile_only,
        args.regular_only,
        args.skip_chart_errors,
    )
    output = args.output or (
        ROOT / "certificates" / f"e6_mw3_modular_slices_gf{args.prime}.json"
    )
    if not args.no_write:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
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
