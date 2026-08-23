#!/usr/bin/env python3
"""Derive the complete P1/P2/P3 incidence over the two-split chart.

The output has three deliberately separate layers:

* exact algebra: the 62-variable incidence presentation, its six local
  syzygies, and an exact five-variable triangular reduction for P3;
* modular evidence: the previously certified P1/P2 Jacobian ranks and an
  independent replay of the P3 remainder test on one F_11 seed;
* claim boundary: no full triple, characteristic-zero point, or rank-31 curve
  is asserted.

SymPy is used only for exact polynomial algebra.  The finite-field replay is
implemented independently with small coefficient arrays.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
from typing import Sequence

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SEARCH_CERTIFICATE = (
    ROOT / "certificates" / "e8_a2_semistable_two_split_target_small_primes.json"
)
DEFAULT_PAIR_CERTIFICATE = (
    ROOT / "certificates"
    / "e8_a2_semistable_two_split_pair_incidence_lifts.json"
)


def canonical_sha256(payload: object) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def polynomial_fingerprint(poly: sp.Poly) -> str:
    terms = [
        [list(monomial), str(coefficient)]
        for monomial, coefficient in poly.terms()
    ]
    return canonical_sha256(terms)


def p3_triangular_reduction() -> dict[str, object]:
    """Eliminate Y3 from degrees 12 down to 6 on the infinity open."""

    t = sp.symbols("t")
    r, x0, x1, x2, x3 = sp.symbols("r x0 x1 x2 x3")
    a0, a1, a2 = sp.symbols("a0 a1 a2")
    b0, b1 = sp.symbols("b0 b1")
    c0, c1, lam = sp.symbols("c0 c1 lambda")
    section_variables = (r, x0, x1, x2, x3)
    all_variables = (
        r, x0, x1, x2, x3,
        a0, a1, a2, b0, b1, c0, c1, lam,
    )

    a = a0 + a1 * t + a2 * t**2
    beta = b0 + b1 * t
    gamma = c0 + c1 * t
    D = t * (t - 1) * (t - lam)
    X = x0 + x1 * t + x2 * t**2 + x3 * t**3 + r**2 * t**4
    right = sp.expand(
        X**3 + 3 * a * X**2 - 6 * D * beta * X + D**2 * gamma
    )
    right_coefficients = [
        sp.Poly(right, t).coeff_monomial(t**degree)
        for degree in range(13)
    ]
    if sp.expand(right_coefficients[12] - r**6) != 0:
        raise AssertionError("P3 leading coefficient is not r^6")

    y: list[sp.Expr | None] = [None] * 7
    y[6] = r**3
    solved: list[dict[str, object]] = []
    for degree in range(11, 5, -1):
        index = degree - 6
        known = sum(
            y[left] * y[degree - left]
            for left in range(index + 1, 6)
            if 0 <= degree - left <= 5
        )
        y[index] = sp.cancel(
            (right_coefficients[degree] - known) / (2 * y[6])
        )
        numerator, denominator = sp.fraction(y[index])
        solved.append({
            "coefficient": f"y{index}",
            "from_degree": degree,
            "denominator": str(sp.factor(denominator)),
            "numerator_sha256": polynomial_fingerprint(
                sp.Poly(numerator, *all_variables, domain=sp.QQ)
            ),
        })

    Y = sum(y[index] * t**index for index in range(7))
    identity = sp.Poly(sp.cancel(Y**2 - right), t)
    for degree in range(6, 13):
        if sp.cancel(identity.coeff_monomial(t**degree)) != 0:
            raise AssertionError(f"triangular P3 degree {degree} did not vanish")

    remainders: list[dict[str, object]] = []
    combined_hash_input: list[str] = []
    for degree in range(6):
        residual = sp.cancel(identity.coeff_monomial(t**degree))
        numerator, denominator = sp.fraction(residual)
        full_poly = sp.Poly(numerator, *all_variables, domain=sp.QQ)
        section_poly = sp.Poly(numerator, *section_variables)
        fingerprint = polynomial_fingerprint(full_poly)
        combined_hash_input.append(fingerprint)
        remainders.append({
            "name": f"N{degree}",
            "coefficient_degree": degree,
            "denominator": str(sp.factor(denominator)),
            "section_variable_total_degree": section_poly.total_degree(),
            "term_count": len(full_poly.terms()),
            "numerator_sha256": fingerprint,
        })

    return {
        "infinity_parameterization": ["x4=r^2", "y6=r^3", "r!=0"],
        "solved_high_coefficients": solved,
        "remainder_equations": remainders,
        "remainder_family_sha256": canonical_sha256(combined_hash_input),
        "exact_equivalence": (
            "On r!=0, the P3 coefficient equations are equivalent to the "
            "six numerator equations N0=...=N5=0 after the displayed "
            "triangular definitions of y5,...,y0."
        ),
    }


def verify_cubic_difference_identity() -> None:
    x_left, x_right, a, D, beta, gamma = sp.symbols(
        "x_left x_right a D beta gamma"
    )

    def cubic(x: sp.Expr) -> sp.Expr:
        return x**3 + 3 * a * x**2 - 6 * D * beta * x + D**2 * gamma

    quotient = (
        x_left**2 + x_left * x_right + x_right**2
        + 3 * a * (x_left + x_right) - 6 * D * beta
    )
    if sp.expand(
        cubic(x_left) - cubic(x_right)
        - (x_left - x_right) * quotient
    ) != 0:
        raise AssertionError("difference-of-cubics identity failed")


def poly_mul_mod(left: Sequence[int], right: Sequence[int], prime: int) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] = (result[i + j] + a * b) % prime
    return result


def poly_add_scaled_mod(
    target: list[int], source: Sequence[int], scale: int, prime: int
) -> None:
    for index, coefficient in enumerate(source):
        target[index] = (target[index] + scale * coefficient) % prime


def poly_evaluate_mod(poly: Sequence[int], value: int, prime: int) -> int:
    answer = 0
    for coefficient in reversed(poly):
        answer = (answer * value + coefficient) % prime
    return answer


def replay_p3_remainders(seed: dict[str, object]) -> dict[str, object]:
    """Exhaust the reduced P3 system for one recorded finite-field surface."""

    prime = int(seed["prime"])
    lam = int(seed["parameters"]["lambda"])
    a = [int(value) for value in seed["a"]]
    beta = [int(value) for value in seed["beta"]]
    gamma = [int(value) for value in seed["gamma"]]
    D = poly_mul_mod(
        poly_mul_mod([0, 1], [-1, 1], prime), [-lam, 1], prime
    )
    D2 = poly_mul_mod(D, D, prime)
    Dbeta = poly_mul_mod(D, beta, prime)
    raw = 0
    node_admissible = 0
    hits: list[dict[str, object]] = []
    audit = hashlib.sha256()

    for r in range(1, prime):
        x4 = r * r % prime
        y6 = r * x4 % prime
        inverse_twice_y6 = pow(2 * y6, -1, prime)
        for lower in itertools.product(range(prime), repeat=4):
            raw += 1
            X = [*lower, x4]
            if (
                poly_evaluate_mod(X, 0, prime) == 0
                or poly_evaluate_mod(X, 1, prime) == 0
                or poly_evaluate_mod(X, lam, prime) == 0
            ):
                continue
            node_admissible += 1
            X2 = poly_mul_mod(X, X, prime)
            X3 = poly_mul_mod(X2, X, prime)
            right = X3 + [0] * (13 - len(X3))
            poly_add_scaled_mod(
                right, poly_mul_mod(a, X2, prime), 3, prime
            )
            poly_add_scaled_mod(
                right, poly_mul_mod(Dbeta, X, prime), -6, prime
            )
            poly_add_scaled_mod(
                right, poly_mul_mod(D2, gamma, prime), 1, prime
            )
            y = [0] * 7
            y[6] = y6
            for degree in range(11, 5, -1):
                index = degree - 6
                known = sum(
                    y[left] * y[degree - left]
                    for left in range(index + 1, 6)
                    if 0 <= degree - left <= 5
                ) % prime
                y[index] = (
                    (right[degree] - known) * inverse_twice_y6
                ) % prime
            low_residuals = [
                (
                    sum(y[left] * y[degree - left]
                        for left in range(degree + 1))
                    - right[degree]
                ) % prime
                for degree in range(6)
            ]
            audit.update(bytes([r, *lower, *low_residuals]))
            if not any(low_residuals):
                hits.append({"r": r, "X_lower": list(lower), "Y": y})

    return {
        "prime": prime,
        "parameters": seed["parameters"],
        "raw_r_Xlower_tuples": raw,
        "node_admissible_tuples": node_admissible,
        "audit_sha256": audit.hexdigest(),
        "P3_hits": hits,
        "P3_hit_count": len(hits),
    }


def compute_certificate(
    search_certificate_path: Path = DEFAULT_SEARCH_CERTIFICATE,
    pair_certificate_path: Path = DEFAULT_PAIR_CERTIFICATE,
) -> dict[str, object]:
    verify_cubic_difference_identity()
    triangular = p3_triangular_reduction()
    search_certificate = json.loads(
        search_certificate_path.read_text(encoding="utf-8")
    )
    pair_certificate = json.loads(
        pair_certificate_path.read_text(encoding="utf-8")
    )
    seeds = search_certificate["representative_pair_seeds"]["seeds"]
    seed_ranks = [
        {
            "seed_index": entry["seed_index"],
            "prime": entry["prime"],
            "P1_P2_jacobian_rank": entry["jacobian_rank"],
            "P1_P2_tangent_dimension": entry["zariski_tangent_dimension"],
        }
        for entry in pair_certificate["seed_analysis"]
    ]
    if any(entry["P1_P2_jacobian_rank"] != 28 for entry in seed_ranks):
        raise AssertionError("unexpected P1/P2 seed rank")
    replay = replay_p3_remainders(seeds[0])
    if replay["P3_hit_count"] != 0:
        raise AssertionError("the independent P3 replay unexpectedly found a hit")

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "e8_a2_semistable_two_split_full_incidence",
        "inputs": {
            "search_certificate": {
                "path": str(search_certificate_path.relative_to(ROOT)),
                "sha256": hashlib.sha256(search_certificate_path.read_bytes()).hexdigest(),
            },
            "pair_incidence_certificate": {
                "path": str(pair_certificate_path.relative_to(ROOT)),
                "sha256": hashlib.sha256(pair_certificate_path.read_bytes()).hexdigest(),
            },
        },
        "exact_full_incidence": {
            "surface_parameters": ["lambda", "y", "z", "W"],
            "section_variables": {
                "P1": ["u1_0..u1_2", "v1_0..v1_4"],
                "P2": ["u2_0..u2_3", "v2_0..v2_5"],
                "P3": ["x3_0..x3_4", "y3_0..y3_6"],
            },
            "section_equations": {
                "P1": {
                    "count": 9,
                    "identity": "V1^2=t(t-lambda)U1^3+3aU1^2-6(t-1)beta U1+(t-1)^2 gamma",
                },
                "P2": {
                    "count": 11,
                    "identity": "V2^2=(t-lambda)U2^3+3aU2^2-6t(t-1)beta U2+t^2(t-1)^2 gamma",
                },
                "P3": {
                    "count": 13,
                    "identity": "Y3^2=X3^3+3aX3^2-6D beta X3+D^2 gamma",
                },
            },
            "quadratic_intersection_witnesses": {
                "P1_P2": {
                    "x": "tU1-U2=h12*Cx12",
                    "y": "tV1-V2=h12*Cy12",
                    "new_variables": 8,
                    "equations": 10,
                },
                "P1_P3": {
                    "x": "t(t-lambda)U1-X3=h13*Cx13",
                    "y": "t(t-lambda)V1-Y3=h13*Cy13",
                    "new_variables": 10,
                    "equations": 12,
                },
                "P2_P3": {
                    "x": "(t-lambda)U2-X3=h23*Cx23",
                    "y": "(t-lambda)V2-Y3=h23*Cy23",
                    "new_variables": 10,
                    "equations": 12,
                },
                "all_h_are_monic_quadratics": True,
                "exact_degree_two_open": "gcd(Cx_ij,Cy_ij)=1",
                "same_y_branch_open": "gcd(h_ij,y_i+y_j)=1",
            },
            "raw_variable_count": 62,
            "raw_equation_count": 67,
            "curve_difference_identity": (
                "(y_i-y_j)(y_i+y_j)=(x_i-x_j)"
                "(x_i^2+x_i*x_j+x_j^2+3a(x_i+x_j)-6D beta)"
            ),
            "local_syzygies": (
                "On each same-y-branch open, the x-factor equation and the "
                "two section equations force h_ij to divide y_i-y_j. Hence "
                "two y-factor coefficient equations per pair are redundant."
            ),
            "syzygetic_equation_count": 6,
            "reduced_local_equation_count": 61,
            "reduced_local_jacobian_shape": [61, 62],
            "rank_bound_at_every_open_triple": 61,
            "smooth_triple_criterion": (
                "A target triple is a smooth point of the expected "
                "one-dimensional incidence locus if this reduced Jacobian "
                "has rank 61."
            ),
            "expected_dimension": 1,
            "observed_full_triple_jacobian_rank": None,
            "reason_no_observed_rank": "No finite-field P1/P2/P3 target triple is known.",
        },
        "P3_triangular_elimination": triangular,
        "Noether_Lefschetz_condition": {
            "base_dimension": 4,
            "P3_incidence_dimension_expected": 3,
            "P1_P2_locus_dimension_certified": 2,
            "P1_P2_P3_locus_dimension_expected": 1,
            "exact_affine_ideal": (
                "I3=<N0,...,N5> saturated by r*X3(0)*X3(1)*X3(lambda) "
                "and by the two-split Kodaira-open factors."
            ),
            "projection_ideal": (
                "J3=I3 intersect k[lambda,y,z,W]. On a finite projection "
                "open, Fitt_0 of the pushed-forward quotient algebra gives "
                "a local equation Phi_P3 for the P3 Noether-Lefschetz divisor."
            ),
            "resultant_candidate": (
                "The saturated primitive non-boundary factor of the Macaulay "
                "resultant of the homogenized N0,...,N5 is a computable "
                "multiple of Phi_P3. The large resultant was not expanded."
            ),
            "restriction_to_pair_locus": (
                "Pulling Phi_P3 back to the smooth two-dimensional P1/P2 "
                "incidence cuts the expected one-dimensional full-triple locus."
            ),
        },
        "modular_evidence": {
            "representative_P1_P2_seed_ranks": seed_ranks,
            "independent_first_seed_P3_remainder_replay": replay,
            "all_source_runs_have_zero_P3_sections": all(
                run["sections"]["P3"] == 0
                for run in search_certificate["runs"]
            ),
            "interpretation": (
                "This proves only the displayed finite-field bounded searches "
                "and the first-seed replay. It does not prove geometric "
                "emptiness over algebraic closures or characteristic zero."
            ),
        },
        "open_conditions": [
            "characteristic is not 2 or 3",
            "the two-split rational-chart denominators are units",
            "the II*+3I3+5I1 Kodaira factors are nonzero",
            "r*x3(0)*x3(1)*x3(lambda) is nonzero",
            "each h_ij is an exact common quadratic factor",
            "gcd(h_ij,y_i+y_j)=1 for the selected same-y branches",
            "the three limiting points at infinity are smooth and distinct",
        ],
        "claim_boundary": {
            "full_modular_triple_found": False,
            "Noether_Lefschetz_polynomial_expanded": False,
            "characteristic_zero_triple_found": False,
            "rank31_curve_found": False,
        },
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--search-certificate", type=Path, default=DEFAULT_SEARCH_CERTIFICATE)
    parser.add_argument("--pair-certificate", type=Path, default=DEFAULT_PAIR_CERTIFICATE)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()
    certificate = compute_certificate(
        arguments.search_certificate.resolve(), arguments.pair_certificate.resolve()
    )
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
