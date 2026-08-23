#!/usr/bin/env python3
"""Reverse-engineer the minimal Riemann--Roch function for ICARM curve #273.

The thirty published independent points P_1,...,P_30 are completed by

    P_31 = -(P_1+...+P_30).

On the integral model y^2+x*y=x^3+a4*x+a6 put Y=2*y+x, so

    Y^2 = D(x) = 4*x^3+x^2+4*a4*x+4*a6.

Riemann--Roch predicts a unique function up to scale

    f=A(x)+Y*B(x),  deg A<=15, deg B=14,

with these 31 points as its zeros.  This script reconstructs that function
by exact rational linear algebra and verifies the complete norm identity

    A(x)^2-D(x)B(x)^2 = c*product_i(x-x(P_i)).

The computation is a concrete audit of the degree-31 representation behind a
rank-30 configuration.  It neither finds a new point nor proves rank 31.
SymPy is used only as an exact QQ linear-algebra engine.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from math import gcd, lcm
from pathlib import Path
from typing import Optional


ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = ROOT / "baseline" / "icarm_curve_273.json"
DEFAULT_CERTIFICATE_PATH = (
    ROOT / "certificates" / "icarm273_rr_reverse_certificate.json"
)

Point = Optional[tuple[Fraction, Fraction]]


def ec_neg(point: Point) -> Point:
    if point is None:
        return None
    x, y = point
    return x, -x - y


def ec_add(point1: Point, point2: Point, a4: int) -> Point:
    """Exact group law on y^2+x*y=x^3+a4*x+a6."""
    if point1 is None:
        return point2
    if point2 is None:
        return point1
    x1, y1 = point1
    x2, y2 = point2
    if x1 == x2 and y1 + y2 + x1 == 0:
        return None
    if point1 == point2:
        denominator = 2 * y1 + x1
        if denominator == 0:
            return None
        slope = (3 * x1 * x1 + a4 - y1) / denominator
    else:
        slope = (y2 - y1) / (x2 - x1)
    intercept = y1 - slope * x1
    x3 = slope * slope + slope - x1 - x2
    y3 = -(slope + 1) * x3 - intercept
    return x3, y3


def point_sum(points: list[tuple[Fraction, Fraction]], a4: int) -> Point:
    total: Point = None
    for point in points:
        total = ec_add(total, point, a4)
    return total


def poly_mul(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    result = [Fraction(0) for _ in range(len(left) + len(right) - 1)]
    for i, coefficient_left in enumerate(left):
        for j, coefficient_right in enumerate(right):
            result[i + j] += coefficient_left * coefficient_right
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def poly_sub(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    width = max(len(left), len(right))
    result = [Fraction(0) for _ in range(width)]
    for index in range(width):
        if index < len(left):
            result[index] += left[index]
        if index < len(right):
            result[index] -= right[index]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def product_roots(roots: list[Fraction]) -> list[Fraction]:
    result = [Fraction(1)]
    for root in roots:
        result = poly_mul(result, [-root, Fraction(1)])
    return result


def fraction_string(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def fraction_digest(values: list[Fraction]) -> str:
    encoded = "\n".join(fraction_string(value) for value in values).encode()
    return hashlib.sha256(encoded).hexdigest()


def canonical_sha256(payload: object) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def primitive_integer_vector(values: list[Fraction]) -> list[int]:
    common_denominator = 1
    for value in values:
        common_denominator = lcm(common_denominator, value.denominator)
    integers = [
        value.numerator * (common_denominator // value.denominator)
        for value in values
    ]
    common_gcd = 0
    for value in integers:
        common_gcd = gcd(common_gcd, abs(value))
    if common_gcd == 0:
        raise AssertionError("zero kernel vector")
    integers = [value // common_gcd for value in integers]
    first_nonzero = next(value for value in integers if value)
    if first_nonzero < 0:
        integers = [-value for value in integers]
    return integers


def reconstruct_function(
    points: list[tuple[Fraction, Fraction]],
) -> tuple[list[Fraction], list[Fraction], int, int]:
    try:
        import sympy as sp
        from sympy.polys.matrices import DomainMatrix
    except ImportError as error:
        raise RuntimeError(
            "SymPy is required for the exact QQ nullspace computation"
        ) from error

    rows = []
    for x, y in points:
        capital_y = 2 * y + x
        powers = [Fraction(1)]
        for _ in range(15):
            powers.append(powers[-1] * x)
        row = powers[:16] + [capital_y * value for value in powers[:15]]
        rows.append(
            [sp.Rational(value.numerator, value.denominator) for value in row]
        )

    matrix = DomainMatrix.from_list_sympy(31, 31, rows)
    nullspace = matrix.nullspace()
    if nullspace.shape != (1, 31):
        raise AssertionError(f"unexpected nullspace shape {nullspace.shape}")
    sympy_values = list(nullspace.to_Matrix().row(0))
    fractions = [Fraction(int(value.p), int(value.q)) for value in sympy_values]
    primitive = primitive_integer_vector(fractions)
    coefficients = [Fraction(value) for value in primitive]
    a_coefficients = coefficients[:16]
    b_coefficients = coefficients[16:]

    for index, (x, y) in enumerate(points, start=1):
        capital_y = 2 * y + x
        a_value = sum(
            coefficient * x**power
            for power, coefficient in enumerate(a_coefficients)
        )
        b_value = sum(
            coefficient * x**power
            for power, coefficient in enumerate(b_coefficients)
        )
        if a_value + capital_y * b_value != 0:
            raise AssertionError(f"function does not vanish at P{index}")

    rank = matrix.rank()
    return a_coefficients, b_coefficients, rank, len(sympy_values) - rank


def compute_certificate() -> dict[str, object]:
    data = json.loads(INPUT_PATH.read_text(encoding="utf-8"))
    a1, a2, a3, a4, a6 = [int(value) for value in data["a_invariants"]]
    if (a1, a2, a3) != (1, 0, 0):
        raise AssertionError("unexpected Weierstrass model")
    original_points = [
        (Fraction(x), Fraction(y)) for x, y in data["points"]
    ]
    if len(original_points) != 30:
        raise AssertionError("expected 30 original points")

    original_sum = point_sum(original_points, a4)
    if original_sum is None:
        raise AssertionError("the sum of the independent points cannot be O")
    dependent_point = ec_neg(original_sum)
    assert dependent_point is not None
    all_points = [*original_points, dependent_point]
    if point_sum(all_points, a4) is not None:
        raise AssertionError("the completed 31-point divisor does not sum to O")
    x_coordinates = [point[0] for point in all_points]
    if len(set(x_coordinates)) != 31:
        raise AssertionError("the completed divisor does not have distinct x-values")

    for index, (x, y) in enumerate(all_points, start=1):
        if y * y + x * y != x**3 + a4 * x + a6:
            raise AssertionError(f"P{index} is not on the curve")

    a_coefficients, b_coefficients, matrix_rank, nullity = reconstruct_function(
        all_points
    )
    if matrix_rank != 30 or nullity != 1:
        raise AssertionError(
            f"evaluation matrix rank/nullity is {matrix_rank}/{nullity}"
        )
    if a_coefficients[-1] == 0 or b_coefficients[-1] == 0:
        raise AssertionError("unexpected degree drop in A or B")

    d_polynomial = [
        Fraction(4 * a6),
        Fraction(4 * a4),
        Fraction(1),
        Fraction(4),
    ]
    norm = poly_sub(
        poly_mul(a_coefficients, a_coefficients),
        poly_mul(d_polynomial, poly_mul(b_coefficients, b_coefficients)),
    )
    roots_polynomial = product_roots(x_coordinates)
    if len(norm) != 32 or len(roots_polynomial) != 32:
        raise AssertionError("expected degree-31 norm and root product")
    scalar = norm[-1]
    if norm != [scalar * coefficient for coefficient in roots_polynomial]:
        raise AssertionError("exact split norm identity failed")

    x31, y31 = dependent_point
    coefficient_values = [*a_coefficients, *b_coefficients]
    max_coefficient_digits = max(
        max(len(str(abs(value.numerator))), len(str(value.denominator)))
        for value in coefficient_values
    )
    payload: dict[str, object] = {
        "schema_version": 1,
        "candidate_id": data["candidate_id"],
        "claim": (
            "The 30 independent ICARM points, completed by the negative of "
            "their sum, are exactly the 31 simple rational zeros of an "
            "explicitly reconstructed f=A(x)+(2y+x)B(x) in L(31O)."
        ),
        "conditional_assumptions": [],
        "curve_model": {
            "equation": "y^2+x*y=x^3+a4*x+a6",
            "a4": str(a4),
            "a6": str(a6),
            "quadratic_model": "Y=2y+x; Y^2=4x^3+x^2+4a4*x+4a6",
        },
        "divisor": {
            "published_point_count": 30,
            "completion": "P31=-(P1+...+P30)",
            "completed_sum": "O",
            "distinct_x_coordinates": 31,
            "P31": [fraction_string(x31), fraction_string(y31)],
            "P31_coordinate_sha256": fraction_digest([x31, y31]),
        },
        "exact_linear_algebra": {
            "evaluation_matrix": "31x31 over Q",
            "basis": "1,x,...,x^15,Y,Yx,...,Yx^14",
            "rank": matrix_rank,
            "nullity": nullity,
        },
        "function": {
            "shape": "A(x)+Y*B(x)",
            "degree_A": len(a_coefficients) - 1,
            "degree_B": len(b_coefficients) - 1,
            "primitive_integer_coefficient_count": len(coefficient_values),
            "max_numerator_or_denominator_digits": max_coefficient_digits,
            "coefficients_A_then_B_sha256": fraction_digest(coefficient_values),
        },
        "norm_identity": {
            "identity": "A^2-(4x^3+x^2+4a4*x+4a6)B^2=c*product_i(x-x_i)",
            "degree": len(norm) - 1,
            "distinct_rational_roots": len(x_coordinates),
            "verified_by_exact_coefficient_comparison": True,
            "scalar_c": fraction_string(scalar),
            "norm_coefficients_sha256": fraction_digest(norm),
            "root_product_coefficients_sha256": fraction_digest(roots_polynomial),
        },
        "interpretation": (
            "This realizes the minimal one-function rank-30 representation "
            "on the new record curve. It validates the reverse-construction "
            "machinery but does not add a 31st independent point."
        ),
        "implementation": {
            "script": "research/icarm273_rr_reverse_certificate.py",
            "arithmetic": "exact Q arithmetic; SymPy DomainMatrix nullspace",
        },
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--compare", type=Path, default=None
    )
    arguments = parser.parse_args()
    certificate = compute_certificate()
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(f"wrote {arguments.output}")
    if arguments.compare:
        committed = json.loads(arguments.compare.read_text(encoding="utf-8"))
        if committed != certificate:
            raise AssertionError("committed certificate does not match")
        print(f"matched {arguments.compare}")
    print("completed divisor: 31 distinct rational points, sum O")
    print("evaluation matrix: rank 30, nullity 1")
    print("function degrees: deg A=15, deg B=14")
    print("norm: degree 31, exactly split over Q")
    print(f"certificate sha256: {certificate['certificate_sha256']}")


if __name__ == "__main__":
    main()
