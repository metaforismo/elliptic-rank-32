#!/usr/bin/env python3
"""Exact GF(31) replay of the finite E6/MW3 canonical-P2 search tranche.

This is an independent standard-Python port of the mathematical search
described in royvanrijn/jacobian-research at commit
e821b0c258bf635f38a1616806d09cb6b8025e387af.  The upstream reconstruction
reduces an E6+A3^2+A1^2 neighbor fibration to seven declared split-root cores
over GF(31).  For every reconstructed surface, the canonical one-denominator
section has

    X2(t) = C(t) + q0*t*(t-1)*(t-lambda)*(t-mu),

so only 27*31 = 837 (pole,q0) cases remain.

The script derives the four triangular P1 eliminations symbolically, verifies
every reconstructed P1 identity, enumerates every declared surface and every
canonical P2 case, and emits an exact JSON certificate.  It does not replay
the upstream 893,730-point proof that the seven input cores are complete, and
it does not by itself construct the characteristic-zero rank-17 fibration.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Iterable


PRIME = 31
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERTIFICATE = ROOT / "certificates" / "e6_mw3_p2_split_cores_gf31.json"
UPSTREAM_COMMIT = "e821b0c258bf635f38a1612ff8f2f769bf4ed04b"
UPSTREAM_URL = (
    "https://github.com/royvanrijn/jacobian-research/blob/"
    + UPSTREAM_COMMIT
    + "/elkies-k3/scripts/search_e6_mw3_p2_split_cores.sage"
)

CORES = [
    (2, 12, 25, 30),
    (4, 16, 6, 23),
    (14, 22, 30, 21),
    (15, 4, 17, 30),
    (17, 23, 4, 19),
    (21, 27, 0, 4),
    (23, 11, 8, 5),
]

Polynomial = list[int]  # coefficients in ascending degree order, modulo 31


def mod(value: int) -> int:
    return value % PRIME


def inverse(value: int) -> int:
    return pow(value % PRIME, -1, PRIME)


def divide(numerator: int, denominator: int) -> int:
    return numerator % PRIME * inverse(denominator) % PRIME


def trim(poly: Polynomial) -> Polynomial:
    result = [value % PRIME for value in poly]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result or [0]


def poly_add(left: Polynomial, right: Polynomial) -> Polynomial:
    width = max(len(left), len(right))
    return trim(
        [
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
            for index in range(width)
        ]
    )


def poly_neg(poly: Polynomial) -> Polynomial:
    return trim([-value for value in poly])


def poly_sub(left: Polynomial, right: Polynomial) -> Polynomial:
    return poly_add(left, poly_neg(right))


def poly_scale(poly: Polynomial, scalar: int) -> Polynomial:
    return trim([scalar * value for value in poly])


def poly_mul(left: Polynomial, right: Polynomial) -> Polynomial:
    result = [0] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            result[i + j] += left_value * right_value
    return trim(result)


def poly_pow(poly: Polynomial, exponent: int) -> Polynomial:
    if exponent < 0:
        raise ValueError("negative polynomial exponent")
    result = [1]
    base = trim(poly)
    power = exponent
    while power:
        if power & 1:
            result = poly_mul(result, base)
        base = poly_mul(base, base)
        power >>= 1
    return result


def poly_derivative(poly: Polynomial) -> Polynomial:
    if len(poly) <= 1:
        return [0]
    return trim([index * poly[index] for index in range(1, len(poly))])


def poly_eval(poly: Polynomial, point: int) -> int:
    value = 0
    for coefficient in reversed(poly):
        value = (value * point + coefficient) % PRIME
    return value


def poly_product(polys: Iterable[Polynomial]) -> Polynomial:
    result = [1]
    for poly in polys:
        result = poly_mul(result, poly)
    return result


def polynomial_sqrt_roots(poly: Polynomial) -> list[Polynomial]:
    """Return every polynomial square root in GF(31)[t]."""
    poly = trim(poly)
    if poly == [0]:
        return [[0]]
    degree = len(poly) - 1
    if degree % 2:
        return []
    root_degree = degree // 2
    leading_roots = [
        value for value in range(PRIME) if value * value % PRIME == poly[-1]
    ]
    roots: list[Polynomial] = []
    for leading in leading_roots:
        coefficients = [0] * (root_degree + 1)
        coefficients[root_degree] = leading
        for target_degree in range(degree - 1, root_degree - 1, -1):
            index = target_degree - root_degree
            known = 0
            for i in range(index + 1, root_degree + 1):
                j = target_degree - i
                if 0 <= j <= root_degree and j != index:
                    known += coefficients[i] * coefficients[j]
            coefficients[index] = divide(
                poly[target_degree] - known, 2 * leading
            )
        candidate = trim(coefficients)
        if poly_pow(candidate, 2) == poly:
            roots.append(candidate)
    return roots


def valuation_at(poly: Polynomial, point: int) -> int:
    multiplicity = 0
    current = trim(poly)
    while current != [0] and poly_eval(current, point) == 0:
        multiplicity += 1
        current = poly_derivative(current)
    return multiplicity


def interpolate_x(
    points: list[int], nodes: list[int], pole: int
) -> Polynomial:
    result = [0]
    for index, (point, node) in enumerate(zip(points, nodes)):
        basis = [1]
        denominator = 1
        for other_index, other in enumerate(points):
            if other_index == index:
                continue
            basis = poly_mul(basis, [-other, 1])
            denominator = denominator * (point - other) % PRIME
        value = node * (point - pole) ** 2 % PRIME
        result = poly_add(result, poly_scale(basis, divide(value, denominator)))
    return result


def derive_triangular_expressions(
    r0_value: int = 4,
    s0_value: int = 18,
    x1_value: int = 27,
):
    """Derive the triangular chart and its residual fixed-slice equations.

    Besides the four eliminated coordinates, return the coefficients of the
    section residual after imposing the public coordinate slice
    ``(r0,s0,x1)=(4,18,27)`` and the resulting affine-linear formula for y1.
    These residual equations are the input for an independent replay of the
    preceding 31^4 core exhaustion.
    """
    try:
        import sympy as sp
    except ImportError as error:
        raise RuntimeError("SymPy is required for symbolic elimination") from error

    t = sp.symbols("t")
    names = "a1 a2 a3 a4 b4 b5 s0 s1 x0 x1 y0 y1 y2"
    symbols = sp.symbols(names)
    (
        a1,
        a2,
        a3,
        a4,
        b4,
        b5,
        s0,
        s1,
        x0,
        x1,
        y0,
        y1,
        y2,
    ) = symbols

    a0 = -3 * s0**2
    b0 = 2 * s0**3
    b1 = -s0 * a1
    b2 = a1**2 / (12 * s0) - s0 * a2
    b3 = (a1**3 + 36 * a1 * a2 * s0**2 - 216 * a3 * s0**4) / (
        216 * s0**3
    )
    a5 = -3 * s1**2 - (a0 + a1 + a2 + a3 + a4)
    b6 = sp.symbols("b6")
    b7 = 2 * s1**3 - (b0 + b1 + b2 + b3 + b4 + b5 + b6) - 1
    derivative_condition = (
        b1 + 2 * b2 + 3 * b3 + 4 * b4 + 5 * b5 + 6 * b6 + 7 * b7 + 8
        + s1 * (a1 + 2 * a2 + 3 * a3 + 4 * a4 + 5 * a5)
    )
    cancel = lambda expression: sp.cancel(expression, modulus=PRIME)
    b6_expression = cancel(
        -derivative_condition.subs(b6, 0) / sp.diff(derivative_condition, b6)
    )
    b7_expression = cancel(b7.subs(b6, b6_expression))

    A = sum(
        coefficient * t**index
        for index, coefficient in enumerate([a0, a1, a2, a3, a4, a5])
    )
    B = sum(
        coefficient * t**index
        for index, coefficient in enumerate(
            [b0, b1, b2, b3, b4, b5, b6_expression, b7_expression, 1]
        )
    )
    x2 = s1 - x0 - x1
    y3 = -1 - y0 - y1 - y2
    X = x0 + x1 * t + x2 * t**2
    Y = y0 + y1 * t + y2 * t**2 + y3 * t**3 + t**4
    residual = sp.Poly(sp.together(Y**2 - X**3 - A * X - B), t)
    coefficients = {degree: residual.coeff_monomial(t**degree) for degree in range(8)}

    e7 = cancel(coefficients[7])
    b4_expression = cancel(-e7.subs(b4, 0) / sp.diff(e7, b4))

    e6 = cancel(coefficients[6].subs(b4, b4_expression))
    b5_expression = cancel(-e6.subs(b5, 0) / sp.diff(e6, b5))

    e5 = cancel(
        coefficients[5].subs(b4, b4_expression).subs(b5, b5_expression)
    )
    y2_expression = cancel(-e5.subs(y2, 0) / sp.diff(e5, y2))

    e4 = cancel(
        coefficients[4]
        .subs(b4, b4_expression)
        .subs(b5, b5_expression)
        .subs(y2, y2_expression)
    )
    a3_expression = cancel(-e4.subs(a3, 0) / sp.diff(e4, a3))

    expressions = {
        "a3": a3_expression,
        "y2": y2_expression,
        "b5": b5_expression,
        "b4": b4_expression,
    }
    if r0_value % PRIME == 0:
        raise ValueError("the coordinate slice requires r0 != 0")
    fixed_values = {
        s0: s0_value,
        x0: r0_value**2 - 2 * s0_value,
        x1: x1_value,
        y0: r0_value * (r0_value**2 - 3 * s0_value),
        y1: (
            a1
            + 3 * (r0_value**2 - s0_value) * x1_value
        )
        / (2 * r0_value),
    }
    triangular_values = {
        b4: b4_expression,
        b5: b5_expression,
        y2: y2_expression,
        a3: a3_expression,
    }
    core_residuals = {}
    for degree, coefficient in coefficients.items():
        specialized = coefficient
        for _ in range(8):
            previous = specialized
            specialized = specialized.subs(triangular_values).subs(fixed_values)
            if specialized == previous:
                break
        specialized = cancel(specialized)
        if specialized != 0:
            core_residuals[degree] = specialized
    if core_residuals:
        raise AssertionError(
            f"fixed-slice P1 residual survived triangular elimination: {core_residuals}"
        )

    discriminant = -16 * (4 * A**3 + 27 * B**2)
    core_equations = {}
    for derivative_order in (2, 3):
        equation = sp.diff(discriminant, t, derivative_order).subs(t, 1)
        for _ in range(8):
            previous = equation
            equation = equation.subs(triangular_values).subs(fixed_values)
            if equation == previous:
                break
        core_equations[f"D{derivative_order}_at_1"] = cancel(equation)
    symbol_map = {str(symbol): symbol for symbol in symbols}
    return expressions, symbol_map, core_equations


def derive_degenerate_branch_equations(
    r0_value: int,
    s0_value: int,
    x1_value: int,
):
    """Derive the exact branch where the P1_5 -> y2 pivot vanishes.

    On the coordinate slice, the third triangular pivot is

        2*(y0+y1).

    This function fixes the unique a1 for which it is zero, retains a3 and y2
    as genuine variables, and returns the four remaining equations P1_5,
    P1_4, Delta''(1), Delta'''(1), together with the still-valid b4/b5
    eliminations.  It is the input to a separate finite-field exhaustion and
    must not be replaced by cancellation of the generic rational formulas.
    """
    try:
        import sympy as sp
    except ImportError as error:
        raise RuntimeError("SymPy is required for symbolic elimination") from error

    if r0_value % PRIME == 0:
        raise ValueError("the coordinate slice requires r0 != 0")
    t = sp.symbols("t")
    names = "a1 a2 a3 a4 b4 b5 s0 s1 x0 x1 y0 y1 y2"
    symbols = sp.symbols(names)
    (
        a1,
        a2,
        a3,
        a4,
        b4,
        b5,
        s0,
        s1,
        x0,
        x1,
        y0,
        y1,
        y2,
    ) = symbols

    a0 = -3 * s0**2
    b0 = 2 * s0**3
    b1 = -s0 * a1
    b2 = a1**2 / (12 * s0) - s0 * a2
    b3 = (a1**3 + 36 * a1 * a2 * s0**2 - 216 * a3 * s0**4) / (
        216 * s0**3
    )
    a5 = -3 * s1**2 - (a0 + a1 + a2 + a3 + a4)
    b6 = sp.symbols("b6")
    b7 = 2 * s1**3 - (b0 + b1 + b2 + b3 + b4 + b5 + b6) - 1
    derivative_condition = (
        b1
        + 2 * b2
        + 3 * b3
        + 4 * b4
        + 5 * b5
        + 6 * b6
        + 7 * b7
        + 8
        + s1 * (a1 + 2 * a2 + 3 * a3 + 4 * a4 + 5 * a5)
    )
    cancel = lambda expression: sp.cancel(expression, modulus=PRIME)
    b6_expression = cancel(
        -derivative_condition.subs(b6, 0) / sp.diff(derivative_condition, b6)
    )
    b7_expression = cancel(b7.subs(b6, b6_expression))

    A = sum(
        coefficient * t**index
        for index, coefficient in enumerate([a0, a1, a2, a3, a4, a5])
    )
    B = sum(
        coefficient * t**index
        for index, coefficient in enumerate(
            [b0, b1, b2, b3, b4, b5, b6_expression, b7_expression, 1]
        )
    )
    X = x0 + x1 * t + (s1 - x0 - x1) * t**2
    Y = y0 + y1 * t + y2 * t**2 + (-1 - y0 - y1 - y2) * t**3 + t**4
    residual = sp.Poly(sp.together(Y**2 - X**3 - A * X - B), t)
    coefficients = {
        degree: residual.coeff_monomial(t**degree) for degree in range(8)
    }

    e7 = cancel(coefficients[7])
    b4_expression = cancel(-e7.subs(b4, 0) / sp.diff(e7, b4))
    e6 = cancel(coefficients[6].subs(b4, b4_expression))
    b5_expression = cancel(-e6.subs(b5, 0) / sp.diff(e6, b5))

    x0_value = (r0_value**2 - 2 * s0_value) % PRIME
    y0_value = r0_value * (r0_value**2 - 3 * s0_value) % PRIME
    constant_in_y1 = 3 * (r0_value**2 - s0_value) * x1_value % PRIME
    a1_value = (-constant_in_y1 - 2 * r0_value * y0_value) % PRIME
    y1_value = (-y0_value) % PRIME
    fixed_values = {
        a1: a1_value,
        s0: s0_value,
        x0: x0_value,
        x1: x1_value,
        y0: y0_value,
        y1: y1_value,
    }

    def specialize(expression):
        expression = expression.subs(b4, b4_expression).subs(b5, b5_expression)
        return cancel(expression.subs(fixed_values))

    discriminant = -16 * (4 * A**3 + 27 * B**2)
    equations = {
        "P1_5": specialize(coefficients[5]),
        "P1_4": specialize(coefficients[4]),
        "D2_at_1": specialize(sp.diff(discriminant, t, 2).subs(t, 1)),
        "D3_at_1": specialize(sp.diff(discriminant, t, 3).subs(t, 1)),
    }
    specialized_expressions = {
        "b5": cancel(b5_expression.subs(fixed_values)),
        "b4": cancel(
            b4_expression.subs(b5, b5_expression).subs(fixed_values)
        ),
    }
    symbol_map = {str(symbol): symbol for symbol in symbols}
    metadata = {
        "a1": a1_value,
        "x0": x0_value,
        "y0": y0_value,
        "y1": y1_value,
    }
    return specialized_expressions, symbol_map, equations, metadata


def evaluate_expression(expression, values: dict[str, int], symbol_map) -> int:
    substitutions = {
        symbol_map[name]: value for name, value in values.items() if name in symbol_map
    }
    result = expression
    for _ in range(8):
        previous = result
        result = result.subs(substitutions)
        if result == previous:
            break
    import sympy as sp

    result = sp.cancel(result, modulus=PRIME)
    numerator, denominator = result.as_numer_denom()
    if numerator.free_symbols or denominator.free_symbols:
        raise AssertionError(f"expression did not specialize: {result}")
    return divide(int(numerator), int(denominator))


def section_residual(
    A: Polynomial, B: Polynomial, X: Polynomial, Y: Polynomial
) -> Polynomial:
    return poly_sub(
        poly_sub(poly_pow(Y, 2), poly_pow(X, 3)),
        poly_add(poly_mul(A, X), B),
    )


def build_p1_candidate(
    core: tuple[int, int, int, int],
    a3: int,
    y2: int,
    b5: int,
    b4: int,
    coordinates: tuple[int, int, int] = (4, 18, 27),
) -> tuple[Polynomial, Polynomial, Polynomial, Polynomial]:
    """Build A,B,P1 from four residual triangular coordinates."""
    r0, s0, x1 = coordinates
    a1, a2, a4, s1 = core
    x0 = mod(r0**2 - 2 * s0)
    y0 = mod(r0 * (r0**2 - 3 * s0))
    y1 = divide(a1 + 3 * (r0**2 - s0) * x1, 2 * r0)
    a0 = mod(-3 * s0**2)
    b0 = mod(2 * s0**3)
    b1 = mod(-s0 * a1)
    b2 = mod(divide(a1**2, 12 * s0) - s0 * a2)
    b3 = divide(
        a1**3 + 36 * a1 * a2 * s0**2 - 216 * a3 * s0**4,
        216 * s0**3,
    )
    a5 = mod(-3 * s1**2 - (a0 + a1 + a2 + a3 + a4))
    A = trim([a0, a1, a2, a3, a4, a5])

    def derivative_condition(b6_value: int) -> tuple[int, int]:
        b7_value = mod(
            2 * s1**3 - (b0 + b1 + b2 + b3 + b4 + b5 + b6_value) - 1
        )
        coefficients = [b0, b1, b2, b3, b4, b5, b6_value, b7_value, 1]
        value = sum(index * coefficient for index, coefficient in enumerate(coefficients))
        value += s1 * sum(index * coefficient for index, coefficient in enumerate(A))
        return mod(value), b7_value

    condition_zero, _ = derivative_condition(0)
    condition_one, _ = derivative_condition(1)
    b6 = divide(-condition_zero, condition_one - condition_zero)
    condition, b7 = derivative_condition(b6)
    if condition != 0:
        raise AssertionError("failed to reconstruct b6")
    B = trim([b0, b1, b2, b3, b4, b5, b6, b7, 1])
    X1 = trim([x0, x1, mod(s1 - x0 - x1)])
    Y1 = trim([y0, y1, y2, mod(-1 - y0 - y1 - y2), 1])
    return A, B, X1, Y1


def brute_force_p1_candidate(
    core: tuple[int, int, int, int],
    coordinates: tuple[int, int, int] = (4, 18, 27),
) -> tuple[int, int, int, int, Polynomial, Polynomial, Polynomial, Polynomial]:
    """Finite-field fallback when a QQ chart cancellation is invalid mod 31.

    For fixed a3,y2,b5 the t^7 coefficient is affine-linear in b4, so the
    complete search needs only 31^3 candidate triples rather than 31^4.
    """
    solutions = []
    for a3 in range(PRIME):
        for y2 in range(PRIME):
            for b5 in range(PRIME):
                candidates = []
                residual_zero = section_residual(
                    *build_p1_candidate(core, a3, y2, b5, 0, coordinates)
                )
                residual_one = section_residual(
                    *build_p1_candidate(core, a3, y2, b5, 1, coordinates)
                )
                coefficient_zero = residual_zero[7] if len(residual_zero) > 7 else 0
                coefficient_one = residual_one[7] if len(residual_one) > 7 else 0
                slope = mod(coefficient_one - coefficient_zero)
                if slope:
                    candidates = [divide(-coefficient_zero, slope)]
                elif coefficient_zero == 0:
                    candidates = list(range(PRIME))
                for b4 in candidates:
                    A, B, X1, Y1 = build_p1_candidate(
                        core, a3, y2, b5, b4, coordinates
                    )
                    if section_residual(A, B, X1, Y1) == [0]:
                        solutions.append((a3, y2, b5, b4, A, B, X1, Y1))
    if len(solutions) != 1:
        raise AssertionError(
            f"core {core} has {len(solutions)} exact P1 candidates, expected 1"
        )
    return solutions[0]


def reconstruct_core(core, expressions, symbol_map) -> dict[str, object]:
    r0, s0, x1 = 4, 18, 27
    a1, a2, a4, s1 = core
    x0 = mod(r0**2 - 2 * s0)
    y0 = mod(r0 * (r0**2 - 3 * s0))
    y1 = divide(a1 + 3 * (r0**2 - s0) * x1, 2 * r0)
    values = {
        "a1": a1,
        "a2": a2,
        "a4": a4,
        "s0": s0,
        "s1": s1,
        "x0": x0,
        "x1": x1,
        "y0": y0,
        "y1": y1,
    }
    for key in ("a3", "y2", "b5", "b4"):
        values[key] = evaluate_expression(expressions[key], values, symbol_map)

    a3, y2, b5, b4 = (values[key] for key in ("a3", "y2", "b5", "b4"))
    A, B, X1, Y1 = build_p1_candidate(core, a3, y2, b5, b4)
    reconstruction_method = "triangular expressions specialized from Q"
    if section_residual(A, B, X1, Y1) != [0]:
        (
            a3,
            y2,
            b5,
            b4,
            A,
            B,
            X1,
            Y1,
        ) = brute_force_p1_candidate(core)
        reconstruction_method = "complete GF(31) triangular fallback"

    discriminant = poly_scale(
        poly_add(poly_scale(poly_pow(A, 3), 4), poly_scale(poly_pow(B, 2), 27)),
        -16,
    )
    node_map: dict[int, list[int]] = {}
    derivative_a = poly_derivative(A)
    derivative_b = poly_derivative(B)
    for point in range(PRIME):
        matches = []
        for node in range(PRIME):
            if (
                poly_eval(A, point) == mod(-3 * node**2)
                and poly_eval(B, point) == mod(2 * node**3)
                and mod(poly_eval(derivative_b, point) + node * poly_eval(derivative_a, point))
                == 0
            ):
                matches.append(node)
        if matches:
            node_map[point] = matches
    repeated = [
        point
        for point in range(PRIME)
        if point not in (0, 1) and valuation_at(discriminant, point) >= 2
    ]
    return {
        "core": list(core),
        "A": A,
        "B": B,
        "X1": X1,
        "Y1": Y1,
        "triangular_coordinates_a3_y2_b5_b4": [a3, y2, b5, b4],
        "reconstruction_method": reconstruction_method,
        "discriminant": discriminant,
        "node_map": node_map,
        "repeated": repeated,
        "multiplicity_at_0": valuation_at(discriminant, 0),
        "multiplicity_at_1": valuation_at(discriminant, 1),
    }


def search_p2(
    A: Polynomial,
    B: Polynomial,
    lam: int,
    mu: int,
    sl: int,
    sm: int,
    s0: int = 18,
) -> tuple[int, list[dict[str, object]]]:
    points = [0, 1, lam, mu]
    nodes_at_one = [
        node
        for node in range(PRIME)
        if poly_eval(A, 1) == mod(-3 * node**2)
        and poly_eval(B, 1) == mod(2 * node**3)
    ]
    if len(nodes_at_one) != 1:
        raise AssertionError(f"unexpected I4 node set {nodes_at_one}")
    nodes = [s0, nodes_at_one[0], sl, sm]
    fiber_product = poly_product([[-point, 1] for point in points])
    tested = 0
    hits: list[dict[str, object]] = []
    for pole in range(PRIME):
        if pole in points:
            continue
        constant = interpolate_x(points, nodes, pole)
        z = [-pole, 1]
        z4 = poly_pow(z, 4)
        z6 = poly_pow(z, 6)
        for q0 in range(PRIME):
            tested += 1
            X2 = poly_add(constant, poly_scale(fiber_product, q0))
            if poly_eval(X2, pole) == 0:
                continue
            H = poly_add(
                poly_pow(X2, 3),
                poly_add(
                    poly_mul(poly_mul(A, X2), z4),
                    poly_mul(B, z6),
                ),
            )
            for Y2 in polynomial_sqrt_roots(H):
                if any(poly_eval(Y2, point) != 0 for point in points):
                    continue
                hits.append({"pole": pole, "q0": q0, "X2": X2, "Y2": Y2})
    return tested, hits


def canonical_sign(poly: Polynomial) -> Polynomial:
    opposite = poly_neg(poly)
    return min(trim(poly), opposite)


def search_p3_polynomial(
    A: Polynomial,
    B: Polynomial,
    s1: int,
    lam: int,
    mu: int,
    sl: int,
    sm: int,
    s0: int = 18,
) -> tuple[int, list[dict[str, object]]]:
    """Enumerate the canonical polynomial P3 profile on a fixed surface.

    P3=(2,3,0,0,0) has P3.O=0.  On the normalized model this makes X a
    quadratic polynomial and the known control section has X(0)=s0.  We do
    not pre-classify its components at the other reducible fibers here.  The
    canonical target is nonidentity at I4@0 and identity at I4@1 and both I2
    fibers; exact post-search component tests are therefore required.
    """
    tested = 0
    candidates: dict[tuple[tuple[int, ...], tuple[int, ...]], dict[str, object]] = {}
    for x1 in range(PRIME):
        for x2 in range(PRIME):
            tested += 1
            X3 = trim([s0, x1, x2])
            H = poly_add(poly_pow(X3, 3), poly_add(poly_mul(A, X3), B))
            for Y3 in polynomial_sqrt_roots(H):
                if poly_eval(Y3, 0) != 0:
                    continue
                signed = canonical_sign(Y3)
                key = (tuple(X3), tuple(signed))
                candidates[key] = {"X3": X3, "Y3_up_to_sign": signed}
    return tested, list(candidates.values())


def p1_plus_p2_equals_negative_p3(
    X1: Polynomial,
    Y1: Polynomial,
    X2: Polynomial,
    Y2: Polynomial,
    pole: int,
    X3: Polynomial,
    Y3: Polynomial,
) -> bool:
    """Verify P1+P2=-P3 in GF(31)(t) by cleared polynomial identities."""
    z = [-pole, 1]
    z2 = poly_pow(z, 2)
    z3 = poly_pow(z, 3)
    slope_numerator = poly_sub(Y2, poly_mul(Y1, z3))
    x_difference_numerator = poly_sub(X2, poly_mul(X1, z2))

    # lambda^2 = x1+x2+x3, with
    # lambda = slope_numerator/(z*x_difference_numerator).
    x_identity = poly_sub(
        poly_pow(slope_numerator, 2),
        poly_mul(
            poly_pow(x_difference_numerator, 2),
            poly_add(poly_mul(poly_add(X1, X3), z2), X2),
        ),
    )
    if x_identity != [0]:
        return False

    # lambda*(x1-x3)-y1 = -y3.
    y_identity = poly_add(
        poly_mul(slope_numerator, poly_sub(X1, X3)),
        poly_mul(
            poly_mul(z, x_difference_numerator),
            poly_sub(Y3, Y1),
        ),
    )
    return y_identity == [0]


def canonical_sha256(payload: object) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def compute_certificate() -> dict[str, object]:
    expressions, symbol_map, _ = derive_triangular_expressions()
    core_records: list[dict[str, object]] = []
    total_surfaces = 0
    total_cases = 0
    total_hits = 0
    total_geometric_hits = 0
    total_dependent_geometric_hits = 0
    for core_index, core in enumerate(CORES, start=1):
        reconstructed = reconstruct_core(core, expressions, symbol_map)
        boundary = not (
            reconstructed["multiplicity_at_0"] == 4
            and reconstructed["multiplicity_at_1"] == 4
        )
        surface_records: list[dict[str, object]] = []
        if not boundary:
            repeated = reconstructed["repeated"]
            node_map = reconstructed["node_map"]
            for first_index, lam in enumerate(repeated):
                for mu in repeated[first_index + 1 :]:
                    for sl in node_map.get(lam, []):
                        for sm in node_map.get(mu, []):
                            tested, hits = search_p2(
                                reconstructed["A"],
                                reconstructed["B"],
                                lam,
                                mu,
                                sl,
                                sm,
                            )
                            total_surfaces += 1
                            total_cases += tested
                            total_hits += len(hits)
                            p3_tested, p3_candidates = search_p3_polynomial(
                                reconstructed["A"],
                                reconstructed["B"],
                                core[3],
                                lam,
                                mu,
                                sl,
                                sm,
                            )
                            geometric_p2: dict[
                                tuple[int, int, tuple[int, ...]], dict[str, object]
                            ] = {}
                            for hit in hits:
                                key = (
                                    hit["pole"],
                                    hit["q0"],
                                    tuple(hit["X2"]),
                                )
                                record = geometric_p2.setdefault(
                                    key,
                                    {
                                        "pole": hit["pole"],
                                        "q0": hit["q0"],
                                        "X2": hit["X2"],
                                        "Y2_up_to_sign": canonical_sign(hit["Y2"]),
                                        "signed_square_roots": 0,
                                    },
                                )
                                record["signed_square_roots"] += 1
                            geometric_records = list(geometric_p2.values())
                            total_geometric_hits += len(geometric_records)
                            for geometric_hit in geometric_records:
                                relations = []
                                for p3_candidate in p3_candidates:
                                    for sign_p2 in (1, -1):
                                        for sign_p3 in (1, -1):
                                            if p1_plus_p2_equals_negative_p3(
                                                reconstructed["X1"],
                                                reconstructed["Y1"],
                                                geometric_hit["X2"],
                                                poly_scale(
                                                    geometric_hit["Y2_up_to_sign"],
                                                    sign_p2,
                                                ),
                                                geometric_hit["pole"],
                                                p3_candidate["X3"],
                                                poly_scale(
                                                    p3_candidate["Y3_up_to_sign"],
                                                    sign_p3,
                                                ),
                                            ):
                                                relations.append(
                                                    {
                                                        "sign_P2": sign_p2,
                                                        "sign_P3": sign_p3,
                                                        "relation": "P1+sign_P2*P2+sign_P3*P3=O",
                                                    }
                                                )
                                geometric_hit["exact_P1_P2_P3_relations"] = relations
                                geometric_hit["independent_rank3_seed"] = not relations
                                if relations:
                                    total_dependent_geometric_hits += 1
                            surface_records.append(
                                {
                                    "lambda": lam,
                                    "mu": mu,
                                    "s_lambda": sl,
                                    "s_mu": sm,
                                    "tested_pole_q0_pairs": tested,
                                    "signed_P2_hits": hits,
                                    "geometric_P2_hits": geometric_records,
                                    "tested_P3_quadratic_x_polynomials": p3_tested,
                                    "canonical_P3_candidates": p3_candidates,
                                }
                            )
        core_records.append(
            {
                "core_index": core_index,
                "parameters_a1_a2_a4_s1": reconstructed["core"],
                "A_coefficients": reconstructed["A"],
                "B_coefficients": reconstructed["B"],
                "P1_X_coefficients": reconstructed["X1"],
                "P1_Y_coefficients": reconstructed["Y1"],
                "triangular_coordinates_a3_y2_b5_b4": reconstructed[
                    "triangular_coordinates_a3_y2_b5_b4"
                ],
                "reconstruction_method": reconstructed["reconstruction_method"],
                "P1_identity_verified": True,
                "delta_multiplicity_at_0": reconstructed["multiplicity_at_0"],
                "delta_multiplicity_at_1": reconstructed["multiplicity_at_1"],
                "repeated_fiber_positions": reconstructed["repeated"],
                "boundary": boundary,
                "surfaces": surface_records,
            }
        )

    payload: dict[str, object] = {
        "schema_version": 1,
        "claim": (
            "Every canonical P2 case on every surface reconstructed from the "
            "seven declared GF(31) split-root E6/P1 cores was tested exactly."
        ),
        "field": "GF(31)",
        "upstream": {
            "repository": "royvanrijn/jacobian-research",
            "commit": UPSTREAM_COMMIT,
            "source_script": UPSTREAM_URL,
            "input_scope": (
                "The seven-core list is taken from the upstream exhaustive "
                "core scan; this replay does not independently re-enumerate "
                "the 893,730 preceding core points."
            ),
        },
        "ansatz": {
            "fiber_configuration": "IV* + I4 + I4 + I2 + I2 + 4 I1",
            "neighbor_root_system": "E6 + A3^2 + A1^2",
            "generic_neighbor_mw_rank_target": 3,
            "P2": "X2=C+q0*t*(t-1)*(t-lambda)*(t-mu), pole r",
            "cases_per_surface": 837,
        },
        "core_records": core_records,
        "summary": {
            "declared_cores": len(CORES),
            "reconstructed_surfaces": total_surfaces,
            "tested_pole_q0_pairs": total_cases,
            "signed_canonical_P2_hits": total_hits,
            "geometric_canonical_P2_hits": total_geometric_hits,
            "dependent_geometric_P2_hits": total_dependent_geometric_hits,
            "independent_rank3_seed_hits": (
                total_geometric_hits - total_dependent_geometric_hits
            ),
        },
        "conditional_assumptions": [],
        "claim_boundary": (
            "A hit is only a finite-field section seed and must still be "
            "checked for Mordell-Weil dependence, lifted, rationally "
            "reconstructed, and verified in characteristic zero. No hit and "
            "no bounded negative result here is an upper bound on ranks over Q."
        ),
        "implementation": {
            "language": "Python plus SymPy for exact symbolic elimination",
            "finite_field_and_polynomial_arithmetic": "independent standard-Python implementation",
            "script": "research/replay_e6_mw3_p2_split_cores.py",
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
            raise AssertionError("committed certificate mismatch")
        print(f"matched {arguments.compare}")
    summary = certificate["summary"]
    print(
        "E6/MW3 GF(31): "
        f"cores={summary['declared_cores']}, "
        f"surfaces={summary['reconstructed_surfaces']}, "
        f"cases={summary['tested_pole_q0_pairs']}, "
        f"geometric_hits={summary['geometric_canonical_P2_hits']}, "
        f"independent_seeds={summary['independent_rank3_seed_hits']}"
    )
    print(f"certificate sha256: {certificate['certificate_sha256']}")


if __name__ == "__main__":
    main()
