#!/usr/bin/env python3
"""Analyze and formally lift the two-split P1/P2 incidence seeds.

The finite-field search records pairs of sections for which

    gcd(t*U1-U2, t*V1-V2)

has degree two.  This script replaces that procedural gcd condition by a
square incidence system with a monic quadratic witness h and quotient
polynomials Cx, Cy:

    t*U1-U2 = h*Cx,       t*V1-V2 = h*Cy.

Together with the coefficient equations for the two section identities this
gives 30 equations in 30 variables.  The Jacobian is evaluated by exact
forward automatic differentiation over F_p.  A nonsingular seed admits a
unique formal p-adic lift on this chart; the digit-by-digit Hensel lift below
is deterministic and deliberately makes no rational-lift claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    ROOT / "certificates" / "e8_a2_semistable_two_split_target_small_primes.json"
)

VARIABLES = [
    "lambda", "y", "z", "W",
    *[f"u1_{i}" for i in range(3)],
    *[f"v1_{i}" for i in range(5)],
    *[f"u2_{i}" for i in range(4)],
    *[f"v2_{i}" for i in range(6)],
    "h_0", "h_1",
    "cx_0", "cx_1",
    *[f"cy_{i}" for i in range(4)],
]

EQUATIONS = [
    *[f"P1_coefficient_{i}" for i in range(9)],
    *[f"P2_coefficient_{i}" for i in range(11)],
    *[f"x_factor_coefficient_{i}" for i in range(4)],
    *[f"y_factor_coefficient_{i}" for i in range(6)],
]


@dataclass(frozen=True)
class Jet:
    """A first-order dual number over a prime field."""

    value: int
    gradient: tuple[int, ...]
    prime: int

    @classmethod
    def constant(cls, value: int, dimension: int, prime: int) -> "Jet":
        return cls(value % prime, (0,) * dimension, prime)

    def _coerce(self, other: object) -> "Jet":
        if isinstance(other, Jet):
            if self.prime != other.prime or len(self.gradient) != len(other.gradient):
                raise ValueError("incompatible jets")
            return other
        if not isinstance(other, int):
            return NotImplemented
        return Jet.constant(other, len(self.gradient), self.prime)

    def __add__(self, other: object) -> "Jet":
        rhs = self._coerce(other)
        if rhs is NotImplemented:
            return NotImplemented
        p = self.prime
        return Jet(
            (self.value + rhs.value) % p,
            tuple((a + b) % p for a, b in zip(self.gradient, rhs.gradient)),
            p,
        )

    __radd__ = __add__

    def __neg__(self) -> "Jet":
        p = self.prime
        return Jet((-self.value) % p, tuple((-a) % p for a in self.gradient), p)

    def __sub__(self, other: object) -> "Jet":
        rhs = self._coerce(other)
        if rhs is NotImplemented:
            return NotImplemented
        return self + (-rhs)

    def __rsub__(self, other: object) -> "Jet":
        lhs = self._coerce(other)
        if lhs is NotImplemented:
            return NotImplemented
        return lhs - self

    def __mul__(self, other: object) -> "Jet":
        rhs = self._coerce(other)
        if rhs is NotImplemented:
            return NotImplemented
        p = self.prime
        return Jet(
            self.value * rhs.value % p,
            tuple(
                (self.value * b + rhs.value * a) % p
                for a, b in zip(self.gradient, rhs.gradient)
            ),
            p,
        )

    __rmul__ = __mul__

    def inverse(self) -> "Jet":
        if self.value == 0:
            raise ZeroDivisionError("jet denominator vanishes")
        p = self.prime
        reciprocal = pow(self.value, -1, p)
        scale = -reciprocal * reciprocal
        return Jet(
            reciprocal,
            tuple(scale * a % p for a in self.gradient),
            p,
        )

    def __truediv__(self, other: object) -> "Jet":
        rhs = self._coerce(other)
        if rhs is NotImplemented:
            return NotImplemented
        return self * rhs.inverse()

    def __rtruediv__(self, other: object) -> "Jet":
        lhs = self._coerce(other)
        if lhs is NotImplemented:
            return NotImplemented
        return lhs * self.inverse()


Scalar = int | Jet


def inverse(value: Scalar, modulus: int) -> Scalar:
    if isinstance(value, Jet):
        return value.inverse()
    return pow(value % modulus, -1, modulus)


def divide(numerator: Scalar, denominator: Scalar, modulus: int) -> Scalar:
    return numerator * inverse(denominator, modulus)


def poly_add(left: Sequence[Scalar], right: Sequence[Scalar]) -> list[Scalar]:
    zero: Scalar = 0
    width = max(len(left), len(right))
    return [
        (left[i] if i < len(left) else zero)
        + (right[i] if i < len(right) else zero)
        for i in range(width)
    ]


def poly_sub(left: Sequence[Scalar], right: Sequence[Scalar]) -> list[Scalar]:
    zero: Scalar = 0
    width = max(len(left), len(right))
    return [
        (left[i] if i < len(left) else zero)
        - (right[i] if i < len(right) else zero)
        for i in range(width)
    ]


def poly_scale(poly: Sequence[Scalar], scalar: Scalar) -> list[Scalar]:
    return [scalar * coefficient for coefficient in poly]


def poly_mul(left: Sequence[Scalar], right: Sequence[Scalar]) -> list[Scalar]:
    if not left or not right:
        return []
    result: list[Scalar] = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] = result[i + j] + a * b
    return result


def padded(poly: Sequence[Scalar], width: int) -> list[Scalar]:
    if len(poly) > width:
        tail = poly[width:]
        if any((entry.value if isinstance(entry, Jet) else entry) != 0 for entry in tail):
            raise AssertionError("unexpected polynomial degree")
    return [poly[i] if i < len(poly) else 0 for i in range(width)]


def surface_polynomials(values: Sequence[Scalar], modulus: int) -> dict[str, list[Scalar]]:
    lam, y_ratio, z_ratio, W = values[:4]
    inv3 = inverse(3, modulus)
    g_ratio = divide(y_ratio * y_ratio + lam - 1, lam, modulus)
    q_ratio = divide(1 - lam + lam * z_ratio, y_ratio, modulus)
    a0 = inv3
    a_at_1 = divide(z_ratio * z_ratio, 3 * g_ratio, modulus)
    a_at_lam = q_ratio * q_ratio * inv3
    delta_1 = a_at_1 - a0
    delta_lam = a_at_lam - a0
    a2 = divide(delta_lam - lam * delta_1, lam * (lam - 1), modulus)
    a1 = delta_1 - a2
    beta = [W * inv3, (z_ratio - 1) * W * inv3]
    gamma = [W * W, (g_ratio - 1) * W * W]
    return {
        "a": [a0, a1, a2],
        "beta": beta,
        "gamma": gamma,
        "g": [g_ratio],
        "q": [q_ratio],
    }


def evaluate_system(values: Sequence[Scalar], modulus: int) -> list[Scalar]:
    if len(values) != len(VARIABLES):
        raise ValueError(f"expected {len(VARIABLES)} variables")
    surface = surface_polynomials(values, modulus)
    a = surface["a"]
    beta = surface["beta"]
    gamma = surface["gamma"]
    lam = values[0]
    u1 = list(values[4:7])
    v1 = list(values[7:12])
    u2 = list(values[12:16])
    v2 = list(values[16:22])
    h = [values[22], values[23], 1]
    cx = list(values[24:26])
    cy = list(values[26:30])

    t = [0, 1]
    tm1 = [-1, 1]
    tlambda = [-lam, 1]
    tz = poly_mul(t, tlambda)
    T = poly_mul(t, tm1)

    p1_right = poly_add(
        poly_add(
            poly_mul(tz, poly_mul(poly_mul(u1, u1), u1)),
            poly_scale(poly_mul(a, poly_mul(u1, u1)), 3),
        ),
        poly_add(
            poly_scale(poly_mul(poly_mul(tm1, beta), u1), -6),
            poly_mul(poly_mul(tm1, tm1), gamma),
        ),
    )
    p2_right = poly_add(
        poly_add(
            poly_mul(tlambda, poly_mul(poly_mul(u2, u2), u2)),
            poly_scale(poly_mul(a, poly_mul(u2, u2)), 3),
        ),
        poly_add(
            poly_scale(poly_mul(poly_mul(T, beta), u2), -6),
            poly_mul(poly_mul(T, T), gamma),
        ),
    )
    p1_residual = padded(poly_sub(poly_mul(v1, v1), p1_right), 9)
    p2_residual = padded(poly_sub(poly_mul(v2, v2), p2_right), 11)

    dx = padded(poly_sub(poly_mul(t, u1), u2), 4)
    dy = padded(poly_sub(poly_mul(t, v1), v2), 6)
    x_factor = padded(poly_sub(dx, poly_mul(h, cx)), 4)
    y_factor = padded(poly_sub(dy, poly_mul(h, cy)), 6)
    result = [*p1_residual, *p2_residual, *x_factor, *y_factor]
    if len(result) != len(EQUATIONS):
        raise AssertionError("incidence system is not square")
    if isinstance(result[0], Jet):
        return result
    return [entry % modulus for entry in result]


def jacobian(values: Sequence[int], prime: int) -> tuple[list[int], list[list[int]]]:
    dimension = len(VARIABLES)
    jets: list[Jet] = []
    for index, value in enumerate(values):
        gradient = [0] * dimension
        gradient[index] = 1
        jets.append(Jet(value % prime, tuple(gradient), prime))
    residuals = evaluate_system(jets, prime)
    return (
        [entry.value for entry in residuals],
        [list(entry.gradient) for entry in residuals],
    )


def matrix_rank(matrix: Sequence[Sequence[int]], prime: int) -> int:
    work = [[entry % prime for entry in row] for row in matrix]
    if not work:
        return 0
    rows = len(work)
    columns = len(work[0])
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, rows) if work[row][column] != 0),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = pow(work[pivot_row][column], -1, prime)
        work[pivot_row] = [scale * entry % prime for entry in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or work[row][column] == 0:
                continue
            factor = work[row][column]
            work[row] = [
                (a - factor * b) % prime
                for a, b in zip(work[row], work[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def determinant_mod(matrix: Sequence[Sequence[int]], prime: int) -> int:
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("determinant requires a square matrix")
    work = [[entry % prime for entry in row] for row in matrix]
    determinant = 1
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column] != 0),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        pivot_value = work[column][column]
        determinant = determinant * pivot_value % prime
        reciprocal = pow(pivot_value, -1, prime)
        for row in range(column + 1, size):
            factor = work[row][column] * reciprocal % prime
            for entry in range(column, size):
                work[row][entry] = (
                    work[row][entry] - factor * work[column][entry]
                ) % prime
    return determinant % prime


def solve_square(matrix: Sequence[Sequence[int]], rhs: Sequence[int], prime: int) -> list[int]:
    size = len(matrix)
    work = [
        [entry % prime for entry in row] + [rhs[index] % prime]
        for index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column] != 0),
            None,
        )
        if pivot is None:
            raise ValueError("singular Jacobian")
        work[column], work[pivot] = work[pivot], work[column]
        scale = pow(work[column][column], -1, prime)
        work[column] = [scale * entry % prime for entry in work[column]]
        for row in range(size):
            if row == column or work[row][column] == 0:
                continue
            factor = work[row][column]
            work[row] = [
                (a - factor * b) % prime
                for a, b in zip(work[row], work[column])
            ]
    return [work[row][-1] for row in range(size)]


def trim(poly: Sequence[int], prime: int) -> list[int]:
    result = [entry % prime for entry in poly]
    while result and result[-1] == 0:
        result.pop()
    return result


def poly_divmod_mod(left: Sequence[int], right: Sequence[int], prime: int) -> tuple[list[int], list[int]]:
    numerator = trim(left, prime)
    denominator = trim(right, prime)
    if not denominator:
        raise ZeroDivisionError("zero polynomial")
    if len(numerator) < len(denominator):
        return [], numerator
    quotient = [0] * (len(numerator) - len(denominator) + 1)
    inverse_lead = pow(denominator[-1], -1, prime)
    while numerator and len(numerator) >= len(denominator):
        degree = len(numerator) - len(denominator)
        coefficient = numerator[-1] * inverse_lead % prime
        quotient[degree] = coefficient
        for index, entry in enumerate(denominator):
            numerator[index + degree] = (
                numerator[index + degree] - coefficient * entry
            ) % prime
        numerator = trim(numerator, prime)
    return trim(quotient, prime), numerator


def poly_gcd_mod(left: Sequence[int], right: Sequence[int], prime: int) -> list[int]:
    a = trim(left, prime)
    b = trim(right, prime)
    while b:
        _, remainder = poly_divmod_mod(a, b, prime)
        a, b = b, remainder
    if not a:
        return []
    scale = pow(a[-1], -1, prime)
    return [entry * scale % prime for entry in a]


def pair_polynomials(seed: dict[str, object], prime: int) -> tuple[list[int], list[int]]:
    p1 = seed["P1"]
    p2 = seed["P2"]
    u1 = p1["U"]
    v1 = p1["V"]
    u2 = p2["U"]
    v2 = p2["V"]
    dx = [
        -u2[0], u1[0] - u2[1], u1[1] - u2[2], u1[2] - u2[3]
    ]
    dy = [
        -v2[0], v1[0] - v2[1], v1[1] - v2[2],
        v1[2] - v2[3], v1[3] - v2[4], v1[4] - v2[5],
    ]
    return trim(dx, prime), trim(dy, prime)


def seed_vector(seed: dict[str, object]) -> tuple[list[int], dict[str, list[int]]]:
    prime = int(seed["prime"])
    dx, dy = pair_polynomials(seed, prime)
    h = poly_gcd_mod(dx, dy, prime)
    if len(h) != 3 or h[-1] != 1:
        raise AssertionError("seed gcd is not monic quadratic")
    cx, rx = poly_divmod_mod(dx, h, prime)
    cy, ry = poly_divmod_mod(dy, h, prime)
    if rx or ry:
        raise AssertionError("gcd witness does not divide both differences")
    if poly_gcd_mod(cx, cy, prime) != [1]:
        raise AssertionError("quadratic witness is not the exact gcd")
    v1 = seed["P1"]["V"]
    v2 = seed["P2"]["V"]
    y_sum = [
        v2[0], v1[0] + v2[1], v1[1] + v2[2],
        v1[2] + v2[3], v1[3] + v2[4], v1[4] + v2[5],
    ]
    if poly_gcd_mod(h, y_sum, prime) != [1]:
        raise AssertionError("same-sign intersection open gcd(h,tV1+V2)=1 failed")
    parameters = seed["parameters"]
    vector = [
        parameters["lambda"], parameters["y"], parameters["z"], parameters["W"],
        *seed["P1"]["U"], *seed["P1"]["V"],
        *seed["P2"]["U"], *seed["P2"]["V"],
        h[0], h[1],
        *padded(cx, 2), *padded(cy, 4),
    ]
    return [int(entry) % prime for entry in vector], {
        "h": h,
        "Cx": padded(cx, 2),
        "Cy": padded(cy, 4),
        "dx": dx,
        "dy": dy,
        "tV1_plus_V2": trim(y_sum, prime),
        "gcd_h_tV1_plus_V2": [1],
    }


def hensel_lift_smooth_slice(
    seed: Sequence[int], prime: int, digits: int
) -> tuple[list[int], list[dict[str, object]], int]:
    """Lift with cy_2 and cy_3 fixed, using the first 28 equations.

    The last two factor equations follow on the certified open set from the
    difference-of-cubic syzygy recorded in the output certificate.  They are
    nevertheless checked after every Hensel digit.
    """
    residual_mod_p, derivative = jacobian(seed, prime)
    if any(residual_mod_p):
        raise AssertionError("Hensel seed is not a solution modulo p")
    active = len(VARIABLES) - 2
    minor = [row[:active] for row in derivative[:active]]
    minor_determinant = determinant_mod(minor, prime)
    if minor_determinant == 0:
        raise ValueError("chosen transverse Jacobian minor is singular")
    current = [entry % prime for entry in seed]
    trace: list[dict[str, object]] = []
    for precision in range(1, digits + 1):
        modulus = prime**precision
        residual = evaluate_system(current, modulus)
        if any(residual):
            raise AssertionError(f"lift failed modulo p^{precision}")
        encoded = ",".join(str(entry) for entry in current).encode("ascii")
        trace.append({
            "precision": precision,
            "modulus": modulus,
            "solution_sha256": hashlib.sha256(encoded).hexdigest(),
        })
        if precision == digits:
            break
        next_modulus = prime ** (precision + 1)
        residual_next = evaluate_system(current, next_modulus)
        step = prime**precision
        if any(entry % step != 0 for entry in residual_next):
            raise AssertionError("current residual lost its certified precision")
        rhs = [-(entry // step) % prime for entry in residual_next[:active]]
        correction = solve_square(minor, rhs, prime)
        for index, delta in enumerate(correction):
            current[index] = (current[index] + step * delta) % next_modulus
        # cy_2 and cy_3 are the two fixed transverse-slice parameters.
        current[-2] %= next_modulus
        current[-1] %= next_modulus
        if any(evaluate_system(current, next_modulus)):
            raise AssertionError(
                "the two syzygetic equations did not follow at the next digit"
            )
    return current, trace, minor_determinant


def canonical_sha256(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def compute_certificate(input_path: Path, digits: int) -> dict[str, object]:
    source = json.loads(input_path.read_text(encoding="utf-8"))
    seeds = source["representative_pair_seeds"]["seeds"]
    records: list[dict[str, object]] = []
    lift_seed: tuple[list[int], int, int] | None = None
    for index, seed in enumerate(seeds):
        prime = int(seed["prime"])
        vector, witness = seed_vector(seed)
        residual, derivative = jacobian(vector, prime)
        if any(residual):
            raise AssertionError(f"seed {index} is not on the incidence scheme")
        rank = matrix_rank(derivative, prime)
        active_minor = [row[:28] for row in derivative[:28]]
        active_determinant = determinant_mod(active_minor, prime)
        record = {
            "seed_index": index,
            "prime": prime,
            "parameters": seed["parameters"],
            "intersection_witness": witness,
            "jacobian_rank": rank,
            "ambient_dimension": len(VARIABLES),
            "zariski_tangent_dimension": len(VARIABLES) - rank,
            "transverse_28_by_28_minor_determinant": active_determinant,
            "smooth_dimension_two_on_same_sign_open": (
                rank == 28 and active_determinant != 0
            ),
        }
        records.append(record)
        if lift_seed is None and prime == 11 and rank == 28 and active_determinant != 0:
            lift_seed = (vector, prime, index)

    lift: dict[str, object] | None = None
    if lift_seed is not None:
        vector, prime, index = lift_seed
        final, trace, minor_determinant = hensel_lift_smooth_slice(
            vector, prime, digits
        )
        modulus = prime**digits
        final_residual = evaluate_system(final, modulus)
        lift = {
            "seed_index": index,
            "prime": prime,
            "digits": digits,
            "modulus": modulus,
            "variables": dict(zip(VARIABLES, final)),
            "precision_trace": trace,
            "transverse_slice": {
                "fixed_variables": ["cy_2", "cy_3"],
                "fixed_values": [final[-2], final[-1]],
                "active_equations": EQUATIONS[:28],
                "active_variables": VARIABLES[:28],
                "jacobian_minor_determinant_mod_p": minor_determinant,
            },
            "all_30_residuals_zero_modulus": all(entry == 0 for entry in final_residual),
            "formal_claim": (
                "After fixing cy_2 and cy_3, the nonsingular transverse slice "
                "has a unique compatible solution over Z_%d.  The full local "
                "incidence scheme is smooth of dimension two on the recorded "
                "same-sign open set." % prime
            ),
            "rational_lift_claim": False,
        }

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "e8_a2_semistable_two_split_pair_incidence_lifts",
        "input": {
            "path": str(input_path.relative_to(ROOT)),
            "sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
            "representative_seeds": len(seeds),
        },
        "incidence_system": {
            "variables": VARIABLES,
            "equations": EQUATIONS,
            "variable_count": len(VARIABLES),
            "equation_count": len(EQUATIONS),
            "P1_identity": (
                "V1^2=t(t-lambda)U1^3+3aU1^2-6(t-1)beta U1+(t-1)^2 gamma"
            ),
            "P2_identity": (
                "V2^2=(t-lambda)U2^3+3aU2^2-6t(t-1)beta U2+t^2(t-1)^2 gamma"
            ),
            "intersection_witness": [
                "tU1-U2=(t^2+h_1*t+h_0)(cx_1*t+cx_0)",
                "tV1-V2=(t^2+h_1*t+h_0)(cy_3*t^3+cy_2*t^2+cy_1*t+cy_0)",
                "gcd(Cx,Cy)=1 is an open exact-degree-two condition",
            ],
            "expected_dimension": 2,
            "dimension_count": (
                "The 20 section coefficient equations leave expected dimension "
                "2. Factoring the cubic tU1-U2 introduces four variables and "
                "four equations. On gcd(h,tV1+V2)=1, the curve-difference "
                "syzygy forces h to divide tV1-V2, so its quotient equations "
                "introduce four variables and only define that quotient."
            ),
            "exact_curve_difference_syzygy": {
                "definitions": [
                    "x1=(t-lambda)*t*U1, x2=(t-lambda)*U2",
                    "y1=(t-lambda)*t*V1, y2=(t-lambda)*V2",
                    "D=(t-lambda)*t*(t-1)",
                    "K=(t-lambda)*((tU1)^2+(tU1)U2+U2^2)+3a(tU1+U2)-6t(t-1)beta",
                ],
                "identity": (
                    "(tV1-V2)(tV1+V2)=(tU1-U2)K"
                ),
                "consequence": (
                    "If tU1-U2=h*Cx and gcd(h,tV1+V2)=1, then h divides "
                    "tV1-V2; the final two displayed y-factor coefficients "
                    "are locally syzygetic."
                ),
            },
        },
        "seed_analysis": records,
        "formal_hensel_lift": lift,
        "P3_good_reduction_obstruction": (
            "The exhaustive source search found no target P3 on any recorded "
            "special fibre. Therefore no integral P3 in the same polynomial "
            "chart with unit node values can exist in these formal residue "
            "tubes: it would reduce to a searched P3. A P3 lift would have to "
            "enter through a bad-reduction or omitted boundary chart."
        ),
        "claim_boundary": {
            "formal_p_adic_pair_lift_found": lift is not None,
            "rational_pair_lift_found": False,
            "P3_lift_found": False,
            "rank31_curve_found": False,
        },
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--digits", type=int, default=8)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()
    if arguments.digits < 2:
        raise SystemExit("--digits must be at least 2")
    certificate = compute_certificate(arguments.input.resolve(), arguments.digits)
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
