#!/usr/bin/env python3
"""Test and construct p-adic lifts of the certified GF(29) MW-rank-3 seed.

The system keeps the normalized coordinate slice

    A(0)=-12, B(0)=16, B_8=1,
    X_1(0)=-3, [t]X_1=1, Y_1(0)=-5, [t^4]Y_1=1,
    X_3(0)=2, Y_3(0)=0,

and treats the remaining surface/section coefficients, the P2 pole, and the
two split-I2 locations/nodes as 40 variables.  Its 47 integral equations are
the three section identities, I4 discriminant vanishing at 0 and 1, and the
two nodal-I2 conditions.  A full-column-rank overdetermined Hensel step is
solved only when all 47 congruences are compatible.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from replay_e6_mw3_p2_split_cores import canonical_sha256


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "certificates" / "e6_mw3_gf29_independence.json"
DEFAULT_OUTPUT = ROOT / "certificates" / "e6_mw3_gf29_hensel_lifts.json"
PRIME = 29


def trim(poly):
    poly = list(poly)
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def poly_add(left, right):
    result = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        result[index] += value
    for index, value in enumerate(right):
        result[index] += value
    return trim(result)


def poly_sub(left, right):
    result = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        result[index] += value
    for index, value in enumerate(right):
        result[index] -= value
    return trim(result)


def poly_mul(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for first, a in enumerate(left):
        for second, b in enumerate(right):
            result[first + second] += a * b
    return trim(result)


def poly_pow(poly, exponent):
    result = [1]
    base = list(poly)
    while exponent:
        if exponent & 1:
            result = poly_mul(result, base)
        base = poly_mul(base, base)
        exponent >>= 1
    return result


def poly_eval(poly, value):
    result = 0
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


def poly_derivative(poly):
    return [index * coefficient for index, coefficient in enumerate(poly)][1:] or [0]


def fixed_length(poly, length):
    if len(poly) > length and any(poly[length:]):
        raise AssertionError(f"unexpected polynomial degree {len(poly) - 1}")
    return poly[:length] + [0] * max(0, length - len(poly))


VARIABLE_NAMES = (
    [f"a{index}" for index in range(1, 6)]
    + [f"b{index}" for index in range(1, 8)]
    + ["p1_x2"]
    + [f"p1_y{index}" for index in range(1, 4)]
    + ["p2_pole"]
    + [f"p2_x{index}" for index in range(5)]
    + [f"p2_y{index}" for index in range(8)]
    + ["p3_x1", "p3_x2"]
    + [f"p3_y{index}" for index in range(1, 5)]
    + ["lambda", "mu", "s_lambda", "s_mu"]
)


def unpack(values):
    if len(values) != len(VARIABLE_NAMES):
        raise ValueError("wrong number of variables")
    index = 0

    def take(count):
        nonlocal index
        result = values[index : index + count]
        index += count
        return result

    A = [-12] + take(5)
    B = [16] + take(7) + [1]
    P1_X = [-3, 1] + take(1)
    P1_Y = [-5] + take(3) + [1]
    pole = take(1)[0]
    P2_X = take(5)
    P2_Y = take(8)
    P3_X = [2] + take(2)
    P3_Y = [0] + take(4)
    lam, mu, s_lam, s_mu = take(4)
    if index != len(values):
        raise AssertionError("variable unpack did not consume the vector")
    return A, B, P1_X, P1_Y, pole, P2_X, P2_Y, P3_X, P3_Y, lam, mu, s_lam, s_mu


def equations(values):
    (
        A,
        B,
        P1_X,
        P1_Y,
        pole,
        P2_X,
        P2_Y,
        P3_X,
        P3_Y,
        lam,
        mu,
        s_lam,
        s_mu,
    ) = unpack(values)
    result = []
    p1_residual = poly_sub(
        poly_sub(poly_pow(P1_Y, 2), poly_pow(P1_X, 3)),
        poly_add(poly_mul(A, P1_X), B),
    )
    result.extend(fixed_length(p1_residual, 9))
    z = [-pole, 1]
    p2_residual = poly_sub(
        poly_sub(poly_pow(P2_Y, 2), poly_pow(P2_X, 3)),
        poly_add(
            poly_mul(poly_mul(A, P2_X), poly_pow(z, 4)),
            poly_mul(B, poly_pow(z, 6)),
        ),
    )
    result.extend(fixed_length(p2_residual, 15))
    p3_residual = poly_sub(
        poly_sub(poly_pow(P3_Y, 2), poly_pow(P3_X, 3)),
        poly_add(poly_mul(A, P3_X), B),
    )
    result.extend(fixed_length(p3_residual, 9))
    discriminant_unit = poly_add(
        [4 * coefficient for coefficient in poly_pow(A, 3)],
        [27 * coefficient for coefficient in poly_pow(B, 2)],
    )
    result.extend(discriminant_unit[:4] + [0] * max(0, 4 - len(discriminant_unit)))
    current = discriminant_unit
    for _order in range(4):
        result.append(poly_eval(current, 1))
        current = poly_derivative(current)
    derivative_A = poly_derivative(A)
    derivative_B = poly_derivative(B)
    for location, node in ((lam, s_lam), (mu, s_mu)):
        result.extend(
            [
                poly_eval(A, location) + 3 * node**2,
                poly_eval(B, location) - 2 * node**3,
                poly_eval(derivative_B, location)
                + node * poly_eval(derivative_A, location),
            ]
        )
    if len(result) != 47:
        raise AssertionError(f"expected 47 equations, got {len(result)}")
    return result


def starting_vector(data, p3_index):
    result = data["results"][p3_index - 1]
    surface_A = data["A"]
    surface_B = data["B"]
    P1 = data["P1"]
    P2 = data["P2"]
    P3 = result["P3"]
    values = (
        surface_A[1:6]
        + surface_B[1:8]
        + [P1["X"][2]]
        + P1["Y"][1:4]
        + [P2["pole"]]
        + P2["X2"]
        + P2["Y2"]
        + P3["X3"][1:3]
        + P3["Y3_up_to_sign"][1:5]
        + [6, 21, 22, 27]
    )
    if len(values) != len(VARIABLE_NAMES):
        raise AssertionError("starting vector has wrong length")
    if any(value % PRIME for value in equations(values)):
        raise AssertionError("starting vector does not solve the system mod 29")
    return values


def jacobian_mod_prime(values):
    base = equations(values)
    matrix = []
    columns = []
    for variable in range(len(values)):
        displaced = list(values)
        displaced[variable] += PRIME
        evaluated = equations(displaced)
        column = [((right - left) // PRIME) % PRIME for left, right in zip(base, evaluated)]
        columns.append(column)
    for equation_index in range(len(base)):
        matrix.append([column[equation_index] for column in columns])
    return matrix


def row_reduce(matrix, target, prime):
    rows = [
        [value % prime for value in row] + [target_value % prime]
        for row, target_value in zip(matrix, target, strict=True)
    ]
    row_count = len(rows)
    column_count = len(matrix[0])
    pivot_rows = []
    pivot_row = 0
    for column in range(column_count):
        selected = next(
            (row for row in range(pivot_row, row_count) if rows[row][column]),
            None,
        )
        if selected is None:
            continue
        rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
        scalar = pow(rows[pivot_row][column], -1, prime)
        rows[pivot_row] = [value * scalar % prime for value in rows[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or rows[row][column] == 0:
                continue
            factor = rows[row][column]
            rows[row] = [
                (left - factor * right) % prime
                for left, right in zip(rows[row], rows[pivot_row], strict=True)
            ]
        pivot_rows.append((pivot_row, column))
        pivot_row += 1
    inconsistent_rows = [
        row for row in rows if not any(row[:column_count]) and row[column_count]
    ]
    return rows, pivot_rows, inconsistent_rows


def solve_full_column_rank(matrix, target, prime):
    rows, pivot_rows, inconsistent_rows = row_reduce(matrix, target, prime)
    column_count = len(matrix[0])
    rank = len(pivot_rows)
    if inconsistent_rows:
        return rank, None
    solution = [0] * column_count
    for row, column in pivot_rows:
        solution[column] = rows[row][column_count]
    return rank, solution


def affine_solution_space(matrix, target, prime):
    """Return rank, one solution, and a basis for the homogeneous kernel."""
    rows, pivot_rows, inconsistent_rows = row_reduce(matrix, target, prime)
    if inconsistent_rows:
        return len(pivot_rows), None, []
    column_count = len(matrix[0])
    pivot_columns = {column for _row, column in pivot_rows}
    free_columns = [
        column for column in range(column_count) if column not in pivot_columns
    ]
    particular = [0] * column_count
    for row, column in pivot_rows:
        particular[column] = rows[row][column_count]
    basis = []
    for free in free_columns:
        vector = [0] * column_count
        vector[free] = 1
        for row, column in pivot_rows:
            vector[column] = (-rows[row][free]) % prime
        basis.append(vector)
    return len(pivot_rows), particular, basis


def lift_one(data, p3_index, digits):
    values = starting_vector(data, p3_index)
    jacobian = jacobian_mod_prime(values)
    rank, homogeneous = solve_full_column_rank(
        jacobian, [0] * len(jacobian), PRIME
    )
    if homogeneous is None:
        raise AssertionError("homogeneous Jacobian system cannot be inconsistent")
    record = {
        "P3_index": p3_index,
        "variables": list(VARIABLE_NAMES),
        "equations": 47,
        "jacobian_rank_mod_29": rank,
        "jacobian_nullity_mod_29": len(values) - rank,
        "requested_29_adic_digits": digits,
        "lifted_digits": 1,
        "status": None,
    }
    modulus = PRIME
    for digit in range(1, digits):
        residual = equations(values)
        if any(value % modulus for value in residual):
            raise AssertionError("current vector is not a solution at its modulus")
        target = [-(value // modulus) % PRIME for value in residual]
        step_rank, correction = solve_full_column_rank(jacobian, target, PRIME)
        if step_rank != rank:
            raise AssertionError("Jacobian rank changed unexpectedly")
        if correction is None:
            record.update(
                {
                    "status": "obstructed",
                    "first_obstructed_modulus": modulus * PRIME,
                    "compatibility_failed_at_digit": digit + 1,
                }
            )
            return record
        values = [
            value + modulus * delta for value, delta in zip(values, correction, strict=True)
        ]
        modulus *= PRIME
        if any(value % modulus for value in equations(values)):
            raise AssertionError("Hensel correction did not lift all equations")
        record["lifted_digits"] = digit + 1
    record.update(
        {
            "status": "lifted",
            "modulus": modulus,
            "residues": [value % modulus for value in values],
        }
    )
    return record


def compute(input_path, digits):
    data = json.loads(input_path.read_text())
    records = [lift_one(data, index, digits) for index in range(1, 5)]
    payload = {
        "schema": "elliptic-rank30/e6-mw3-overdetermined-hensel/v1",
        "field_prime": PRIME,
        "source_certificate": str(input_path),
        "source_certificate_sha256": data["certificate_sha256"],
        "normalization": {
            "A0": -12,
            "B0": 16,
            "B8": 1,
            "P1_X0": -3,
            "P1_X1": 1,
            "P1_Y0": -5,
            "P1_Y4": 1,
            "P3_X0": 2,
            "P3_Y0": 0,
        },
        "records": records,
        "claim_boundary": (
            "A lifted record is a compatible p-adic solution of the stated "
            "normalized surface/section system, not yet a rational model. An "
            "obstruction excludes only this normalized coefficient lift above "
            "the listed modular point."
        ),
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--digits", type=int, default=4)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    payload = compute(args.input, args.digits)
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "records": [
                    {
                        "P3_index": record["P3_index"],
                        "jacobian_rank_mod_29": record["jacobian_rank_mod_29"],
                        "status": record["status"],
                        "lifted_digits": record["lifted_digits"],
                        "first_obstructed_modulus": record.get(
                            "first_obstructed_modulus"
                        ),
                    }
                    for record in payload["records"]
                ],
                "certificate_sha256": payload["certificate_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
