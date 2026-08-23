#!/usr/bin/env python3
"""Hensel system allowing the GF(29) slice coordinates to deform."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import hensel_lift_e6_mw3_gf29 as fixed
from replay_e6_mw3_p2_split_cores import canonical_sha256


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "certificates" / "e6_mw3_gf29_independence.json"
DEFAULT_OUTPUT = ROOT / "certificates" / "e6_mw3_gf29_flexible_slice_lifts.json"
PRIME = 29


VARIABLE_NAMES = (
    ["r0", "s0", "x1"]
    + [f"a{index}" for index in range(1, 6)]
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

    r0, s0, x1 = take(3)
    A = [-3 * s0**2] + take(5)
    B = [2 * s0**3] + take(7) + [1]
    P1_X = [r0**2 - 2 * s0, x1] + take(1)
    P1_Y = [r0 * (r0**2 - 3 * s0)] + take(3) + [1]
    pole = take(1)[0]
    P2_X = take(5)
    P2_Y = take(8)
    P3_X = [s0] + take(2)
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
    p1_residual = fixed.poly_sub(
        fixed.poly_sub(fixed.poly_pow(P1_Y, 2), fixed.poly_pow(P1_X, 3)),
        fixed.poly_add(fixed.poly_mul(A, P1_X), B),
    )
    result.extend(fixed.fixed_length(p1_residual, 9))
    z = [-pole, 1]
    p2_residual = fixed.poly_sub(
        fixed.poly_sub(fixed.poly_pow(P2_Y, 2), fixed.poly_pow(P2_X, 3)),
        fixed.poly_add(
            fixed.poly_mul(fixed.poly_mul(A, P2_X), fixed.poly_pow(z, 4)),
            fixed.poly_mul(B, fixed.poly_pow(z, 6)),
        ),
    )
    result.extend(fixed.fixed_length(p2_residual, 15))
    p3_residual = fixed.poly_sub(
        fixed.poly_sub(fixed.poly_pow(P3_Y, 2), fixed.poly_pow(P3_X, 3)),
        fixed.poly_add(fixed.poly_mul(A, P3_X), B),
    )
    result.extend(fixed.fixed_length(p3_residual, 9))
    discriminant_unit = fixed.poly_add(
        [4 * coefficient for coefficient in fixed.poly_pow(A, 3)],
        [27 * coefficient for coefficient in fixed.poly_pow(B, 2)],
    )
    result.extend(
        discriminant_unit[:4] + [0] * max(0, 4 - len(discriminant_unit))
    )
    current = discriminant_unit
    for _order in range(4):
        result.append(fixed.poly_eval(current, 1))
        current = fixed.poly_derivative(current)
    derivative_A = fixed.poly_derivative(A)
    derivative_B = fixed.poly_derivative(B)
    for location, node in ((lam, s_lam), (mu, s_mu)):
        result.extend(
            [
                fixed.poly_eval(A, location) + 3 * node**2,
                fixed.poly_eval(B, location) - 2 * node**3,
                fixed.poly_eval(derivative_B, location)
                + node * fixed.poly_eval(derivative_A, location),
            ]
        )
    if len(result) != 47:
        raise AssertionError("wrong equation count")
    return result


def starting_vector(data, p3_index):
    return [1, 2, 1] + fixed.starting_vector(data, p3_index)


def jacobian_mod_prime(values):
    base = equations(values)
    columns = []
    for variable in range(len(values)):
        displaced = list(values)
        displaced[variable] += PRIME
        evaluated = equations(displaced)
        columns.append(
            [
                ((right - left) // PRIME) % PRIME
                for left, right in zip(base, evaluated, strict=True)
            ]
        )
    return [
        [column[row] for column in columns] for row in range(len(base))
    ]


def lift_greedy(data, p3_index, digits):
    values = starting_vector(data, p3_index)
    if any(value % PRIME for value in equations(values)):
        raise AssertionError("flexible seed does not solve the system mod 29")
    jacobian = jacobian_mod_prime(values)
    rank, _zero = fixed.solve_full_column_rank(
        jacobian, [0] * len(jacobian), PRIME
    )
    modulus = PRIME
    record = {
        "P3_index": p3_index,
        "variables": len(values),
        "equations": 47,
        "jacobian_rank_mod_29": rank,
        "jacobian_nullity_mod_29": len(values) - rank,
        "lifted_digits": 1,
    }
    while record["lifted_digits"] < digits:
        residual = equations(values)
        target = [-(value // modulus) % PRIME for value in residual]
        _rank, correction = fixed.solve_full_column_rank(jacobian, target, PRIME)
        if correction is None:
            record.update(
                {
                    "status": "chosen_path_obstructed",
                    "first_failed_modulus": modulus * PRIME,
                }
            )
            return record
        values = [
            value + modulus * delta
            for value, delta in zip(values, correction, strict=True)
        ]
        modulus *= PRIME
        if any(value % modulus for value in equations(values)):
            raise AssertionError("greedy correction failed")
        record["lifted_digits"] += 1
    record.update(
        {
            "status": "lifted_on_greedy_path",
            "modulus": modulus,
            "residues": [value % modulus for value in values],
        }
    )
    return record


def compute(input_path, digits):
    data = json.loads(input_path.read_text())
    records = [lift_greedy(data, index, digits) for index in range(1, 5)]
    payload = {
        "schema": "elliptic-rank30/e6-mw3-flexible-slice-hensel/v1",
        "field_prime": PRIME,
        "source_certificate": str(input_path),
        "source_certificate_sha256": data["certificate_sha256"],
        "records": records,
        "claim_boundary": (
            "The slice coordinates r0,s0,x1 are allowed to deform. A greedy "
            "path fixes all free correction parameters to zero; failure of that "
            "path is not an exhaustive obstruction, and success is only a "
            "finite p-adic congruence."
        ),
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--digits", type=int, default=6)
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
                        key: record.get(key)
                        for key in (
                            "P3_index",
                            "jacobian_rank_mod_29",
                            "jacobian_nullity_mod_29",
                            "status",
                            "lifted_digits",
                            "first_failed_modulus",
                        )
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
