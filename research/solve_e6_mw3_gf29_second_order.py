#!/usr/bin/env python3
"""Derive second-order Hensel obstruction equations for the GF(29) seed."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import sympy as sp

import hensel_lift_e6_mw3_gf29 as lift
from replay_e6_mw3_p2_split_cores import canonical_sha256


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "certificates" / "e6_mw3_gf29_independence.json"
DEFAULT_OUTPUT = ROOT / "certificates" / "e6_mw3_gf29_second_order.json"
PRIME = 29


def transpose(matrix):
    return [list(column) for column in zip(*matrix, strict=True)]


def vector_add(left, right):
    return [(a + b) % PRIME for a, b in zip(left, right, strict=True)]


def vector_scale(scalar, vector):
    return [scalar * value % PRIME for value in vector]


def affine_point(particular, basis, parameters):
    result = list(particular)
    for parameter, vector in zip(parameters, basis, strict=True):
        result = vector_add(result, vector_scale(parameter, vector))
    return result


def dot(left, right):
    return sum(a * b for a, b in zip(left, right, strict=True)) % PRIME


def independent_polynomial_rows(rows):
    """Keep a row basis over GF(29), preserving explicit coefficient vectors."""
    if not rows:
        return []
    matrix = [[value % PRIME for value in row] for row in rows]
    kept = []
    pivot_row = 0
    for column in range(len(matrix[0])):
        selected = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if selected is None:
            continue
        matrix[pivot_row], matrix[selected] = matrix[selected], matrix[pivot_row]
        kept.append(selected)
        inverse = pow(matrix[pivot_row][column], -1, PRIME)
        matrix[pivot_row] = [value * inverse % PRIME for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (left - factor * right) % PRIME
                for left, right in zip(matrix[row], matrix[pivot_row], strict=True)
            ]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return matrix[:pivot_row]


def obstruction_model(data, p3_index, system=lift):
    start = system.starting_vector(data, p3_index)
    jacobian = system.jacobian_mod_prime(start)
    first_target = [-(value // PRIME) % PRIME for value in system.equations(start)]
    rank, particular, kernel_basis = lift.affine_solution_space(
        jacobian, first_target, PRIME
    )
    if particular is None:
        raise AssertionError("seed has no lift modulo 29^2")
    left_rank, zero, left_kernel = lift.affine_solution_space(
        transpose(jacobian), [0] * len(jacobian[0]), PRIME
    )
    if zero is None or left_rank != rank:
        raise AssertionError("left/right Jacobian ranks disagree")

    def obstruction(parameters):
        correction = affine_point(particular, kernel_basis, parameters)
        lifted = [
            value + PRIME * delta
            for value, delta in zip(start, correction, strict=True)
        ]
        residual = system.equations(lifted)
        if any(value % (PRIME**2) for value in residual):
            raise AssertionError("affine correction failed modulo 29^2")
        second_target = [-(value // (PRIME**2)) % PRIME for value in residual]
        return [dot(vector, second_target) for vector in left_kernel]

    dimension = len(kernel_basis)
    zero_parameters = (0,) * dimension
    constant = obstruction(zero_parameters)
    unit_values = []
    double_values = []
    for index in range(dimension):
        unit = [0] * dimension
        unit[index] = 1
        double = [0] * dimension
        double[index] = 2
        unit_values.append(obstruction(tuple(unit)))
        double_values.append(obstruction(tuple(double)))
    inverse_two = pow(2, -1, PRIME)
    diagonal = []
    linear = []
    for index in range(dimension):
        quadratic = [
            (double_values[index][row] - 2 * unit_values[index][row] + constant[row])
            * inverse_two
            % PRIME
            for row in range(len(left_kernel))
        ]
        diagonal.append(quadratic)
        linear.append(
            [
                (unit_values[index][row] - constant[row] - quadratic[row]) % PRIME
                for row in range(len(left_kernel))
            ]
        )
    cross = {}
    for first in range(dimension):
        for second in range(first + 1, dimension):
            parameters = [0] * dimension
            parameters[first] = parameters[second] = 1
            value = obstruction(tuple(parameters))
            cross[(first, second)] = [
                (
                    value[row]
                    - unit_values[first][row]
                    - unit_values[second][row]
                    + constant[row]
                )
                % PRIME
                for row in range(len(left_kernel))
            ]

    monomials = [(0, ())]
    monomials.extend((1, (index,)) for index in range(dimension))
    monomials.extend((2, (index, index)) for index in range(dimension))
    monomials.extend(
        (2, (first, second))
        for first in range(dimension)
        for second in range(first + 1, dimension)
    )
    coefficient_rows = []
    nonzero_row_indices = []
    for row in range(len(left_kernel)):
        coefficients = [constant[row]]
        coefficients.extend(linear[index][row] for index in range(dimension))
        coefficients.extend(diagonal[index][row] for index in range(dimension))
        coefficients.extend(
            cross[(first, second)][row]
            for first in range(dimension)
            for second in range(first + 1, dimension)
        )
        if any(coefficients):
            coefficient_rows.append(coefficients)
            nonzero_row_indices.append(row)
    independent_rows = independent_polynomial_rows(coefficient_rows)

    def evaluate_model(row, parameters):
        cursor = 0
        value = row[cursor]
        cursor += 1
        for parameter in parameters:
            value += row[cursor] * parameter
            cursor += 1
        for parameter in parameters:
            value += row[cursor] * parameter * parameter
            cursor += 1
        for first in range(dimension):
            for second in range(first + 1, dimension):
                value += row[cursor] * parameters[first] * parameters[second]
                cursor += 1
        return value % PRIME

    verification_points = [
        tuple((3 * index + 5 * variable + 1) % PRIME for variable in range(dimension))
        for index in range(5)
    ]
    for parameters in verification_points:
        actual = obstruction(parameters)
        modeled = [
            evaluate_model(row, parameters) for row in coefficient_rows
        ]
        expected = [actual[index] for index in nonzero_row_indices]
        if modeled != expected or any(
            actual[index]
            for index in range(len(actual))
            if index not in nonzero_row_indices
        ):
            raise AssertionError("quadratic obstruction interpolation failed")

    return {
        "P3_index": p3_index,
        "jacobian_rank": rank,
        "tangent_parameter_dimension": dimension,
        "cokernel_dimension": len(left_kernel),
        "nonzero_obstruction_rows": len(coefficient_rows),
        "linearly_independent_obstruction_rows": len(independent_rows),
        "monomial_order": [
            {"degree": degree, "indices": list(indices)} for degree, indices in monomials
        ],
        "obstruction_coefficient_rows": independent_rows,
        "first_correction_particular": particular,
        "first_correction_kernel_basis": kernel_basis,
    }


def sympy_polynomials(record):
    dimension = record["tangent_parameter_dimension"]
    variables = sp.symbols(" ".join(f"u{index}" for index in range(dimension)))
    if dimension == 1:
        variables = (variables,)
    polynomials = []
    for row in record["obstruction_coefficient_rows"]:
        cursor = 0
        expression = sp.Integer(row[cursor])
        cursor += 1
        for variable in variables:
            expression += row[cursor] * variable
            cursor += 1
        for variable in variables:
            expression += row[cursor] * variable**2
            cursor += 1
        for first in range(dimension):
            for second in range(first + 1, dimension):
                expression += row[cursor] * variables[first] * variables[second]
                cursor += 1
        polynomials.append(expression)
    return variables, polynomials


def add_groebner(record, order):
    variables, polynomials = sympy_polynomials(record)
    basis = sp.groebner(polynomials, *variables, modulus=PRIME, order=order)
    encoded = [str(poly.as_expr()) for poly in basis.polys]
    record["groebner_order"] = order
    record["groebner_basis"] = encoded
    record["second_order_solution_exists"] = not any(
        poly.as_expr() == 1 for poly in basis.polys
    )


def compute(input_path, p3_indices, groebner_order):
    data = json.loads(input_path.read_text())
    records = []
    for p3_index in p3_indices:
        record = obstruction_model(data, p3_index)
        if groebner_order:
            add_groebner(record, groebner_order)
        records.append(record)
    payload = {
        "schema": "elliptic-rank30/e6-mw3-second-order-obstruction/v1",
        "field_prime": PRIME,
        "source_certificate": str(input_path),
        "source_certificate_sha256": data["certificate_sha256"],
        "records": records,
        "claim_boundary": (
            "The quadratic systems are exactly the compatibility conditions "
            "for lifting the normalized modular seed from 29^2 to 29^3. A "
            "solution is not yet a full p-adic or rational lift."
        ),
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--p3", default="1,2,3,4")
    parser.add_argument("--groebner-order", choices=("grevlex", "lex"))
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    indices = [int(value) for value in args.p3.split(",")]
    payload = compute(args.input, indices, args.groebner_order)
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "records": [
                    {
                        "P3_index": record["P3_index"],
                        "tangent_dimension": record["tangent_parameter_dimension"],
                        "cokernel_dimension": record["cokernel_dimension"],
                        "independent_quadrics": record[
                            "linearly_independent_obstruction_rows"
                        ],
                        "second_order_solution_exists": record.get(
                            "second_order_solution_exists"
                        ),
                        "groebner_basis_size": len(record.get("groebner_basis", [])),
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
