#!/usr/bin/env python3
"""Resolve the next Hensel obstruction on the P3 #1 branch."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import hensel_lift_e6_mw3_gf29 as lift
import lift_e6_mw3_gf29_p3_1_branch as branch
from replay_e6_mw3_p2_split_cores import canonical_sha256


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "certificates" / "e6_mw3_gf29_independence.json"
DEFAULT_OUTPUT = ROOT / "certificates" / "e6_mw3_gf29_p3_1_third_order.json"
PRIME = 29


def correction_space_with_next(
    values, modulus, jacobian, left_kernel, system=lift
):
    """Corrections solving one digit and making the following digit solvable."""
    rank, particular, basis = lift.affine_solution_space(
        jacobian, branch.target_at(values, modulus, system), PRIME
    )
    if particular is None:
        return None

    def signature(parameters):
        correction = branch.affine_point(particular, basis, parameters)
        _candidate, result = branch.next_compatibility(
            values, modulus, correction, left_kernel, system
        )
        return result

    constant = signature([0] * len(basis))
    columns = []
    for index in range(len(basis)):
        parameters = [0] * len(basis)
        parameters[index] = 1
        value = signature(parameters)
        columns.append(
            [
                (right - left) % PRIME
                for left, right in zip(constant, value, strict=True)
            ]
        )
    matrix = [
        [column[row] for column in columns] for row in range(len(left_kernel))
    ]
    lookahead_rank, parameter_particular, parameter_basis = (
        lift.affine_solution_space(
            matrix, [(-value) % PRIME for value in constant], PRIME
        )
    )
    if parameter_particular is None:
        return None
    correction_particular = branch.affine_point(
        particular, basis, parameter_particular
    )
    correction_basis = []
    zero_correction = branch.affine_point(particular, basis, [0] * len(basis))
    for parameter_vector in parameter_basis:
        moved = branch.affine_point(particular, basis, parameter_vector)
        correction_basis.append(
            [
                (right - left) % PRIME
                for left, right in zip(zero_correction, moved, strict=True)
            ]
        )
    return {
        "jacobian_rank": rank,
        "initial_free_dimension": len(basis),
        "lookahead_rank": lookahead_rank,
        "particular": correction_particular,
        "basis": correction_basis,
    }


def initial_mod_29_squared(data, v, jacobian, left_kernel):
    values = lift.starting_vector(data, 1)
    _rank, particular, basis = lift.affine_solution_space(
        jacobian, branch.target_at(values, PRIME), PRIME
    )
    parameters = [
        (23 - 8 * v) % PRIME,
        24,
        v,
        (14 + 10 * v) % PRIME,
        (17 + 4 * v) % PRIME,
        25,
        21,
    ]
    correction = branch.affine_point(particular, basis, parameters)
    values, signature = branch.next_compatibility(
        values, PRIME, correction, left_kernel
    )
    if any(signature):
        raise AssertionError("invalid quadratic-branch parameter")
    return values


def interpolate_quadratic(function, dimension):
    zero = (0,) * dimension
    constant = function(zero)
    units = []
    doubles = []
    for index in range(dimension):
        unit = [0] * dimension
        unit[index] = 1
        double = [0] * dimension
        double[index] = 2
        units.append(function(tuple(unit)))
        doubles.append(function(tuple(double)))
    inverse_two = pow(2, -1, PRIME)
    diagonal = []
    linear = []
    for index in range(dimension):
        quadratic = [
            (doubles[index][row] - 2 * units[index][row] + constant[row])
            * inverse_two
            % PRIME
            for row in range(len(constant))
        ]
        diagonal.append(quadratic)
        linear.append(
            [
                (units[index][row] - constant[row] - quadratic[row]) % PRIME
                for row in range(len(constant))
            ]
        )
    cross = {}
    for first in range(dimension):
        for second in range(first + 1, dimension):
            parameters = [0] * dimension
            parameters[first] = parameters[second] = 1
            value = function(tuple(parameters))
            cross[(first, second)] = [
                (
                    value[row]
                    - units[first][row]
                    - units[second][row]
                    + constant[row]
                )
                % PRIME
                for row in range(len(constant))
            ]
    rows = []
    for row in range(len(constant)):
        coefficients = [constant[row]]
        coefficients.extend(linear[index][row] for index in range(dimension))
        coefficients.extend(diagonal[index][row] for index in range(dimension))
        coefficients.extend(
            cross[(first, second)][row]
            for first in range(dimension)
            for second in range(first + 1, dimension)
        )
        if any(coefficients):
            rows.append(coefficients)
    return rows


def evaluate_quadratic(coefficients, parameters):
    cursor = 0
    value = coefficients[cursor]
    cursor += 1
    for parameter in parameters:
        value += coefficients[cursor] * parameter
        cursor += 1
    for parameter in parameters:
        value += coefficients[cursor] * parameter * parameter
        cursor += 1
    for first in range(len(parameters)):
        for second in range(first + 1, len(parameters)):
            value += coefficients[cursor] * parameters[first] * parameters[second]
            cursor += 1
    return value % PRIME


def compute(input_path, branch_parameter):
    data = json.loads(input_path.read_text())
    start = lift.starting_vector(data, 1)
    jacobian = lift.jacobian_mod_prime(start)
    rank, _zero, left_kernel = lift.affine_solution_space(
        branch.transpose(jacobian), [0] * len(jacobian[0]), PRIME
    )
    values2 = initial_mod_29_squared(
        data, branch_parameter % PRIME, jacobian, left_kernel
    )
    second_space = correction_space_with_next(
        values2, PRIME**2, jacobian, left_kernel
    )
    if second_space is None:
        raise AssertionError("quadratic branch cannot reach 29^4")
    dimension = len(second_space["basis"])

    def values3(parameters):
        correction = branch.affine_point(
            second_space["particular"], second_space["basis"], parameters
        )
        candidate = [
            value + PRIME**2 * delta
            for value, delta in zip(values2, correction, strict=True)
        ]
        branch.target_at(candidate, PRIME**3)
        return candidate

    # At digit 3, the matrix describing how the seven free corrections affect
    # digit-5 compatibility is constant on this digit-3 parameter space.
    base3 = values3((0,) * dimension)
    rank3, eta_particular, eta_basis = lift.affine_solution_space(
        jacobian, branch.target_at(base3, PRIME**3), PRIME
    )
    if eta_particular is None:
        raise AssertionError("base digit-3 point unexpectedly cannot reach digit 4")

    def digit5_signature(values3_point, eta_parameters):
        _rank, particular, basis = lift.affine_solution_space(
            jacobian, branch.target_at(values3_point, PRIME**3), PRIME
        )
        eta = branch.affine_point(particular, basis, eta_parameters)
        _values4, signature = branch.next_compatibility(
            values3_point, PRIME**3, eta, left_kernel
        )
        return signature

    constant0 = digit5_signature(base3, [0] * len(eta_basis))
    eta_columns = []
    for index in range(len(eta_basis)):
        parameters = [0] * len(eta_basis)
        parameters[index] = 1
        value = digit5_signature(base3, parameters)
        eta_columns.append(
            [
                (right - left) % PRIME
                for left, right in zip(constant0, value, strict=True)
            ]
        )
    eta_matrix = [
        [column[row] for column in eta_columns] for row in range(len(left_kernel))
    ]
    eta_rank, _eta_zero, eta_left_kernel = lift.affine_solution_space(
        branch.transpose(eta_matrix), [0] * len(eta_matrix[0]), PRIME
    )

    def obstruction(parameters):
        values3_point = values3(parameters)
        constant = digit5_signature(values3_point, [0] * len(eta_basis))
        return [branch.dot(vector, constant) for vector in eta_left_kernel]

    obstruction_rows = interpolate_quadratic(obstruction, dimension)
    solutions = []
    for parameters in itertools.product(range(PRIME), repeat=dimension):
        if all(
            evaluate_quadratic(row, parameters) == 0
            for row in obstruction_rows
        ):
            solutions.append(list(parameters))
    if not solutions:
        status = "third_order_obstructed"
        chosen = None
        values4 = None
        eta_parameters = None
    else:
        status = "reaches_29_to_the_4_and_prepares_digit_5"
        chosen = solutions[0]
        chosen3 = values3(tuple(chosen))
        constant = digit5_signature(chosen3, [0] * len(eta_basis))
        _eta_rank, eta_parameters = lift.solve_full_column_rank(
            eta_matrix, [(-value) % PRIME for value in constant], PRIME
        )
        if eta_parameters is None:
            raise AssertionError("chosen obstruction root did not solve eta system")
        _rank, particular, basis = lift.affine_solution_space(
            jacobian, branch.target_at(chosen3, PRIME**3), PRIME
        )
        eta = branch.affine_point(particular, basis, eta_parameters)
        values4, signature = branch.next_compatibility(
            chosen3, PRIME**3, eta, left_kernel
        )
        if any(signature):
            raise AssertionError("chosen third-order root failed digit-5 compatibility")
        branch.target_at(values4, PRIME**4)
    payload = {
        "schema": "elliptic-rank30/e6-mw3-p3-1-third-order/v1",
        "field_prime": PRIME,
        "source_certificate": str(input_path),
        "source_certificate_sha256": data["certificate_sha256"],
        "quadratic_branch_parameter": branch_parameter % PRIME,
        "jacobian_rank": rank,
        "digit_3_parameter_dimension": dimension,
        "digit_5_eta_matrix_rank": eta_rank,
        "digit_5_eta_cokernel_dimension": len(eta_left_kernel),
        "third_order_obstruction_rows": obstruction_rows,
        "solutions": solutions,
        "solution_count": len(solutions),
        "chosen_solution": chosen,
        "chosen_eta_parameters": eta_parameters,
        "status": status,
        "modulus": PRIME**4 if values4 is not None else PRIME**3,
        "residues": [value % (PRIME**4) for value in values4]
        if values4 is not None
        else None,
        "claim_boundary": (
            "The exhaustive finite-field solve covers all residual digit-3 "
            "parameters above the selected quadratic branch. Reaching 29^4 "
            "with digit-5 compatibility is not yet an infinite p-adic or "
            "rational lift."
        ),
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--branch-parameter", type=int, default=0)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    payload = compute(args.input, args.branch_parameter)
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": payload["status"],
                "digit_3_parameter_dimension": payload[
                    "digit_3_parameter_dimension"
                ],
                "third_order_equations": len(
                    payload["third_order_obstruction_rows"]
                ),
                "solution_count": payload["solution_count"],
                "chosen_solution": payload["chosen_solution"],
                "certificate_sha256": payload["certificate_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
