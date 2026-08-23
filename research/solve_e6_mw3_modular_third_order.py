#!/usr/bin/env python3
"""Exact next p-adic obstruction above one modular second-order point."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import hensel_lift_e6_mw3_gf29 as linear
import hensel_lift_e6_mw3_gf29_flexible_slice as system
import hensel_lift_e6_mw3_modular_flexible as modular
import lift_e6_mw3_modular_branch as branch
import solve_e6_mw3_modular_second_order as quadratic
from replay_e6_mw3_p2_split_cores import canonical_sha256


ROOT = Path(__file__).resolve().parents[1]


def correction_space_with_next(values, modulus, jacobian, left_kernel, prime):
    """Corrections solving one digit and making the next digit solvable."""
    rank, particular, basis = linear.affine_solution_space(
        jacobian, branch.target_at(values, modulus, prime), prime
    )
    if particular is None:
        return None

    def signature(parameters):
        correction = branch.affine_point(particular, basis, parameters, prime)
        _candidate, result = branch.next_compatibility(
            values, modulus, correction, left_kernel, prime
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
                (right - left) % prime
                for left, right in zip(constant, value, strict=True)
            ]
        )
    matrix = [
        [column[row] for column in columns] for row in range(len(left_kernel))
    ]
    lookahead_rank, parameter_particular, parameter_basis = (
        linear.affine_solution_space(
            matrix, [(-value) % prime for value in constant], prime
        )
    )
    if parameter_particular is None:
        return None
    correction_particular = branch.affine_point(
        particular, basis, parameter_particular, prime
    )
    zero_correction = branch.affine_point(
        particular, basis, [0] * len(basis), prime
    )
    correction_basis = []
    for parameter_vector in parameter_basis:
        moved = branch.affine_point(
            particular, basis, parameter_vector, prime
        )
        correction_basis.append(
            [
                (right - left) % prime
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


def interpolate_quadratic(function, dimension, prime):
    zero = (0,) * dimension
    constant = function(zero)
    units = []
    doubles = []
    for index in range(dimension):
        unit = [0] * dimension
        double = [0] * dimension
        unit[index] = 1
        double[index] = 2
        units.append(function(tuple(unit)))
        doubles.append(function(tuple(double)))
    inverse_two = pow(2, -1, prime)
    diagonal = []
    linear_rows = []
    for index in range(dimension):
        square = [
            (doubles[index][row] - 2 * units[index][row] + constant[row])
            * inverse_two
            % prime
            for row in range(len(constant))
        ]
        diagonal.append(square)
        linear_rows.append(
            [
                (units[index][row] - constant[row] - square[row]) % prime
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
                % prime
                for row in range(len(constant))
            ]
    rows = []
    for row in range(len(constant)):
        coefficients = [constant[row]]
        coefficients.extend(linear_rows[index][row] for index in range(dimension))
        coefficients.extend(diagonal[index][row] for index in range(dimension))
        coefficients.extend(
            cross[(first, second)][row]
            for first in range(dimension)
            for second in range(first + 1, dimension)
        )
        if any(coefficients):
            rows.append(coefficients)
    return quadratic.independent_polynomial_rows(rows, prime)


def select_record(records, joint_index, p2_index, p3_index):
    matches = [
        record
        for record in records
        if (
            record["joint_index"],
            record["P2_index"],
            record["P3_index"],
        )
        == (joint_index, p2_index, p3_index)
    ]
    if len(matches) != 1:
        raise ValueError("target does not identify one record")
    return matches[0]


def compute(
    slices_path,
    second_order_path,
    joint_index,
    p2_index,
    p3_index,
    first_parameters,
):
    slices = json.loads(slices_path.read_text())
    second = json.loads(second_order_path.read_text())
    prime = slices["field_prime"]
    if second["source_slices_sha256"] != slices["certificate_sha256"]:
        raise ValueError("second-order certificate does not match slices")
    source_record = select_record(
        second["records"], joint_index, p2_index, p3_index
    )
    first_parameters = [value % prime for value in first_parameters]
    if len(first_parameters) != source_record["tangent_parameter_dimension"]:
        raise ValueError("wrong number of first parameters")
    if any(
        quadratic.evaluate_quadratic(row, first_parameters, prime)
        for row in source_record["obstruction_coefficient_rows"]
    ):
        raise ValueError("first parameters do not solve the p^3 obstruction")

    joint = slices["joint_candidates"][joint_index - 1]
    surface = joint["surface"]
    start = modular.starting_vector(
        surface,
        joint["coordinate_slice"],
        surface["geometric_P2_hits"][p2_index - 1],
        surface["P3_candidates"][p3_index - 1],
        prime,
    )
    jacobian = modular.jacobian_mod_prime(start, prime)
    rank, first_particular, first_basis = linear.affine_solution_space(
        jacobian, branch.target_at(start, prime, prime), prime
    )
    left_rank, zero, left_kernel = linear.affine_solution_space(
        branch.transpose(jacobian), [0] * len(jacobian[0]), prime
    )
    if first_particular is None or zero is None or left_rank != rank:
        raise AssertionError("invalid Jacobian spaces")
    first_correction = branch.affine_point(
        first_particular, first_basis, first_parameters, prime
    )
    values2, signature = branch.next_compatibility(
        start, prime, first_correction, left_kernel, prime
    )
    if any(signature):
        raise AssertionError("first obstruction root failed exact replay")

    second_space = correction_space_with_next(
        values2, prime**2, jacobian, left_kernel, prime
    )
    if second_space is None:
        raise AssertionError("first branch cannot reach p^4")
    dimension = len(second_space["basis"])

    def values3(parameters):
        correction = branch.affine_point(
            second_space["particular"],
            second_space["basis"],
            parameters,
            prime,
        )
        candidate = [
            value + prime**2 * delta
            for value, delta in zip(values2, correction, strict=True)
        ]
        branch.target_at(candidate, prime**3, prime)
        return candidate

    base3 = values3((0,) * dimension)
    _rank3, eta_particular, eta_basis = linear.affine_solution_space(
        jacobian, branch.target_at(base3, prime**3, prime), prime
    )
    if eta_particular is None:
        raise AssertionError("base digit-3 point cannot reach digit 4")

    def digit5_signature(values3_point, eta_parameters):
        _rank, particular, basis = linear.affine_solution_space(
            jacobian,
            branch.target_at(values3_point, prime**3, prime),
            prime,
        )
        if particular is None:
            raise AssertionError("digit-3 point cannot reach digit 4")
        eta = branch.affine_point(particular, basis, eta_parameters, prime)
        _values4, result = branch.next_compatibility(
            values3_point, prime**3, eta, left_kernel, prime
        )
        return result

    constant0 = digit5_signature(base3, [0] * len(eta_basis))
    eta_columns = []
    for index in range(len(eta_basis)):
        parameters = [0] * len(eta_basis)
        parameters[index] = 1
        value = digit5_signature(base3, parameters)
        eta_columns.append(
            [
                (right - left) % prime
                for left, right in zip(constant0, value, strict=True)
            ]
        )
    eta_matrix = [
        [column[row] for column in eta_columns]
        for row in range(len(left_kernel))
    ]
    eta_rank, eta_zero, eta_left_kernel = linear.affine_solution_space(
        branch.transpose(eta_matrix), [0] * len(eta_matrix[0]), prime
    )
    if eta_zero is None:
        raise AssertionError("eta left kernel computation failed")

    def obstruction(parameters):
        point = values3(parameters)
        constant = digit5_signature(point, [0] * len(eta_basis))
        return [
            branch.dot(vector, constant, prime) for vector in eta_left_kernel
        ]

    rows = interpolate_quadratic(obstruction, dimension, prime)
    verification_points = [
        tuple(
            (2 + 7 * sample + 5 * index) % prime
            for index in range(dimension)
        )
        for sample in range(7)
    ]
    for parameters in verification_points:
        actual = obstruction(parameters)
        modeled_nonzero = any(
            quadratic.evaluate_quadratic(row, parameters, prime)
            for row in rows
        )
        if modeled_nonzero != any(actual):
            raise AssertionError("higher obstruction interpolation failed")

    record = {
        "joint_index": joint_index,
        "P2_index": p2_index,
        "P3_index": p3_index,
        "first_parameter_solution": first_parameters,
        "jacobian_rank": rank,
        "tangent_parameter_dimension": dimension,
        "digit3_initial_free_dimension": second_space["initial_free_dimension"],
        "digit3_lookahead_rank": second_space["lookahead_rank"],
        "digit5_eta_matrix_rank": eta_rank,
        "digit5_eta_cokernel_dimension": len(eta_left_kernel),
        "linearly_independent_obstruction_rows": len(rows),
        "obstruction_coefficient_rows": rows,
        "digit3_correction_particular": second_space["particular"],
        "digit3_correction_basis": second_space["basis"],
        "eta_matrix": eta_matrix,
        "verification_points": [list(point) for point in verification_points],
    }
    payload = {
        "schema": "elliptic-rank30/e6-mw3-modular-third-order/v1",
        "field_prime": prime,
        "source_slices": str(slices_path),
        "source_slices_sha256": slices["certificate_sha256"],
        "source_second_order": str(second_order_path),
        "source_second_order_sha256": second["certificate_sha256"],
        "records": [record],
        "claim_boundary": (
            "The displayed quadratic ideal is the complete digit-5 "
            "compatibility condition over every digit-3 correction above "
            "one selected p^3-compatible branch. A root prepares a lift to "
            "p^5; it is not yet an infinite p-adic or rational lift."
        ),
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--slices", type=Path, required=True)
    parser.add_argument("--second-order", type=Path, required=True)
    parser.add_argument("--joint", type=int, required=True)
    parser.add_argument("--p2", type=int, required=True)
    parser.add_argument("--p3", type=int, required=True)
    parser.add_argument("--parameters", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--singular-output", type=Path)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    payload = compute(
        args.slices,
        args.second_order,
        args.joint,
        args.p2,
        args.p3,
        [int(value) for value in args.parameters.split(",")],
    )
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if args.singular_output:
        args.singular_output.parent.mkdir(parents=True, exist_ok=True)
        quadratic.write_singular(
            payload["records"][0],
            payload["field_prime"],
            args.singular_output,
            False,
        )
    record = payload["records"][0]
    print(
        json.dumps(
            {
                "digit3_parameter_dimension": record[
                    "tangent_parameter_dimension"
                ],
                "digit3_lookahead_rank": record["digit3_lookahead_rank"],
                "digit5_eta_matrix_rank": record["digit5_eta_matrix_rank"],
                "digit5_eta_cokernel_dimension": record[
                    "digit5_eta_cokernel_dimension"
                ],
                "independent_quadrics": record[
                    "linearly_independent_obstruction_rows"
                ],
                "certificate_sha256": payload["certificate_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
