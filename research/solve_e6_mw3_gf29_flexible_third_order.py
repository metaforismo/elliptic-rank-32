#!/usr/bin/env python3
"""Third-order obstruction system for the flexible GF(29) P3 #1 branch."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

import hensel_lift_e6_mw3_gf29 as linear
import hensel_lift_e6_mw3_gf29_flexible_slice as system
import lift_e6_mw3_gf29_flexible_branch as flexible_branch
import lift_e6_mw3_gf29_p3_1_branch as branch
import solve_e6_mw3_gf29_p3_1_third_order as third
from replay_e6_mw3_p2_split_cores import canonical_sha256


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "certificates" / "e6_mw3_gf29_independence.json"
DEFAULT_OUTPUT = ROOT / "certificates" / "e6_mw3_gf29_flexible_third_order.json"
PRIME = 29


def compute(
    input_path,
    groebner_order,
    line_parameter=18,
    first_parameters=None,
    p3_index=1,
):
    data = json.loads(input_path.read_text())
    start = system.starting_vector(data, p3_index)
    jacobian = system.jacobian_mod_prime(start)
    rank, first_particular, first_basis = linear.affine_solution_space(
        jacobian, branch.target_at(start, PRIME, system), PRIME
    )
    left_rank, _zero, left_kernel = linear.affine_solution_space(
        branch.transpose(jacobian), [0] * len(jacobian[0]), PRIME
    )
    if rank != left_rank:
        raise AssertionError("left/right Jacobian ranks disagree")
    if first_parameters is None:
        first_parameters = flexible_branch.first_parameters_on_line(line_parameter)
    else:
        first_parameters = [value % PRIME for value in first_parameters]
        if len(first_parameters) != len(first_basis):
            raise ValueError("first_parameters has the wrong dimension")
    first_correction = branch.affine_point(
        first_particular, first_basis, first_parameters
    )
    values2, signature = branch.next_compatibility(
        start, PRIME, first_correction, left_kernel, system
    )
    if any(signature):
        raise AssertionError("flexible first correction failed")
    second_space = third.correction_space_with_next(
        values2, PRIME**2, jacobian, left_kernel, system
    )
    if second_space is None:
        raise AssertionError("flexible branch cannot prepare digit 4")
    dimension = len(second_space["basis"])

    def values3(parameters):
        correction = branch.affine_point(
            second_space["particular"], second_space["basis"], parameters
        )
        candidate = [
            value + PRIME**2 * delta
            for value, delta in zip(values2, correction, strict=True)
        ]
        branch.target_at(candidate, PRIME**3, system)
        return candidate

    base3 = values3((0,) * dimension)
    _rank3, _eta_particular, eta_basis = linear.affine_solution_space(
        jacobian, branch.target_at(base3, PRIME**3, system), PRIME
    )

    def digit5_signature(values3_point, eta_parameters):
        _rank, particular, basis = linear.affine_solution_space(
            jacobian, branch.target_at(values3_point, PRIME**3, system), PRIME
        )
        if particular is None:
            raise AssertionError("digit-3 point cannot reach digit 4")
        eta = branch.affine_point(particular, basis, eta_parameters)
        _values4, result = branch.next_compatibility(
            values3_point, PRIME**3, eta, left_kernel, system
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
                (right - left) % PRIME
                for left, right in zip(constant0, value, strict=True)
            ]
        )
    eta_matrix = [
        [column[row] for column in eta_columns] for row in range(len(left_kernel))
    ]
    eta_rank, _eta_zero, eta_left_kernel = linear.affine_solution_space(
        branch.transpose(eta_matrix), [0] * len(eta_matrix[0]), PRIME
    )

    def obstruction(parameters):
        point = values3(parameters)
        constant = digit5_signature(point, [0] * len(eta_basis))
        return [branch.dot(vector, constant) for vector in eta_left_kernel]

    rows = third.interpolate_quadratic(obstruction, dimension)
    # Verify interpolation at deterministic nonsparse points.
    for sample in range(5):
        parameters = tuple(
            (3 + 5 * sample + 7 * index) % PRIME for index in range(dimension)
        )
        actual = obstruction(parameters)
        modeled_nonzero = any(
            third.evaluate_quadratic(row, parameters) for row in rows
        )
        if modeled_nonzero != any(actual):
            raise AssertionError("third-order quadratic interpolation failed")

    variables = sp.symbols(" ".join(f"z{index}" for index in range(dimension)))
    polynomials = []
    for row in rows:
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
    groebner = None
    solution_exists = None
    basis_strings = None
    if groebner_order:
        groebner = sp.groebner(
            polynomials, *variables, modulus=PRIME, order=groebner_order
        )
        solution_exists = not any(poly.as_expr() == 1 for poly in groebner.polys)
        basis_strings = [str(poly.as_expr()) for poly in groebner.polys]
    payload = {
        "schema": "elliptic-rank30/e6-mw3-flexible-third-order/v1",
        "field_prime": PRIME,
        "source_certificate": str(input_path),
        "source_certificate_sha256": data["certificate_sha256"],
        "P3_index": p3_index,
        "second_order_line_parameter": line_parameter % PRIME,
        "first_parameter_solution": first_parameters,
        "jacobian_rank": rank,
        "digit_3_parameter_dimension": dimension,
        "digit_5_eta_matrix_rank": eta_rank,
        "digit_5_eta_cokernel_dimension": len(eta_left_kernel),
        "third_order_obstruction_rows": rows,
        "groebner_order": groebner_order,
        "groebner_basis": basis_strings,
        "solution_exists": solution_exists,
        "claim_boundary": (
            "These are the exact digit-5 compatibility equations on every "
            "digit-3 parameter above one explicit flexible second-order branch."
        ),
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--groebner-order", choices=("grevlex", "lex"))
    parser.add_argument("--line-parameter", type=int, default=18)
    parser.add_argument("--p3", type=int, default=1)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    payload = compute(
        args.input,
        args.groebner_order,
        args.line_parameter,
        p3_index=args.p3,
    )
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "digit_3_parameter_dimension": payload[
                    "digit_3_parameter_dimension"
                ],
                "P3_index": payload["P3_index"],
                "second_order_line_parameter": payload[
                    "second_order_line_parameter"
                ],
                "third_order_equations": len(
                    payload["third_order_obstruction_rows"]
                ),
                "eta_cokernel_dimension": payload[
                    "digit_5_eta_cokernel_dimension"
                ],
                "solution_exists": payload["solution_exists"],
                "groebner_basis_size": len(payload["groebner_basis"] or []),
                "certificate_sha256": payload["certificate_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
