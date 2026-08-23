#!/usr/bin/env python3
"""Exact p^3 compatibility equations for modular E6/MW3 seeds.

For every independently certified P1/P2/P3 triple, interpolate and verify the
quadratic obstruction on the full affine space of corrections modulo p.  The
result is the complete condition for a lift from p^2 to p^3 inside the 43
variable flexible-slice system.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import hensel_lift_e6_mw3_gf29 as linear
import hensel_lift_e6_mw3_gf29_flexible_slice as system
import hensel_lift_e6_mw3_modular_flexible as modular
from replay_e6_mw3_p2_split_cores import canonical_sha256


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SLICES = ROOT / "certificates" / "e6_mw3_modular_slices_gf23.json"
DEFAULT_INDEPENDENCE = (
    ROOT / "certificates" / "e6_mw3_modular_independence_gf23.json"
)
DEFAULT_OUTPUT = ROOT / "certificates" / "e6_mw3_modular_second_order_gf23.json"


def transpose(matrix):
    return [list(column) for column in zip(*matrix, strict=True)]


def dot(left, right, prime):
    return sum(a * b for a, b in zip(left, right, strict=True)) % prime


def affine_point(particular, basis, parameters, prime):
    result = list(particular)
    for parameter, vector in zip(parameters, basis, strict=True):
        result = [
            (left + parameter * right) % prime
            for left, right in zip(result, vector, strict=True)
        ]
    return result


def independent_polynomial_rows(rows, prime):
    if not rows:
        return []
    matrix = [[value % prime for value in row] for row in rows]
    pivot_row = 0
    for column in range(len(matrix[0])):
        selected = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if selected is None:
            continue
        matrix[pivot_row], matrix[selected] = matrix[selected], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], -1, prime)
        matrix[pivot_row] = [value * inverse % prime for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (left - factor * right) % prime
                for left, right in zip(matrix[row], matrix[pivot_row], strict=True)
            ]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return matrix[:pivot_row]


def evaluate_quadratic(row, parameters, prime):
    dimension = len(parameters)
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
    if cursor != len(row):
        raise AssertionError("quadratic row has wrong length")
    return value % prime


def obstruction_model(surface, coordinate_slice, p2_data, p3_data, prime):
    start = modular.starting_vector(
        surface, coordinate_slice, p2_data, p3_data, prime
    )
    jacobian = modular.jacobian_mod_prime(start, prime)
    first_target = [
        -(value // prime) % prime for value in system.equations(start)
    ]
    rank, particular, kernel_basis = linear.affine_solution_space(
        jacobian, first_target, prime
    )
    if particular is None:
        raise AssertionError("seed has no lift modulo p^2")
    left_rank, zero, left_kernel = linear.affine_solution_space(
        transpose(jacobian), [0] * len(jacobian[0]), prime
    )
    if zero is None or left_rank != rank:
        raise AssertionError("left/right Jacobian ranks disagree")

    def obstruction(parameters):
        correction = affine_point(particular, kernel_basis, parameters, prime)
        lifted = [
            value + prime * delta
            for value, delta in zip(start, correction, strict=True)
        ]
        residual = system.equations(lifted)
        if any(value % (prime**2) for value in residual):
            raise AssertionError("affine correction failed modulo p^2")
        target = [-(value // (prime**2)) % prime for value in residual]
        return [dot(vector, target, prime) for vector in left_kernel]

    dimension = len(kernel_basis)
    zero_parameters = (0,) * dimension
    constant = obstruction(zero_parameters)
    unit_values = []
    double_values = []
    for index in range(dimension):
        unit = [0] * dimension
        double = [0] * dimension
        unit[index] = 1
        double[index] = 2
        unit_values.append(obstruction(tuple(unit)))
        double_values.append(obstruction(tuple(double)))
    inverse_two = pow(2, -1, prime)
    diagonal = []
    linear_rows = []
    for index in range(dimension):
        quadratic = [
            (
                double_values[index][row]
                - 2 * unit_values[index][row]
                + constant[row]
            )
            * inverse_two
            % prime
            for row in range(len(left_kernel))
        ]
        diagonal.append(quadratic)
        linear_rows.append(
            [
                (unit_values[index][row] - constant[row] - quadratic[row])
                % prime
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
                % prime
                for row in range(len(left_kernel))
            ]

    coefficient_rows = []
    nonzero_row_indices = []
    for row in range(len(left_kernel)):
        coefficients = [constant[row]]
        coefficients.extend(linear_rows[index][row] for index in range(dimension))
        coefficients.extend(diagonal[index][row] for index in range(dimension))
        coefficients.extend(
            cross[(first, second)][row]
            for first in range(dimension)
            for second in range(first + 1, dimension)
        )
        if any(coefficients):
            coefficient_rows.append(coefficients)
            nonzero_row_indices.append(row)
    independent_rows = independent_polynomial_rows(coefficient_rows, prime)

    verification_points = [
        tuple(
            (3 * sample + 5 * variable + 1) % prime
            for variable in range(dimension)
        )
        for sample in range(7)
    ]
    for parameters in verification_points:
        actual = obstruction(parameters)
        modeled = [
            evaluate_quadratic(row, parameters, prime)
            for row in coefficient_rows
        ]
        expected = [actual[index] for index in nonzero_row_indices]
        if modeled != expected or any(
            actual[index]
            for index in range(len(actual))
            if index not in nonzero_row_indices
        ):
            raise AssertionError("quadratic obstruction interpolation failed")

    return {
        "jacobian_rank": rank,
        "tangent_parameter_dimension": dimension,
        "cokernel_dimension": len(left_kernel),
        "nonzero_obstruction_rows": len(coefficient_rows),
        "linearly_independent_obstruction_rows": len(independent_rows),
        "obstruction_coefficient_rows": independent_rows,
        "first_correction_particular": particular,
        "first_correction_kernel_basis": kernel_basis,
        "verification_points": [list(point) for point in verification_points],
    }


def singular_expression(row, dimension):
    terms = []
    cursor = 0
    if row[cursor]:
        terms.append(str(row[cursor]))
    cursor += 1
    for index in range(dimension):
        if row[cursor]:
            terms.append(f"{row[cursor]}*u{index}")
        cursor += 1
    for index in range(dimension):
        if row[cursor]:
            terms.append(f"{row[cursor]}*u{index}^2")
        cursor += 1
    for first in range(dimension):
        for second in range(first + 1, dimension):
            if row[cursor]:
                terms.append(f"{row[cursor]}*u{first}*u{second}")
            cursor += 1
    return "+".join(terms) or "0"


def write_singular(record, prime, path, field_equations):
    dimension = record["tangent_parameter_dimension"]
    variables = [f"u{index}" for index in range(dimension)]
    equations = [
        singular_expression(row, dimension)
        for row in record["obstruction_coefficient_rows"]
    ]
    if field_equations:
        equations.extend(f"{name}^{prime}-{name}" for name in variables)
    lines = [
        f"ring r={prime},({','.join(variables)}),dp;",
        "option(redSB);",
        "ideal I=" + ",\n".join(equations) + ";",
        "int start_timer=timer;",
        "ideal G=std(I);",
        'print("elapsed_ms="); print(timer-start_timer);',
        'print("basis_size="); print(size(G));',
        'print("dimension="); print(dim(G));',
        'print("independent_set="); print(indepSet(G));',
        'print("remainder_of_one="); print(reduce(1,G));',
    ]
    if field_equations:
        lines.append('print("vector_space_dimension="); print(vdim(G));')
    lines.append("quit;")
    path.write_text("\n".join(lines) + "\n")


def compute(slices_path, independence_path):
    slices = json.loads(slices_path.read_text())
    independence = json.loads(independence_path.read_text())
    prime = slices["field_prime"]
    if independence["source_certificate_sha256"] != slices["certificate_sha256"]:
        raise ValueError("independence certificate does not match slices")
    records = []
    for certified in independence["records"]:
        if not certified["independent"]:
            continue
        joint_index = certified["joint_index"]
        joint = slices["joint_candidates"][joint_index - 1]
        p2_index = certified["P2_index"]
        p3_index = certified["P3_index"]
        result = obstruction_model(
            joint["surface"],
            joint["coordinate_slice"],
            joint["surface"]["geometric_P2_hits"][p2_index - 1],
            joint["surface"]["P3_candidates"][p3_index - 1],
            prime,
        )
        records.append(
            {
                "joint_index": joint_index,
                "P2_index": p2_index,
                "P3_index": p3_index,
                "coordinate_slice": joint["coordinate_slice"],
                "core": joint["core"],
                **result,
            }
        )
    payload = {
        "schema": "elliptic-rank30/e6-mw3-modular-second-order/v1",
        "field_prime": prime,
        "source_slices": str(slices_path),
        "source_slices_sha256": slices["certificate_sha256"],
        "source_independence": str(independence_path),
        "source_independence_sha256": independence["certificate_sha256"],
        "records": records,
        "claim_boundary": (
            "Each quadratic ideal is the complete compatibility condition "
            "for lifting the corresponding flexible modular seed from p^2 "
            "to p^3. A solution is not yet an infinite p-adic or rational "
            "characteristic-zero lift."
        ),
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--slices", type=Path, default=DEFAULT_SLICES)
    parser.add_argument(
        "--independence", type=Path, default=DEFAULT_INDEPENDENCE
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--singular-dir", type=Path)
    parser.add_argument("--field-equations", action="store_true")
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    payload = compute(args.slices, args.independence)
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if args.singular_dir:
        args.singular_dir.mkdir(parents=True, exist_ok=True)
        for record in payload["records"]:
            label = (
                f"j{record['joint_index']}_p2{record['P2_index']}_"
                f"p3{record['P3_index']}"
            )
            write_singular(
                record,
                payload["field_prime"],
                args.singular_dir / f"{label}.sing",
                args.field_equations,
            )
    print(
        json.dumps(
            {
                "records": [
                    {
                        key: record[key]
                        for key in (
                            "joint_index",
                            "P2_index",
                            "P3_index",
                            "jacobian_rank",
                            "tangent_parameter_dimension",
                            "cokernel_dimension",
                            "linearly_independent_obstruction_rows",
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
