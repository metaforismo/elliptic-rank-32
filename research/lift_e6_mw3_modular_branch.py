#!/usr/bin/env python3
"""Follow an explicit flexible E6/MW3 modular branch to p-adic precision."""

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
DEFAULT_SECOND_ORDER = (
    ROOT / "certificates" / "e6_mw3_modular_second_order_gf23.json"
)


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


def target_at(values, modulus, prime):
    residual = system.equations(values)
    if any(value % modulus for value in residual):
        raise AssertionError("current vector is not a solution at its modulus")
    return [-(value // modulus) % prime for value in residual]


def next_compatibility(values, modulus, correction, left_kernel, prime):
    candidate = [
        value + modulus * delta
        for value, delta in zip(values, correction, strict=True)
    ]
    next_modulus = modulus * prime
    next_target = target_at(candidate, next_modulus, prime)
    return candidate, [
        dot(vector, next_target, prime) for vector in left_kernel
    ]


def choose_linear_lookahead(values, modulus, jacobian, left_kernel, prime):
    rank, particular, basis = linear.affine_solution_space(
        jacobian, target_at(values, modulus, prime), prime
    )
    if particular is None:
        return None, {"reason": "current Hensel equation is inconsistent"}
    candidate, constant = next_compatibility(
        values, modulus, particular, left_kernel, prime
    )
    columns = []
    for index in range(len(basis)):
        parameters = [0] * len(basis)
        parameters[index] = 1
        correction = affine_point(particular, basis, parameters, prime)
        _candidate, signature = next_compatibility(
            values, modulus, correction, left_kernel, prime
        )
        columns.append(
            [
                (right - left) % prime
                for left, right in zip(constant, signature, strict=True)
            ]
        )
    lookahead_matrix = [
        [column[row] for column in columns]
        for row in range(len(left_kernel))
    ]
    lookahead_rank, parameters = linear.solve_full_column_rank(
        lookahead_matrix,
        [(-value) % prime for value in constant],
        prime,
    )
    if parameters is None:
        return None, {
            "reason": "next Hensel equation cannot be made compatible",
            "jacobian_rank": rank,
            "free_parameters": len(basis),
            "lookahead_rank": lookahead_rank,
        }
    correction = affine_point(particular, basis, parameters, prime)
    candidate, signature = next_compatibility(
        values, modulus, correction, left_kernel, prime
    )
    if any(signature):
        raise AssertionError("linear look-ahead correction failed")
    return candidate, {
        "jacobian_rank": rank,
        "free_parameters": len(basis),
        "lookahead_rank": lookahead_rank,
        "lookahead_nullity": len(basis) - lookahead_rank,
        "chosen_parameters": parameters,
    }


def balanced(value, modulus):
    value %= modulus
    return value - modulus if value > modulus // 2 else value


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
        raise ValueError("target does not identify one second-order record")
    return matches[0]


def compute(
    slices_path,
    second_order_path,
    joint_index,
    p2_index,
    p3_index,
    first_parameters,
    digits,
):
    slices = json.loads(slices_path.read_text())
    second_order = json.loads(second_order_path.read_text())
    prime = slices["field_prime"]
    if second_order["source_slices_sha256"] != slices["certificate_sha256"]:
        raise ValueError("second-order certificate does not match slices")
    record = select_record(
        second_order["records"], joint_index, p2_index, p3_index
    )
    dimension = record["tangent_parameter_dimension"]
    first_parameters = [value % prime for value in first_parameters]
    if len(first_parameters) != dimension:
        raise ValueError(f"expected {dimension} first parameters")
    if any(
        evaluate_quadratic(row, first_parameters, prime)
        for row in record["obstruction_coefficient_rows"]
    ):
        raise ValueError("first parameters do not solve the p^3 obstruction")

    joint = slices["joint_candidates"][joint_index - 1]
    surface = joint["surface"]
    values = modular.starting_vector(
        surface,
        joint["coordinate_slice"],
        surface["geometric_P2_hits"][p2_index - 1],
        surface["P3_candidates"][p3_index - 1],
        prime,
    )
    jacobian = modular.jacobian_mod_prime(values, prime)
    rank, particular, basis = linear.affine_solution_space(
        jacobian, target_at(values, prime, prime), prime
    )
    if rank != record["jacobian_rank"] or particular is None:
        raise AssertionError("recomputed Jacobian data disagree")
    left_rank, zero, left_kernel = linear.affine_solution_space(
        transpose(jacobian), [0] * len(jacobian[0]), prime
    )
    if zero is None or left_rank != rank:
        raise AssertionError("left/right Jacobian ranks disagree")
    first_correction = affine_point(particular, basis, first_parameters, prime)
    values, signature = next_compatibility(
        values, prime, first_correction, left_kernel, prime
    )
    if any(signature):
        raise AssertionError("declared p^3 obstruction point is not compatible")
    modulus = prime**2
    steps = [
        {
            "from_digit": 1,
            "to_digit": 2,
            "chosen_parameters": first_parameters,
            "next_compatibility_zero": True,
        }
    ]
    status = "lifted"
    while modulus < prime**digits:
        candidate, metadata = choose_linear_lookahead(
            values, modulus, jacobian, left_kernel, prime
        )
        if candidate is None:
            status = "chosen_branch_obstructed"
            steps.append(metadata)
            break
        values = candidate
        previous = modulus
        modulus *= prime
        steps.append(
            {"modulus_before": previous, "modulus_after": modulus, **metadata}
        )
    if any(value % modulus for value in system.equations(values)):
        raise AssertionError("final p-adic congruence failed")
    lifted_digits = next(
        exponent for exponent in range(1, digits + 2) if prime**exponent == modulus
    )
    payload = {
        "schema": "elliptic-rank30/e6-mw3-modular-adic-branch/v1",
        "field_prime": prime,
        "source_slices": str(slices_path),
        "source_slices_sha256": slices["certificate_sha256"],
        "source_second_order": str(second_order_path),
        "source_second_order_sha256": second_order["certificate_sha256"],
        "joint_index": joint_index,
        "P2_index": p2_index,
        "P3_index": p3_index,
        "first_parameter_solution": first_parameters,
        "jacobian_rank": rank,
        "jacobian_nullity": len(values) - rank,
        "status": status,
        "lifted_digits": lifted_digits,
        "modulus": modulus,
        "variable_names": list(system.VARIABLE_NAMES),
        "residues": [value % modulus for value in values],
        "balanced_residues": [balanced(value, modulus) for value in values],
        "steps": steps,
        "claim_boundary": (
            "This is an exact compatible solution modulo the recorded power "
            "of p. Repeated deterministic linear look-ahead is evidence for "
            "a formal branch, not yet a proof of an infinite p-adic branch, "
            "an algebraic characteristic-zero point, or rational coefficients."
        ),
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


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
    return value % prime


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--slices", type=Path, default=DEFAULT_SLICES)
    parser.add_argument("--second-order", type=Path, default=DEFAULT_SECOND_ORDER)
    parser.add_argument("--joint", type=int, required=True)
    parser.add_argument("--p2", type=int, required=True)
    parser.add_argument("--p3", type=int, required=True)
    parser.add_argument("--parameters", required=True)
    parser.add_argument("--digits", type=int, default=20)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    parameters = [int(value) for value in args.parameters.split(",")]
    payload = compute(
        args.slices,
        args.second_order,
        args.joint,
        args.p2,
        args.p3,
        parameters,
        args.digits,
    )
    output = args.output or (
        ROOT
        / "certificates"
        / (
            f"e6_mw3_modular_adic_gf{payload['field_prime']}_"
            f"j{args.joint}_p2{args.p2}_p3{args.p3}.json"
        )
    )
    if not args.no_write:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "target": [args.joint, args.p2, args.p3],
                "status": payload["status"],
                "lifted_digits": payload["lifted_digits"],
                "modulus": payload["modulus"],
                "step_nullities": [
                    step.get("lookahead_nullity") for step in payload["steps"][1:]
                ],
                "certificate_sha256": payload["certificate_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
