#!/usr/bin/env python3
"""Follow the one-dimensional P3 #1 branch to high 29-adic precision."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import hensel_lift_e6_mw3_gf29 as lift
from replay_e6_mw3_p2_split_cores import canonical_sha256


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "certificates" / "e6_mw3_gf29_independence.json"
DEFAULT_OUTPUT = ROOT / "certificates" / "e6_mw3_gf29_p3_1_adic_branch.json"
PRIME = 29


def transpose(matrix):
    return [list(column) for column in zip(*matrix, strict=True)]


def dot(left, right):
    return sum(a * b for a, b in zip(left, right, strict=True)) % PRIME


def affine_point(particular, basis, parameters):
    result = list(particular)
    for parameter, vector in zip(parameters, basis, strict=True):
        result = [
            (left + parameter * right) % PRIME
            for left, right in zip(result, vector, strict=True)
        ]
    return result


def target_at(values, modulus, system=lift):
    residual = system.equations(values)
    if any(value % modulus for value in residual):
        raise AssertionError("current vector is not a solution at its modulus")
    return [-(value // modulus) % PRIME for value in residual]


def next_compatibility(values, modulus, correction, left_kernel, system=lift):
    candidate = [
        value + modulus * delta
        for value, delta in zip(values, correction, strict=True)
    ]
    next_modulus = modulus * PRIME
    next_target = target_at(candidate, next_modulus, system)
    return candidate, [dot(vector, next_target) for vector in left_kernel]


def choose_linear_lookahead(values, modulus, jacobian, left_kernel, system=lift):
    rank, particular, basis = lift.affine_solution_space(
        jacobian, target_at(values, modulus, system), PRIME
    )
    if particular is None:
        return None, {"reason": "current Hensel equation is inconsistent"}
    candidate, constant = next_compatibility(
        values, modulus, particular, left_kernel, system
    )
    columns = []
    for index in range(len(basis)):
        parameters = [0] * len(basis)
        parameters[index] = 1
        correction = affine_point(particular, basis, parameters)
        _candidate, signature = next_compatibility(
            values, modulus, correction, left_kernel, system
        )
        columns.append(
            [
                (right - left) % PRIME
                for left, right in zip(constant, signature, strict=True)
            ]
        )
    lookahead_matrix = [
        [column[row] for column in columns] for row in range(len(left_kernel))
    ]
    lookahead_rank, parameters = lift.solve_full_column_rank(
        lookahead_matrix, [(-value) % PRIME for value in constant], PRIME
    )
    if parameters is None:
        return None, {
            "reason": "next Hensel equation cannot be made compatible",
            "jacobian_rank": rank,
            "free_parameters": len(basis),
            "lookahead_rank": lookahead_rank,
        }
    correction = affine_point(particular, basis, parameters)
    candidate, signature = next_compatibility(
        values, modulus, correction, left_kernel, system
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


def compute(input_path, branch_parameter, digits):
    data = json.loads(input_path.read_text())
    values = lift.starting_vector(data, 1)
    jacobian = lift.jacobian_mod_prime(values)
    rank, first_particular, first_basis = lift.affine_solution_space(
        jacobian, target_at(values, PRIME), PRIME
    )
    if first_particular is None or len(first_basis) != 7:
        raise AssertionError("unexpected first-order lift space")
    left_rank, _zero, left_kernel = lift.affine_solution_space(
        transpose(jacobian), [0] * len(jacobian[0]), PRIME
    )
    if left_rank != rank:
        raise AssertionError("left/right Jacobian ranks disagree")

    v = branch_parameter % PRIME
    first_parameters = [
        (23 - 8 * v) % PRIME,
        24,
        v,
        (14 + 10 * v) % PRIME,
        (17 + 4 * v) % PRIME,
        25,
        21,
    ]
    first_correction = affine_point(
        first_particular, first_basis, first_parameters
    )
    values, signature = next_compatibility(
        values, PRIME, first_correction, left_kernel
    )
    if any(signature):
        raise AssertionError("declared second-order branch is not compatible")
    modulus = PRIME**2
    steps = [
        {
            "from_digit": 1,
            "to_digit": 2,
            "quadratic_branch_parameter": v,
            "chosen_parameters": first_parameters,
            "next_compatibility_zero": True,
        }
    ]
    status = "lifted"
    while modulus < PRIME**digits:
        candidate, metadata = choose_linear_lookahead(
            values, modulus, jacobian, left_kernel
        )
        if candidate is None:
            status = "obstructed"
            steps.append(
                {
                    "from_digit": len(steps) + 1,
                    "to_digit": len(steps) + 2,
                    **metadata,
                }
            )
            break
        values = candidate
        previous_modulus = modulus
        modulus *= PRIME
        steps.append(
            {
                "from_digit": len(steps) + 1,
                "to_digit": len(steps) + 2,
                "modulus_before": previous_modulus,
                "modulus_after": modulus,
                **metadata,
            }
        )
    if any(value % modulus for value in lift.equations(values)):
        raise AssertionError("final p-adic vector failed exact congruence check")
    payload = {
        "schema": "elliptic-rank30/e6-mw3-p3-1-adic-branch/v1",
        "field_prime": PRIME,
        "source_certificate": str(input_path),
        "source_certificate_sha256": data["certificate_sha256"],
        "branch_parameter_mod_29": v,
        "status": status,
        "lifted_digits": round(len(str(modulus)) / len(str(PRIME))),
        "modulus": modulus,
        "variable_names": list(lift.VARIABLE_NAMES),
        "residues": [value % modulus for value in values],
        "balanced_residues": [balanced(value, modulus) for value in values],
        "steps": steps,
        "claim_boundary": (
            "This is an exact compatible solution modulo the recorded power "
            "of 29. It is evidence for a formal branch, but is not a proof of "
            "a rational characteristic-zero point or a full infinite p-adic lift."
        ),
    }
    payload["lifted_digits"] = next(
        exponent for exponent in range(1, digits + 2) if PRIME**exponent == modulus
    )
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--branch-parameter", type=int, default=0)
    parser.add_argument("--digits", type=int, default=12)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    payload = compute(args.input, args.branch_parameter, args.digits)
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": payload["status"],
                "branch_parameter_mod_29": payload["branch_parameter_mod_29"],
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
