#!/usr/bin/env python3
"""Follow an explicit flexible-slice branch above the GF(29) P3 #1 seed."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import hensel_lift_e6_mw3_gf29 as linear
import hensel_lift_e6_mw3_gf29_flexible_slice as system
import lift_e6_mw3_gf29_p3_1_branch as branch
from replay_e6_mw3_p2_split_cores import canonical_sha256


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "certificates" / "e6_mw3_gf29_independence.json"
DEFAULT_OUTPUT = ROOT / "certificates" / "e6_mw3_gf29_flexible_branch.json"
PRIME = 29
FIRST_PARAMETERS = [4, 18, 22, 2, 13, 0, 17, 0, 0, 0]


def first_parameters_on_line(parameter):
    value = parameter % PRIME
    return [
        4,
        value,
        22,
        (26 - 11 * value) % PRIME,
        (21 + 6 * value) % PRIME,
        (9 + 14 * value) % PRIME,
        17,
        0,
        0,
        0,
    ]


def compute(input_path, digits):
    data = json.loads(input_path.read_text())
    values = system.starting_vector(data, 1)
    jacobian = system.jacobian_mod_prime(values)
    rank, first_particular, first_basis = linear.affine_solution_space(
        jacobian, branch.target_at(values, PRIME, system), PRIME
    )
    if first_particular is None or len(first_basis) != len(FIRST_PARAMETERS):
        raise AssertionError("unexpected flexible first-order lift space")
    left_rank, _zero, left_kernel = linear.affine_solution_space(
        branch.transpose(jacobian), [0] * len(jacobian[0]), PRIME
    )
    if left_rank != rank:
        raise AssertionError("left/right Jacobian ranks disagree")
    first_correction = branch.affine_point(
        first_particular, first_basis, FIRST_PARAMETERS
    )
    values, signature = branch.next_compatibility(
        values, PRIME, first_correction, left_kernel, system
    )
    if any(signature):
        raise AssertionError("explicit flexible branch fails second order")
    modulus = PRIME**2
    steps = [
        {
            "from_digit": 1,
            "to_digit": 2,
            "chosen_parameters": FIRST_PARAMETERS,
            "next_compatibility_zero": True,
        }
    ]
    status = "lifted"
    while modulus < PRIME**digits:
        candidate, metadata = branch.choose_linear_lookahead(
            values, modulus, jacobian, left_kernel, system
        )
        if candidate is None:
            status = "greedy_lookahead_obstructed"
            steps.append(metadata)
            break
        values = candidate
        previous = modulus
        modulus *= PRIME
        steps.append(
            {"modulus_before": previous, "modulus_after": modulus, **metadata}
        )
    if any(value % modulus for value in system.equations(values)):
        raise AssertionError("final flexible branch congruence failed")
    lifted_digits = next(
        exponent for exponent in range(1, digits + 2) if PRIME**exponent == modulus
    )
    payload = {
        "schema": "elliptic-rank30/e6-mw3-flexible-adic-branch/v1",
        "field_prime": PRIME,
        "source_certificate": str(input_path),
        "source_certificate_sha256": data["certificate_sha256"],
        "P3_index": 1,
        "first_parameter_solution": FIRST_PARAMETERS,
        "jacobian_rank_mod_29": rank,
        "jacobian_nullity_mod_29": len(values) - rank,
        "status": status,
        "lifted_digits": lifted_digits,
        "modulus": modulus,
        "variable_names": list(system.VARIABLE_NAMES),
        "residues": [value % modulus for value in values],
        "steps": steps,
        "claim_boundary": (
            "The slice coordinates deform and the first nonlinear obstruction "
            "is solved exactly. Later look-ahead steps choose one point in each "
            "remaining affine correction space; a later failure of that greedy "
            "choice is not an exhaustive branch obstruction."
        ),
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--digits", type=int, default=12)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    payload = compute(args.input, args.digits)
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
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
