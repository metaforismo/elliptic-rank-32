#!/usr/bin/env python3
"""Probe flexible p-adic lifts of certified modular E6/MW3 seeds.

The input consists of a modular-slice search certificate and its matching
independence certificate.  For every certified triple P1/P2/P3, this script
forms the same 47 integral equations used in the GF(29) experiment, while
allowing the three slice coordinates r0, s0, x1 to deform.  It reports the
Jacobian rank and follows the deterministic affine correction with all free
parameters set to zero.

A failed greedy branch is not an obstruction to every p-adic branch.  A
successful finite lift is not yet an algebraic characteristic-zero point.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import hensel_lift_e6_mw3_gf29 as linear
import hensel_lift_e6_mw3_gf29_flexible_slice as system
from replay_e6_mw3_p2_split_cores import canonical_sha256


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SLICES = ROOT / "certificates" / "e6_mw3_modular_slices_gf23.json"
DEFAULT_INDEPENDENCE = (
    ROOT / "certificates" / "e6_mw3_modular_independence_gf23.json"
)
DEFAULT_OUTPUT = ROOT / "certificates" / "e6_mw3_modular_flexible_lifts_gf23.json"


def starting_vector(surface, coordinate_slice, p2_data, p3_data, prime):
    def padded(values, length):
        if len(values) > length:
            raise ValueError("polynomial exceeds its expected degree")
        return list(values) + [0] * (length - len(values))

    r0, s0, x1 = coordinate_slice
    all_roots = list(surface["valid_I2_roots"])
    if len(all_roots) < 2:
        raise ValueError("expected at least two split I2 roots")
    # The integral system records two movable I2 fibers.  A modular surface
    # may accidentally have more; retaining the first two gives a valid seed
    # while allowing every additional fiber to disappear in characteristic 0.
    roots = all_roots[:2]
    nodes = []
    for root in roots:
        choices = surface["singular_nodes"][str(root)]
        if len(choices) != 1:
            raise ValueError("expected one singular node at each I2 root")
        nodes.append(choices[0])
    surface_A = padded(surface["A"], 6)
    surface_B = padded(surface["B"], 9)
    p1_X = padded(surface["X1"], 3)
    p1_Y = padded(surface["Y1"], 5)
    p2_X = padded(p2_data["X2"], 5)
    p2_Y = padded(p2_data["Y2"], 8)
    p3_X = padded(p3_data["X3"], 3)
    p3_Y = padded(p3_data["Y3_up_to_sign"], 5)
    values = (
        [r0, s0, x1]
        + surface_A[1:6]
        + surface_B[1:8]
        + [p1_X[2]]
        + p1_Y[1:4]
        + [p2_data["pole"]]
        + p2_X
        + p2_Y
        + p3_X[1:3]
        + p3_Y[1:5]
        + roots
        + nodes
    )
    if len(values) != len(system.VARIABLE_NAMES):
        raise AssertionError("starting vector has wrong length")
    if any(value % prime for value in system.equations(values)):
        raise AssertionError("starting vector does not solve the system mod p")
    return values


def jacobian_mod_prime(values, prime):
    base = system.equations(values)
    columns = []
    for variable in range(len(values)):
        displaced = list(values)
        displaced[variable] += prime
        evaluated = system.equations(displaced)
        columns.append(
            [
                ((right - left) // prime) % prime
                for left, right in zip(base, evaluated, strict=True)
            ]
        )
    return [
        [column[row] for column in columns] for row in range(len(base))
    ]


def lift_greedy(surface, coordinate_slice, p2_data, p3_data, prime, digits):
    values = starting_vector(
        surface, coordinate_slice, p2_data, p3_data, prime
    )
    jacobian = jacobian_mod_prime(values, prime)
    rank, _zero = linear.solve_full_column_rank(
        jacobian, [0] * len(jacobian), prime
    )
    modulus = prime
    record = {
        "variables": len(values),
        "equations": len(system.equations(values)),
        "jacobian_rank_mod_p": rank,
        "jacobian_nullity_mod_p": len(values) - rank,
        "lifted_digits": 1,
    }
    while record["lifted_digits"] < digits:
        residual = system.equations(values)
        if any(value % modulus for value in residual):
            raise AssertionError("current vector is not a solution modulo p^n")
        target = [-(value // modulus) % prime for value in residual]
        _rank, correction = linear.solve_full_column_rank(
            jacobian, target, prime
        )
        if correction is None:
            record.update(
                {
                    "status": "chosen_path_obstructed",
                    "first_failed_modulus": modulus * prime,
                }
            )
            return record
        values = [
            value + modulus * delta
            for value, delta in zip(values, correction, strict=True)
        ]
        modulus *= prime
        if any(value % modulus for value in system.equations(values)):
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


def compute(slices_path, independence_path, digits):
    slices = json.loads(slices_path.read_text())
    independence = json.loads(independence_path.read_text())
    prime = slices["field_prime"]
    if independence["field"] != f"GF({prime})(t)":
        raise ValueError("slice and independence fields disagree")
    if independence["source_certificate_sha256"] != slices["certificate_sha256"]:
        raise ValueError("independence certificate does not match slice input")
    records = []
    for certified in independence["records"]:
        if not certified["independent"]:
            continue
        joint_index = certified["joint_index"]
        joint = slices["joint_candidates"][joint_index - 1]
        p2_index = certified["P2_index"]
        p3_index = certified["P3_index"]
        p2_data = joint["surface"]["geometric_P2_hits"][p2_index - 1]
        p3_data = joint["surface"]["P3_candidates"][p3_index - 1]
        result = lift_greedy(
            joint["surface"],
            joint["coordinate_slice"],
            p2_data,
            p3_data,
            prime,
            digits,
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
        "schema": "elliptic-rank30/e6-mw3-modular-flexible-hensel/v1",
        "field_prime": prime,
        "source_slices": str(slices_path),
        "source_slices_sha256": slices["certificate_sha256"],
        "source_independence": str(independence_path),
        "source_independence_sha256": independence["certificate_sha256"],
        "records": records,
        "claim_boundary": (
            "The slice coordinates r0,s0,x1 may deform. Each record follows "
            "the deterministic affine correction with all free variables set "
            "to zero. Failure is not an exhaustive p-adic obstruction; "
            "success is only a finite congruence, not a characteristic-zero "
            "algebraic solution."
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
    parser.add_argument("--digits", type=int, default=8)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    payload = compute(args.slices, args.independence, args.digits)
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
                            "joint_index",
                            "P2_index",
                            "P3_index",
                            "jacobian_rank_mod_p",
                            "jacobian_nullity_mod_p",
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
