#!/usr/bin/env python3
"""Classify finite reducible-fiber components of modular E6/MW3 seeds.

The public E6 construction uses the component order

    (IV*, I4@0, I4@1, I2@lambda, I2@mu)

and the target profiles

    P1 = (1,0,1,0,0), P2 = (1,1,2,1,1), P3 = (2,3,0,0,0).

The finite entries can be tested without resolving the surface explicitly.
At an I_n fiber a section is on the identity component exactly when its
specialization is a smooth point (including O).  At I4, a nodal section has
class 2 exactly when its double is smooth; if its double is still nodal, its
class is odd (1 or 3).  The latter ambiguity is exchanged by P -> -P and is
all that remains before choosing section signs and the IV* labels.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import check_e6_mw3_candidate_relations_gf31 as group
import replay_e6_mw3_p2_split_cores as sections


ROOT = Path(__file__).resolve().parents[1]


def evaluate_polynomial(coefficients, value, prime):
    result = 0
    for coefficient in reversed(coefficients):
        result = (result * value + coefficient) % prime
    return result


def evaluate_rational(function, value, prime):
    numerator, denominator = function
    denominator_value = evaluate_polynomial(denominator, value, prime)
    if denominator_value == 0:
        return None
    return (
        evaluate_polynomial(numerator, value, prime)
        * pow(denominator_value, -1, prime)
        % prime
    )


def reduction_status(point, fiber_point, node, prime):
    """Return O, node, smooth, or inconsistent for a specialized section."""
    if point is None:
        return "O"
    x_value = evaluate_rational(point[0], fiber_point, prime)
    y_value = evaluate_rational(point[1], fiber_point, prime)
    if x_value is None or y_value is None:
        return "O" if x_value is None and y_value is None else "inconsistent"
    if x_value == node % prime and y_value == 0:
        return "node"
    return "smooth"


def i4_component_parity(point, A, fiber_point, node, prime):
    status = reduction_status(point, fiber_point, node, prime)
    if status in ("O", "smooth"):
        return {
            "section_reduction": status,
            "double_reduction": None,
            "component_class": "0",
        }
    if status != "node":
        return {
            "section_reduction": status,
            "double_reduction": None,
            "component_class": "inconsistent",
        }
    double = group.multiply(2, point, A)
    double_status = reduction_status(double, fiber_point, node, prime)
    if double_status in ("O", "smooth"):
        component = "2"
    elif double_status == "node":
        component = "odd"
    else:
        component = "inconsistent"
    return {
        "section_reduction": status,
        "double_reduction": double_status,
        "component_class": component,
    }


def i2_component(point, fiber_point, node, prime):
    status = reduction_status(point, fiber_point, node, prime)
    if status in ("O", "smooth"):
        component = "0"
    elif status == "node":
        component = "1"
    else:
        component = "inconsistent"
    return {"section_reduction": status, "component_class": component}


def polynomial_point(x_coefficients, y_coefficients):
    return group.polynomial(x_coefficients), group.polynomial(y_coefficients)


def p2_point(data):
    z = [-data["pole"], 1]
    y_coefficients = data.get("Y2", data.get("Y2_up_to_sign"))
    if y_coefficients is None:
        raise KeyError("P2 record lacks Y2/Y2_up_to_sign")
    return (
        group.rational(data["X2"], sections.poly_pow(z, 2)),
        group.rational(y_coefficients, sections.poly_pow(z, 3)),
    )


def classify_point(point, A, finite_fibers, prime):
    result = {}
    for name in ("I4_0", "I4_1"):
        fiber = finite_fibers[name]
        result[name] = i4_component_parity(
            point, A, fiber["point"], fiber["node"], prime
        )
    for name in ("I2_lambda", "I2_mu"):
        fiber = finite_fibers[name]
        result[name] = i2_component(
            point, fiber["point"], fiber["node"], prime
        )
    return result


def component_vector(classification):
    return [
        classification[name]["component_class"]
        for name in ("I4_0", "I4_1", "I2_lambda", "I2_mu")
    ]


TARGET_FINITE = {
    "P1": ["0", "odd", "0", "0"],
    "P2": ["odd", "2", "1", "1"],
    "P3": ["odd", "0", "0", "0"],
}


def classify_joint(joint, prime, joint_index):
    surface = joint["surface"]
    coordinate_slice = joint["coordinate_slice"]
    core = joint["core"]
    roots = surface["valid_I2_roots"]
    if len(roots) < 2:
        raise AssertionError("joint surface has fewer than two I2 roots")
    nodes = {int(key): value for key, value in surface["singular_nodes"].items()}
    A = group.polynomial(surface["A"])
    B = group.polynomial(surface["B"])
    P1 = polynomial_point(surface["X1"], surface["Y1"])
    if not group.point_on_curve(P1, A, B):
        raise AssertionError("P1 failed the function-field curve equation")

    records = []
    for left, right in itertools.combinations(roots, 2):
        finite_fibers = {
            "I4_0": {"point": 0, "node": coordinate_slice[1]},
            "I4_1": {"point": 1, "node": core[3]},
            "I2_lambda": {"point": left, "node": nodes[left][0]},
            "I2_mu": {"point": right, "node": nodes[right][0]},
        }
        p1_classification = classify_point(P1, A, finite_fibers, prime)
        p1_vector = component_vector(p1_classification)
        for p2_index, p2_data in enumerate(surface["geometric_P2_hits"], start=1):
            P2 = p2_point(p2_data)
            if not group.point_on_curve(P2, A, B):
                raise AssertionError("P2 failed the function-field curve equation")
            p2_classification = classify_point(P2, A, finite_fibers, prime)
            p2_vector = component_vector(p2_classification)
            for p3_index, p3_data in enumerate(surface["P3_candidates"], start=1):
                P3 = polynomial_point(
                    p3_data["X3"], p3_data["Y3_up_to_sign"]
                )
                if not group.point_on_curve(P3, A, B):
                    raise AssertionError("P3 failed the function-field curve equation")
                p3_classification = classify_point(P3, A, finite_fibers, prime)
                p3_vector = component_vector(p3_classification)
                matches = {
                    "P1": p1_vector == TARGET_FINITE["P1"],
                    "P2": p2_vector == TARGET_FINITE["P2"],
                    "P3": p3_vector == TARGET_FINITE["P3"],
                }
                records.append(
                    {
                        "joint_index": joint_index,
                        "branch": joint["branch"],
                        "coordinate_slice": coordinate_slice,
                        "core": core,
                        "I2_roots": [left, right],
                        "P2_index": p2_index,
                        "P3_index": p3_index,
                        "P1": {
                            "finite_component_vector": p1_vector,
                            "fibers": p1_classification,
                        },
                        "P2": {
                            "finite_component_vector": p2_vector,
                            "fibers": p2_classification,
                        },
                        "P3": {
                            "finite_component_vector": p3_vector,
                            "fibers": p3_classification,
                        },
                        "matches": matches,
                        "matches_target_finite_profile": all(matches.values()),
                    }
                )
    return records


def compute(input_path):
    source = json.loads(input_path.read_text())
    prime = source["field_prime"]
    group.configure_prime(prime)
    records = []
    for joint_index, joint in enumerate(source["joint_candidates"], start=1):
        records.extend(classify_joint(joint, prime, joint_index))
    matching = [record for record in records if record["matches_target_finite_profile"]]
    payload = {
        "schema": "elliptic-rank30/e6-mw3-finite-component-profile/v1",
        "field_prime": prime,
        "source_certificate": str(input_path),
        "source_certificate_sha256": source["certificate_sha256"],
        "target_component_order": ["I4@0", "I4@1", "I2@lambda", "I2@mu"],
        "target_finite_profiles_up_to_sign": TARGET_FINITE,
        "triples_classified": len(records),
        "triples_matching_target_finite_profile": len(matching),
        "records": records,
        "proof": (
            "Identity versus nonidentity is decided by smooth versus nodal "
            "specialization. For I4, the component homomorphism sends doubling "
            "to multiplication by 2 in Z/4Z: a nodal point has class 2 iff its "
            "double is smooth, and odd class iff its double remains nodal."
        ),
        "claim_boundary": (
            "This certifies the four finite component entries up to the sign "
            "ambiguity 1 versus 3. It does not yet choose IV* labels, certify "
            "mutual section intersections, or prove a characteristic-zero lift."
        ),
    }
    payload["certificate_sha256"] = sections.canonical_sha256(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    payload = compute(args.input)
    output = args.output or (
        ROOT
        / "certificates"
        / f"e6_mw3_component_profiles_gf{payload['field_prime']}.json"
    )
    if not args.no_write:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "field_prime": payload["field_prime"],
                "triples_classified": payload["triples_classified"],
                "triples_matching_target_finite_profile": payload[
                    "triples_matching_target_finite_profile"
                ],
                "certificate_sha256": payload["certificate_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
