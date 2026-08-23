#!/usr/bin/env python3
"""Certify three E6/MW3 sections by simultaneous good specialization.

For a prime ell different from the ground-field characteristic, a good
specialization whose group order is prime to ell proves that the generic
ell-torsion is trivial.  If the images of P1, P2, P3 are independent in a
product of finite quotients E_t(F_p)/ell E_t(F_p), an integral relation would
therefore have all coefficients divisible by ell, could be divided by ell,
and iterated to a contradiction.

This gives an exact Mordell--Weil independence certificate without numerical
heights or a full component-by-component resolution of the elliptic surface.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

from replay_e6_mw3_p2_split_cores import canonical_sha256


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "certificates" / "e6_mw3_modular_slice_gf29_r1_s2_x1.json"
DEFAULT_OUTPUT = ROOT / "certificates" / "e6_mw3_gf29_independence.json"


def poly_eval(coefficients, value, prime):
    result = 0
    for coefficient in reversed(coefficients):
        result = (result * value + coefficient) % prime
    return result


def inverse(value, prime):
    if value % prime == 0:
        raise ZeroDivisionError("finite-field inverse of zero")
    return pow(value, -1, prime)


def point_add(left, right, a, prime):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if (y1 + y2) % prime == 0:
            return None
        slope = (3 * x1 * x1 + a) * inverse(2 * y1, prime) % prime
    else:
        slope = (y2 - y1) * inverse(x2 - x1, prime) % prime
    x3 = (slope * slope - x1 - x2) % prime
    y3 = (slope * (x1 - x3) - y1) % prime
    return x3, y3


def point_multiply(coefficient, point, a, prime):
    if coefficient < 0:
        return point_multiply(-coefficient, negate(point, prime), a, prime)
    result = None
    addend = point
    while coefficient:
        if coefficient & 1:
            result = point_add(result, addend, a, prime)
        addend = point_add(addend, addend, a, prime)
        coefficient >>= 1
    return result


def negate(point, prime):
    if point is None:
        return None
    return point[0], (-point[1]) % prime


def on_curve(point, a, b, prime):
    if point is None:
        return True
    x, y = point
    return (y * y - x * x * x - a * x - b) % prime == 0


def enumerate_points(a, b, prime):
    points = [None]
    for x in range(prime):
        rhs = (x * x * x + a * x + b) % prime
        points.extend((x, y) for y in range(prime) if y * y % prime == rhs)
    return points


def discriminant(a, b, prime):
    return (-16 * (4 * a * a * a + 27 * b * b)) % prime


def specialize_polynomial_point(x_coefficients, y_coefficients, t, prime):
    return (
        poly_eval(x_coefficients, t, prime),
        poly_eval(y_coefficients, t, prime),
    )


def specialize_p2(data, t, prime):
    denominator = (t - data["pole"]) % prime
    if denominator == 0:
        raise ZeroDivisionError("P2 pole")
    return (
        poly_eval(data["X2"], t, prime) * inverse(denominator**2, prime) % prime,
        poly_eval(data["Y2"], t, prime) * inverse(denominator**3, prime) % prime,
    )


def linear_combination(coefficients, points, a, prime):
    result = None
    for coefficient, point in zip(coefficients, points, strict=True):
        result = point_add(
            result, point_multiply(coefficient, point, a, prime), a, prime
        )
    return result


def specialization_records(surface, p2_data, p3_data, prime):
    records = []
    for t in range(prime):
        if t == p2_data["pole"]:
            continue
        a = poly_eval(surface["A"], t, prime)
        b = poly_eval(surface["B"], t, prime)
        if discriminant(a, b, prime) == 0:
            continue
        points = [
            specialize_polynomial_point(surface["X1"], surface["Y1"], t, prime),
            specialize_p2(p2_data, t, prime),
            specialize_polynomial_point(
                p3_data["X3"], p3_data["Y3_up_to_sign"], t, prime
            ),
        ]
        if not all(on_curve(point, a, b, prime) for point in points):
            raise AssertionError(f"section failed after specialization t={t}")
        group = enumerate_points(a, b, prime)
        records.append(
            {
                "t": t,
                "A_t": a,
                "B_t": b,
                "group_order": len(group),
                "points": points,
                "group": group,
            }
        )
    return records


def quotient_kernel(record, ell, prime):
    ell_multiples = {
        point_multiply(ell, point, record["A_t"], prime)
        for point in record["group"]
    }
    kernel = []
    for coefficients in itertools.product(range(ell), repeat=3):
        value = linear_combination(
            coefficients, record["points"], record["A_t"], prime
        )
        if value in ell_multiples:
            kernel.append(coefficients)
    return set(kernel), len(ell_multiples)


def find_certificate_for_p3(
    surface, p2_data, p3_data, prime, ell_candidates
):
    records = specialization_records(surface, p2_data, p3_data, prime)
    for ell in ell_candidates:
        torsion_witness = next(
            (record for record in records if record["group_order"] % ell), None
        )
        if torsion_witness is None:
            continue
        full_kernel = set(itertools.product(range(ell), repeat=3))
        useful = []
        for record in records:
            local_kernel, ell_image_size = quotient_kernel(record, ell, prime)
            intersection = full_kernel & local_kernel
            if len(intersection) < len(full_kernel):
                useful.append(
                    {
                        "t": record["t"],
                        "A_t": record["A_t"],
                        "B_t": record["B_t"],
                        "group_order": record["group_order"],
                        "points": [list(point) for point in record["points"]],
                        "ell_multiple_image_size": ell_image_size,
                        "quotient_size": record["group_order"] // ell_image_size,
                        "kernel_size_before": len(full_kernel),
                        "kernel_size_after": len(intersection),
                        "kernel_after": [
                            list(coefficients) for coefficients in sorted(intersection)
                        ],
                    }
                )
                full_kernel = intersection
            if full_kernel == {(0, 0, 0)}:
                return {
                    "ell": ell,
                    "torsion_exclusion_specialization": {
                        "t": torsion_witness["t"],
                        "A_t": torsion_witness["A_t"],
                        "B_t": torsion_witness["B_t"],
                        "group_order": torsion_witness["group_order"],
                    },
                    "quotient_specializations": useful,
                    "simultaneous_kernel": [[0, 0, 0]],
                    "independent": True,
                }
    return {
        "ell_candidates": list(ell_candidates),
        "independent": False,
        "reason": "no tested prime produced a trivial simultaneous kernel",
    }


def compute(input_path, ell_candidates):
    source = json.loads(input_path.read_text())
    prime = source["field_prime"]
    if prime != 29:
        raise ValueError("the current certificate is pinned to the GF(29) seed")
    if len(source["joint_candidates"]) != 1:
        raise AssertionError("expected exactly one joint candidate")
    joint = source["joint_candidates"][0]
    surface = joint["surface"]
    results = []
    for index, p3_data in enumerate(surface["P3_candidates"], start=1):
        result = find_certificate_for_p3(
            surface,
            surface["geometric_P2_hits"][0],
            p3_data,
            prime,
            tuple(ell_candidates),
        )
        results.append({"P3_index": index, "P3": p3_data, **result})
    payload = {
        "schema": "elliptic-rank30/e6-mw3-specialization-independence/v1",
        "field": "GF(29)(t)",
        "source_certificate": str(input_path),
        "source_certificate_sha256": source["certificate_sha256"],
        "coordinate_slice": joint["coordinate_slice"],
        "core": joint["core"],
        "A": surface["A"],
        "B": surface["B"],
        "P1": {"X": surface["X1"], "Y": surface["Y1"]},
        "P2": surface["geometric_P2_hits"][0],
        "results": results,
        "all_four_triples_independent": all(
            result["independent"] for result in results
        ),
        "proof": (
            "At every listed t the surface and sections specialize exactly to "
            "the recorded finite group. The simultaneous kernel in the product "
            "of E_t(F_29)/ell E_t(F_29) is zero. The torsion witness has group "
            "order prime to ell, so good-reduction injectivity excludes generic "
            "ell-torsion. Any integral relation is therefore ell-divisible, can "
            "be divided by ell, and infinite descent forces all coefficients zero."
        ),
        "claim_boundary": (
            "This proves Mordell--Weil rank at least 3 over GF(29)(t) for each "
            "listed triple. It does not prove that the modular surface or its "
            "sections lift to characteristic zero."
        ),
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--ell", default="2,3,5,7,11,13")
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    ell_candidates = [int(value) for value in args.ell.split(",")]
    payload = compute(args.input, ell_candidates)
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "all_four_triples_independent": payload[
                    "all_four_triples_independent"
                ],
                "results": [
                    {
                        "P3_index": result["P3_index"],
                        "ell": result.get("ell"),
                        "independent": result["independent"],
                        "specializations": len(
                            result.get("quotient_specializations", [])
                        ),
                    }
                    for result in payload["results"]
                ],
                "certificate_sha256": payload["certificate_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
