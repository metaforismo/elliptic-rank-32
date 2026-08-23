#!/usr/bin/env sage-python
"""Independent Sage audit of the GF(29)(t) specialization certificate."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "certificates" / "e6_mw3_gf29_independence.json"


def polynomial(coefficients, variable):
    return sum(coefficient * variable**index for index, coefficient in enumerate(coefficients))


def audit(path):
    data = json.loads(path.read_text())
    field = GF(29)
    ring = PolynomialRing(field, "t")
    t = ring.gen()
    function_field = ring.fraction_field()
    A = function_field(polynomial(data["A"], t))
    B = function_field(polynomial(data["B"], t))
    curve = EllipticCurve(function_field, [A, B])
    P1 = curve(
        function_field(polynomial(data["P1"]["X"], t)),
        function_field(polynomial(data["P1"]["Y"], t)),
    )
    p2_data = data["P2"]
    denominator = function_field(t - p2_data["pole"])
    P2 = curve(
        function_field(polynomial(p2_data["X2"], t)) / denominator**2,
        function_field(polynomial(p2_data["Y2"], t)) / denominator**3,
    )
    audited = []
    for result in data["results"]:
        p3_data = result["P3"]
        P3 = curve(
            function_field(polynomial(p3_data["X3"], t)),
            function_field(polynomial(p3_data["Y3_up_to_sign"], t)),
        )
        ell = result["ell"]
        simultaneous = set(itertools.product(range(ell), repeat=3))
        for record in result["quotient_specializations"]:
            fiber = EllipticCurve(
                field, [field(record["A_t"]), field(record["B_t"])]
            )
            if fiber.order() != record["group_order"]:
                raise AssertionError("Sage group order disagrees with certificate")
            points = [fiber(*coordinates) for coordinates in record["points"]]
            ell_image = {ell * point for point in fiber.points()}
            if len(ell_image) != record["ell_multiple_image_size"]:
                raise AssertionError("Sage ell-multiple image size disagrees")
            local_kernel = {
                coefficients
                for coefficients in itertools.product(range(ell), repeat=3)
                if sum(
                    (coefficient * point for coefficient, point in zip(coefficients, points)),
                    fiber(0),
                )
                in ell_image
            }
            simultaneous &= local_kernel
            stored = {tuple(item) for item in record["kernel_after"]}
            if simultaneous != stored:
                raise AssertionError("Sage quotient kernel disagrees")
        if simultaneous != {(0, 0, 0)}:
            raise AssertionError("nontrivial simultaneous kernel")
        witness = result["torsion_exclusion_specialization"]
        witness_fiber = EllipticCurve(
            field, [field(witness["A_t"]), field(witness["B_t"])]
        )
        if witness_fiber.order() != witness["group_order"]:
            raise AssertionError("Sage torsion-witness order disagrees")
        if witness_fiber.order() % ell == 0:
            raise AssertionError("torsion witness is not prime to ell")
        audited.append(
            {
                "P3_index": result["P3_index"],
                "ell": ell,
                "generic_points_constructed": all(point in curve for point in (P1, P2, P3)),
                "simultaneous_kernel": [[0, 0, 0]],
                "torsion_witness_order": int(witness_fiber.order()),
            }
        )
    return audited


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    args = parser.parse_args()
    print(json.dumps({"sage_audit": audit(args.input)}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
