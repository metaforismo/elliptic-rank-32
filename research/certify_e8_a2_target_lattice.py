#!/usr/bin/env python3
"""Certify the target height profiles, intersections, torsion, and saturation.

The calculation assumes an elliptic K3 surface in characteristic zero with
reducible fibres exactly II*+3I3 and three sections with the displayed target
height Gram.  It is a lattice consequence, not an existence proof for those
sections.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path


GRAM_NUMERATOR = (
    (8, -1, 0),
    (-1, 10, 0),
    (0, 0, 12),
)


def determinant_3(matrix: tuple[tuple[int, ...], ...]) -> int:
    return (
        matrix[0][0]
        * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1]
        * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2]
        * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def diagonal_profiles(height: Fraction) -> list[tuple[int, int]]:
    profiles: list[tuple[int, int]] = []
    for intersection_with_zero in range(5):
        for nonidentity_i3_count in range(4):
            observed = (
                Fraction(4)
                + 2 * intersection_with_zero
                - Fraction(2 * nonidentity_i3_count, 3)
            )
            if observed == height:
                profiles.append((intersection_with_zero, nonidentity_i3_count))
    return profiles


def local_mutual_correction(label_left: int, label_right: int) -> Fraction:
    if label_left == 0 or label_right == 0:
        return Fraction(0)
    if label_left == label_right:
        return Fraction(2, 3)
    return Fraction(1, 3)


def section_intersection(
    pairing: Fraction,
    zero_left: int,
    zero_right: int,
    labels_left: tuple[int, int, int],
    labels_right: tuple[int, int, int],
) -> Fraction:
    correction = sum(
        (
            local_mutual_correction(left, right)
            for left, right in zip(labels_left, labels_right, strict=True)
        ),
        Fraction(0),
    )
    return (
        Fraction(2 + zero_left + zero_right)
        - correction
        - pairing
    )


def canonical_sha256(payload: object) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def fraction_string(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def compute_certificate() -> dict[str, object]:
    target_heights = [Fraction(8, 3), Fraction(10, 3), Fraction(4)]
    profiles = [diagonal_profiles(height) for height in target_heights]
    if profiles != [[(0, 2)], [(0, 1)], [(0, 0), (1, 3)]]:
        raise AssertionError("unexpected diagonal height profiles")

    # Exclude the second P3 profile using integrality of P2.P3.
    alternative_p3_intersections = sorted(
        {
            section_intersection(
                Fraction(0),
                0,
                1,
                (1, 0, 0),
                labels_p3,
            )
            for labels_p3 in itertools.product((1, 2), repeat=3)
        }
    )
    if alternative_p3_intersections != [Fraction(7, 3), Fraction(8, 3)]:
        raise AssertionError("P3 alternative-profile calculation failed")

    # Fix fibre order so P2 is nonidentity at the first I3.  Exhaust the P1
    # supports and labels; integrality and the target off-diagonal select the
    # unique orbit P1=(1,1,0), P2=(2,0,0), up to permutations and reversal.
    p2_labels = (2, 0, 0)
    compatible_p1: list[tuple[tuple[int, int, int], str]] = []
    all_p1 = []
    for support in itertools.combinations(range(3), 2):
        for nonzero_labels in itertools.product((1, 2), repeat=2):
            labels = [0, 0, 0]
            for index, label in zip(support, nonzero_labels, strict=True):
                labels[index] = label
            labels_tuple = tuple(labels)
            intersection = section_intersection(
                Fraction(-1, 3), 0, 0, labels_tuple, p2_labels
            )
            all_p1.append((labels_tuple, intersection))
            if intersection.denominator == 1 and intersection == 2:
                compatible_p1.append((labels_tuple, fraction_string(intersection)))
    expected_compatible = [
        ((1, 1, 0), "2"),
        ((1, 2, 0), "2"),
        ((1, 0, 1), "2"),
        ((1, 0, 2), "2"),
    ]
    if compatible_p1 != expected_compatible:
        raise AssertionError("component-profile orbit calculation failed")

    canonical_labels = {
        "P1": (1, 1, 0),
        "P2": (2, 0, 0),
        "P3": (0, 0, 0),
    }
    pairings = {
        ("P1", "P2"): Fraction(-1, 3),
        ("P1", "P3"): Fraction(0),
        ("P2", "P3"): Fraction(0),
    }
    intersections: dict[str, str] = {}
    for (left, right), pairing in pairings.items():
        value = section_intersection(
            pairing,
            0,
            0,
            canonical_labels[left],
            canonical_labels[right],
        )
        if value != 2:
            raise AssertionError("canonical intersection is not two")
        intersections[f"{left}.{right}"] = fraction_string(value)

    determinant_scaled = determinant_3(GRAM_NUMERATOR)
    determinant_height = Fraction(determinant_scaled, 27)
    if determinant_scaled != 948 or determinant_height != Fraction(316, 9):
        raise AssertionError("Gram determinant calculation failed")

    # B=3<,> is an even integral pairing on the full MW group: the only local
    # denominators are thirds and 3<P,P>=12+6(P.O)-2n.  Classify index-two
    # candidates v=(aP1+bP2+cP3)/2 by Mv being integral.
    order_two_dual_classes: list[dict[str, object]] = []
    for vector in itertools.product((0, 1), repeat=3):
        if vector == (0, 0, 0):
            continue
        matrix_times = tuple(
            sum(GRAM_NUMERATOR[row][column] * vector[column] for column in range(3))
            for row in range(3)
        )
        if all(value % 2 == 0 for value in matrix_times):
            norm_numerator = sum(
                vector[row]
                * GRAM_NUMERATOR[row][column]
                * vector[column]
                for row in range(3)
                for column in range(3)
            )
            norm = Fraction(norm_numerator, 4)
            order_two_dual_classes.append(
                {
                    "vector_mod_2": list(vector),
                    "B_norm": fraction_string(norm),
                    "is_even_isotropic": norm.denominator == 1 and norm % 2 == 0,
                }
            )
    if order_two_dual_classes != [
        {
            "vector_mod_2": [0, 0, 1],
            "B_norm": "3",
            "is_even_isotropic": False,
        }
    ]:
        raise AssertionError("index-two discriminant-form calculation failed")

    possible_indices = [
        index
        for index in range(1, 100)
        if determinant_scaled % (index * index) == 0
    ]
    if possible_indices != [1, 2]:
        raise AssertionError("possible-overlattice index calculation failed")

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "e8_a2_target_gram_lattice_primitivity",
        "assumptions": [
            "elliptic Jacobian K3 surface in characteristic zero",
            "reducible fibres are exactly II*+3I3",
            "three sections have the target height Gram",
        ],
        "target_height_gram": {
            "matrix": [
                [fraction_string(Fraction(value, 3)) for value in row]
                for row in GRAM_NUMERATOR
            ],
            "three_times_gram": [list(row) for row in GRAM_NUMERATOR],
            "determinant": fraction_string(determinant_height),
            "scaled_determinant": determinant_scaled,
        },
        "height_formula": {
            "formula": "h(P)=4+2*(P.O)-2*n(P)/3",
            "P1_profiles_(P.O,n)": [list(item) for item in profiles[0]],
            "P2_profiles_(P.O,n)": [list(item) for item in profiles[1]],
            "P3_initial_profiles_(P.O,n)": [list(item) for item in profiles[2]],
            "P3_alternative_intersections_with_P2": [
                fraction_string(value) for value in alternative_p3_intersections
            ],
            "P3_forced_profile": [0, 0],
        },
        "unique_component_orbit": {
            "canonical_labels_at_three_I3": {
                name: list(labels) for name, labels in canonical_labels.items()
            },
            "equivalences": [
                "permute the three I3 fibres",
                "simultaneously reverse labels 1 and 2 at any I3",
            ],
            "interpretation": (
                "P2 support is contained in P1 support and their labels are "
                "opposite at the overlap"
            ),
            "pairwise_section_intersections": intersections,
        },
        "torsion": {
            "status": "trivial",
            "proof": (
                "Every section has h(P)>=4-3*(2/3)=2 unless it is zero; "
                "a nonzero torsion section would have height zero."
            ),
        },
        "saturation": {
            "integral_even_pairing": "B=3 times the Shioda height pairing",
            "possible_indices_from_det_B": possible_indices,
            "order_two_dual_classes": order_two_dual_classes,
            "conclusion": (
                "The only nonzero order-two dual class is P3/2; its B-norm "
                "is 3, so it is not isotropic for the even discriminant "
                "form. No proper even integral overlattice exists."
            ),
            "generated_rank_three_lattice_is_primitive": True,
        },
        "neron_severi_consequence": {
            "trivial_lattice_absolute_discriminant": 27,
            "if_MW_rank_is_exactly_three": (
                "abs(disc(NS))=27*(316/9)=948"
            ),
            "if_more_independent_sections_exist": (
                "the target rank-three sublattice remains primitive, but "
                "the full NS discriminant must include the larger MW lattice"
            ),
        },
        "claim_boundary": {
            "sections_exist": False,
            "picard_rank_is_19": False,
            "rank31_curve_found": False,
        },
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()
    certificate = compute_certificate()
    encoded = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if arguments.compare is not None:
        expected = json.loads(arguments.compare.read_text(encoding="utf-8"))
        if certificate != expected:
            raise SystemExit("certificate mismatch")
    if arguments.output is not None:
        arguments.output.write_text(encoded, encoding="utf-8")
    elif arguments.compare is None:
        print(encoded, end="")


if __name__ == "__main__":
    main()
