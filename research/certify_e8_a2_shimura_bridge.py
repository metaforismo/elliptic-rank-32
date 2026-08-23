#!/usr/bin/env python3
"""Certify the lattice bridge from the E8+A2^3 target to X(6,79).

The target three-section configuration gives an integral positive-definite
essential lattice of rank 17.  This script constructs that lattice directly
from the fibre components and section intersections, then compares its finite
quadratic form with both

* the transparent A11 + K6 determinant-948 neighbor seed; and
* the rank-three period lattice used by the X(6,79) neighbor workflow.

All computations use exact standard-library arithmetic.  The result is a
lattice-theoretic bridge.  It is not an explicit moduli map, Weierstrass
surface, rational P3 section, or rank-31 curve.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Iterable


Matrix = list[list[int]]
RationalMatrix = list[list[Fraction]]


PERIOD_LATTICE: Matrix = [
    [-316, 0, 288],
    [0, 474, -15],
    [288, -15, -262],
]

K6: Matrix = [
    [2, -1, 0, 0, 0, 0],
    [-1, 2, -1, 0, 0, 0],
    [0, -1, 2, -1, 0, 0],
    [0, 0, -1, 8, -1, 0],
    [0, 0, 0, -1, 2, -1],
    [0, 0, 0, 0, -1, 2],
]


def zero_matrix(size: int) -> Matrix:
    return [[0 for _ in range(size)] for _ in range(size)]


def a_cartan(rank: int) -> Matrix:
    result = zero_matrix(rank)
    for index in range(rank):
        result[index][index] = 2
        if index + 1 < rank:
            result[index][index + 1] = -1
            result[index + 1][index] = -1
    return result


def e8_cartan() -> Matrix:
    result = zero_matrix(8)
    for index in range(8):
        result[index][index] = 2
    # Arms have lengths 2, 4, 1 from the trivalent vertex 2.
    for left, right in [
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)
    ]:
        result[left][right] = -1
        result[right][left] = -1
    return result


def block_diagonal(*blocks: Matrix) -> Matrix:
    size = sum(len(block) for block in blocks)
    result = zero_matrix(size)
    offset = 0
    for block in blocks:
        for row in range(len(block)):
            for column in range(len(block)):
                result[offset + row][offset + column] = block[row][column]
        offset += len(block)
    return result


def target_essential_lattice() -> Matrix:
    """Return the positive essential Gram for E8+A2^3 plus three sections.

    In NS use Q_i=P_i-O-2F, which is orthogonal to the hyperbolic plane.
    Then -Q_i^2=4 and -Q_i.Q_j=0.  A section meeting a simple fibre
    component once gives pairing -1 in the positive essential lattice.
    The component orbit is P1=(1,1,0), P2=(2,0,0), P3=(0,0,0).
    """

    result = block_diagonal(
        e8_cartan(), a_cartan(2), a_cartan(2), a_cartan(2), [[4]], [[4]], [[4]]
    )
    # A2 blocks start at 8, 10, 12; adjusted sections at 14, 15, 16.
    result[14][8] = result[8][14] = -1
    result[14][10] = result[10][14] = -1
    result[15][9] = result[9][15] = -1
    return result


def transparent_seed_lattice() -> Matrix:
    return block_diagonal(a_cartan(11), K6)


def bareiss_determinant(matrix: Matrix) -> int:
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("determinant requires a square matrix")
    if size == 0:
        return 1
    work = [row[:] for row in matrix]
    sign = 1
    previous_pivot = 1
    for column in range(size - 1):
        if work[column][column] == 0:
            swap = next(
                (
                    row
                    for row in range(column + 1, size)
                    if work[row][column] != 0
                ),
                None,
            )
            if swap is None:
                return 0
            work[column], work[swap] = work[swap], work[column]
            sign = -sign
        pivot = work[column][column]
        for row in range(column + 1, size):
            for target in range(column + 1, size):
                numerator = (
                    work[row][target] * pivot
                    - work[row][column] * work[column][target]
                )
                if numerator % previous_pivot:
                    raise AssertionError("non-exact Bareiss division")
                work[row][target] = numerator // previous_pivot
        previous_pivot = pivot
        for row in range(column + 1, size):
            work[row][column] = 0
    return sign * work[-1][-1]


def solve_exact(matrix: Matrix, rhs: list[int]) -> list[Fraction]:
    size = len(matrix)
    if len(rhs) != size or any(len(row) != size for row in matrix):
        raise ValueError("invalid exact linear system")
    augmented: RationalMatrix = [
        [Fraction(value) for value in row] + [Fraction(rhs[index])]
        for index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot_row = next(
            (
                row
                for row in range(column, size)
                if augmented[row][column] != 0
            ),
            None,
        )
        if pivot_row is None:
            raise ArithmeticError("singular matrix")
        augmented[column], augmented[pivot_row] = (
            augmented[pivot_row], augmented[column]
        )
        pivot = augmented[column][column]
        augmented[column] = [value / pivot for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [
                    value - factor * pivot_value
                    for value, pivot_value in zip(
                        augmented[row], augmented[column], strict=True
                    )
                ]
    return [augmented[row][-1] for row in range(size)]


def ldl(matrix: Matrix) -> tuple[RationalMatrix, list[Fraction]]:
    size = len(matrix)
    lower: RationalMatrix = [
        [Fraction(int(row == column)) for column in range(size)]
        for row in range(size)
    ]
    diagonal = [Fraction(0) for _ in range(size)]
    for index in range(size):
        diagonal[index] = Fraction(matrix[index][index]) - sum(
            lower[index][prior] ** 2 * diagonal[prior]
            for prior in range(index)
        )
        if diagonal[index] == 0:
            raise ArithmeticError("zero LDL pivot")
        for row in range(index + 1, size):
            lower[row][index] = (
                Fraction(matrix[row][index])
                - sum(
                    lower[row][prior]
                    * lower[index][prior]
                    * diagonal[prior]
                    for prior in range(index)
                )
            ) / diagonal[index]
    return lower, diagonal


def signature(matrix: Matrix) -> tuple[int, int]:
    _lower, diagonal = ldl(matrix)
    return (
        sum(value > 0 for value in diagonal),
        sum(value < 0 for value in diagonal),
    )


def fraction_string(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def discriminant_generator(
    matrix: Matrix, witness: list[int]
) -> dict[str, object]:
    dual = solve_exact(matrix, witness)
    order = math.lcm(*(value.denominator for value in dual))
    norm = sum(
        Fraction(witness[index]) * dual[index]
        for index in range(len(witness))
    )
    determinant = abs(bareiss_determinant(matrix))
    if order != determinant:
        raise AssertionError("witness does not generate the discriminant group")
    scaled_mod_two = int(norm * determinant) % (2 * determinant)
    return {
        "witness": witness,
        "order": order,
        "norm": fraction_string(norm),
        "scaled_norm_mod_2det": scaled_mod_two,
    }


def enumerate_norm_two(matrix: Matrix) -> list[tuple[int, ...]]:
    """Exact Fincke-Pohst enumeration using an LDL triangularization."""

    lower, diagonal = ldl(matrix)
    if any(value <= 0 for value in diagonal):
        raise ValueError("root enumeration needs a positive-definite matrix")
    size = len(matrix)
    vector = [0 for _ in range(size)]
    roots: list[tuple[int, ...]] = []

    def recurse(index: int, partial: Fraction) -> None:
        if index < 0:
            if partial == 2:
                roots.append(tuple(vector))
            return
        remaining = Fraction(2) - partial
        if remaining < 0:
            return
        center = sum(
            lower[row][index] * vector[row]
            for row in range(index + 1, size)
        )
        radius = math.sqrt(float(remaining / diagonal[index]))
        low = math.ceil(float(-center) - radius - 1e-12)
        high = math.floor(float(-center) + radius + 1e-12)
        for coordinate in range(low, high + 1):
            contribution = diagonal[index] * (Fraction(coordinate) + center) ** 2
            if contribution <= remaining:
                vector[index] = coordinate
                recurse(index - 1, partial + contribution)
        vector[index] = 0

    recurse(size - 1, Fraction(0))
    for root in roots:
        norm = sum(
            root[row] * matrix[row][column] * root[column]
            for row in range(size)
            for column in range(size)
        )
        if norm != 2:
            raise AssertionError("root enumerator returned a non-root")
    return roots


def valuation_unit(value: Fraction, prime: int) -> tuple[int, Fraction]:
    exponent = 0
    numerator = value.numerator
    denominator = value.denominator
    while numerator % prime == 0:
        numerator //= prime
        exponent += 1
    while denominator % prime == 0:
        denominator //= prime
        exponent -= 1
    return exponent, Fraction(numerator, denominator)


def legendre_fraction(value: Fraction, prime: int) -> int:
    residue = (
        value.numerator * pow(value.denominator, -1, prime)
    ) % prime
    symbol = pow(residue, (prime - 1) // 2, prime)
    return -1 if symbol == prime - 1 else 1


def hilbert_symbol(left: Fraction, right: Fraction, prime: int | str) -> int:
    if prime == "infinity":
        return -1 if left < 0 and right < 0 else 1
    if not isinstance(prime, int):
        raise TypeError("finite place must be an integer prime")
    alpha, unit_left = valuation_unit(left, prime)
    beta, unit_right = valuation_unit(right, prime)
    if prime == 2:
        left_mod_8 = (
            unit_left.numerator * pow(unit_left.denominator, -1, 8)
        ) % 8
        right_mod_8 = (
            unit_right.numerator * pow(unit_right.denominator, -1, 8)
        ) % 8
        exponent = (
            ((left_mod_8 - 1) // 2) * ((right_mod_8 - 1) // 2)
            + alpha * ((right_mod_8 * right_mod_8 - 1) // 8)
            + beta * ((left_mod_8 * left_mod_8 - 1) // 8)
        ) % 2
        return -1 if exponent else 1
    sign = -1 if (alpha * beta * ((prime - 1) // 2)) % 2 else 1
    if beta % 2:
        sign *= legendre_fraction(unit_left, prime)
    if alpha % 2:
        sign *= legendre_fraction(unit_right, prime)
    return sign


def canonical_sha256(payload: object) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def compute_certificate() -> dict[str, object]:
    target = target_essential_lattice()
    seed = transparent_seed_lattice()
    period = PERIOD_LATTICE

    if bareiss_determinant(target) != 948:
        raise AssertionError("unexpected target determinant")
    if bareiss_determinant(seed) != 948:
        raise AssertionError("unexpected seed determinant")
    if bareiss_determinant(period) != -948:
        raise AssertionError("unexpected period determinant")
    if signature(target) != (17, 0) or signature(seed) != (17, 0):
        raise AssertionError("essential lattice is not positive definite")
    if signature(period) != (2, 1):
        raise AssertionError("period lattice has the wrong signature")
    if any(matrix[index][index] % 2 for matrix in (target, seed, period) for index in range(len(matrix))):
        raise AssertionError("one of the lattices is not even")

    target_witness = [0] * 17
    for index in (8, 12, 16):
        target_witness[index] = 1
    seed_witness = [0] * 17
    for index in (0, 11):
        seed_witness[index] = 1
    period_witness = [1, 1, 0]

    target_generator = discriminant_generator(target, target_witness)
    seed_generator = discriminant_generator(seed, seed_witness)
    period_generator = discriminant_generator(period, period_witness)
    modulus = 2 * 948
    target_multiplier = 67
    seed_multiplier = 101
    if (
        target_generator["scaled_norm_mod_2det"]
        * target_multiplier
        * target_multiplier
        - period_generator["scaled_norm_mod_2det"]
    ) % modulus:
        raise AssertionError("target and period discriminant forms differ")
    if (
        seed_generator["scaled_norm_mod_2det"]
        * seed_multiplier
        * seed_multiplier
        - period_generator["scaled_norm_mod_2det"]
    ) % modulus:
        raise AssertionError("seed and period discriminant forms differ")

    roots = enumerate_norm_two(target)
    mixed_roots = [root for root in roots if any(root[index] for index in range(14, 17))]
    if len(roots) != 240 + 3 * 6 or mixed_roots:
        raise AssertionError("target has unexpected roots")

    _period_lower, period_diagonal = ldl(period)
    if period_diagonal != [Fraction(-316), Fraction(474), Fraction(1, 158)]:
        raise AssertionError("unexpected period diagonalization")
    # For <a,b,c>, C^+(q)=(-ab,-ac).  Here this is
    # (149784,2)=(6*158^2,2), hence the quaternion (6,2).
    quaternion = (Fraction(6), Fraction(2))
    local_symbols = {
        str(place): hilbert_symbol(*quaternion, place)
        for place in (2, 3, 79, "infinity")
    }
    if local_symbols != {"2": -1, "3": -1, "79": 1, "infinity": 1}:
        raise AssertionError("unexpected Clifford ramification")

    t_value = Fraction(14, 13)
    u_value = Fraction(16064, 2197)
    shimura_rhs = (
        16 * t_value**6 - 19 * t_value**4 + 88 * t_value**2 - 48
    )
    if shimura_rhs != u_value**2:
        raise AssertionError("published non-CM orbit representative is off C")

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "e8_a2_target_to_x6_79_lattice_bridge",
        "exact_claim": (
            "The positive essential lattice forced by the primitive E8+A2^3 "
            "three-section target is even of rank 17 and determinant 948, has "
            "root system exactly E8+A2^3, and has the same cyclic finite "
            "quadratic form as both the transparent determinant-948 neighbor "
            "seed and the X(6,79) period lattice used by the rank-17 workflow."
        ),
        "target_essential_lattice": {
            "construction": (
                "E8+A2^3 fibre roots plus Q_i=P_i-O-2F for component orbit "
                "P1=(1,1,0), P2=(2,0,0), P3=(0,0,0)"
            ),
            "gram_construction": {
                "root_blocks": ["E8", "A2", "A2", "A2"],
                "adjusted_section_diagonal": [4, 4, 4],
                "nonzero_root_section_pairings": [
                    "Q1-A2(1).component1=-1",
                    "Q1-A2(2).component1=-1",
                    "Q2-A2(1).component2=-1",
                ],
            },
            "rank": 17,
            "signature": [17, 0],
            "determinant": 948,
            "root_system": "E8+A2^3",
            "root_count": len(roots),
            "roots_using_adjusted_section_coordinates": len(mixed_roots),
            "mordell_weil_rank_from_root_complement": 3,
            "discriminant_group": "Z/948Z",
            "discriminant_generator": target_generator,
            "matrix_sha256": canonical_sha256(target),
        },
        "transparent_neighbor_seed": {
            "construction": "A11 direct-sum K6(det=79)",
            "determinant": 948,
            "discriminant_group": "Z/948Z",
            "discriminant_generator": seed_generator,
            "matrix_sha256": canonical_sha256(seed),
        },
        "x6_79_period_lattice": {
            "gram_matrix": period,
            "signature": [2, 1],
            "determinant": -948,
            "discriminant_group": "Z/948Z",
            "discriminant_generator": period_generator,
            "ldl_diagonal": [fraction_string(value) for value in period_diagonal],
            "even_clifford_quaternion": "(6,2)",
            "hilbert_symbols": local_symbols,
            "quaternion_finite_ramification": [2, 3],
            "remaining_squarefree_level_factor": 79,
        },
        "finite_quadratic_form_isometries": {
            "modulus_2det": modulus,
            "target_generator_multiplier_to_period": target_multiplier,
            "seed_generator_multiplier_to_period": seed_multiplier,
            "target_congruence": "1709*67^2 == 485 (mod 1896)",
            "seed_congruence": "1589*101^2 == 485 (mod 1896)",
        },
        "shimura_anchor": {
            "curve": "u^2=16*t^6-19*t^4+88*t^2-48",
            "verified_rational_point": ["14/13", "16064/2197"],
            "N": 474,
            "factorization": "6*79",
            "primary_sources": [
                "https://arxiv.org/abs/0709.2908",
                "https://arxiv.org/abs/0802.1301",
            ],
            "source_claim": (
                "Elkies identifies this rational non-CM orbit with the "
                "rank-17 elliptic K3 construction on X(6,79)/<w_474>."
            ),
        },
        "interpretation": {
            "proved": (
                "The E8+A2^3 target is lattice-compatible with the exact "
                "X(6,79) period datum and with the existing neighbor seed; "
                "it is not merely a numerical determinant coincidence."
            ),
            "strong_inference": (
                "The target should be an alternate elliptic fibration in the "
                "same determinant-948 K3 neighbor landscape."
            ),
            "missing_bridge": (
                "Construct an explicit stable isometry/neighbor chain and "
                "transport the published rational Shimura point to the "
                "two-split Weierstrass chart with all three sections."
            ),
        },
        "claim_boundary": {
            "finite_quadratic_form_bridge_proved": True,
            "explicit_integral_neighbor_chain_from_target_proved": False,
            "explicit_moduli_map_proved": False,
            "rational_target_surface_found": False,
            "P3_found": False,
            "rank31_curve_found": False,
        },
        "implementation": {
            "language": "Python standard library",
            "script": "research/certify_e8_a2_shimura_bridge.py",
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
    if arguments.output:
        arguments.output.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(f"wrote {arguments.output}")
    if arguments.compare:
        committed = json.loads(arguments.compare.read_text(encoding="utf-8"))
        if committed != certificate:
            raise AssertionError(f"certificate mismatch: {arguments.compare}")
        print(f"matched {arguments.compare}")
    print("target essential lattice: det 948, roots E8+A2^3")
    print("discriminant form matches X(6,79) period and neighbor seed")
    print("Clifford quaternion ramifies exactly at 2 and 3; level factor 79")
    print(f"certificate sha256: {certificate['certificate_sha256']}")
    print("NO explicit moduli map, P3, or rank-31 curve")


if __name__ == "__main__":
    main()
