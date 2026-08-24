#!/usr/bin/env python3
"""Independent exact lower-bound certificate for ICARM curve #302.

The proof uses only rational arithmetic and explicitly enumerated finite
groups.  It verifies the 31 published points on the generalized Weierstrass
model, transports them to an integral short model, and maps them into a
deterministically selected product of quotients E(F_p)/2E(F_p).  Full binary
row rank proves that the points are Z-linearly independent.

No analytic-rank estimate, BSD, GRH, floating point arithmetic, SageMath,
Magma, or PARI/GP is used for the lower bound.  The separately reported
exact-rank statement under GRH+BSD is retained only as a conditional note.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = ROOT / "baseline" / "icarm_curve_302.json"
DEFAULT_CERTIFICATE_PATH = (
    ROOT / "baseline" / "rank31_icarm302_mod2_certificate.json"
)
RANK29_VERIFIER_PATH = ROOT / "baseline" / "verify_rank29_mod2.py"

EXPECTED_POINT_COUNT = 31
EXPECTED_MATRIX_RANK = 31
MIN_SEARCH_PRIME = 5
MAX_SEARCH_PRIME = 5000

RationalPoint = tuple[Fraction, Fraction]


def load_rank29_verifier():
    """Load the existing exact short-Weierstrass finite-group arithmetic."""

    spec = importlib.util.spec_from_file_location(
        "rank29_mod2_verifier_for_rank31", RANK29_VERIFIER_PATH
    )
    if spec is None or spec.loader is None:
        raise ImportError(RANK29_VERIFIER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_inputs() -> tuple[dict[str, object], list[int], list[RationalPoint]]:
    data = json.loads(INPUT_PATH.read_text(encoding="utf-8"))
    ainvs = [int(value) for value in data["a_invariants"]]
    if len(ainvs) != 5:
        raise AssertionError(f"expected five a-invariants, found {len(ainvs)}")
    points = [(Fraction(x), Fraction(y)) for x, y in data["points"]]
    if len(points) != EXPECTED_POINT_COUNT:
        raise AssertionError(
            f"expected {EXPECTED_POINT_COUNT} points, found {len(points)}"
        )
    boundary = data.get("reported_exact_rank_boundary")
    if boundary != {
        "claim": "rank equals 31 under GRH+BSD",
        "conditional_assumptions": ["GRH", "BSD"],
    }:
        raise AssertionError("unexpected conditional exact-rank boundary")
    return data, ainvs, points


def verify_generalized_rational_points(
    ainvs: list[int], points: Iterable[RationalPoint]
) -> None:
    a1, a2, a3, a4, a6 = ainvs
    for index, (x, y) in enumerate(points, start=1):
        residual = (
            y * y
            + a1 * x * y
            + a3 * y
            - (x * x * x + a2 * x * x + a4 * x + a6)
        )
        if residual != 0:
            raise AssertionError(f"P{index} is not on the curve: residual={residual}")


def short_model_coefficients(invariants: dict[str, int]) -> tuple[int, int]:
    """Return A,B for Y^2=X^3+A*X+B under the standard integral map."""

    return -27 * invariants["c4"], -54 * invariants["c6"]


def to_short_point(
    point: RationalPoint, a1: int, a3: int, b2: int
) -> RationalPoint:
    x, y = point
    return 36 * x + 3 * b2, 108 * (2 * y + a1 * x + a3)


def verify_short_model_map(
    ainvs: list[int],
    invariants: dict[str, int],
    short_a: int,
    short_b: int,
    points: Iterable[RationalPoint],
) -> None:
    a1, _a2, a3, _a4, _a6 = ainvs
    for index, point in enumerate(points, start=1):
        x_short, y_short = to_short_point(point, a1, a3, invariants["b2"])
        residual = y_short * y_short - (
            x_short * x_short * x_short + short_a * x_short + short_b
        )
        if residual != 0:
            raise AssertionError(
                f"P{index} does not map to the short model: residual={residual}"
            )


def reduce_short_point(
    verifier,
    point: RationalPoint,
    prime: int,
    a1: int,
    a3: int,
    b2: int,
):
    x, y = point
    if x.denominator % prime == 0 or y.denominator % prime == 0:
        raise AssertionError(f"point is nonintegral at p={prime}")
    x_mod = x.numerator * verifier.inverse_mod(x.denominator, prime) % prime
    y_mod = y.numerator * verifier.inverse_mod(y.denominator, prime) % prime
    short_x = (36 * x_mod + 3 * b2) % prime
    short_y = (108 * (2 * y_mod + a1 * x_mod + a3)) % prime
    return short_x, short_y


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor <= math.isqrt(value):
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def search_primes() -> Iterable[int]:
    for value in range(MIN_SEARCH_PRIME, MAX_SEARCH_PRIME + 1):
        if is_prime(value):
            yield value


def points_are_integral_at(points: Iterable[RationalPoint], prime: int) -> bool:
    return all(
        x.denominator % prime != 0 and y.denominator % prime != 0
        for x, y in points
    )


def find_local_independence_certificate(
    verifier,
    ainvs: list[int],
    invariants: dict[str, int],
    short_a: int,
    short_b: int,
    rational_points: list[RationalPoint],
) -> tuple[list[dict[str, object]], list[list[int]], int, list[int], dict[str, int]]:
    """Scan primes increasingly and stop at the first rank-31 product."""

    a1, _a2, a3, _a4, _a6 = ainvs
    matrix_rows: list[list[int]] = [[] for _ in rational_points]
    local_records: list[dict[str, object]] = []
    counters = {
        "primes_considered": 0,
        "bad_reduction_skipped": 0,
        "nonintegral_point_reduction_skipped": 0,
        "quotient_dimension_not_two_skipped": 0,
        "last_prime_considered": 0,
    }

    for prime in search_primes():
        counters["primes_considered"] += 1
        counters["last_prime_considered"] = prime
        if invariants["discriminant"] % prime == 0:
            counters["bad_reduction_skipped"] += 1
            continue
        if not points_are_integral_at(rational_points, prime):
            counters["nonintegral_point_reduction_skipped"] += 1
            continue

        finite_a = short_a % prime
        finite_b = short_b % prime
        group = verifier.enumerate_curve_points(prime, finite_a, finite_b)
        vector_map, basis, cosets = verifier.quotient_by_doubling(
            group, prime, finite_a
        )
        if len(cosets) != 4:
            counters["quotient_dimension_not_two_skipped"] += 1
            continue

        vectors: list[str] = []
        for index, rational_point in enumerate(rational_points):
            reduced = reduce_short_point(
                verifier,
                rational_point,
                prime,
                a1,
                a3,
                invariants["b2"],
            )
            if reduced not in vector_map:
                raise AssertionError(
                    f"P{index + 1} does not reduce to the short curve at p={prime}"
                )
            vector = vector_map[reduced]
            if len(vector) != 2:
                raise AssertionError("selected local quotient is not two-dimensional")
            matrix_rows[index].extend(vector)
            vectors.append("".join(str(bit) for bit in vector))

        rank, pivot_columns, _rref = verifier.gf2_rref(matrix_rows)
        local_records.append(
            {
                "prime": prime,
                "group_order": len(group),
                "quotient_dimension": 2,
                "quotient_basis_representatives": [
                    verifier.point_to_json(point) for point in basis
                ],
                "point_vectors": vectors,
                "matrix_rank_after_prime": rank,
            }
        )
        if rank == EXPECTED_MATRIX_RANK:
            return local_records, matrix_rows, rank, pivot_columns, counters

    rank, pivot_columns, _rref = verifier.gf2_rref(matrix_rows)
    raise AssertionError(
        f"prime scan through {MAX_SEARCH_PRIME} reached rank {rank}, not 31; "
        f"pivots={pivot_columns}"
    )


def find_torsion_reduction_pair(
    verifier,
    invariants: dict[str, int],
    short_a: int,
    short_b: int,
) -> tuple[list[dict[str, int]], dict[str, int]]:
    """Find the first safe coprime pair encountered in the prime scan."""

    good_records: list[dict[str, int]] = []
    primes_considered = 0
    for prime in search_primes():
        primes_considered += 1
        if invariants["discriminant"] % prime == 0:
            continue
        group_order = len(
            verifier.enumerate_curve_points(prime, short_a % prime, short_b % prime)
        )
        current = {"prime": prime, "group_order": group_order}
        for earlier in good_records:
            earlier_prime = earlier["prime"]
            earlier_order = earlier["group_order"]
            if math.gcd(earlier_order, group_order) != 1:
                continue
            # A torsion point of order equal to a residue characteristic must
            # be excluded using the other reduction prime.
            if group_order % earlier_prime == 0:
                continue
            if earlier_order % prime == 0:
                continue
            return [earlier, current], {
                "primes_considered": primes_considered,
                "last_prime_considered": prime,
            }
        good_records.append(current)
    raise AssertionError(
        f"no safe coprime torsion pair found through {MAX_SEARCH_PRIME}"
    )


def compute_certificate() -> dict[str, object]:
    verifier = load_rank29_verifier()
    data, ainvs, rational_points = load_inputs()
    a1, a2, a3, a4, a6 = ainvs

    invariants = verifier.generalized_invariants(a1, a2, a3, a4, a6)
    if invariants["discriminant"] == 0:
        raise AssertionError("singular curve")
    if str(invariants["discriminant"]) != data["reported_discriminant"]:
        raise AssertionError("reported discriminant mismatch")
    verify_generalized_rational_points(ainvs, rational_points)

    short_a, short_b = short_model_coefficients(invariants)
    verify_short_model_map(
        ainvs, invariants, short_a, short_b, rational_points
    )

    (
        local_records,
        matrix_rows,
        rank,
        pivot_columns,
        local_search,
    ) = find_local_independence_certificate(
        verifier,
        ainvs,
        invariants,
        short_a,
        short_b,
        rational_points,
    )
    if rank != EXPECTED_MATRIX_RANK:
        raise AssertionError(f"local image matrix has rank {rank}, expected 31")

    torsion_records, torsion_search = find_torsion_reduction_pair(
        verifier, invariants, short_a, short_b
    )
    first, second = torsion_records
    if math.gcd(first["group_order"], second["group_order"]) != 1:
        raise AssertionError("torsion reduction orders are not coprime")

    invariant_record = {
        name: str(invariants[name])
        for name in ("b2", "b4", "b6", "b8", "c4", "c6", "discriminant")
    }
    column_count = len(matrix_rows[0])
    payload: dict[str, object] = {
        "schema_version": 1,
        "candidate_id": data["candidate_id"],
        "claim": (
            "The 31 listed rational points are Z-linearly independent; "
            "rank E(Q) >= 31."
        ),
        "proof_type": "exact reduction to a product of E(F_p)/2E(F_p)",
        "source": data["source"],
        "input_payload_sha256": verifier.canonical_sha256(data),
        "curve": {
            "a_invariants": [str(value) for value in ainvs],
            "invariants": invariant_record,
            "reported_conductor": data["reported_conductor"],
            "short_model": {
                "equation": "Y^2 = X^3 + A*X + B",
                "A": str(short_a),
                "B": str(short_b),
                "map": {
                    "X": f"36*x+3*({invariants['b2']})",
                    "Y": f"108*(2*y+({a1})*x+({a3}))",
                },
                "formula": {
                    "X": "36*x+3*b2",
                    "Y": "108*(2*y+a1*x+a3)",
                    "A": "-27*c4",
                    "B": "-54*c6",
                },
            },
        },
        "point_count": len(rational_points),
        "point_membership": (
            "verified exactly over Q on both the generalized and short models"
        ),
        "torsion_certificate": {
            "reductions": torsion_records,
            "search": {
                "order": "primes p >= 5 in increasing order",
                "maximum_prime": MAX_SEARCH_PRIME,
                **torsion_search,
            },
            "group_order_gcd": math.gcd(
                first["group_order"], second["group_order"]
            ),
            "cross_characteristic_exclusions_verified": True,
            "conclusion": "E(Q)_tors is trivial",
            "argument": (
                f"At the good primes {first['prime']} and {second['prime']} the "
                f"group orders are {first['group_order']} and "
                f"{second['group_order']}, with gcd 1. A rational torsion point "
                "of prime order ell injects under good reduction away from ell. "
                "Orders ell different from both residue characteristics would "
                "divide the gcd; the two residue-characteristic cases are "
                "excluded by the opposite reduction, as checked explicitly."
            ),
        },
        "local_prime_search": {
            "order": "primes p >= 5 in increasing order",
            "maximum_prime": MAX_SEARCH_PRIME,
            "selection": (
                "good reduction, every listed point affine-integral at p, and "
                "dim_F2 E(F_p)/2E(F_p) = 2"
            ),
            "termination": "first selected-prime prefix with binary row rank 31",
            **local_search,
        },
        "local_quotients": local_records,
        "binary_matrix": {
            "orientation": (
                f"31 point rows by {column_count} local-quotient coordinates"
            ),
            "rows": ["".join(str(bit) for bit in row) for row in matrix_rows],
            "row_count": len(matrix_rows),
            "column_count": column_count,
            "rank": rank,
            "pivot_columns_zero_based": pivot_columns,
        },
        "independence_argument": (
            "Any integral relation reduces in every quotient E(F_p)/2E(F_p). "
            "The binary matrix has row rank 31, so every coefficient is even. "
            "Halving the coefficients gives a point killed by 2; trivial "
            "rational torsion makes it O. Infinite descent forces every "
            "coefficient to be zero."
        ),
        "conditional_assumptions": [],
        "conditional_exact_rank_note": (
            "The separate reported assertion rank E(Q)=31 uses GRH+BSD and is "
            "not part of this unconditional lower-bound certificate."
        ),
        "implementation": {
            "language": "Python standard library",
            "script": "baseline/verify_rank31_icarm302.py",
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "shared_exact_group_arithmetic": "baseline/verify_rank29_mod2.py",
            "shared_exact_group_arithmetic_sha256": hashlib.sha256(
                RANK29_VERIFIER_PATH.read_bytes()
            ).hexdigest(),
        },
    }
    payload["certificate_sha256"] = verifier.canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-certificate", type=Path)
    parser.add_argument(
        "--certificate", type=Path, default=DEFAULT_CERTIFICATE_PATH
    )
    args = parser.parse_args()

    computed = compute_certificate()
    if args.write_certificate is not None:
        args.write_certificate.parent.mkdir(parents=True, exist_ok=True)
        args.write_certificate.write_text(
            json.dumps(computed, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(f"wrote {args.write_certificate}")
    else:
        committed = json.loads(args.certificate.read_text(encoding="utf-8"))
        if committed != computed:
            raise AssertionError("committed certificate does not match recomputation")

    matrix = computed["binary_matrix"]
    torsion = computed["torsion_certificate"]["reductions"]
    local_primes = [record["prime"] for record in computed["local_quotients"]]
    print("exact generalized-model point checks: 31/31")
    print("exact short-model image checks: 31/31")
    print(
        "torsion: trivial "
        f"(#E(F_{torsion[0]['prime']})={torsion[0]['group_order']}, "
        f"#E(F_{torsion[1]['prime']})={torsion[1]['group_order']})"
    )
    print(f"local quotient primes: {local_primes}")
    print(
        f"local quotient matrix: {matrix['row_count']}x{matrix['column_count']}, "
        f"rank={matrix['rank']}"
    )
    print(f"certificate sha256: {computed['certificate_sha256']}")
    print("UNCONDITIONAL RESULT: rank E(Q) >= 31")
    print("CONDITIONAL SEPARATION: exact rank 31 requires GRH+BSD")


if __name__ == "__main__":
    main()
