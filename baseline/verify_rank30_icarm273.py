#!/usr/bin/env python3
"""Independent exact certificate for the 30 points on ICARM curve #273.

The proof uses only rational arithmetic and finite groups.  It verifies the
30 points exactly and maps them into a product of quotients
E(F_p)/2E(F_p).  The resulting 30 by 30 binary matrix has full rank.

No analytic-rank estimate, BSD, GRH, floating point arithmetic, SageMath,
Magma, or PARI/GP is used for the lower bound.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = ROOT / "baseline" / "icarm_curve_273.json"
DEFAULT_CERTIFICATE_PATH = ROOT / "baseline" / "rank30_icarm273_mod2_certificate.json"
RANK29_VERIFIER_PATH = ROOT / "baseline" / "verify_rank29_mod2.py"

LOCAL_PRIMES = (43, 61, 101, 211, 223, 241, 263, 271, 283, 311, 313, 457, 521, 569, 571)
EXPECTED_LOCAL_ORDERS = {
    43: 56,
    61: 72,
    101: 120,
    211: 212,
    223: 248,
    241: 260,
    263: 292,
    271: 296,
    283: 304,
    311: 344,
    313: 336,
    457: 464,
    521: 516,
    569: 600,
    571: 604,
}
TORSION_PRIMES = (11, 59)
EXPECTED_TORSION_ORDERS = {11: 18, 59: 73}
EXPECTED_MATRIX_RANK = 30


def load_rank29_verifier():
    spec = importlib.util.spec_from_file_location(
        "rank29_mod2_verifier", RANK29_VERIFIER_PATH
    )
    if spec is None or spec.loader is None:
        raise ImportError(RANK29_VERIFIER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_inputs() -> tuple[dict[str, object], list[int], list[tuple[Fraction, Fraction]]]:
    data = json.loads(INPUT_PATH.read_text(encoding="utf-8"))
    ainvs = [int(value) for value in data["a_invariants"]]
    if ainvs[:3] != [1, 0, 0]:
        raise AssertionError(f"unexpected a1,a2,a3: {ainvs[:3]}")
    points = [(Fraction(x), Fraction(y)) for x, y in data["points"]]
    if len(points) != 30:
        raise AssertionError(f"expected 30 points, found {len(points)}")
    return data, ainvs, points


def compute_certificate() -> dict[str, object]:
    verifier = load_rank29_verifier()
    data, ainvs, rational_points = load_inputs()
    a1, a2, a3, a4, a6 = ainvs

    invariants = verifier.generalized_invariants(a1, a2, a3, a4, a6)
    if invariants["discriminant"] == 0:
        raise AssertionError("singular curve")
    if str(invariants["discriminant"]) != data["reported_discriminant"]:
        raise AssertionError("reported discriminant mismatch")
    verifier.verify_rational_points(a4, a6, rational_points)

    short_a, short_b = verifier.short_model_coefficients(a4, a6)
    local_records: list[dict[str, object]] = []
    matrix_rows: list[list[int]] = [[] for _ in rational_points]

    for prime in LOCAL_PRIMES:
        if invariants["discriminant"] % prime == 0:
            raise AssertionError(f"p={prime} is not a good-reduction prime")
        finite_a = short_a % prime
        finite_b = short_b % prime
        group = verifier.enumerate_curve_points(prime, finite_a, finite_b)
        if len(group) != EXPECTED_LOCAL_ORDERS[prime]:
            raise AssertionError(
                f"#E(F_{prime})={len(group)}, expected {EXPECTED_LOCAL_ORDERS[prime]}"
            )
        vector_map, basis, cosets = verifier.quotient_by_doubling(
            group, prime, finite_a
        )
        if len(cosets) != 4:
            raise AssertionError(
                f"E(F_{prime})/2E(F_{prime}) has size {len(cosets)}, expected 4"
            )

        vectors: list[str] = []
        for index, rational_point in enumerate(rational_points):
            reduced = verifier.reduce_rational_point(rational_point, prime)
            if reduced not in vector_map:
                raise AssertionError(
                    f"P{index + 1} does not reduce to E(F_{prime})"
                )
            vector = vector_map[reduced]
            matrix_rows[index].extend(vector)
            vectors.append("".join(str(bit) for bit in vector))

        local_records.append(
            {
                "prime": prime,
                "group_order": len(group),
                "quotient_dimension": 2,
                "quotient_basis_representatives": [
                    verifier.point_to_json(point) for point in basis
                ],
                "point_vectors": vectors,
            }
        )

    rank, pivot_columns, _rref = verifier.gf2_rref(matrix_rows)
    if rank != EXPECTED_MATRIX_RANK:
        raise AssertionError(f"local image matrix has rank {rank}, expected 30")
    if pivot_columns != list(range(30)):
        raise AssertionError(f"unexpected pivot columns: {pivot_columns}")

    torsion_records: list[dict[str, int]] = []
    for prime in TORSION_PRIMES:
        if invariants["discriminant"] % prime == 0:
            raise AssertionError(f"torsion prime p={prime} has bad reduction")
        group_order = len(
            verifier.enumerate_curve_points(
                prime, short_a % prime, short_b % prime
            )
        )
        if group_order != EXPECTED_TORSION_ORDERS[prime]:
            raise AssertionError(
                f"#E(F_{prime})={group_order}, expected {EXPECTED_TORSION_ORDERS[prime]}"
            )
        torsion_records.append({"prime": prime, "group_order": group_order})

    payload: dict[str, object] = {
        "schema_version": 1,
        "candidate_id": data["candidate_id"],
        "claim": (
            "The 30 listed rational points are Z-linearly independent; "
            "rank E(Q) >= 30."
        ),
        "proof_type": "exact reduction to a product of E(F_p)/2E(F_p)",
        "source": data["source"],
        "curve": {
            "a_invariants": [str(value) for value in ainvs],
            "discriminant": str(invariants["discriminant"]),
            "reported_conductor": data["reported_conductor"],
            "short_model": {
                "equation": "Y^2 = X^3 + A*X + B",
                "A": str(short_a),
                "B": str(short_b),
                "map": {"X": "36*x+3", "Y": "108*(2*y+x)"},
            },
        },
        "point_count": len(rational_points),
        "point_membership": "verified exactly over Q",
        "torsion_certificate": {
            "reductions": torsion_records,
            "conclusion": "E(Q)_tors is trivial",
            "argument": (
                "At the good primes 11 and 59 the group orders are 18 and 73. "
                "A rational torsion point of prime order ell injects at good "
                "reduction primes away from ell. If ell is neither 11 nor 59, "
                "then ell divides gcd(18,73)=1. The cases ell=11 and ell=59 "
                "are excluded by reduction at 59 and 11 respectively."
            ),
        },
        "local_quotients": local_records,
        "binary_matrix": {
            "orientation": "30 point rows by 30 local-quotient coordinates",
            "rows": ["".join(str(bit) for bit in row) for row in matrix_rows],
            "rank": rank,
            "pivot_columns_zero_based": pivot_columns,
        },
        "independence_argument": (
            "Any integral relation reduces in every quotient E(F_p)/2E(F_p). "
            "The binary matrix has row rank 30, so every coefficient is even. "
            "Halving the coefficients gives a point killed by 2; trivial "
            "rational torsion makes it O. Infinite descent forces every "
            "coefficient to be zero."
        ),
        "conditional_assumptions": [],
        "conditional_exact_rank_note": (
            "The separate assertion rank E(Q)=30 uses GRH+BSD and is not part "
            "of this lower-bound certificate."
        ),
        "implementation": {
            "language": "Python standard library",
            "script": "baseline/verify_rank30_icarm273.py",
            "shared_exact_group_arithmetic": "baseline/verify_rank29_mod2.py",
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
    print("exact point checks: 30/30")
    print("torsion: trivial (#E(F_11)=18, #E(F_59)=73)")
    print(f"local quotient matrix: 30x30, rank={matrix['rank']}")
    print(f"certificate sha256: {computed['certificate_sha256']}")
    print("UNCONDITIONAL RESULT: rank E(Q) >= 30")


if __name__ == "__main__":
    main()
