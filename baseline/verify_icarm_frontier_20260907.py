#!/usr/bin/env python3
"""Replay the four rank >= 30 entries in the 2026-09-07 ICARM projection.

Only exact arithmetic is used. Each lower bound is proved independently by
point substitution, trivial torsion, and full row rank in finite mod-2
quotients. No analytic upper bound or nonexistence of rank 32 is proved.
The optional full-database input additionally audits the source projection.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "baseline/icarm_frontier_20260907.json"
CERTIFICATE = ROOT / "certificates/icarm_frontier_20260907.json"
HELPER_PATH = ROOT / "baseline/verify_rank31_icarm302.py"
SOURCE_URL = "https://elliptic-rank.icarm.cloud/database.json"


def load_helpers():
    spec = importlib.util.spec_from_file_location("frontier_rank31_helpers", HELPER_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(HELPER_PATH)
    helpers = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helpers)
    return helpers, helpers.load_rank29_verifier()


def validate_snapshot(data: dict) -> None:
    if data["schema_version"] != 1 or data["source"]["url"] != SOURCE_URL:
        raise AssertionError("unexpected snapshot schema or source")
    if data["retrieved_date_utc"] != "2026-09-07":
        raise AssertionError("unexpected snapshot date")
    if [c["id"] for c in data["curves"]] != [273, 302, 398, 582]:
        raise AssertionError("unexpected or duplicate frontier curve IDs")
    histogram = {int(k): v for k, v in data["source_rank_histogram"].items()}
    if any(type(v) is not int or v < 1 for v in histogram.values()):
        raise AssertionError("invalid source rank histogram")
    if sum(histogram.values()) != data["source_curve_count"]:
        raise AssertionError("source census mismatch")
    selected = Counter(c["rank_lower_bound"] for c in data["curves"])
    if selected != Counter({k: v for k, v in histogram.items() if k >= 30}):
        raise AssertionError("incomplete high-rank projection")
    baseline = json.loads((ROOT / "baseline/icarm_curve_302.json").read_text())
    record = next(c for c in data["curves"] if c["id"] == 302)
    if (record["ainvs"] != baseline["a_invariants"]
            or record["points"] != baseline["points"]
            or record["discriminant"] != baseline["reported_discriminant"]):
        raise AssertionError("rank-31 baseline changed")


def audit_source_database(data: dict, raw: bytes) -> None:
    if hashlib.sha256(raw).hexdigest() != data["raw_database_sha256"]:
        raise AssertionError("raw source SHA-256 mismatch")
    if len(raw) != data["raw_database_bytes"]:
        raise AssertionError("raw source length mismatch")
    source = json.loads(raw)
    curves = source["curves"]
    if source["count"] != len(curves) or len({c["id"] for c in curves}) != len(curves):
        raise AssertionError("incomplete or duplicate source records")
    if source["count"] != data["source_curve_count"]:
        raise AssertionError("source count changed")
    if {k: source[k] for k in data["source"]} != data["source"]:
        raise AssertionError("source attribution changed")
    histogram = Counter(c["rank_lower_bound"] for c in curves)
    if histogram != Counter({int(k): v for k, v in data["source_rank_histogram"].items()}):
        raise AssertionError("source histogram changed")
    selected = sorted((c for c in curves if c["rank_lower_bound"] >= 30), key=lambda c: c["id"])
    if len(selected) != len(data["curves"]):
        raise AssertionError("source frontier count changed")
    for actual, frozen in zip(selected, data["curves"], strict=True):
        if {k: actual[k] for k in frozen} != frozen:
            raise AssertionError(f"source projection changed for #{frozen['id']}")


def certify_curve(record: dict, helpers, group) -> dict:
    ainvs = [int(value) for value in record["ainvs"]]
    if len(ainvs) != 5:
        raise AssertionError("five a-invariants required")
    points = [(Fraction(x), Fraction(y)) for x, y in record["points"]]
    count = record["rank_lower_bound"]
    if type(count) is not int or count < 1 or len(points) != count:
        raise AssertionError("point count does not match the claimed lower bound")
    if len({x for x, _ in points}) != count:
        raise AssertionError("equal x-coordinates cannot give independent points")
    inv = group.generalized_invariants(*ainvs)
    if not inv["discriminant"] or str(inv["discriminant"]) != record["discriminant"]:
        raise AssertionError("singular curve or incorrect discriminant")
    helpers.verify_generalized_rational_points(ainvs, points)
    short_a, short_b = helpers.short_model_coefficients(inv)
    helpers.verify_short_model_map(ainvs, inv, short_a, short_b, points)
    torsion, _ = helpers.find_torsion_reduction_pair(group, inv, short_a, short_b)
    rows: list[list[int]] = [[] for _ in points]
    local = []
    rank = 0
    for prime in helpers.search_primes():
        if inv["discriminant"] % prime == 0 or not helpers.points_are_integral_at(points, prime):
            continue
        finite_group = group.enumerate_curve_points(prime, short_a % prime, short_b % prime)
        vectors, basis, cosets = group.quotient_by_doubling(finite_group, prime, short_a % prime)
        if len(cosets) != 4:
            continue
        images = []
        for row, point in zip(rows, points, strict=True):
            reduced = helpers.reduce_short_point(group, point, prime, ainvs[0], ainvs[2], inv["b2"])
            vector = vectors[reduced]
            if len(vector) != 2:
                raise AssertionError("unexpected local quotient dimension")
            row.extend(vector)
            images.append("".join(map(str, vector)))
        rank, pivots, _ = group.gf2_rref(rows)
        local.append({
            "prime": prime,
            "group_order": len(finite_group),
            "basis": [group.point_to_json(point) for point in basis],
            "point_vectors": images,
            "matrix_rank_after_prime": rank,
        })
        if rank == count:
            break
    if rank != count:
        raise AssertionError(f"#{record['id']}: local image rank {rank}, expected {count}")
    return {
        "id": record["id"],
        "source_url": f"https://elliptic-rank.icarm.cloud/curve/{record['id']}",
        "source_record_sha256": group.canonical_sha256(record),
        "unconditional_lower_bound": rank,
        "point_count": count,
        "point_membership_both_models": True,
        "discriminant": str(inv["discriminant"]),
        "j_invariant": str(Fraction(inv["c4"] ** 3, inv["discriminant"])),
        "torsion_reductions": torsion,
        "torsion_trivial": True,
        "local_quotients": local,
        "binary_matrix": {
            "rows": ["".join(map(str, row)) for row in rows],
            "row_count": count,
            "column_count": len(rows[0]),
            "rank": rank,
            "pivot_columns_zero_based": pivots,
        },
        "exact_rank_upper_bound_proved": False,
    }


def compute_certificate(data: dict) -> dict:
    validate_snapshot(data)
    helpers, group = load_helpers()
    records = [certify_curve(c, helpers, group) for c in data["curves"]]
    payload = {
        "schema_version": 1,
        "snapshot": str(INPUT.relative_to(ROOT)),
        "snapshot_payload_sha256": group.canonical_sha256(data),
        "raw_database_sha256": data["raw_database_sha256"],
        "retrieved_date_utc": data["retrieved_date_utc"],
        "source_curve_count": data["source_curve_count"],
        "source_maximum_rank_lower_bound": max(map(int, data["source_rank_histogram"])),
        "curves": records,
        "all_j_invariants_pairwise_distinct": len({c["j_invariant"] for c in records}) == len(records),
        "proof": "Full binary row rank forces relation coefficients to be even; trivial rational torsion allows infinite descent. Hence each listed set is Z-linearly independent.",
        "claim_boundary": "The source census is a dated observation, not nonexistence of rank 32. Distinct j-invariants exclude isomorphisms and twists between these four curves, but do not exclude a common parameter family. No upper bound or saturation claim is made.",
        "rank32_curve_found": False,
        "implementation_sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (Path(__file__).resolve(), HELPER_PATH, helpers.RANK29_VERIFIER_PATH)
        },
    }
    payload["certificate_sha256"] = group.canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-database", type=Path)
    parser.add_argument("--write-certificate", type=Path)
    args = parser.parse_args()
    data = json.loads(INPUT.read_text())
    validate_snapshot(data)
    if args.source_database:
        audit_source_database(data, args.source_database.read_bytes())
    computed = compute_certificate(data)
    if args.write_certificate:
        with args.write_certificate.open("x", encoding="utf-8") as output:
            output.write(json.dumps(computed, indent=2, ensure_ascii=False) + "\n")
    elif computed != json.loads(CERTIFICATE.read_text()):
        raise AssertionError("stored certificate differs from exact replay")
    print(json.dumps({
        "status": "PASS",
        "full_source_projection_verified": args.source_database is not None,
        "curves": [{"id": c["id"], "rank_lower_bound": c["unconditional_lower_bound"],
                    "matrix_columns": c["binary_matrix"]["column_count"]} for c in computed["curves"]],
        "certificate_sha256": computed["certificate_sha256"],
        "rank32_curve_found": False,
    }, indent=2))


if __name__ == "__main__":
    main()
