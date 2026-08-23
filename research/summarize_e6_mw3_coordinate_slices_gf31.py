#!/usr/bin/env python3
"""Summarize the exact P2/P3 search across five GF(31) coordinate slices."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

from replay_e6_mw3_p2_split_cores import (
    canonical_sha256,
    canonical_sign,
    search_p2,
    search_p3_polynomial,
)


ROOT = Path(__file__).resolve().parents[1]
SLICE_CERTIFICATES = [
    ROOT / "certificates" / "e6_mw3_core_slice_r1_s1_x1_gf31.json",
    ROOT / "certificates" / "e6_mw3_core_slice_r2_s1_x1_gf31.json",
    ROOT / "certificates" / "e6_mw3_core_slice_r3_s1_x1_gf31.json",
    ROOT / "certificates" / "e6_mw3_core_slice_r1_s2_x1_gf31.json",
    ROOT / "certificates" / "e6_mw3_core_slice_r1_s1_x2_gf31.json",
]
RELATION_CERTIFICATE = (
    ROOT
    / "certificates"
    / "e6_mw3_candidate_r3_s1_x1_core_1_0_14_23_relations.json"
)
DEGENERATE_CERTIFICATE = (
    ROOT / "certificates" / "e6_mw3_degenerate_pivots_gf31.json"
)
DEFAULT_OUTPUT = ROOT / "certificates" / "e6_mw3_coordinate_slices_gf31_summary.json"


def geometric_p2_key(hit):
    return (
        hit["pole"],
        hit["q0"],
        tuple(hit["X2"]),
        tuple(canonical_sign(hit["Y2"])),
    )


def compute():
    relation = json.loads(RELATION_CERTIFICATE.read_text())
    degenerate = json.loads(DEGENERATE_CERTIFICATE.read_text())
    slice_records = []
    total_quadruples = 0
    total_core_solutions = 0
    total_surfaces = 0
    total_ordered_p2_cases = 0
    total_unique_p2_cases = 0
    total_signed_ordered_hits = 0
    total_geometric_hits = 0
    total_p3_cases = 0
    total_p3_candidates = 0
    joint_candidates = []

    for path in SLICE_CERTIFICATES:
        source = json.loads(path.read_text())
        coordinates = source["coordinate_slice"]
        counts = source["counts"]
        total_quadruples += counts["quadruples_examined"]
        total_core_solutions += counts["simultaneous_core_solutions"]
        surfaces = []
        for record in source["valid_two_I2_surfaces"]:
            total_surfaces += 1
            roots = record["valid_I2_roots"]
            if len(roots) != 2:
                raise AssertionError(f"expected two I2 roots, got {roots}")
            nodes = {
                int(key): value for key, value in record["singular_nodes"].items()
            }
            A, B = record["A"], record["B"]
            s0 = coordinates["s0"]
            ordered_hits = []
            ordered_cases = 0
            for lam, mu in (tuple(roots), tuple(reversed(roots))):
                tested, hits = search_p2(
                    A,
                    B,
                    lam,
                    mu,
                    nodes[lam][0],
                    nodes[mu][0],
                    s0=s0,
                )
                ordered_cases += tested
                ordered_hits.extend(hits)
            unique_hits = {
                geometric_p2_key(hit): hit for hit in ordered_hits
            }
            p3_cases, p3_hits = search_p3_polynomial(
                A,
                B,
                record["core"][3],
                roots[0],
                roots[1],
                nodes[roots[0]][0],
                nodes[roots[1]][0],
                s0=s0,
            )
            total_ordered_p2_cases += ordered_cases
            total_unique_p2_cases += ordered_cases // 2
            total_signed_ordered_hits += len(ordered_hits)
            total_geometric_hits += len(unique_hits)
            total_p3_cases += p3_cases
            total_p3_candidates += len(p3_hits)
            surface = {
                "core": record["core"],
                "I2_roots": roots,
                "ordered_P2_cases": ordered_cases,
                "geometric_P2_hits": list(unique_hits.values()),
                "P3_cases": p3_cases,
                "P3_candidates": p3_hits,
            }
            if unique_hits and p3_hits:
                candidate = {
                    "coordinate_slice": coordinates,
                    "core": record["core"],
                    "geometric_P2_hits": list(unique_hits.values()),
                    "P3_candidates": p3_hits,
                }
                if (
                    relation["coordinate_slice"] == coordinates
                    and relation["core_a1_a2_a4_s1"] == record["core"]
                ):
                    candidate["bounded_relation_certificate_sha256"] = relation[
                        "certificate_sha256"
                    ]
                    candidate["relations"] = relation["relations"]
                joint_candidates.append(candidate)
            surfaces.append(surface)
        slice_records.append(
            {
                "path": str(path.relative_to(ROOT)),
                "certificate_sha256": source["certificate_sha256"],
                "coordinate_slice": coordinates,
                "counts": counts,
                "surfaces": surfaces,
            }
        )

    if len(joint_candidates) != 1:
        raise AssertionError(f"expected one joint P2/P3 candidate, got {len(joint_candidates)}")
    if [1, 0, 2] not in joint_candidates[0].get("relations", []):
        raise AssertionError("expected exact dependency P1+2*P3=O")

    payload = {
        "schema": "elliptic-rank30/e6-mw3-coordinate-slices-gf31/v1",
        "field": "GF(31)",
        "slice_records": slice_records,
        "relation_certificate": str(RELATION_CERTIFICATE.relative_to(ROOT)),
        "relation_certificate_sha256": relation["certificate_sha256"],
        "degenerate_pivot_certificate": str(
            DEGENERATE_CERTIFICATE.relative_to(ROOT)
        ),
        "degenerate_pivot_certificate_sha256": degenerate[
            "certificate_sha256"
        ],
        "degenerate_pivot_totals": degenerate["totals"],
        "totals": {
            "coordinate_slices": len(SLICE_CERTIFICATES),
            "core_quadruples_examined": total_quadruples,
            "simultaneous_core_solutions": total_core_solutions,
            "valid_two_I2_surfaces": total_surfaces,
            "ordered_P2_cases": total_ordered_p2_cases,
            "unique_P2_cases_up_to_I2_label_swap": total_unique_p2_cases,
            "signed_ordered_P2_hits": total_signed_ordered_hits,
            "geometric_P2_hits": total_geometric_hits,
            "P3_quadratic_X_cases": total_p3_cases,
            "P3_candidates": total_p3_candidates,
            "joint_P2_P3_candidates": len(joint_candidates),
            "independent_rank3_seeds_after_exact_relation_filter": 0,
        },
        "joint_candidates": joint_candidates,
        "conclusion": (
            "The only joint P2/P3 candidate satisfies P1+2*P3=O; these five "
            "direct coordinate slices and their exceptional pivot planes contain "
            "no surviving MW-rank-3 seed."
        ),
        "claim_boundary": (
            "Complete for the regular charts in the five listed GF(31) coordinate "
            "slices and, through the linked certificate, for the exceptional "
            "y0+y1=0 planes in those five slices plus the public slice. Other "
            "slices, primes, and K3 neighbors remain open."
        ),
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    payload = compute()
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"totals": payload["totals"], "certificate_sha256": payload["certificate_sha256"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
