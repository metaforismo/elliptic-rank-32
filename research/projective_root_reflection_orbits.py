#!/usr/bin/env python3
"""Exact F2 line orbits under explicit integral coordinate-root reflections.

This finite symmetry audit does not construct neighbors or count lattice
isometry classes. In particular the bad-prime p=2 boundary is still open.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from research import sample_projective_lines as sampler  # noqa: E402
from research import search_e8_a2_target_neighbor_bridge as source  # noqa: E402


def quadratic_masks(gram):
    return [
        sum((gram[i][j] % 2) << j for j in range(i))
        | (((gram[i][i] // 2) % 2) << i)
        for i in range(len(gram))
    ]


def q_mod2(mask, rows):
    result = 0
    remaining = mask
    while remaining:
        bit = remaining & -remaining
        result ^= (rows[bit.bit_length() - 1] & mask).bit_count() % 2
        remaining ^= bit
    return result


def coordinate_reflection(gram, root_index):
    """s_i(v)=v-(v,e_i)e_i, for the norm-two integral root e_i."""
    size = len(gram)
    if gram[root_index][root_index] != 2:
        raise ValueError("coordinate vector is not a norm-two root")
    matrix = [[int(i == j) for j in range(size)] for i in range(size)]
    for j in range(size):
        matrix[root_index][j] -= gram[root_index][j]
    product = [[sum(gram[i][k] * matrix[k][j] for k in range(size)) for j in range(size)] for i in range(size)]
    transported = [[sum(matrix[k][i] * product[k][j] for k in range(size)) for j in range(size)] for i in range(size)]
    if transported != gram:
        raise AssertionError("root reflection failed its exact Gram identity")
    square = [[sum(matrix[i][k] * matrix[k][j] for k in range(size)) for j in range(size)] for i in range(size)]
    if square != [[int(i == j) for j in range(size)] for i in range(size)]:
        raise AssertionError("root reflection is not an integral involution")
    return matrix


def enumerate_orbits(gram):
    sampler.validate_gram(gram)
    dimension = len(gram)
    if dimension > 20:
        raise ValueError("exact F2 enumeration is bounded to dimension at most 20")
    qrows = quadratic_masks(gram)
    isotropic = {mask for mask in range(1, 1 << dimension) if q_mod2(mask, qrows) == 0}
    gram_rows = [sum((entry % 2) << j for j, entry in enumerate(row)) for row in gram]
    root_indices = [i for i in range(dimension) if gram[i][i] == 2]
    reflections = [coordinate_reflection(gram, i) for i in root_indices]
    generators = [(1 << i, gram_rows[i]) for i in root_indices]
    labels = {}
    orbits = []
    for representative in sorted(isotropic):
        if representative in labels:
            continue
        orbit_index = len(orbits)
        labels[representative] = orbit_index
        queue = [representative]
        for current in queue:
            for bit, row in generators:
                neighbor = current ^ (bit if (row & current).bit_count() % 2 else 0)
                if neighbor not in isotropic:
                    raise AssertionError("reflection left the nonzero isotropic set")
                if neighbor not in labels:
                    labels[neighbor] = orbit_index
                    queue.append(neighbor)
                elif labels[neighbor] != orbit_index:
                    raise AssertionError("orbit partition is not disjoint")
        orbits.append({
            "representative_mask": representative,
            "representative_vector": [(representative >> i) & 1 for i in range(dimension)],
            "size": len(queue),
            "is_bilinear_radical_line": all((row & representative).bit_count() % 2 == 0 for row in gram_rows),
        })
    if set(labels) != isotropic or sum(item["size"] for item in orbits) != len(isotropic):
        raise AssertionError("incomplete isotropic-line orbit partition")
    return {
        "dimension": dimension,
        "gram_sha256": sampler.payload_sha256(gram),
        "nonzero_isotropic_lines": len(isotropic),
        "coordinate_root_indices_zero_based": root_indices,
        "integral_reflection_matrices": reflections,
        "orbit_count": len(orbits),
        "orbit_size_histogram": {str(size): count for size, count in sorted(Counter(item["size"] for item in orbits).items())},
        "all_nonzero_isotropic_masks_sha256": sampler.payload_sha256(sorted(isotropic)),
        "complete_mask_to_orbit_index_sha256": sampler.payload_sha256(sorted(labels.items())),
        "explicit_orbits": orbits if generators else None,
        "singleton_orbit_rule": None if generators else "No coordinate-root generators: every nonzero isotropic mask is its own orbit.",
        "all_generators_are_exact_integral_involutory_isometries": True,
        "all_orbits_closed_and_partition_all_nonzero_isotropic_lines": True,
    }, labels


def build_certificate():
    summaries = {}
    seed = "rank32-20260907-v1"
    for name, gram in (
        ("origin", source.target_essential_lattice()),
        ("transparent", source.transparent_seed_lattice()),
        ("rootless", source.FROZEN_ROOTLESS_GRAM),
    ):
        record, labels = enumerate_orbits(gram)
        sampled = sampler.sample_window(gram, 2, seed, 0, 32, 4096)
        sampled_masks = [sum(value << i for i, value in enumerate(line["vector"])) for line in sampled["selected_lines"]]
        record["initial_hash_window"] = {
            "seed": seed,
            "window_sha256": sampled["window_sha256"],
            "selected_line_count": 32,
            "distinct_reflection_subgroup_orbits_touched": len({labels[mask] for mask in sampled_masks}),
        }
        summaries[name] = record
    paths = [
        "research/projective_root_reflection_orbits.py",
        "research/sample_projective_lines.py",
        "research/search_e8_a2_target_neighbor_bridge.py",
        "research/certify_e8_a2_shimura_bridge.py",
    ]
    record = {
        "schema_version": 1,
        "certificate_kind": "exact-initial-F2-coordinate-root-reflection-orbits",
        "prime": 2,
        "vector_encoding": "bit i is coordinate i, both zero-based; F2 projective lines are nonzero vectors",
        "initial_lattices": summaries,
        "source_sha256": {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest() for path in paths},
        "claim_boundary": {
            "only_the_subgroup_generated_by_listed_coordinate_root_reflections": True,
            "full_integral_automorphism_group_computed": False,
            "lattice_isometry_classes_counted": False,
            "neighbors_constructed": False,
            "p2_neighbor_enumeration_complete": False,
            "rank32_curve_found": False,
            "next_step": "Test orbit representatives with exact neighbor and bad-prime gates; line symmetry alone does not establish complete neighbor coverage.",
        },
    }
    return {**record, "record_sha256": sampler.payload_sha256(record)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-certificate", type=Path)
    parser.add_argument("--compare", type=Path)
    args = parser.parse_args()
    record = build_certificate()
    if args.compare and json.loads(args.compare.read_text()) != record:
        raise SystemExit("certificate mismatch")
    if args.write_certificate:
        with args.write_certificate.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "record_sha256": record["record_sha256"],
        "summary": {name: {key: item[key] for key in ("nonzero_isotropic_lines", "orbit_count", "initial_hash_window")} for name, item in record["initial_lattices"].items()},
        "claim_boundary": record["claim_boundary"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
