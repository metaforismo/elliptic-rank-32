#!/usr/bin/env python3
"""Certify all joint P1/P2/P3 triples in a modular-slice certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from certify_e6_mw3_gf29_independence import find_certificate_for_p3
from replay_e6_mw3_p2_split_cores import canonical_sha256


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "certificates" / "e6_mw3_modular_slices_gf23.json"


def compute(input_path, ell_candidates):
    source = json.loads(input_path.read_text())
    prime = source["field_prime"]
    records = []
    for joint_index, joint in enumerate(source["joint_candidates"], start=1):
        surface = joint["surface"]
        for p2_index, p2_data in enumerate(
            surface["geometric_P2_hits"], start=1
        ):
            for p3_index, p3_data in enumerate(
                surface["P3_candidates"], start=1
            ):
                result = find_certificate_for_p3(
                    surface,
                    p2_data,
                    p3_data,
                    prime,
                    tuple(ell_candidates),
                )
                records.append(
                    {
                        "joint_index": joint_index,
                        "branch": joint["branch"],
                        "coordinate_slice": joint["coordinate_slice"],
                        "core": joint["core"],
                        "P2_index": p2_index,
                        "P2": p2_data,
                        "P3_index": p3_index,
                        "P3": p3_data,
                        **result,
                    }
                )
    certified = [record for record in records if record["independent"]]
    payload = {
        "schema": "elliptic-rank30/e6-mw3-modular-independence/v1",
        "field": f"GF({prime})(t)",
        "source_certificate": str(input_path),
        "source_certificate_sha256": source["certificate_sha256"],
        "ell_candidates": list(ell_candidates),
        "triples_tested": len(records),
        "triples_certified_independent": len(certified),
        "records": records,
        "proof": (
            "For each certified record, the diagonal specialization map has "
            "zero kernel on the three displayed classes modulo ell. A separate "
            "good fiber has order prime to ell, excluding generic ell-torsion. "
            "Infinite ell-divisibility therefore rules out every nonzero "
            "integral relation among P1, P2, P3."
        ),
        "claim_boundary": (
            f"This proves the displayed modular triples independent over "
            f"GF({prime})(t). It does not prove a characteristic-zero lift."
        ),
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--ell", default="2,3,5,7,11,13")
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    ell_candidates = [int(value) for value in args.ell.split(",")]
    payload = compute(args.input, ell_candidates)
    output = args.output or (
        ROOT
        / "certificates"
        / f"e6_mw3_modular_independence_gf{json.loads(args.input.read_text())['field_prime']}.json"
    )
    if not args.no_write:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "triples_tested": payload["triples_tested"],
                "triples_certified_independent": payload[
                    "triples_certified_independent"
                ],
                "certified": [
                    {
                        "joint_index": record["joint_index"],
                        "P2_index": record["P2_index"],
                        "P3_index": record["P3_index"],
                        "ell": record.get("ell"),
                        "specializations": len(
                            record.get("quotient_specializations", [])
                        ),
                    }
                    for record in payload["records"]
                    if record["independent"]
                ],
                "certificate_sha256": payload["certificate_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
