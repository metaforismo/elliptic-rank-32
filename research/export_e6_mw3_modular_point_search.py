#!/usr/bin/env python3
"""Export one modular obstruction ideal to the compact C search format."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--joint", type=int, required=True)
    parser.add_argument("--p2", type=int, required=True)
    parser.add_argument("--p3", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    records = [
        record
        for record in data["records"]
        if (
            record["joint_index"],
            record["P2_index"],
            record["P3_index"],
        )
        == (args.joint, args.p2, args.p3)
    ]
    if len(records) != 1:
        raise ValueError("target does not identify exactly one record")
    record = records[0]
    rows = record["obstruction_coefficient_rows"]
    dimension = record["tangent_parameter_dimension"]
    monomials = 1 + 2 * dimension + dimension * (dimension - 1) // 2
    if any(len(row) != monomials for row in rows):
        raise AssertionError("unexpected obstruction row length")
    lines = [f"{data['field_prime']} {dimension} {len(rows)} {monomials}"]
    lines.extend(" ".join(map(str, row)) for row in rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n")
    print(
        json.dumps(
            {
                "prime": data["field_prime"],
                "dimension": dimension,
                "rows": len(rows),
                "monomials": monomials,
                "output": str(args.output),
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
