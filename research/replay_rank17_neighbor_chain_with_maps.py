#!/usr/bin/env sage-python
"""Replay the frozen rank-17 chain and serialize exact basis transformations.

The historical replay retained each prime, neighbor vector, and LLL-reduced
child Gram matrix, but not the rational basis matrix returned by Sage.  This
companion replay uses ``return_matrix=True`` and records

    T = B * U,       T^t G_parent T = G_child,

where ``B`` is Sage's p-neighbor basis and ``U`` is the exact integral LLL
change of basis.  The matrices close the algebraic parent-to-child lattice
transitions.  They are not an explicit Neron-Severi marking, K3 neighbor map,
rational section, or rank-31 certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

from sage.all import Matrix, QQ, ZZ, identity_matrix

from replay_rank17_neighbor_chain import (
    DEFAULT_SOURCE,
    EXPECTED_TARGET_HASH,
    extract_discovery_script,
    patch_script,
)


EXPECTED_CHAIN_SHA256 = (
    "f549651c06190e97f53947475ef9ee3149ddb175b2db2a289dd7db02d60b9d4b"
)


def replace_once(script: str, old: str, new: str, label: str) -> str:
    if script.count(old) != 1:
        raise RuntimeError(f"expected one {label} patch site, found {script.count(old)}")
    return script.replace(old, new, 1)


def patch_basis_capture(script: str) -> str:
    old_state = """def lll_reduce_gram(G):
    U = G.LLL_gram()
    R = U.transpose() * G * U
    return Matrix(ZZ, R)

def form_from_gram(G):
    return QuadraticForm(ZZ, G)

def theta_data(Q):
    # For an even lattice with Hessian/Gram G, Sage's quadratic
    # form is q(x)=x^T G x/2.  Thus q^1 counts norm-2 roots and
    # q^2 counts norm-4 vectors.
    series = Q.theta_series(3)
    return int(series[1]), int(series[2])

def state_from_form(Q, genus, parent=None, move=None, source=None):
    G = lll_reduce_gram(Matrix(ZZ, Q.Hessian_matrix()))
    if G.det() != 948 or any(G[i, i] % 2 for i in range(17)):
        raise AssertionError(('invalid neighbor Gram', G.det(), G.diagonal()))
    if Genus(G) != genus:
        raise AssertionError(('neighbor left target genus', Genus(G), genus))
    RQ = form_from_gram(G)
    roots, norm4 = theta_data(RQ)
    return {
        'gram': G,
        'form': RQ,
        'hash': matrix_hash(G),
        'root_count': roots,
        'norm4_count': norm4,
        'parent': parent,
        'move': move,
        'source': source,
    }
"""
    new_state = """def lll_reduce_gram(G):
    U = G.LLL_gram()
    R = U.transpose() * G * U
    return Matrix(ZZ, R), Matrix(ZZ, U)

def serial_rational_matrix(M):
    return [
        [
            [int(M[i, j].numerator()), int(M[i, j].denominator())]
            for j in range(M.ncols())
        ]
        for i in range(M.nrows())
    ]

def form_from_gram(G):
    return QuadraticForm(ZZ, G)

def theta_data(Q):
    # For an even lattice with Hessian/Gram G, Sage's quadratic
    # form is q(x)=x^T G x/2.  Thus q^1 counts norm-2 roots and
    # q^2 counts norm-4 vectors.
    series = Q.theta_series(3)
    return int(series[1]), int(series[2])

def state_from_form(
    Q,
    genus,
    parent=None,
    move=None,
    source=None,
    parent_gram=None,
    basis_from_parent=None,
):
    G, U = lll_reduce_gram(Matrix(ZZ, Q.Hessian_matrix()))
    parent_transform = None
    if basis_from_parent is not None:
        if parent_gram is None:
            raise AssertionError('basis transform supplied without parent Gram')
        parent_transform = basis_from_parent * U
        if parent_transform.transpose() * parent_gram * parent_transform != G:
            raise AssertionError('B*U does not transport parent Gram to child Gram')
        if abs(parent_transform.det()) != 1:
            raise AssertionError(('unexpected neighbor basis determinant', parent_transform.det()))
    if G.det() != 948 or any(G[i, i] % 2 for i in range(17)):
        raise AssertionError(('invalid neighbor Gram', G.det(), G.diagonal()))
    if Genus(G) != genus:
        raise AssertionError(('neighbor left target genus', Genus(G), genus))
    RQ = form_from_gram(G)
    roots, norm4 = theta_data(RQ)
    return {
        'gram': G,
        'form': RQ,
        'hash': matrix_hash(G),
        'root_count': roots,
        'norm4_count': norm4,
        'parent': parent,
        'move': move,
        'source': source,
        'parent_transform': parent_transform,
    }
"""
    script = replace_once(script, old_state, new_state, "LLL/state")

    old_public_tail = """        'gram': serial_matrix(state['gram']),
    }

try:
"""
    new_public_tail = """        'gram': serial_matrix(state['gram']),
    }

def public_map_record(state):
    T = state['parent_transform']
    if state['parent'] is None or state['move'] is None or T is None:
        raise AssertionError('cannot serialize a seed as a neighbor move')
    determinant = T.det()
    return {
        'parent_hash': state['parent'],
        'child_hash': state['hash'],
        'prime': state['move']['prime'],
        'vector': state['move']['vector'],
        'basis_transform_child_in_parent': serial_rational_matrix(T),
        'basis_identity_verified': True,
        'determinant': [int(determinant.numerator()), int(determinant.denominator())],
        'sage_api': 'QuadraticForm.find_p_neighbor_from_vec(return_matrix=True)',
        'lll_composition': 'T = B * U',
    }

try:
"""
    script = replace_once(script, old_public_tail, new_public_tail, "map serializer")

    old_neighbor = """                    neighbor = parent['form'].find_p_neighbor_from_vec(p, vector)
                    generated += 1
                    state = state_from_form(
                        neighbor,
                        target_genus,
                        parent=parent['hash'],
                        move={'prime': int(p), 'vector': [int(v) for v in vector]},
                    )
"""
    new_neighbor = """                    neighbor_basis = parent['form'].find_p_neighbor_from_vec(
                        p, vector, return_matrix=True
                    )
                    neighbor = parent['form'](neighbor_basis)
                    generated += 1
                    state = state_from_form(
                        neighbor,
                        target_genus,
                        parent=parent['hash'],
                        move={'prime': int(p), 'vector': [int(v) for v in vector]},
                        parent_gram=parent['gram'],
                        basis_from_parent=neighbor_basis,
                    )
"""
    script = replace_once(script, old_neighbor, new_neighbor, "p-neighbor basis")

    old_chain = """    state_by_hash = {state['hash']: state for state in archive}
    target_neighbor_chains = []
    for target_state in targets:
        chain = []
        cursor = target_state
        while cursor is not None:
            chain.append(public_record(cursor, index_by_hash))
            cursor = state_by_hash.get(cursor['parent']) if cursor['parent'] is not None else None
        chain.reverse()
        target_neighbor_chains.append(chain)

    result.update({
"""
    new_chain = """    state_by_hash = {state['hash']: state for state in archive}
    target_neighbor_chains = []
    target_neighbor_maps = []
    for target_state in targets:
        state_chain = []
        cursor = target_state
        while cursor is not None:
            state_chain.append(cursor)
            cursor = state_by_hash.get(cursor['parent']) if cursor['parent'] is not None else None
        state_chain.reverse()
        target_neighbor_chains.append([
            public_record(state, index_by_hash) for state in state_chain
        ])
        target_neighbor_maps.append([
            public_map_record(state) for state in state_chain[1:]
        ])

    result.update({
"""
    script = replace_once(script, old_chain, new_chain, "target-chain traversal")
    old_field = """        'target_neighbor_chains': target_neighbor_chains,
        'target_found': bool(verified_targets),
"""
    new_field = """        'target_neighbor_chains': target_neighbor_chains,
        'target_neighbor_maps': target_neighbor_maps,
        'target_found': bool(verified_targets),
"""
    return replace_once(script, old_field, new_field, "target map result")


def rational_matrix(data: list) -> Matrix:
    if len(data) != 17 or any(len(row) != 17 for row in data):
        raise ValueError("basis transform is not 17 by 17")
    return Matrix(
        QQ,
        [
            [QQ(entry[0]) / QQ(entry[1]) for entry in row]
            for row in data
        ],
    )


def serial_rational_matrix(matrix: Matrix) -> list[list[list[int]]]:
    return [
        [
            [int(matrix[i, j].numerator()), int(matrix[i, j].denominator())]
            for j in range(matrix.ncols())
        ]
        for i in range(matrix.nrows())
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    base = patch_script(extract_discovery_script(args.source), output)
    script = patch_basis_capture(base)
    replay_script = output / "map-replayed-discovery.py"
    replay_script.write_text(script)
    completed = subprocess.run(
        [sys.executable, str(replay_script)],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    (output / "map-replay.log").write_text(completed.stdout)
    print(completed.stdout, end="")
    if completed.returncode:
        raise SystemExit(f"map replay failed with exit code {completed.returncode}")

    discovery = json.loads((output / "rank17-root-elimination.json").read_text())
    if discovery.get("status") != "target_found":
        raise SystemExit(f"map replay did not reach target: {discovery.get('status')}")
    chains = discovery.get("target_neighbor_chains") or []
    matches = [chain for chain in chains if chain and chain[-1]["hash"] == EXPECTED_TARGET_HASH]
    if len(matches) != 1:
        raise SystemExit(f"expected one target chain, found {len(matches)}")
    chain = matches[0]
    chain_raw = (json.dumps(chain, indent=2, sort_keys=True) + "\n").encode()
    chain_sha256 = hashlib.sha256(chain_raw).hexdigest()
    if chain_sha256 != EXPECTED_CHAIN_SHA256:
        raise SystemExit(f"map replay changed frozen chain: {chain_sha256}")

    map_sets = discovery.get("target_neighbor_maps") or []
    matching_maps = [
        maps
        for maps in map_sets
        if maps and maps[-1].get("child_hash") == EXPECTED_TARGET_HASH
    ]
    if len(matching_maps) != 1:
        raise SystemExit(f"expected one target map chain, found {len(matching_maps)}")
    transitions = matching_maps[0]
    if len(transitions) != 7:
        raise SystemExit(f"expected seven basis transformations, found {len(transitions)}")

    composed = identity_matrix(QQ, 17)
    for index, transition in enumerate(transitions):
        parent = chain[index]
        child = chain[index + 1]
        if transition["parent_hash"] != parent["hash"]:
            raise SystemExit(f"parent hash mismatch at transition {index}")
        if transition["child_hash"] != child["hash"]:
            raise SystemExit(f"child hash mismatch at transition {index}")
        if transition["prime"] != child["move"]["prime"]:
            raise SystemExit(f"prime mismatch at transition {index}")
        if transition["vector"] != child["move"]["vector"]:
            raise SystemExit(f"neighbor vector mismatch at transition {index}")
        transform = rational_matrix(transition["basis_transform_child_in_parent"])
        parent_gram = Matrix(ZZ, parent["gram"])
        child_gram = Matrix(ZZ, child["gram"])
        if transform.transpose() * parent_gram * transform != child_gram:
            raise SystemExit(f"Gram transport failed at transition {index}")
        if abs(transform.det()) != 1:
            raise SystemExit(f"determinant gate failed at transition {index}")
        composed *= transform

    seed_gram = Matrix(ZZ, chain[0]["gram"])
    target_gram = Matrix(ZZ, chain[-1]["gram"])
    if composed.transpose() * seed_gram * composed != target_gram:
        raise SystemExit("composed seed-to-target Gram transport failed")
    if abs(composed.det()) != 1:
        raise SystemExit("composed basis determinant is not a unit")

    payload = {
        "schema_version": 1,
        "status": "completed",
        "claim_status": (
            "exact rational parent-to-child basis transformations for the "
            "frozen determinant-948 lattice chain"
        ),
        "chain_sha256": chain_sha256,
        "source_workflow": str(args.source),
        "source_workflow_sha256": hashlib.sha256(args.source.read_bytes()).hexdigest(),
        "map_replayed_discovery_sha256": hashlib.sha256(
            replay_script.read_bytes()
        ).hexdigest(),
        "target_hash": EXPECTED_TARGET_HASH,
        "transition_count": len(transitions),
        "transitions": transitions,
        "target_basis_in_seed_coordinates": serial_rational_matrix(composed),
        "composed_basis_identity_verified": True,
        "composed_determinant": [
            int(composed.det().numerator()),
            int(composed.det().denominator()),
        ],
        "truth_note": (
            "Each T=B*U is checked by T^t G_parent T=G_child, with B returned "
            "by Sage's p-neighbor API and U the integral LLL basis change. "
            "This closes the algebraic lattice transitions only. It does not "
            "construct an integral chain from the E8+A2^3 target lattice, a "
            "geometric Neron-Severi marking, a K3 neighbor map, a rational P3 "
            "section, or a rank-31 curve over Q."
        ),
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["record_sha256"] = hashlib.sha256(raw).hexdigest()
    maps_path = output / "exact-neighbor-maps.json"
    maps_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    summary = {
        "status": "completed",
        "chain_sha256": chain_sha256,
        "map_record_sha256": payload["record_sha256"],
        "map_file_sha256": hashlib.sha256(maps_path.read_bytes()).hexdigest(),
        "map_replayed_discovery_sha256": hashlib.sha256(
            replay_script.read_bytes()
        ).hexdigest(),
        "transition_count": 7,
        "target_hash": EXPECTED_TARGET_HASH,
        "composed_basis_identity_verified": True,
        "claim_boundary": (
            "exact lattice basis maps only; no E8+A2^3 integral bridge, "
            "geometric K3 map, rational P3, or rank-31 curve"
        ),
    }
    (output / "map-summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
