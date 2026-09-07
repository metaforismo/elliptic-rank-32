"""Reproducible, bounded full-coordinate sampling of isotropic projective lines.

This is a hash-based sample, not an exhaustive enumeration or a proof of
uniform coverage. Every draw fills every coordinate; byte rejection avoids
modulo bias, then zero/non-isotropic/repeated projective lines are rejected.
The last nonzero coordinate is normalized to one. No Sage dependency is used,
so the complete draw stream and its costs can be replayed by the importer.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from typing import Any


ALGORITHM = "sha256-shake256-full-coordinate-projective-v1"
FIELD_BYTES_PER_DRAW = 4096
MAX_DRAWS = 1_000_000


class SamplingBudgetError(RuntimeError):
    """No search conclusion follows when a sampling budget is exhausted."""


def payload_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def validate_parameters(prime: int, seed: str, offset: int, limit: int, max_draws: int) -> None:
    if type(prime) is not int or not 2 <= prime <= 251 or any(
        prime % divisor == 0 for divisor in range(2, math.isqrt(prime) + 1)
    ):
        raise ValueError("hash sampling requires a prime between 2 and 251")
    if not isinstance(seed, str) or re.fullmatch(r"[A-Za-z0-9_.:-]{1,128}", seed) is None:
        raise ValueError("sampling seed must contain 1 to 128 portable ASCII characters")
    if type(offset) is not int or offset < 0 or type(limit) is not int or limit < 1:
        raise ValueError("sampling offset must be nonnegative and limit positive")
    if type(max_draws) is not int or not 1 <= max_draws <= MAX_DRAWS:
        raise ValueError("sampling draw budget must be between 1 and 1000000")
    if offset + limit > max_draws:
        raise ValueError("requested unique lines exceed the raw draw budget")


def validate_gram(gram: list[list[int]]) -> None:
    if not isinstance(gram, list) or not 1 <= len(gram) <= 128:
        raise ValueError("sampling Gram dimension must be between 1 and 128")
    size = len(gram)
    if any(
        not isinstance(row, list) or len(row) != size or any(type(x) is not int for x in row)
        for row in gram
    ):
        raise ValueError("sampling Gram must be a square integer matrix")
    if any(gram[i][i] % 2 for i in range(size)) or any(
        gram[i][j] != gram[j][i] for i in range(size) for j in range(i)
    ):
        raise ValueError("sampling Gram must be symmetric with even diagonal")


def q_value(gram: list[list[int]], vector: list[int]) -> int:
    """Evaluate v^t G v / 2 without division modulo two."""
    return sum(
        gram[i][i] // 2 * vector[i] ** 2
        + sum(gram[i][j] * vector[i] * vector[j] for j in range(i))
        for i in range(len(gram))
    )


def normalize_projective(vector: list[int], prime: int) -> tuple[int, ...] | None:
    reduced = [entry % prime for entry in vector]
    pivot = next((entry for entry in reversed(reduced) if entry), None)
    if pivot is None:
        return None
    inverse = pow(pivot, -1, prime)
    return tuple(entry * inverse % prime for entry in reduced)


def draw_vector(domain: bytes, draw_index: int, prime: int, dimension: int) -> tuple[list[int], int]:
    stream = hashlib.shake_256(domain + draw_index.to_bytes(8, "big")).digest(FIELD_BYTES_PER_DRAW)
    cutoff = 256 - 256 % prime
    values: list[int] = []
    rejected = 0
    for byte in stream:
        if byte >= cutoff:
            rejected += 1
            continue
        values.append(byte % prime)
        if len(values) == dimension:
            return values, rejected
    raise SamplingBudgetError("per-draw field-byte budget exhausted; no mathematical conclusion")


def sample_window(
    gram: list[list[int]], prime: int, seed: str, offset: int, limit: int, max_draws: int
) -> dict[str, Any]:
    validate_parameters(prime, seed, offset, limit, max_draws)
    validate_gram(gram)
    domain_payload = {
        "algorithm": ALGORITHM,
        "gram_sha256": payload_sha256(gram),
        "prime": prime,
        "seed": seed,
    }
    domain = bytes.fromhex(payload_sha256(domain_payload))
    counters = {
        "raw_draws": 0,
        "rejected_field_bytes": 0,
        "zero_vectors": 0,
        "nonisotropic_vectors": 0,
        "duplicate_isotropic_lines": 0,
        "accepted_unique_isotropic_lines": 0,
    }
    seen: set[tuple[int, ...]] = set()
    selected: list[dict[str, Any]] = []
    for draw_index in range(1, max_draws + 1):
        raw, rejected = draw_vector(domain, draw_index, prime, len(gram))
        counters["raw_draws"] += 1
        counters["rejected_field_bytes"] += rejected
        vector = normalize_projective(raw, prime)
        if vector is None:
            counters["zero_vectors"] += 1
            continue
        if q_value(gram, list(vector)) % prime:
            counters["nonisotropic_vectors"] += 1
            continue
        if vector in seen:
            counters["duplicate_isotropic_lines"] += 1
            continue
        seen.add(vector)
        counters["accepted_unique_isotropic_lines"] += 1
        ordinal = len(seen)
        if ordinal > offset:
            selected.append({
                "ordinal_one_based": ordinal,
                "raw_draw_index_one_based": draw_index,
                "vector": list(vector),
            })
        if len(selected) == limit:
            record = {
                **domain_payload,
                "dimension": len(gram),
                "offset": offset,
                "limit": limit,
                "max_draws": max_draws,
                "field_bytes_per_draw_budget": FIELD_BYTES_PER_DRAW,
                "counters": counters,
                "selected_lines": selected,
            }
            return {**record, "window_sha256": payload_sha256(record)}
    raise SamplingBudgetError(
        f"raw sampling budget exhausted: {counters}; requested offset={offset}, limit={limit}; "
        "no mathematical conclusion"
    )
