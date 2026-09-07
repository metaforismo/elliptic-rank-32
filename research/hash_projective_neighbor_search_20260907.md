# Full-coordinate projective neighbor sampling, 2026-09-07

## Hypothesis and controlled change

The two earlier 32-line lexicographic windows only exercised a small subset
of coordinates in each parent's current Gram basis. LLL changes that basis,
so this was not a proof of confinement to a fixed original subspace. It was
nevertheless a concrete reason to test a different selection rule.

Keep the three exact initial lattices, genus checks, neighbor construction,
rational basis maps, endpoint-balanced beam, and isometry promotion gates.
Use the same six rounds, beam size eight, primes 2 and 5, 32 lines per
parent/prime, 10,000 successful-move ceiling, and 500 exact-isometry ceiling.
Change only the selection of projective isotropic lines: use the fixed seed
`rank32-20260907-v1` and all 17 coordinates. There is no selected new curve,
rational P3, or elliptic-rank improvement at implementation time.

## Exact sampling and charged work

`research/sample_projective_lines.py` uses a SHA-256 domain containing the
algorithm version, exact Gram hash, prime, and seed. Each one-based raw-draw
index drives a SHAKE-256 byte stream. Byte rejection avoids modulo bias when
forming field coordinates. Zero vectors are rejected; every remaining line
is normalized by setting its last nonzero coordinate to one. The sampler
evaluates `q(v)=v^t G v/2` over the integers, checks `q(v) mod p=0`, and rejects
already accepted normalized isotropic lines. The offset counts accepted
unique lines, not raw draws. A window ends after the requested number of
post-offset lines or fails closed at 4,096 raw draws. A separate fixed byte
budget also bounds each raw draw.

The result is a finite hash-based sample, not an exhaustive enumeration or
a statistical proof of uniform graph coverage. Raw draws, field-byte
rejections, zero/non-isotropic vectors, duplicates, accepted lines, and
prepared-but-unattempted terminal lines are distinct costs.

The Sage 10.9
[neighbor implementation](https://github.com/sagemath/sage/blob/10.9/src/sage/quadratic_forms/quadratic_form__neighbors.py)
accepts a vector with `q(v)` divisible by `p` and can return its exact basis
matrix. We retain the existing exact determinant, intersection-index,
Gram-identity, evenness, positivity, and genus checks. The already documented
non-maximal/even `p=2` boundary remains incomplete; any other error aborts
without a mathematical negative conclusion.

## Evidence contract

Producer schema 2 explicitly distinguishes `hash-projective` and historical
`lexicographic` selection. The new `sampling-windows.jsonl` is included in
both the closed artifact manifest and the result's evidence hashes. Each
attempt links to its window hash. The sampler source is hash-bound alongside
the search and lattice-formula sources.

The standard-library importer replays the full draw stream from each
recorded parent Gram, checks window contents and counters, verifies the
configured round/side/parent/prime order, and compares every attempt with its
expected sampled vector and ordinal. Only a final bridge or move-budget
stop may leave a window partially attempted. Missing windows, relabeled
seeds, rehashed tampering, and partial schema downgrades are rejected.
Legacy schema-1 artifacts retain their own evidence and claim boundaries.

This replay checks sample selection, not all Sage neighbor computations.
As before, the importer independently recomputes every rational matrix on
any promoted bridge; it does not independently recompute every unpromoted
search move. A lattice bridge alone would still not provide a marked K3
fibration switch, a moduli map, rational P3, or 32 independent rational points.

## Initial local controls

The initial 32-line windows have the following charged draw counts. Each
window contains a nonzero entry in every one of its 17 coordinates across
the selected vectors; this is an observed check, not a coverage theorem.

| Initial lattice | Prime | Selected unique lines | Raw draws | Non-isotropic draws |
| --- | ---: | ---: | ---: | ---: |
| Origin E8+A2^3 | 2 | 32 | 79 | 47 |
| Origin E8+A2^3 | 5 | 32 | 155 | 123 |
| Transparent seed | 2 | 32 | 60 | 28 |
| Transparent seed | 5 | 32 | 169 | 137 |
| Frozen rootless | 2 | 32 | 75 | 43 |
| Frozen rootless | 5 | 32 | 192 | 160 |

These six windows have no zero vectors or duplicate accepted lines. Their
field-byte rejection counts are respectively 0, 13, 0, 11, 0, and 13.
Toy F2 and F5 isotropic projective spaces match independent exhaustive
enumeration. Tests cover reproducibility, offset slicing, normalization,
budget failure, invalid inputs, and importer tampering.

## Runtime validation

The local full suite completed 181 tests in 102.020 seconds: 180 passed and
one Sage-only API test was explicitly skipped. Both edited workflows passed
`actionlint`; `git diff --check` was clean. The new importer tests include
12 separately rehashed tampering cases and a partial-final-window check.

The pinned-Sage workflow first exercises real neighbor construction on two
sampled lines per initial lattice and prime, then runs the bounded search.
This Sage-specific test is explicitly skipped on local CPython; no Sage
success is implied by the local sampler tests. Full search evidence will be
imported only after successful runtime and archive/source verification.
