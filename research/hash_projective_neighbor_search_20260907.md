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

### Controlled run: comparison budget exhausted

[Run 34164060337](https://github.com/metaforismo/elliptic-rank-31/actions/runs/34164060337)
at producer `9136eef09b53da17445171548fa8cf3882141802` passed all 11 sampler
tests under Sage, including the real neighbor API. The research run then
stopped in its third round at the configured 500-qfisom ceiling:

- 2,091 attempted and successfully constructed moves, with no construction
  failures or repeated serialized child Grams;
- 2,094 stored Gram presentations;
- 500 exact negative comparisons and one required comparison left unchecked;
- status `inconclusive_isometry_budget_exhausted`, explicitly non-importable.

This is **not** an imported finite-negative search certificate. The existing
importer correctly rejects this status. The download's archive hash matches
GitHub's digest:
`592ccd238653c23f6ccb59c85d65b992884e533911dec1e1104ad2db71f030be`
(artifact `10033632129`). The closed manifest and source hashes passed;
all 66 prepared sampling windows were separately replayed exactly. They
account for 7,511 raw draws, 5,399 non-isotropic rejections, 344 rejected
field bytes, and 2,112 prepared unique lines, of which the final 21 were not
attempted after the stop. No zero vectors or duplicate isotropic lines
occurred in those windows.

The origin side has 776 theta profiles among 940 states, versus the earlier
small lexicographic windows' much narrower profile sets. This is an observed
change in exploration, not evidence of elliptic-rank progress or an
isometry. A strict follow-up uses the same producer, seed, six rounds,
beam and move limits, but raises only the exact-comparison ceiling to
10,000: [run 34164507132](https://github.com/metaforismo/elliptic-rank-31/actions/runs/34164507132).
It found and independently verified an exact bridge, as detailed below.

Hosted standard-Python validation at the same producer also passed:
[run 34164091625](https://github.com/metaforismo/elliptic-rank-31/actions/runs/34164091625)
completed 181 tests in 74.968 seconds, with only the Sage-specific test
skipped there; that test passed in the separate Sage job above.

### Structural follow-up

The [exact root-reflection orbit audit](initial_f2_root_reflection_orbits.md)
reduces the initial 65,279 isotropic F2 lines to 95 and 34 subgroup orbits
on the two rootful endpoints. This exposes remaining redundancy and rare
orbits in the hash sample. The orbit representatives are certified, but
their complete even-neighbor coverage and any geometric transport remain
unproved.

### Positive follow-up: seven-edge rootless bridge

Run `34164507132` completed successfully in its fourth round, with 4,595
successful moves, 4,598 Gram presentations, no construction errors, and
5,439 exact isometry comparisons (5,438 negative, then one positive).
The origin path has primes `5,2,5`; the rootless-to-meeting path has primes
`2,5,2,5`. Thus traversal from origin to the frozen rootless endpoint uses
`5,2,5,5,2,5,2`. The meeting lattices have theta profile `(14,2508)`;
their actual integral isometry, not this fingerprint, establishes the meeting.

Artifact `10033864952` has verified ZIP SHA-256
`eb22ad6fd4f5584375122b32686be66ef4b112e25e874e9f3c93bb08dd584e4f`.
The source-locked importer recomputed all seven path moves, all intersection
indices, the meeting isometry, and the full rational basis composite.
The resulting compact certificate is
`certificates/e8_a2_target_neighbor_bridge_run_34164507132.json`, record hash
`4db1da644af9a647f1264463ab8f809cdac826563a78e8aaf47931bd5beadd13`.
The original `exact-bridge.json` is retained byte-for-byte in
`certificates/e8_a2_rootless_bridge_34164507132/` for local replay.

All 144 prepared sampling windows were replayed: 15,976 raw draws,
11,368 non-isotropic rejections, 721 field-byte rejections, no zero or
duplicate vectors, 4,608 prepared unique lines and 4,595 attempted lines.
The final 13 lines were not attempted after the bridge. The first 2,091
moves and attempts, 2,094 states, 66 windows, and first 500 completed
isometry comparisons match the earlier stopped run exactly. This verifies
the claimed continuation with only the comparison limit changed.

The [separate integral hyperbolic-extension certificate](stable_hyperbolic_neighbor_transport.md)
goes one algebraic step further. Neither certificate proves a geometric
K3 transport or an elliptic curve with 32 independent rational points.
