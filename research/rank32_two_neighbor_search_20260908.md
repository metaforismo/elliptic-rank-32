# Rank-32 continuation: degree-two routes and smaller abstract markings

Date: 2026-09-08, Europe/Rome. GitHub records run timestamps in UTC.

## Scope and repository identity

The repository is now `metaforismo/elliptic-rank-32`; its active development
branch is `codex/rank32-research`. The name denotes the open goal. The
certified public baseline remains rank at least 31, and this project has
not produced a rank-32 curve. Historical source hashes, run identities and
rank-31 certificates must not be rewritten to match the new name.

## Why seek a degree-two route?

The certified seven-edge bridge uses primes `5,2,5,5,2,5,2`. Its integral
U-stabilization is exact but does not provide nef representatives or
rational functions on a K3 surface.

[Elkies--Kumar, Section 5](https://arxiv.org/pdf/1209.3527#page=15)
construct geometric neighbor steps from divisor classes and Riemann--Roch
pencils. [Kumar, Appendix A](https://arxiv.org/pdf/1105.1715#page=44)
gives the degree-two cases and explicit changes of variables. These sources
motivate, but do not guarantee, the usefulness of a lattice route using
only prime two. Such a route would still require an actual marked source
surface, effectivity/nefness, and exact birational maps.

## Frozen experiment

- unchanged exact neighbor and isometry producer;
- bidirectional search from the same three frozen initial lattices;
- prime list `[2]`, eight rounds, beam size eight;
- full-coordinate hash-projective sampling, seed `rank32-20260908-p2-v1`;
- 32 accepted lines per parent, offset zero, at most 4,096 raw draws/window;
- at most 10,000 successful moves and 20,000 exact isometry comparisons.

Without an earlier stop, this schedule attempts at most
`3*32 + 7*3*8*32 = 5472` lines. The larger comparison budget is deliberate:
the earlier full-coordinate control stopped at 500 comparisons, while its
successful continuation first found a bridge at comparison 5,439.
Primes, seed, rounds and budget all differ from that earlier experiment;
this is not advertised as a one-variable controlled comparison.

The pinned Sage API check, raw-draw replay, exact path replay and bad-prime
boundaries are unchanged. An exhausted budget is inconclusive. A completed
finite-negative search excludes only its successfully constructed moves,
not the entire p=2 graph or the existence of a geometric path.

## Observed all-two bridge

[Actions run 34174072896](https://github.com/metaforismo/elliptic-rank-32/actions/runs/34174072896)
ran the frozen configuration on source
`f9ab37d1ecf9d570167f9006a798060251cc4437`. It succeeded after 4,986
attempted/successful moves, zero failed moves, and 6,509 exact isometry
comparisons. The positive path has seven forward origin edges and eight
inverted endpoint edges, all with prime two. Its endpoint is the
transparent `A11 + K6` lattice, **not** the rootless endpoint of the older
seven-edge bridge.

The archive SHA256 is
`f398e90268d9a667e2c0c60ac24c978568940feb9c43b93576bbb323352647f6`.
The independently imported compact record is
`certificates/e8_a2_target_neighbor_bridge_run_34174072896.json`; the raw
path is frozen under `certificates/e8_a2_transparent_p2_bridge_34174072896/`.
Replay covers all path moves, the meeting isometry, the rational
17-dimensional composite and 156 sampling windows: 10,175 raw draws,
5,183 nonisotropic vectors, 4,992 prepared lines, 4,986 attempts, and six
unused lines after the positive stop. This does not enumerate the p=2 graph.

`research/stable_degree_two_bridge.py` constructs integral 19-dimensional
extensions for all 15 edges and composes them, inverting the endpoint
chain. The final determinant is one; adjacent abstract fiber intersections
are two, but the end-to-end intersection is 1,798,927,074,466. The latter is
not a degree-two geometric map. Its certificate is
`certificates/stable_degree_two_bridge_34174072896.json`, record SHA256
`0a82f6899ae3345afb98086cb35e5ee5cd1271293baf563a0612a4aa860b1610`.

## Smaller integral markings on the original rootless bridge

These reductions concern the **older seven-edge rootless bridge**, not the
15-edge transparent bridge. The original certificate remains frozen.
With column-basis convention `M^T S_origin M = S_endpoint`,
`S = U + (-G)` and `U = [[0,1],[1,0]]`, the old/new fiber intersection is
`M[1][0]` (zero-based indexing).

| Marking | Maximum absolute matrix entry | Fiber intersection |
| --- | ---: | ---: |
| Original integral stabilization | 15,171,339,103,561,459,844 | 5,715,299,020,407,871 |
| Shorter lift/witness selection | 61,781,047 | 8,701,436 |
| Bounded Eichler descent | 10,437,320 | 6,050,128 |
| Bounded root reflections and cleanup | 1,320,133 | 710,298 |

All four markings are exact integral isometries with determinant one.
The three new certificates are, in table order,
`short_hyperbolic_rootless_marking_20260908.json`,
`reduced_hyperbolic_fiber_degree_20260908.json`, and
`root_reflection_hyperbolic_marking_20260908.json`, under `certificates/`.
The reproducing scripts are `research/reduce_stable_hyperbolic_marking.py`,
`research/reduce_hyperbolic_fiber_degree.py`, and
`research/reduce_hyperbolic_root_reflections.py`, respectively.

### Short lifts and pairing-one witnesses

Let `q(y) = y^T G y / 2`. Enumerate centered nonzero scalar multiples of
the supplied projective line modulo p. For each coordinate with
`(G y)_j != 0 mod p`, try both recorded representatives of the correction
`y -> y + p*a*e_j` that makes `q(y) = 0 mod p^2`. Test the first four
distinct lifts ordered by `(q, maxabs, vector)`.

The isotropic fiber is `F = (q(y)/p, p, y)`. Pairing-one witnesses include
the original Bezout witness and sparse vectors `(a,0,c*e_j)` solving
`p*a - (G y)_j*c = 1`. Replace a witness Z1 by
`Z = Z1 - (Z1^2/2)*F`, so `Z^2 = F^2 = 0` and `F.Z = 1`.
The lifted frame formula from the original proof then yields the full
integral matrix. Explicit child-basis membership and Gram/determinant
checks reject inadmissible lifts. Selection is finite and deterministic,
not a shortest-vector or global marking oracle.

### Exact Eichler descent

For an integral vector a, use the isometry

```text
        [ 1   q(a)   a^T G ]
E(a) =  [ 0     1       0 ] .
        [ 0     a       I ]
```

Direct multiplication gives `E(a)^T S E(a) = S` and determinant one.
For an isotropic `F=(x,b,z)` with `b>0`, translating `z -> z+b*a`
changes x to `q(z+b*a)/b`. A hyperbolic-coordinate swap strictly decreases
b when `0 < q(z+b*a)/b < b`. On the right, the same calculation uses
`G_endpoint^-1` times the tail of row one of M. These coordinates need
not themselves be integral; the chosen shift and matrix action are.

Exact Gram--Schmidt, Babai rounding and a bounded nearest-plane search
seek a shift with squared norm below two. Each search is capped at 100,000
nodes and stops on its first sufficient witness; at most 2,048 descent
iterations are allowed. The recorded rootless run used 1,904 nodes across
ten actions without exhausting a node budget. Search termination is not a
global optimality theorem. Zero-coordinate cusp boundaries fail closed.

### Bounded square-minus-two reflections

For positive k and an integral u with `q(u) = 1 mod k`, put
`rho=((q(u)-1)/k,k,u)`. Then `rho^2=-2`, and
`R_rho = I + rho*(rho^T S)` is an integral determinant-minus-one involution.
For isotropic `F=(x,b,z)` with `b>0`, direct expansion gives

```text
F.rho = b/(2*k) * ( (u-k*z/b)^T G (u-k*z/b) - 2 ).
```

Thus an exact sphere search with strict radius-squared two finds negative
pairings. Accept only `0 < b + k*(F.rho) < b`. The right action applies
the same test to `M^-1 e`, using the endpoint Gram matrix. Every accepted
action is independently replayed; no floating-point root selection is used.

The recorded run allowed root heights 1 through 64, at most 128 reflections
and 500,000 search nodes. It accepted 22 reflections after 125,793 nodes;
no node budget was exhausted. A final Eichler cleanup supplies the marking
in the table. The finite-height stopping statement belongs to the recorded
pre-cleanup marking and searched heights, not all roots or a global minimum.

## Verification and open geometric boundary

The full local suite loaded 214 tests: 213 passed, zero errors, and one
Sage-only test skipped (606.623 seconds). That unchanged API test passed
in the separately pinned Sage run above: 11/11 tests, no skips (1.328 seconds).

The 17 new reduction/transport tests recompute all four new certificates.
SymPy independently verifies all 22 new edge matrices, the four final
Gram identities and integral inverses, all 22 reflection involutions, and
the Eichler trace. Toy exhaustive CVP/root searches test strict inequalities,
while input/hash/trace tampering and budget/cusp boundaries are rejected or
reported as inconclusive. Reproduce with:

```bash
python3 -m unittest discover -s tests -p 'test_hyperbolic_marking_reductions.py' -v
python3 research/stable_degree_two_bridge.py \
  --compare certificates/stable_degree_two_bridge_34174072896.json
python3 research/reduce_hyperbolic_root_reflections.py \
  --compare certificates/root_reflection_hyperbolic_marking_20260908.json
```

The next geometric input remains an explicit marked source K3 surface over
Q. For a usable neighbor step, we must identify actual divisor classes,
prove nefness (or give the necessary effective corrections), compute its
Riemann--Roch pencil, and verify forward/inverse birational maps and section
transport. Small matrix coefficients and adjacent intersection two do not
provide those objects. No rational P3 or rank-32 curve has been found here.
