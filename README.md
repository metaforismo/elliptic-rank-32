# Elliptic-rank frontier: rank 31 certified, rank 32 open

Canonical repository: [`metaforismo/elliptic-rank-31`](https://github.com/metaforismo/elliptic-rank-31).
The repository name records the milestone it was created to pursue; the next
open record target is now rank 32.

This repository now studies the open constructive target

\[
\operatorname{rank} E(\mathbb Q) \ge 32
\]

by optimizing for a final proof certificate, not for a large analytic-rank score.

## Current truth status

**rank-31 baseline certified; rank 32 unsolved**

ICARM curve #302, submitted by Ava Howell on 2026-08-23 and credited to
Claude, Levent Alpöge, and Ava Howell, supplies 31 explicitly certified
independent rational points.  This repository independently replays that
unconditional lower bound with exact rational and finite-group arithmetic.
The separate statement that its rank is exactly 31 uses GRH+BSD.  This project
did not discover #302 and has not found a curve with 32 certified independent
rational points.

## Rank-31 record baseline

For

```text
y^2 + x*y + y = x^3 + x^2 + A*x + B
A = -1284727764113567728281797636015784768866707681415849262157224232063
B = 560368321454261339256859338901915312332769858684945406858043869199456710681989058863306170127006181
```

the files under `baseline/` verify all 31 public points exactly, prove trivial
rational torsion, and give a full-rank binary image matrix in a
product of finite quotients \(E(\mathbf F_p)/2E(\mathbf F_p)\). Therefore,
unconditionally,

\[
\operatorname{rank} E(\mathbb Q)\ge31.
\]

The reported exact-rank-31 statement uses GRH+BSD and is kept separate from
this certificate.

## Historical rank-30 baseline

ICARM curve #273, submitted publicly by `ranksunbounded` on 2026-08-20, was
the first public curve with 30 explicitly certified independent rational
points.  The repository retains its independent exact certificate for

```text
y^2 + x*y = x^3 + A*x + B
A = -201769035260418549083594900060734240952308696994802735114305555
B = 1151107939141058565733479426024323225135665982951300586808823640527729578307228357301072889377
```

The 30-point lower bound is unconditional; the separate exact-rank-30
statement uses GRH+BSD.

## Historical rank-29 baseline

For

```text
y^2 + x*y = x^3 - A*x + B
A = 27006183241630922218434652145297453784768054621836357954737385
B = 55258058551342376475736699591118191821521067032535079608372404779149413277716173425636721497
```

this repository verifies:

- all 29 complete rational points in `points.json` exactly;
- nonzero discriminant and a global minimal model;
- trivial rational torsion;
- an exact rank-29 local mod-2 independence certificate;
- canonical-height matrices in SageMath 10.9 and Magma V2.29-9;
- nonzero regulator and positive-definite numerical height matrix;
- exact conductor, local reduction data, and global root number `-1`;
- p-saturation of the listed subgroup for every prime below `4096`, independently in Sage and Magma.

Therefore, unconditionally,

\[
\operatorname{rank} E(\mathbb Q) \ge 29.
\]

The package does **not** promote the published conditional upper bound to an
unconditional exact-rank statement.

## Repository layout

```text
curve.json, points.json    historical exact rank-29 inputs
baseline/                  exact rank-29, rank-30, and rank-31 verifiers
certificates/              compact machine-readable proof summaries
research/                  construction theorems and bounded experiments
tests/                     exact regression and certificate tests
verify_exact.py            dependency-free historical verifier
```

## Reproduce the exact certificate

```bash
python3 verify_exact.py
python3 baseline/verify_rank30_icarm273.py
python3 baseline/verify_rank31_icarm302.py
python3 -m unittest discover -s tests -v
```

The optional exact research checks use SymPy when available. The lower-bound
certificates themselves use only the Python standard library.

## A long and useful failed route

Treating 13 extra points over a rank-17 fibration as 13 independent quadratic
splitting characters forces a multiquadratic auxiliary curve of genus at least
`20481`. Low genus is lost already at the fourth independent character. This is
a restricted obstruction to that construction mechanism, not a universal rank bound.

## What the failure reveals

One character need not carry only one point. The package verifies an explicit
family in which one nontrivial quadratic character carries three independent
Mordell-Weil directions. Product-twist channels can also contribute sections
that were not used to construct the original covers.

## The decisive change of setting

Search over complete **Galois-character height packets**. For every low-genus
branch code, compute all twist Mordell-Weil lattices, same-character rank
multiplicities, product twists, successive minima, and global-solubility data.
Only then sieve specializations.

## Current rank-32 frontier

The leading inherited route is to recover the hidden rank-17 K3 family and
search for a specialization with 15 exceptional directions, i.e.
\(32=17+15\).  The calculations below were begun for the earlier
\(31=17+14\) target and remain reusable construction evidence.  The exact
GF(31) E6/MW3 replay in
`research/replay_e6_mw3_p2_split_cores.py` exhausts one published seven-core
tranche: its only geometric P2 hit satisfies
\(P_1+P_2+P_3=O\), so that tranche contains no independent rank-3 modular
seed.

The independent coordinate-slice engine then re-derived the two core
equations and exhaustively scanned five GF(31) slices: 4,468,650 core
quadruples, 9,224 simultaneous core solutions, 12 exact surfaces with the
required two I2 fibers, 10,044 canonical P2 cases up to fiber-label swap, and
11,532 polynomial-P3 cases. One new geometric P2 hit survived, but the only
joint P2/P3 surface satisfies the exact relation

\[
P_1+2P_3=O.
\]

Thus these five nondegenerate slices contain no MW-rank-3 seed. This is a new
bounded negative result, not a global obstruction. See
`research/rank32_program.md`, the historical
`research/rank31_program.md`, and
`certificates/e6_mw3_coordinate_slices_gf31_summary.json` for exact claim
boundaries.

The exceptional planes where the generic triangular pivot vanishes have also
been exhausted separately on the public slice and the five new slices.  The
replacement chain (P_{1,5}\to a_3\), (P_{1,4}\to y_2) reduces them to
172,980 exact triples.  They contain 354 simultaneous core solutions but no
surface with both required split (I_2) fibers, so no section search or lift
survives from those planes.

The proposed all-`IV` `E8+A2^3` shortcut has also been audited and rejected:
its section ansatz omitted two high coefficients, and its constant `j=0`
complex multiplication forces even geometric Mordell--Weil rank.  A certified
replacement is now available.  The full marked semistable chart

\[
A=-3(a^2+2D\beta),\qquad
B=2(a^3+3aD\beta)+D^2\gamma,\qquad
a\gamma-3\beta^2=dD
\]

has dimension four and generically fibres `II*+3I3+5I1`.  The target
rank-three Mordell--Weil lattice has Gram matrix

\[
\frac13\begin{pmatrix}8&-1&0\\-1&10&0\\0&0&12\end{pmatrix},
\]

is primitive, and would give Neron--Severi discriminant `948` when the full
Mordell--Weil rank is three.  This is a recognition theorem, not an existence
claim.

The arithmetic search space is now sharper than the original all-split chart:
the target forces splitting only at two of the three `I3` fibres.  An exact
four-parameter two-split chart and its three intrinsic boundary charts cover
the marked target locus before the residual Kodaira and simple-contact open
conditions.  Exhaustive searches over GF(5), GF(7), GF(11), and GF(13) tested
`71,302,414` section candidates on the dense chart and those boundaries.

The intersection audit found one genuine missing stratum: two sections can
meet at infinity.  Using the exact condition
`I_affine+I_infinity=2` gives 30 corrected `P1/P2` target pairs, ten of them on
this infinity-collision stratum.  All `9,456,854` corresponding `P3` candidates
were exhausted, with zero `P3` and zero complete Gram triples.  Deep nodal
contacts and polynomial degree drops are now proved impossible on the declared
Kodaira open, so they do not hide another target component in these charts.

A deterministic pair-first extension sampled 20,000,000 first-stage
incidences at each of `p=17,19,23` (60,000,000 total).  It performed
`88,931,043` section tests, found three corrected total-intersection-two
`P1/P2` pairs, and exhaustively tested all `5,315,813` downstream `P3`
quartics on those pairs, again finding none.  The downstream searches are
exhaustive conditional on the sampled `P1` hits; the first stage is not
exhaustive at these larger primes.

All eight retained dense-chart representatives are smooth points of a
two-dimensional `P1/P2` incidence locus.  One has been lifted formally through
`11^8`, but only over `Z_11`: no rational lift, `P3`, characteristic-zero target
triple or characteristic-zero curve has been obtained from this K3 route; in
particular it has not produced a rank-32 candidate.

The complete three-section system is now explicit: 62 variables, 67 raw
equations, six exact local syzygies, and a reduced `61 x 62` Jacobian.  On the
infinity open, the `P3` equation triangularizes to six numerator equations of
degrees `21,19,17,15,14,12`.  Their saturated elimination is the concrete
Noether--Lefschetz divisor still to compute.

There is also an exact discriminant-form and Clifford identification with
Elkies' `X(6,79)` landscape.  The positive `E8+A2^3` target lattice has rank
17, determinant 948, exactly 258 roots, and the same cyclic finite quadratic
form as both the published period lattice and the transparent determinant-948
neighbor seed.  Its Clifford invariants recover quaternion discriminant 6 and
level 79.  This strongly identifies the correct K3 landscape, but an integral
target-to-endpoint bridge, explicit neighbor/moduli map, and transported
rational Weierstrass model are still missing.

Two Sage 10.9 Actions computations have now completed.  Run
[`32673843229`](https://github.com/metaforismo/elliptic-rank-31/actions/runs/32673843229)
exactly replayed the public seven-edge Kneser chain with primes
`5,5,2,2,2,2,2`, retained every rational basis transition, verified all
identities `T^t G_parent T=G_child`, and verified the composed basis identity.
It reproduces chain SHA-256
`f549651c06190e97f53947475ef9ee3149ddb175b2db2a289dd7db02d60b9d4b`
and terminal Gram hash
`620a5e06473684d3e8015c0172f63c09c901e742ec02e77ba0aa35a923aa0295`.

The audited bounded search
[`32674002260`](https://github.com/metaforismo/elliptic-rank-31/actions/runs/32674002260)
attempted 4,992 deterministic projective lines and serialized 4,990 exact
neighbor moves, 24 repeated child Grams, and 4,969 Gram presentations.  It
found no common cross-side theta fingerprint, required no `qfisom` comparison,
and found no exact target-to-endpoint bridge.  Two enumerated `p=2` vectors
were rejected at Sage's non-maximal/even boundary; because `2 | 948`, this is
not a complete enumeration of the 2-neighbor graph.  The negative scope is
therefore the 4,990 successfully constructed moves in the recorded beams, not
all 2-neighbors.  The historical engine was audited post hoc and is
superseded for future bounded-negative runs by fail-closed error and
isometry-budget handling.

The next controlled experiment repeats the offset-zero window with an
endpoint-balanced beam of size eight: three origin states directed toward the
transparent endpoint, three toward the rootless endpoint, and two diversity
states.  These computations are exact positive-lattice evidence; they do not
provide a marked K3 fibration switch, rational `P3`, or rank-32 curve.

A primary-source audit also sharpens the missing geometric layer.  A Kneser
`p`-neighbor between positive lattices is not automatically an elliptic
neighbor step on one K3 surface.  The latter additionally requires a marked
Neron--Severi lattice, effective old and new fiber divisors, a two-dimensional
Riemann--Roch pencil, a section witness, and mutually inverse birational maps.
The available Elkies--Kumar algorithms make this explicit for geometric
2- and 3-neighbors; the two leading prime-5 moves in the public lattice chain
remain routing heuristics until separately geometrized.  See
`research/kneser_to_weierstrass_transport_requirements.md`.

Exact counts, proofs, and claim boundaries are in
`research/rank32_program.md`, `research/rank31_program.md`,
`research/e8_a2_semistable_two_split_full_incidence.md`, and
`research/e8_a2_shimura_bridge.md`.  Compact run audits are in
`certificates/rank17_exact_neighbor_chain_run_32673843229.json` and
`certificates/e8_a2_target_neighbor_bridge_run_32674002260.json`.

## Certificate policy

A future rank-32 promotion requires a specified Weierstrass equation, 32 full
rational points, exact substitutions, torsion, a canonical-height matrix,
positive determinant, an independent exact independence certificate, saturation,
software versions, provenance, hashes, and verification in a second system.

See `STATUS.frontier.json` and `research/rank32_program.md`.
