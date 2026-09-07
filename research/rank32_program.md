# Rank 32 programme after ICARM curve #302

## 2026-09-07 checkpoint

The [new source audit](frontier_audit_20260907.md) confirms that the 631-entry
ICARM snapshot still has maximum lower bound 31, only at #302. We independently
certified the two newly posted rank-30 entries #398 and #582 as well; all four
rank-at-least-30 entries have distinct exact j-invariants. No new construction
family for #302 has been reproduced.

The previously pending endpoint-balanced lattice run is now imported: 7,868
successful moves, no cross-side theta fingerprint or bridge, four `p=2`
boundary rejections excluded from its negative scope. The new offset-32 run
also completed and was imported: 7,872 successful moves, no construction
errors, and 32 negative exact isometry checks at shared theta profile
`(40,2570)`. There is still no bridge. The current small lexicographic windows
are coordinate-biased; reproducible sampling across full projective
coordinates is the next unimplemented experiment. See the audit for exact
run IDs, hashes, and the separate hosted-CI dependency failure.

## New baseline

On 2026-08-23, [ICARM curve #302](https://elliptic-rank.icarm.cloud/curve/302)
was submitted by Ava Howell and credited to Claude, Levent Alpöge, and Ava
Howell.  The leaderboard records 31 independent rational points on

\[
y^2+xy+y=x^3+x^2+Ax+B,
\]

where

```text
A = -1284727764113567728281797636015784768866707681415849262157224232063
B = 560368321454261339256859338901915312332769858684945406858043869199456710681989058863306170127006181
```

The lower bound `rank E(Q) >= 31` is unconditional: it comes from explicit
independent points.  The leaderboard commentary's exact-rank-31 statement is
conditional on GRH+BSD and must remain separate.  This repository did not
discover the curve; `baseline/verify_rank31_icarm302.py` independently replays
the lower-bound certificate.

The next record problem is therefore

\[
\boxed{\operatorname{rank}E(\mathbf Q)\ge 32}.
\]

## What changes mathematically

The inherited elliptic-K3 accounting changes from `31=17+14` to
`32=17+15`.  Every finite-field incidence calculation already completed for
the rank-17 source remains valid, but a successful specialization must now
supply one additional independent exceptional direction.

The norm representation changes in parallel.  On a short model
`y^2=R(x)`, 32 independent points can be completed by the negative of their
sum to a 33-point packet cut out by

\[
A(x)^2-R(x)B(x)^2=c\prod_{i=1}^{33}(x-r_i),
\qquad \deg A\le16,\quad\deg B\le15.
\]

This is an exact representation, not a practical search-space reduction.
Only symmetric packets, overlapping low-degree packets, or genuinely
low-dimensional coefficient families should be searched.

## Ranked research routes

1. **Recover the construction provenance of #302.**  The public entry gives
   the curve and witnesses but no family parameter or search method.  A
   reproducible family would reveal which nearby specialization directions
   were almost rank 32 and is more valuable than blind coefficient search.
2. **Use #302 as a certificate baseline, not as the main point-search target.**
   Searching for a 32nd point on the same curve conflicts with the reported
   GRH+BSD exact-rank-31 calculation.  It remains useful for auditing descent,
   saturation, and discovery heuristics.
3. **Complete the rank-17 K3 reconstruction, then demand `17+15`.**  The exact
   determinant-948 lattice work, public Kneser-chain replay, and bounded
   target-to-endpoint search are construction prerequisites only.  A usable
   result still requires a marked Neron--Severi transport, Riemann--Roch
   pencil, section witness, and forward/inverse birational maps.
4. **Compute the `P3` Noether--Lefschetz divisor.**  Saturate the six exact
   remainder equations over the `P1/P2` incidence surface and lift only a
   finite-field point already carrying the complete same-surface section
   packet.
5. **Search structured degree-33 norm packets and new finite-field charts.**
   Group-law independence must be tested before any expensive lift; analytic
   scores and small residuals are triage signals, never certificates.

## Promotion gate

Nothing is a rank-32 result until the repository contains all of:

1. a nonsingular Weierstrass equation over `Q`;
2. 32 complete rational points satisfying it exactly;
3. an exact independence certificate, preferably a full-rank image in a
   product of finite quotients `E(F_p)/ell E(F_p)`;
4. a torsion argument and saturation evidence appropriate to the claim;
5. source provenance, deterministic inputs, software versions, and hashes;
6. a clean separation between the unconditional lower bound and every
   GRH/BSD-dependent upper bound.

The current status is therefore: **rank 31 independently reproducible; rank
32 open**.
