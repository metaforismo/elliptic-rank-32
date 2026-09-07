# Public frontier and resumed search: 2026-09-07

## Public result

The verified public lower-bound record remains **31**, supplied by
[ICARM #302](https://elliptic-rank.icarm.cloud/curve/302) on 2026-08-23 and
credited to Claude, Levent Alpöge, and Ava Howell. We have not found a
verifiable rank-32 announcement in the sources checked on 2026-09-07.
The exact-rank-31 assertion remains conditional on GRH+BSD; the independent
31-point lower bound is unconditional.

The [complete ICARM database](https://elliptic-rank.icarm.cloud/database.json)
downloaded on 2026-09-07 contains 631 distinct curve IDs. Its maximum
`rank_lower_bound` is 31, attained only by #302. There are three entries with
lower bound 30 and none with lower bound at least 32. The full response is
1,834,264 bytes, SHA-256
`8c736cf4eb5cd72306f523d9f663c05bd154da623dc47ff91ac6bd679f686d3b`.
The committed input is a **projection**: the full source histogram, source
metadata, and all four records above the rank-30 threshold. It is not a
mirror of the complete database. The projection was compared to the actual
full download, including unique IDs, count, histogram, and selected fields.

[Dujella's rank history](https://web.math.pmf.unizg.hr/~duje/tors/rankhist.html)
now includes 31, and his [rank-31 page](https://web.math.pmf.unizg.hr/~duje/tors/rk31.html)
lists the same equation and 31 points. The introductory text of the history
page still mentions 30; the explicit record table and dedicated page are
the relevant updated portions.
[Epoch's solution update](https://epoch.ai/frontiermath/open-problems/elliptic-curve-rank)
now marks its rank-at-least-30 problem solved by AI and records the later
rank-31 curve. The old problem statement lower down is historical, not a
current frontier assessment.

This is evidence about these checked sources on this date, **not** a proof
that no unpublished or unindexed rank-32 example exists.

## New exact baseline checks

The frozen inputs are in `baseline/icarm_frontier_20260907.json`, with
verification by `baseline/verify_icarm_frontier_20260907.py` and output in
`certificates/icarm_frontier_20260907.json`.

| ICARM ID | Public submission | Certified lower bound | Exact binary image matrix | Torsion reduction orders |
| --- | --- | ---: | --- | --- |
| [273](https://elliptic-rank.icarm.cloud/curve/273) | 2026-08-20 | 30 | 30 x 30, rank 30 | #E(F19)=28; #E(F23)=33 |
| [302](https://elliptic-rank.icarm.cloud/curve/302) | 2026-08-23 | 31 | 31 x 34, rank 31 | #E(F17)=26; #E(F31)=43 |
| [398](https://elliptic-rank.icarm.cloud/curve/398) | 2026-08-28 | 30 | 30 x 32, rank 30 | #E(F11)=18; #E(F23)=31 |
| [582](https://elliptic-rank.icarm.cloud/curve/582) | 2026-09-04 | 30 | 30 x 32, rank 30 | #E(F19)=28; #E(F43)=53 |

All 121 point substitutions hold exactly on both generalized and short
models. Every listed subgroup is Z-linearly independent: full binary row
rank forces all coefficients of a putative relation to be even; trivial
rational torsion permits infinite descent. Each torsion argument checks
coprime finite group orders and the cross-characteristic exclusions. These
checks use only Python's standard library and the existing exact finite-group
arithmetic. They do not establish saturation or an upper rank bound.

The #302 equation, ordered points, and discriminant match our August
baseline exactly. Its binary matrix is unchanged. The four exact
j-invariants are pairwise distinct. Hence none of these four curves is
isomorphic to, or a twist of, another over Q. This does **not** exclude their
arising in a common parameter family.

The #398 commentary is by Elkies and Klagsbrun and dates their discovery
to September 2025; the database submission occurred on 2026-08-28. Keep
these dates separate. The linked 2026-08-29 NMBRTHRY announcement could not
be retrieved during this audit (web access failed; direct GET timed out).
No claim about its construction method or conditional upper bound is made
from secondary summaries. The #582 entry has no construction commentary.
We still have no reproduced family parameter for #302.

Reproduce the arithmetic from the committed projection:

```sh
python3 baseline/verify_icarm_frontier_20260907.py
python3 baseline/verify_rank31_icarm302.py
python3 -m unittest discover -s tests -v
```

To additionally recheck the projection against a retained **identical** full
download, add `--source-database /path/to/database.json`. A later changed
download must fail the pinned source hash; it requires a separately dated
snapshot, not silently rewriting this one.

## Recovered bounded lattice experiment

[Actions run 32677113429](https://github.com/metaforismo/elliptic-rank-31/actions/runs/32677113429)
completed successfully on 2026-08-24, but its result was not imported at the
previous checkpoint. It is now verified against commit
`02d06d68168fed0b586c4cfa1337a9545f3e3092` and artifact 9503176418, whose
downloaded ZIP SHA-256 matches GitHub's digest:
`825b17a1ef89b91ed1b43a010fc240bd9bac96f28974f284e9d42eb43d239a82`.
The compact certificate is
`certificates/e8_a2_target_neighbor_bridge_run_32677113429.json`.

- Configuration: bidirectional, 6 rounds, beam 8, primes 2 and 5, line
  offset 0, 32 lines per parent/prime, at most 10,000 successful moves and
  500 isometry checks.
- Origin beam: three states directed toward each endpoint and two
  deterministic diversity states.
- 7,872 attempted lines; 7,868 successful moves; four classified `p=2`
  boundary rejections; 21 repeated child Grams; 7,850 retained presentations
  across all discovered states (not 7,850 isometry classes).
- No equal cross-side theta fingerprint; zero `qfisom` calls; no bridge.

The frozen state's first two theta coefficients give 29 distinct profiles
on the origin side, 31 on the transparent side, and 252 on the rootless
side. A close origin/transparent pair is `(50,2598)` versus `(50,2618)`.
Closeness in this heuristic has no theorem-level implication for a lattice
meeting, much less a rational curve.

The archived vectors also expose a sampling limitation: all 3,936 successful
`p=5` moves have coordinates 5--16 zero in their **parent's current basis**;
all 3,932 successful `p=2` moves have coordinates 8--16 zero. Every one has
coordinate 17 equal to one. Thus the small lexicographic windows are highly
coordinate-biased. Since LLL changes the basis after each move, this does
not prove confinement to a single fixed subspace of the original lattice.
If the offset-32 control is again negative, a reproducible sampler covering
the full projective coordinates is a better next experiment than merely
repeating many adjacent tiny windows. It must retain exact isotropy tests,
normalization/deduplication, attempted-vector counts, and all basis maps.

The importer checks the complete SHA manifest, source/runtime bindings,
record hashes, counters, rounds, error classification, and promotion gates.
It does not independently reconstruct every Sage neighbor. Its bounded
negative statement excludes only the successfully constructed moves in
these recorded finite beams. The four failed `p=2` constructions are outside
that scope; the full 2-neighbor graph was not enumerated.

## Resumption

The same producer commit was dispatched as
[run 34161893470](https://github.com/metaforismo/elliptic-rank-31/actions/runs/34161893470)
with `line_offset=32` and otherwise identical bounds. Its deterministic
line ordinals are 33--64 rather than 1--32 for each retained Gram and prime.
The line windows differ; discovered lattices may overlap or be isometric.
This run completed successfully at 21:18 UTC and its artifact was imported
against the same clean producer checkout. Artifact 10033092263 has ZIP
SHA-256 `300930964746109f45c7e0e5892bc53ea96eab49f25e1dbcf4f0f10e12d96806`,
matching the GitHub digest. The compact audit is
`certificates/e8_a2_target_neighbor_bridge_run_34161893470.json`.

It attempted and successfully constructed **7,872** moves, with zero
construction failures, 24 repeated child Grams, and 7,851 presentations.
It found 36 origin theta profiles, 40 transparent profiles, and 268 rootless
profiles. The only shared origin/endpoint profile is `(40,2570)` on the
rootless side. All 32 corresponding pairs were checked with exact PARI
`qfisom` and returned non-isometric; no comparison was skipped or exhausted
the budget. **No bridge was found in this bounded experiment.** These 32
negative decisions are supplied by the pinned Sage/PARI run, not independently
re-executed by the standard-library artifact importer.

The sampling limitation persists: the 3,936 successful `p=5` vectors still
have coordinates 5--16 zero in their parent bases. For `p=2`, coordinates
9--16 are zero. The next proposed experiment is therefore full-coordinate
deterministic projective sampling, not another adjacent small window. No
implementation or run of that new sampler is claimed here.

Even a positive lattice bridge would leave the marked K3 transport,
Riemann--Roch pencil, rational section, and rank-32 specialization unproved.
The constructive target remains `32=17+15`; neither these source checks nor
the lattice computations constitute a discovery of a new rank-record curve.

## Local validation

All 167 regression tests passed on Python 3.12 in 227.372 seconds, including
11 new source-projection and exact-certificate tests. The unchanged rank-31
certificate also replays independently. Python compilation, JSON parsing,
and `git diff --check` passed. These are local results, separate from the
Sage Actions experiment and any hosted CI run.

After adding the final run certificate, the focused lattice-certificate
tests rechecked all four committed run audits and the unchanged rank-31 and
frontier verifiers were replayed again. No mathematical producer code changed
after the full-suite run.

The separately dispatched
[exact-verification run 34162566731](https://github.com/metaforismo/elliptic-rank-31/actions/runs/34162566731),
on commit `9cd1c51524ed1b32286ffe66d17035934231d565`, **failed**:
five pre-existing symbolic test modules cannot import SymPy in the Python
3.13 runner. The log reports `ModuleNotFoundError: No module named 'sympy'`,
157 loaded tests, five errors, and ten skips. All 11 new frontier tests pass
on that runner. The existing workflow does not install SymPy; locally the
167-test pass used SymPy 1.14.0. The proposed correction is to install a
version-pinned symbolic-test dependency before discovery and rerun the
complete suite, without deleting or skipping tests. No workflow correction
has yet been applied, and the hosted suite must not be reported green.

## Acknowledgement

This research made use of the Elliptic Curve Rank Leaderboard maintained by
the Institute for Computer-Aided Reasoning in Mathematics (ICARM), supported
by NSF Grant DMS 2425401. Original submitters and mathematical authors are
credited in the frozen source records and linked primary sources.
