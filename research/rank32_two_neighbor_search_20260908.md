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

## Frozen next experiment

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

## Local algebraic follow-up

The existing stabilized composite has largest absolute entry
15,171,339,103,561,459,844 and old/new fiber intersection
5,715,299,020,407,871. We will test shorter congruent isotropic lifts and
smaller pairing-one witnesses while preserving exact Gram identities and
unimodularity. Any improved marking is a separate certificate; the original
proof remains frozen. Small coefficients alone do not certify a nef class,
compatible periods, a K3 equation, a rational P3, or rank 32.
