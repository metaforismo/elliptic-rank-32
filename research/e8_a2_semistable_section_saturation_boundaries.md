# Section-saturation boundaries in the two-split `E8+A2^3` chart

Date: 2026-08-24

## Scope

This note classifies the section degenerations deliberately omitted by the
simple-contact search engines for the semistable

```text
II* + I3 + I3 + I3 + 5 I1
```

chart.  It separates impossible contacts, wrong component profiles, genuine
intersection strata at infinity, and open characteristic-zero questions.

Throughout,

```text
D=t(t-1)(t-lambda),
y^2=X^2(X+3a)-6D beta X+D^2 gamma,
a gamma-3 beta^2=dD,
```

and

```text
H=4a^2d+12a beta gamma-32beta^3+D gamma^2.
```

The Kodaira open has `gamma1 != 0`, `gcd(D,H)=1`, and squarefree `H`.

## Deep nodal contacts are impossible on the Kodaira open

Let `r` be a root of `D`, write

```text
D=(t-r) delta(t),  delta0=D'(r),
rho^2=3a(r),       s=3beta(r)/rho,
```

so that `gamma(r)=s^2`.  A section through the node has

```text
X=(t-r)U,  y=(t-r)V
```

and reduced square equation

```text
V^2=S=(t-r)U^3+3aU^2-6delta beta U+delta^2 gamma.
```

At the node,

```text
S(r)=(rho U(r)-delta0 s)^2.
```

If the displayed branch term vanishes, the derivative of the surface
relation gives the exact identity

```text
S'(r)
 = delta0^3 (3d rho+s^3)/rho^3
 = 27 delta0^3 H(r)/(4 rho^6).
```

This is nonzero because `D` and `H` are coprime.  But `V(r)=0` would force
`(V^2)'(r)=0`, a contradiction.  Therefore the simple nonidentity-branch
factors used by the engine can never vanish for an actual section on this
open.  Their saturation removes no target component.

By contrast, a section intended to meet the identity component can pass
through the node.  Such a section meets a nonidentity component after
resolution and changes its diagonal height.  For polynomial sections
`P.O=0`, the target diagonals force exactly the counts

```text
P1: two nonidentity I3 components,
P2: one nonidentity I3 component,
P3: zero nonidentity I3 components.
```

Thus the engine's nonvanishing tests at the identity-profile nodes exclude
genuine sections, but only sections with the wrong target Gram.

## Degree drops and singular individual limits are impossible

For a polynomial section put `n=deg(X)<=4`.  Since `gamma1 != 0`, the term
`D^2 gamma` has degree seven.

- If `n<=2`, the right side has degree seven.
- If `n=3`, the unique highest term is `X^3`, of degree nine.
- Only `n=4` can give a square; then the right side has degree twelve and
  `deg(y)=6`.

Consequently every polynomial target section has exact degrees

```text
P1: deg(U,V)=(2,4),
P2: deg(U,V)=(3,5),
P3: deg(X,Y)=(4,6).
```

Its limiting point at infinity is automatically a nonzero smooth point of
`Y^2=X^3`.  Individual degree drops, the cuspidal limiting point `(0,0)`, and
intersection with the zero section are not omitted components.

The height equation also shows that the only non-polynomial diagonal
alternative is `P3.O=1` with `P3` nonidentity at all three `I3` fibres.
It cannot realize the target off-diagonal pairing with `P2`: its sole local
overlap contributes `1/3` or `2/3`, so the required intersection would be
`3-1/3` or `3-2/3`, not an integer.  Hence the polynomial branch is complete
for the displayed Gram.

## Coincident limits at infinity are a genuine target stratum

For two full section coordinate pairs `(X_i,Y_i)`, let `u=1/t` and set

```text
I_infinity = min(
  ord_u(u^4 (X_1-X_2)(1/u)),
  ord_u(u^6 (Y_1-Y_2)(1/u))
).
```

On the smooth limiting cubic this is the exact local intersection
multiplicity.  Therefore the target condition is

```text
I_affine + I_infinity = 2,
```

not merely `I_affine=2` with distinct limiting points.

For `P1/P2`, if

```text
U1=u0+u1*t+u2*t^2,   V1=v0+...+v4*t^4,
U2=w0+...+w3*t^3,    V2=b0+...+b5*t^5,
```

then `I_infinity>=1` is exactly

```text
u2=w3,  v4=b5.
```

It is at least two precisely when, in addition,

```text
u1-lambda*u2 = w2-lambda*w3,
v3-lambda*v4 = b4-lambda*b5.
```

The analogous equations for `P1/P3` and `P2/P3` are obtained by equating
the degree `(4,6)` leading coefficients, and then the next reversed
coefficients.  No new surface chart is required: these are boundary strata
inside the section incidence variety.

## Exact finite-field audit

The total-intersection engine exhausts the dense two-split chart at
`p=5,7,11,13`.

| prime | affine 1 + infinity 1 | affine 2 + infinity 0 | corrected total-two pairs | P3 tests | P3 sections |
|---:|---:|---:|---:|---:|---:|
| 5 | 0 | 0 | 0 | 0 | 0 |
| 7 | 0 | 0 | 0 | 0 | 0 |
| 11 | 0 | 4 | 4 | 644,204 | 0 |
| 13 | 6 | 10 | 16 | 5,940,688 | 0 |

At `p=11`, two previously retained affine-degree-two pairs have
`I_infinity=1`, hence total intersection three; they are correctly removed.
At `p=13`, six previously omitted pairs have one affine and one infinite
intersection.  Exhausting all corrected dense-chart surfaces performs
`6,584,892` `P3` tests and finds no `P3` section.

An independent replay of all printed intrinsic-boundary pairs found four
additional `s_lambda_zero` pairs over `GF(11)` with split `(1,1)`.  A strict-C
targeted replay performed `644,204` further `P3` tests and found no `P3`.
The other boundary-chart target pairs through `GF(13)` have
`I_infinity=0`.

Across the dense chart and the three intrinsic surface-boundary charts for
the four primes, the corrected figures are therefore

```text
raw surface parameter tuples:             37,712
rational-chart tuples:                    32,224
Kodaira-open surfaces:                    17,380
retained P1/P2/P3 sections:                1,558
corrected total-two P1/P2 pairs:              30
  of which infinity-collision pairs:          10
P3 tests on corrected eligible surfaces: 9,456,854
P3 sections:                                   0
complete target triples:                       0
```

Relative to the old simple-infinity runs, `2,871,962` tests were newly run
on genuine infinity-collision pairs.  The corrected aggregate increases by
`2,549,860`, because `322,102` old tests belonged to the two `GF(11)` pairs
whose true intersection is three.

Durable dense-chart artifacts:

```text
research/search_e8_a2_semistable_two_split_total_intersections.c
certificates/e8_a2_semistable_two_split_total_intersections_small_primes.json
tests/test_e8_a2_semistable_two_split_total_intersections.py
```

The strict build and `GF(5)` replay test pass.

## Proved / candidate / open boundary

### Proved

- Deep nonidentity nodal contacts are impossible on the exact Kodaira open.
- Identity-node zeros give the wrong diagonal component counts.
- Polynomial degree drops and singular individual infinity limits are
  impossible.
- Pairwise infinity coincidence is a genuine incidence stratum, measured by
  the displayed exact local length.
- The corrected searches over `GF(5)`, `GF(7)`, `GF(11)`, and `GF(13)` find
  no complete target triple.

### Candidate

- Infinity-collision strata can contribute target pairs: ten exact modular
  pairs occur in the audited domains.  None yet carries `P3`.
- A larger prime or a different characteristic-zero component may contain a
  full triple.

### Open

- No characteristic-zero target triple has been constructed or excluded.
- The finite-field negatives do not prove nonexistence over `Q`.
- Primes above 13 and other lattice-neighbor charts remain unsearched with
  the corrected total-intersection filter.

## Next decisive test

Integrate `I_affine+I_infinity` into both production engines, then search a
new good prime using the cascade on total-two pairs.  The highest-value
variant is to enumerate `P3` through an elimination or square-sieve indexed
by the corrected pair surfaces, rather than scanning all quartics on every
surface.  Any surviving modular triple must then be checked for exact local
profiles, lifted, and certified over `Q`; absence at another prime remains a
bounded negative result.
