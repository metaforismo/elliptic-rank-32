# Rank 31 research programme after ICARM curve #273

## Claim boundary

The current record curve #273 has 30 explicitly certified independent
rational points.  That lower bound is unconditional.  The separate statement
that its rank is exactly 30 uses GRH+BSD and is not part of the lower-bound
certificate.

Consequently, the primary rank-31 question is not “find one more point on
curve #273.”  A point there would collide with a strong conditional upper
bound.  The productive target is a **new curve** with 31 independently
certified points.

## Three representations of the problem

### 1. Specialization accounting

The preceding public record programme starts with a rank-17 elliptic K3
fibration and searches exceptional rational fibers.  The accounting is

\[
29=17+12.
\]

The new rank-30 curve is structurally consistent with, but not proved to come
from, the same programme.  If it does, the likely accounting is

\[
30=17+13.
\]

The next target is therefore

\[
\boxed{31=17+14}.
\]

This route has the best empirical precedent, but it cannot be searched
reproducibly until an explicit rank-17 family, its generic sections, and its
specialization parameter are known.

### 2. Minimal Riemann--Roch norm

On a short Weierstrass curve \(y^2=R(x)\), a single function with 32 rational
zeros has the form

\[
f=A(x)+yB(x),\qquad \deg A=16,\quad\deg B\le14,
\]

and its norm is

\[
A(x)^2-R(x)B(x)^2
=c\prod_{i=1}^{32}(x-r_i).
\]

The points satisfy one forced relation, so their span has rank at most 31.
Conversely, any 31 independent points can be completed by the negative of
their sum to obtain exactly such a degree-32 representation.  This is an exact
equivalence, not by itself a lower-dimensional search.  It becomes useful
only after adding structure such as symmetry, a low-dimensional coefficient
family, or overlapping point packets.

The reverse experiment on curve #273 reconstructs and verifies the analogous
degree-31 norm from its 30 generators.  This validates the machinery but does
not produce an additional independent point.

### 3. Neighbor-fibration reconstruction

The recovered rank-17 lattice has determinant

\[
948=2^2\cdot3\cdot79
\]

and is matched to the \(X(6,79)\) Shimura datum.  A promising neighbor has

\[
E_6+A_3^2+A_1^2,
\]

generic Mordell--Weil rank 3, and fiber configuration

\[
IV^*+I_4+I_4+I_2+I_2+4I_1.
\]

Placing \(IV^*\) at infinity reduces the Weierstrass degrees to
\(\deg A\le5\), \(\deg B\le8\).  This is the most concrete route to recovering
the hidden rank-17 family.

## Exact GF(31) tranches: closed negatively

The public E6/MW3 reduction supplied seven declared split-root cores over
\(\mathbf F_{31}\).  The canonical one-denominator second section reduces to

\[
X_2=C+q_0t(t-1)(t-\lambda)(t-\mu),
\]

leaving exactly \(27\cdot31=837\) pole/coefficient pairs per valid surface.

`research/replay_e6_mw3_p2_split_cores.py` independently reconstructs the
triangular P1 equations, verifies the Weierstrass identities, and exhausts the
finite P2 search.  Relative to the declared seven-core input it finds:

- seven cores examined;
- two nonboundary reconstructed surfaces;
- 1,674 exact \((r,q_0)\) cases;
- one geometric P2 hit (two square-root signs);
- one canonical polynomial P3 on the same surface;
- the exact function-field relation

\[
\boxed{P_1+P_2+P_3=O}.
\]

Thus the only hit is dependent and the tranche contains **no independent
rank-3 modular seed**.  This is a bounded negative theorem, not a rank bound.
The replay does not independently repeat the upstream 893,730-point proof
that its preliminary filtering procedure selects exactly those seven inputs.

The new script `research/replay_e6_mw3_core_exhaustion_gf31.py` does,
however, independently derive the two residual equations

\[
\Delta''(1)=\Delta'''(1)=0
\]

on the fixed coordinate slice and exhaust all 893,730 quadruples with
\(s_1\ne0\).  It exactly reproduces 1,853 simultaneous solutions of those two
equations.  The subsequent geometric reconstruction also explains why a
mere repeated discriminant root is insufficient: a candidate with
\(A=B=0\) is additive rather than a split \(I_2\) fiber.

### Five independently searched coordinate slices

The same exact engine was applied to

\[
(r_0,s_0,x_1)=(1,1,1),(2,1,1),(3,1,1),(1,2,1),(1,1,2).
\]

Across these slices the certificates record:

- 4,468,650 core quadruples exhausted;
- 9,224 simultaneous solutions of the two core equations;
- 12 exact nonboundary surfaces with two verified split \(I_2\) fibers;
- 10,044 canonical \(P_2\) cases up to swapping the two \(I_2\) labels;
- 11,532 polynomial-\(P_3\) cases;
- one geometric \(P_2\) hit and one joint \(P_2/P_3\) surface.

The joint candidate lies on slice \((3,1,1)\), core
\((a_1,a_2,a_4,s_1)=(1,0,14,23)\), but exact group arithmetic in
\(\mathbf F_{31}(t)\) gives

\[
\boxed{P_1+2P_3=O}.
\]

Consequently none of these five slices supplies the desired independent
MW-rank-3 seed.  The result is complete only on their nondegenerate
triangular charts.

The exceptional plane (y_0+y_1=0) has now been treated separately.  There
the failed pivot is replaced exactly by

\[
P_{1,5}\longrightarrow a_3,\qquad P_{1,4}\longrightarrow y_2.
\]

Across the public slice and the five new slices this leaves 172,980 triples.
The exhaustive certificate finds 354 simultaneous core solutions and 325
exact nonboundary (I_4/I_4) surfaces, but zero surfaces carrying both split
\(I_2\) fibers.  Thus the regular and exceptional branches of these slices are
both closed.  All other slices and primes, and other K3 neighbors, remain open.

## Correction and replacement of the `E8+A2^3` chart

The proposed all-additive family `II*+3IV+II` does not realize the desired
rank-three neighbor.  First, its published polynomial-section ansatz omitted
two nonzero high coefficients.  Second, and independently, the family has
constant `j=0`; complex multiplication makes its geometric Mordell--Weil rank
even.  A Picard-rank-19 / Mordell--Weil-rank-3 member is therefore impossible
in that family.

Two exact replacement charts are now available.

The smaller non-isotrivial chart

\[
y^2=x^3-3k^2T^2x+T^2\bigl(c(t-\lambda)^3+2k^3T\bigr),
\qquad T=t(t-1),
\]

has generic fibre configuration `II*+2IV+I3+3I1`.  Its discriminant, split
`I3` normalization, local component labels, and profiled polynomial-section
systems are certified exactly.  Exhaustive finite-field searches find each of
the three required diagonal-height profiles separately, but no complete target
triple has yet survived the exact common-surface filters.

More importantly, this mixed chart is only the `beta=0` specialization of the
full four-dimensional semistable chart

\[
\begin{aligned}
D&=t(t-1)(t-\lambda),\\
A&=-3(a^2+2D\beta),\\
B&=2(a^3+3aD\beta)+D^2\gamma,
\end{aligned}
\]

where `deg(a)<=2`, `deg(beta),deg(gamma)<=1`, and

\[
a\gamma-3\beta^2=dD.
\]

On this identity,

\[
\Delta=-432D^3
\left(4a^2d+12a\beta\gamma-32\beta^3+D\gamma^2\right),
\]

giving `II*+3I3+5I1` on an explicit open set.  The parameter count is four,
exactly `20-rank(U+E8+A2^3)`, so this is the correct ambient chart for the
discriminant-948 Noether--Lefschetz curve rather than a dimensionally narrow
probe.

For three polynomial sections, the target Mordell--Weil Gram is equivalent to
nonidentity-component counts `(2,1,0)`, opposite labels on the unique overlap
of the first two sections, and all three pairwise section intersections equal
to `2`.  Its determinant is `316/9`; multiplying by `disc(A2^3)=27` gives
exactly `948`.  This is an exact recognition criterion, not yet an existence
proof.

### The target lattice is primitive

The height formula

\[
\hat h(P)=4+2(P\mathbin{.}O)-\frac{2}{3}n(P)
\]

forces, up to permutations and simultaneous component-label reversal, the
unique component orbit

```text
P1 = (1,1,0),   P2 = (2,0,0),   P3 = (0,0,0).
```

It also forces `P1.P2=P1.P3=P2.P3=2` and trivial torsion.  After multiplying
the height pairing by three, the integral even Gram matrix is

\[
M=\begin{pmatrix}8&-1&0\\-1&10&0\\0&0&12\end{pmatrix},
\qquad \det M=948.
\]

The only possible nonzero order-two dual class is represented by `P3/2`; its
`M`-norm is `3`, so it is not even-isotropic.  Hence there is no proper even
integral overlattice: if the three sections exist, their rank-three lattice is
primitive.  The calculation is replayed by
`research/certify_e8_a2_target_lattice.py`.

### Necessary two-split chart and its boundaries

The target profiles do not require the `I3` fibre at `t=1` to split.  Requiring
all three fibres to split therefore imposed an unnecessary square condition.
On the dense necessary chart, with parameters `(lambda,y,z,W)`, put

\[
g=\frac{y^2+\lambda-1}{\lambda},\qquad
q=\frac{1-\lambda+\lambda z}{y},\qquad
a_2=\frac{(g-z)^2}{3gy^2},
\]

and interpolate

```text
3a(0)=1,       3a(1)=z^2/g,      3a(lambda)=q^2,
3beta(0)=W,    3beta(1)=zW,
gamma(0)=W^2, gamma(1)=gW^2.
```

These formulas satisfy `a*gamma-3*beta^2=dD` identically and give a dense
birational four-parameter chart in which precisely the two target-forced
fibres are split.  The three missing intrinsic divisors are covered by the
exact charts `s0=0`, `s_lambda=0`, and `beta(1)=gamma(1)=0`.  Together they
cover the marked two-split locus before imposing that the residual quintic is
squarefree and coprime to `D`.  The identities are independently recomputed by
the two chart certificates under `certificates/e8_a2_semistable_two_split_*`.

### Corrected exhaustive small-field result

The profiled section cascade was run on the dense chart and all three intrinsic
boundary charts over GF(5), GF(7), GF(11), and GF(13):

The original search counted only affine intersections between sections.  The
exact local intersection at infinity is

\[
I_\infty(P,Q)=\min\!\left(\operatorname{ord}_u u^4\Delta X(1/u),
\operatorname{ord}_u u^6\Delta Y(1/u)\right),
\]

so the target equation is `I_affine+I_infinity=2`.  Replaying the dense chart
and all three intrinsic boundaries with this correction gives:

| search locus | Kodaira-open surfaces | section tests | corrected target `P1/P2` pairs | subsequent `P3` tests | full triples |
|---|---:|---:|---:|---:|---:|
| dense two-split chart | 11,744 | 50,540,496 | 20 | 6,584,892 | 0 |
| three boundary charts | 5,636 | 20,761,918 | 10 | 2,871,962 | 0 |
| **total** | **17,380** | **71,302,414** | **30** | **9,456,854** | **0** |

Ten target pairs lie on the previously omitted infinity-collision stratum:
six in the dense chart over `GF(13)` and four on the `s_lambda=0` boundary over
`GF(11)`.  Conversely, two old `GF(11)` pairs have two affine intersections
and one at infinity, hence total intersection three, and are no longer target
pairs.  The corrected 30 pairs are genuine finite-field hits with the required
profiles and total intersection two.  None is a rank-three seed because no
`P3` survives on the same surface.

The other apparent saturation boundaries are now closed exactly.  On a
nonidentity branch at a root `r` of `D`, vanishing of the simple-contact term
would force

\[
S'(r)=\frac{27D'(r)^3H(r)}{4\rho^6}\ne0,
\]

contradicting that `S=V^2` has a double zero.  Moreover `gamma_1 != 0` forces
every polynomial target section to have exact degrees `(deg X,deg Y)=(4,6)`;
degree drops and cuspidal individual limits at infinity are impossible.  Thus
pairwise collision at infinity is the only genuine target stratum omitted by
the earlier simple-contact search.  See
`research/e8_a2_semistable_section_saturation_boundaries.md`.

The negative conclusion remains bounded by the four listed fields, the
displayed charts, and the Kodaira open; it is not a nonexistence theorem over
`Q`.

For the eight deduplicated dense-chart representatives, the exact 30-by-30
incidence Jacobian has rank 28.  Thus the `P1/P2` incidence is smooth of
dimension two at every representative seed.  One transverse slice lifts
deterministically over `Z_11` through modulus `11^8`, with all 30 residuals
zero at each precision.  It is only a formal 11-adic `P1/P2` lift.  The special
fibre contains no `P3`, so the same good-reduction tube cannot acquire the
target third section without crossing a denominator, degree drop, deeper
contact, or another degeneration.  See
`research/e8_a2_semistable_two_split_pair_incidence.md`.

### Larger-prime pair-first search

A second engine samples the first `(surface,P1)` incidence without replacement
and then exhausts every `P2` and `P3` conditional on each sampled `P1` hit.
With deterministic seed `20260824` it processed 20,000,000 first-stage
incidences at each of `p=17,19,23`:

| prime | first-stage coverage | `P1` hits | `P2` sections | target pairs | exhaustive downstream `P3` tests | `P3` |
|---:|---:|---:|---:|---:|---:|---:|
| 17 | 6.63% | 172 | 9 | 2 | 2,839,714 | 0 |
| 19 | 2.94% | 115 | 5 | 1 | 2,476,099 | 0 |
| 23 | 0.735% | 69 | 6 | 0 | 0 | 0 |

This is `88,931,043` section tests in total and `5,315,813` exhaustive `P3`
tests after the three sampled target pairs.  Only the first stage is sampled;
the result is not exhaustive at these larger primes.

### Exact full incidence and `P3` elimination

The full three-section incidence now has an exact presentation with 62
variables and 67 raw equations.  The difference-of-cubics identity gives two
local syzygies for each of the three section pairs, leaving a reduced
`61 x 62` Jacobian.  A full triple of rank 61 would therefore be a smooth point
of the expected one-dimensional Noether--Lefschetz locus.  No such point is
known, so no observed rank is reported.

On the infinity open `x4=r^2`, `y6=r^3`, `r!=0`, the `P3` equations solve
triangularly for `y5,...,y0` and reduce to six exact numerator equations
`N0=...=N5=0` in `(r,x0,x1,x2,x3)`.  Their total degrees are
`21,19,17,15,14,12`.  Saturating this ideal by the declared open factors and
eliminating the section variables defines the desired divisor in the
four-dimensional surface chart.  The global saturated Macaulay
resultant/Fitting equation has not yet been expanded.  This is the concrete
symbolic bottleneck, rather than an unspecified search for a third section.

### Exact bridge to Elkies' `X(6,79)` landscape

The positive essential lattice forced by the `E8+A2^3` three-section target is
even of rank 17 and determinant 948.  Exact short-vector enumeration gives
precisely 258 roots, hence root system `E8+A2^3`, and its discriminant group is
cyclic of order 948 with generator norm `1709/948 mod 2Z`.

The rank-three period lattice used by the existing `X(6,79)` workflow is

\[
T=\begin{pmatrix}-316&0&288\\0&474&-15\\288&-15&-262\end{pmatrix}.
\]

It has signature `(2,1)`, determinant `-948`, and cyclic discriminant form
generated by `485/948`.  The congruence

\[
1709\cdot67^2\equiv485\pmod{1896}
\]

proves that the finite quadratic forms agree.  The transparent
`A11+K6(det=79)` neighbor seed is joined independently by the analogous
multiplier 101.  Finally, the even Clifford algebra of `T` has quaternion
discriminant 6 and the remaining level factor is 79, reproducing the
`X(6,79)` Shimura datum.

This is a genuine lattice bridge, not only equality of determinants.  It is
still not an explicit stable isometry, neighbor chain, moduli map, or
Weierstrass model.  Elkies' published Shimura curve

\[
u^2=16t^6-19t^4+88t^2-48
\]

has the rational non-CM point `t=14/13`, `u=16064/2197`; transporting that
point through an explicit neighbor map into the two-split chart is now a
concrete construction problem.  See
`research/e8_a2_shimura_bridge.md`, [Elkies 2007](https://arxiv.org/abs/0709.2908),
and [Elkies 2008](https://arxiv.org/abs/0802.1301).

## Route ranking after the experiment

1. **Build the explicit `X(6,79)` neighbor/moduli map.**  Connect the certified
   `E8+A2^3` essential lattice to the frozen determinant-948 chain and
   transport Elkies' rational Shimura point into the two-split chart.  This
   replaces blind parameter search by a specified rational input, but a
   transported three-section model must still be checked directly.
2. **Compute the `P3` Noether--Lefschetz divisor over the `P1/P2` surface.**
   Use the six exact remainder equations, saturate their boundaries, and
   compute a Fitting/elimination equation or a square-sieve at new good
   primes.  Lift only a special fibre already containing the full same-surface
   triple; the present 11-adic pair tube is locally excluded from doing so.
3. **Recover the E6/MW3 family on another finite-field chart or prime.**  Six
   GF(31) tranches (the published slice plus five new slices), including their
   exceptional pivot planes, are closed, but
   the neighbor construction remains the shortest route to an explicit
   rank-17 model.  Move to new coordinate slices and good primes, and
   enforce group-law independence before Hensel lifting.
4. **Obtain the discoverers' specialization certificate.**  An explicit
   family parameter for curve #273 would immediately turn the problem into a
   reproducible \(17+14\) specialization search.
5. **Use structured norm-32 subfamilies.**  Search only after a symmetry or
   packet design reduces dimension and after local finite-quotient conditions
   are integrated into candidate generation.
6. **Do not prioritize direct point search on #273.**  It is useful for
   auditing descent and the conditional upper-bound boundary, not as the main
   path to a new record.

## Certification gate for a genuine record

A candidate is not rank 31 until all of the following are exact:

1. a nonsingular curve over \(\mathbf Q\);
2. 31 explicit rational points satisfying its equation;
3. an exact independence certificate, preferably via a product of finite
   quotients \(E(\mathbf F_p)/\ell E(\mathbf F_p)\);
4. a clear separation between the unconditional lower bound and every
   analytic or BSD-dependent upper bound.

Scores, numerical regulators, root numbers, finite-field seeds, and small
residuals are discovery diagnostics only.
