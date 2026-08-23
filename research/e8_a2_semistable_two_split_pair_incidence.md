# Two-split P1/P2 incidence and formal lifting

## Claim boundary

This note analyzes the finite-field P1/P2 seeds already present in
`certificates/e8_a2_semistable_two_split_target_small_primes.json`.  It proves
that the representative seeds are smooth points of a two-dimensional
incidence scheme and gives one deterministic formal lift over `Z_11`.

It does **not** produce a rational lift, a P3 section, or a rank-31 elliptic
curve over `Q`.

## Exact chart and section equations

Use the four parameters `(lambda,y,z,W)` and put

```text
g = (y^2+lambda-1)/lambda,
q = (1-lambda+lambda*z)/y.
```

The coefficients of `a=a0+a1*t+a2*t^2`, `beta=b0+b1*t`, and
`gamma=c0+c1*t` are determined by

```text
a(0)=1/3,       a(1)=z^2/(3g),  a(lambda)=q^2/3,
beta(0)=W/3,    beta(1)=zW/3,
gamma(0)=W^2,  gamma(1)=gW^2.
```

Equivalently, if

```text
delta1      = z^2/(3g)-1/3,
deltalambda = q^2/3-1/3,
```

then

```text
a2 = (deltalambda-lambda*delta1)/(lambda*(lambda-1)),
a1 = delta1-a2,
b0 = W/3, b1=(z-1)W/3,
c0 = W^2, c1=(g-1)W^2.
```

For

```text
U1=u1_0+u1_1*t+u1_2*t^2,
V1=v1_0+...+v1_4*t^4,
U2=u2_0+...+u2_3*t^3,
V2=v2_0+...+v2_5*t^5,
```

the two section conditions are the coefficient equations of

```text
V1^2 = t(t-lambda)U1^3 + 3aU1^2
       - 6(t-1)beta U1 + (t-1)^2 gamma,                 (P1)

V2^2 = (t-lambda)U2^3 + 3aU2^2
       - 6t(t-1)beta U2 + t^2(t-1)^2 gamma.             (P2)
```

They give respectively 9 and 11 scalar coefficient equations.

To make the resolved intersection-two condition an algebraic incidence
rather than a procedural gcd test, introduce

```text
h  = t^2+h1*t+h0,
Cx = cx1*t+cx0,
Cy = cy3*t^3+cy2*t^2+cy1*t+cy0
```

and impose all coefficients of

```text
tU1-U2-h*Cx = 0,                                        (X-factor)
tV1-V2-h*Cy = 0.                                        (Y-factor)
```

The open conditions `gcd(Cx,Cy)=1` and `gcd(h,tV1+V2)=1` say that `h` is the
exact common quadratic factor and selects the same-y branch.

Thus the displayed presentation has 30 variables and 30 equations:

```text
4 surface + 3 U1 + 5 V1 + 4 U2 + 6 V2 + 2 h + 2 Cx + 4 Cy = 30,
9 P1 + 11 P2 + 4 X-factor + 6 Y-factor = 30.
```

## Why the dimension is two, not zero

Let

```text
x1=(t-lambda)tU1,  x2=(t-lambda)U2,
y1=(t-lambda)tV1,  y2=(t-lambda)V2,
D =(t-lambda)t(t-1).
```

Subtracting the two shifted Weierstrass equations and factoring the difference
of the cubics gives the exact identity

```text
(tV1-V2)(tV1+V2) = (tU1-U2) K,

K=(t-lambda)((tU1)^2+(tU1)U2+U2^2)
  +3a(tU1+U2)-6t(t-1)beta.
```

If `tU1-U2=h*Cx` and `gcd(h,tV1+V2)=1`, this identity forces `h` to divide
`tV1-V2`.  The six Y-factor equations therefore only define the four
coefficients of `Cy`; locally two of them are syzygetic.  Geometrically, the
monic quadratic factor records a discrete choice of two intersection points;
it does not cut two additional moduli.

This agrees with Noether--Lefschetz dimension: the fixed P1/P2 lattice raises
the generic Neron--Severi rank from 16 to 18, leaving a two-dimensional locus.

## Exact finite-field Jacobian result

The script
`research/analyze_e8_a2_semistable_two_split_pair_lifts.py` uses forward
automatic differentiation over the exact prime field.  For all eight
representative seeds, the full 30-by-30 Jacobian has rank 28 and tangent
dimension 2.  The first 28 equations against the first 28 variables form an
invertible transverse minor; the two free variables are `cy2,cy3`.

| seed | prime | rank | tangent dimension | transverse determinant |
|---:|---:|---:|---:|---:|
| 0 | 11 | 28 | 2 | 4 |
| 1 | 11 | 28 | 2 | 7 |
| 2 | 11 | 28 | 2 | 4 |
| 3 | 13 | 28 | 2 | 6 |
| 4 | 13 | 28 | 2 | 2 |
| 5 | 13 | 28 | 2 | 4 |
| 6 | 13 | 28 | 2 | 9 |
| 7 | 13 | 28 | 2 | 3 |

Every determinant is nonzero in its stated field, and every seed satisfies
both exact-degree-two open gcd conditions.

## Formal lift

For seed 0 over `F_11`, fix `cy2=3` and `cy3=10`.  The transverse determinant
is `4 mod 11`, so multivariate Hensel gives a unique point of this transverse
slice over `Z_11`.  Digit-by-digit lifting was carried through modulus

```text
11^8 = 214358881.
```

All 30 residuals, including the two syzygetic equations omitted from the
Newton solve, vanish modulo `11^k` after every step for `1 <= k <= 8`.  The
complete residues and per-precision hashes are in
`certificates/e8_a2_semistable_two_split_pair_incidence_lifts.json`.

This is a formal 11-adic lift of a P1/P2 pair.  Nothing in the calculation
identifies the 11-adic coordinates with rational numbers, and the prior
exhaustive P3 search on the special fibre remains negative.

More sharply, an integral P3 in the same polynomial chart, with the target
node values remaining units, would reduce to a P3 on this `F_11` special
fibre.  The exhaustive special-fibre search excludes that.  Hence this formal
P1/P2 tube cannot acquire the target P3 without leaving the good-reduction
chart (for example through a denominator, degree drop, deeper nodal contact,
or another boundary chart).  A useful next lifting seed must therefore first
contain a P3 modulo `p`, or explicitly model one of those degeneration modes.
