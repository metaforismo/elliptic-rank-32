# Complete two-split P1/P2/P3 incidence

## Result and boundary

The complete target-section incidence has an exact local presentation with
62 variables and 61 independent equations.  Its expected dimension is one.
This is the Noether--Lefschetz curve obtained by asking for P3 inside the
already certified two-dimensional P1/P2 locus.

No point of this complete incidence is currently known.  In particular, this
document does not claim a modular triple, a rational triple, or a rank-31
elliptic curve over `Q`.

## Surface and section identities

On the two-split chart, let

```text
D=t(t-1)(t-lambda),
Y^2=X^3+3aX^2-6D beta X+D^2 gamma.
```

The section charts are

```text
P1: X=t(t-lambda)U1, Y=t(t-lambda)V1, deg(U1,V1)<=(2,4),
P2: X=(t-lambda)U2,   Y=(t-lambda)V2,   deg(U2,V2)<=(3,5),
P3: X=X3,             Y=Y3,             deg(X3,Y3)<=(4,6).
```

Coefficient comparison gives 9, 11, and 13 equations respectively.

## Three resolved-intersection witnesses

For monic quadratics `h12,h13,h23`, impose

```text
tU1-U2                 = h12 Cx12,
tV1-V2                 = h12 Cy12,

t(t-lambda)U1-X3       = h13 Cx13,
t(t-lambda)V1-Y3       = h13 Cy13,

(t-lambda)U2-X3        = h23 Cx23,
(t-lambda)V2-Y3        = h23 Cy23.
```

The quotient degree bounds are `(1,3)` for `(Cx12,Cy12)` and `(2,4)` for
the other two pairs.  The exact-degree-two open is

```text
gcd(Cx_ij,Cy_ij)=1.
```

The selected same-y branch also requires

```text
gcd(h_ij,Y_i+Y_j)=1.
```

The difference-of-cubics identity

```text
(Y_i-Y_j)(Y_i+Y_j)
=(X_i-X_j)(X_i^2+X_i X_j+X_j^2+3a(X_i+X_j)-6D beta)
```

then proves that the x-factor equation forces `h_ij | (Y_i-Y_j)`.  Thus two
of the displayed y-factor coefficient equations are locally redundant for
each pair.

The count is therefore

```text
variables:
  4 surface + 8 P1 + 10 P2 + 12 P3
  + 8 witness12 + 10 witness13 + 10 witness23 = 62;

raw equations:
  9 + 11 + 13 + 10 + 12 + 12 = 67;

independent local equations:
  67 - 3*2 = 61.
```

Consequently the reduced Jacobian has shape `61 x 62`.  At a future target
triple, rank 61 is the exact smoothness test for a one-dimensional incidence
locus.  No full-triple Jacobian rank is reported before such a point exists.

## Exact P3 elimination

On the safe infinity open, write

```text
X3=x0+x1*t+x2*t^2+x3*t^3+r^2*t^4,
Y3=y0+y1*t+...+y5*t^5+r^3*t^6,
r != 0.
```

The degree-12 coefficient is then automatic.  Coefficients 11 down to 6
solve uniquely and successively for `y5,...,y0`.  Substitution leaves six
numerator equations

```text
N0(r,x0,x1,x2,x3)=...=N5(r,x0,x1,x2,x3)=0.
```

Their total degrees in the five section variables and generic expanded term
counts are:

| equation | degree | terms |
|---|---:|---:|
| N0 | 21 | 566 |
| N1 | 19 | 395 |
| N2 | 17 | 257 |
| N3 | 15 | 150 |
| N4 | 14 | 96 |
| N5 | 12 | 64 |

The six equations in five section variables are the compact algebraic form
of the P3 Noether--Lefschetz condition.  Let

```text
I3=<N0,...,N5> : open_factors^infinity.
```

Eliminating the five section variables gives

```text
J3=I3 intersect k[lambda,y,z,W].
```

On an open where the projection is finite, the zeroth Fitting ideal of the
pushed-forward quotient algebra gives a local equation `Phi_P3` for the
Noether--Lefschetz divisor.  Equivalently, the saturated non-boundary factor
of the Macaulay resultant of the homogenized `N0,...,N5` is a computable
multiple of `Phi_P3`.  The resultant itself is not expanded here; doing so
without saturation would mix the desired divisor with `r=0`, node, and
Kodaira-boundary factors.

Pulling `Phi_P3` back to the smooth two-dimensional P1/P2 incidence gives the
expected one-dimensional full-triple locus.

## Completeness domain

The presentation is equivalent to the original target equations only on the
following explicitly recorded open:

- characteristic different from 2 and 3;
- all denominators of the two-split rational chart are units;
- the exact `II*+3I3+5I1` Kodaira-open factors are nonzero;
- `r X3(0) X3(1) X3(lambda) != 0`;
- every common quadratic is exact and selects the same-y branch;
- the three limiting points at infinity are smooth and pairwise distinct.

Boundary charts, degree drops, deeper nodal contacts, and bad reductions are
not covered by this equivalence and cannot be excluded by its negative
finite-field evidence.

## Modular evidence

All eight existing representative P1/P2 seeds have certified partial
Jacobian rank 28 and tangent dimension 2.  The source search found no P3 on
any of them over their displayed prime fields.

As an independent replay, the triangular P3 system was exhausted on the first
`F_11` seed:

```text
raw (r,x0,x1,x2,x3) tuples: 146410
node-admissible tuples:       110000
P3 hits:                      0
audit SHA-256:
33baf9e529d2e66b4cd347f168822844649c766b12970437d55b555977c40e9a
```

This is a bounded rational-point computation over `F_11`.  It is not a proof
that the geometric P3 fibre is empty over an algebraic closure, and it is not
a characteristic-zero obstruction.
