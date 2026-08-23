# Correction to the proposed `E8 + A2^3` rank-jump system

## Result

The public one-section ansatz

```text
x = q^2 + r,
y = q^3 + s,
deg(q)=2, deg(r),deg(s)<=1
```

does not give an eight-equation system.  For

```text
G(t) = [t(t-1)(t-lambda)]^2(t-mu)
```

the omitted high coefficients of `y^2-x^3-G` are

```text
[t^9] = -3 q2^4 r1,
[t^8] after r1=0 = -3 q2^4 r0.
```

Thus, on the declared chart `q2 != 0` in characteristic different from 3,
the old ansatz forces `r=0`.  Retaining only coefficients `t^0,...,t^7`
silently drops two necessary equations.

## Correct degree-balanced ansatz

The missing term is the first binomial correction to `(q^2+r)^(3/2)`:

```text
x = q^2 + r,
y = q^3 + (3/2) q r + h,
deg(q)=2, deg(r),deg(h)<=1.
```

Indeed,

```text
y^2-x^3
  = 2 q^3 h - (3/4) q^2 r^2 + 3 q r h + h^2 - r^3,
```

which has degree at most seven identically.  There are now genuinely eight
coefficient equations.  On `q2 != 0`, the top two give

```text
h1 = 1/(2 q2^3),
mu = -2 h0 q2^3 - 6 h1 q1 q2^2 - 2 lambda
     + (3/4) q2^2 r1^2 - 2.
```

At `t=0`, write the point on the cuspidal cubic as `(u^2,u^3)`.  Then

```text
r0 = u^2-q0^2,
h0 = (u-q0)^2(2u+q0)/2.
```

For `u != 0`, the coefficient of `t` forces `u != q0` and

```text
r1 = q1(u-q0) + 1/[3 q2^3(u-q0)].
```

For `u=0`, that coefficient vanishes identically and `r1` remains free.
Both branches therefore reduce to four residual equations in five variables.

The corrected identity has exact modular solutions.  The C enumerator
`search_e8_a2_rank3_jump_corrected.c` verifies all thirteen possible
coefficients, rather than trusting only the residual system.

## Independent obstruction to the proposed target family

There is a second, conceptual obstruction.  Every surface in

```text
y^2 = x^3 + G(t)
```

has constant `j=0`.  Over the algebraic closure it has the automorphism

```text
sigma: (x,y) -> (zeta_3 x,y),
sigma^2+sigma+1 = 0.
```

Consequently the free geometric Mordell--Weil group tensored with `Q` is a
vector space over `Q(zeta_3)`, so its integral rank is even.  With fibers
`II* + 3 IV + II`, the reducible root rank is 14 and Shioda--Tate gives

```text
rho = 16 + rank(MW).
```

Hence `rho` is even in this family.  It cannot contain the intended generic
Picard-rank-19 / Mordell--Weil-rank-3 fibration.  Corrected modular sections
are valid points on other `j=0` rank-jump loci, but they do not repair this
parity obstruction.

## Consequence for the rank-31 program

The isotrivial two-parameter candidate is not the desired bridge from the
rank-3 neighbor to the rootless rank-17 fibration.

There is, however, a non-isotrivial replacement chart.  Put

```text
T = t(t-1),  z = t-lambda,
A = -3 k^2 T^2,
B = T^2 (c z^3 + 2 k^3 T).
```

Then

```text
Delta = -432 c T^4 z^3 (c z^3 + 4 k^3 T).
```

Under the explicit nonvanishing and squarefreeness conditions certified by
`derive_e8_a2_mixed_family.py`, this has fibres

```text
II* + IV + IV + I3 + 3 I1.
```

The exact local component calculation shows that the target Gram is equivalent
to three polynomial sections with nonidentity-fibre counts `(2,1,0)`, opposite
component labels on the unique overlap of the first two sections, and all three
pairwise section intersections equal to `2`.  Its Mordell--Weil determinant is
`316/9`; multiplying by `disc(A2^3)=27` recovers exactly `948`.

The mixed chart is still only a two-dimensional subfamily of the full
four-dimensional `E8+A2^3` moduli chart.  The full semistable chart is

```text
D = t(t-1)(t-lambda),
a = a0+a1*t+a2*t^2,
beta = beta0+beta1*t,
gamma = gamma0+gamma1*t,

A = -3(a^2+2D beta),
B = 2(a^3+3aD beta)+D^2 gamma,
```

subject to the polynomial identity

```text
a gamma - 3 beta^2 = d D.
```

On that identity,

```text
Delta = -432 D^3
        (4a^2 d + 12a beta gamma - 32 beta^3 + D gamma^2),
```

which generically gives `II* + 3 I3 + 5 I1`.  The mixed chart is exactly the
specialization

```text
beta=0,  a=kT,  gamma=cz,  d=kc.
```

Therefore a negative search in the mixed chart is not a global obstruction.
The full four-dimensional semistable chart and the non-isotrivial E6 neighbor
remain viable.  In either chart, exact component profiles and height-pairing
intersections must be enforced before independence tests or p-adic lifting.
