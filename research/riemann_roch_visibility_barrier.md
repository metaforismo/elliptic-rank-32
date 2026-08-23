# The one-function Riemann--Roch visibility barrier

## The theorem

Let \(E/k\) be an elliptic curve with identity \(O\), and let \(f\) be a
nonzero rational function whose pole divisor is exactly \(N(O)\). Suppose
that the zero divisor of \(f\) is

\[
(f)_0=P_1+\cdots+P_N,
\]

where the \(P_i\in E(k)\) are counted with multiplicity. Then

\[
(f)=P_1+\cdots+P_N-N(O).
\]

The divisor is principal, so its class in \(\operatorname{Pic}^0(E)\) is
zero. Under the canonical identification
\(\operatorname{Pic}^0(E)\cong E\), this says

\[
\boxed{P_1+\cdots+P_N=O.}
\]

Thus the rational points exposed as all zeros of a single function always
satisfy at least one exact integral relation. In particular, their
Mordell--Weil span has rank at most \(N-1\).

This statement is unconditional and is independent of analytic rank, BSD,
GRH, height approximations, or a choice of software.

## Consequence for split near-squares

Work on a nonsingular short Weierstrass curve

\[
E:\quad y^2=R(x),\qquad \deg R=3,
\]

over a field of characteristic different from \(2\) and \(3\). Let
\(Q(x)\) have degree \(n\ge2\), and take

\[
f=y-Q(x).
\]

At \(O\), the pole orders of \(x\) and \(y\) are \(2\) and \(3\).
Consequently \(f\) has exact pole order \(2n\). Its norm to \(k(x)\) is

\[
\operatorname{Norm}(y-Q)
=(y-Q)(-y-Q)
=Q(x)^2-R(x).
\]

Suppose

\[
Q(x)^2-R(x)=c\prod_{i=1}^{2n}(x-r_i)
\]

has \(2n\) distinct roots \(r_i\in k\). Then

\[
P_i=(r_i,Q(r_i))
\]

is a zero of \(f\), and the theorem gives

\[
P_1+\cdots+P_{2n}=O.
\]

For \(n=15\), a completely split degree-30 near-square therefore produces
thirty rational points with a forced relation:

\[
\boxed{\operatorname{rank}\langle P_1,\ldots,P_{30}\rangle\le29.}
\]

Searching for thirty independent points solely as the thirty zeros of
\(y-Q(x)\) is structurally impossible. This rules out a representation, not
rank 30 itself.

## Why pole order 31 is the first one-function target

The affine coordinate ring of \(E\) is
\(k[x,y]/(y^2-R(x))\). Every function with no finite poles is uniquely of the
form

\[
f=A(x)+yB(x)
\]

with \(A,B\in k[x]\). The two candidate leading pole orders at \(O\) are

\[
2\deg A
\quad\text{and}\quad
2\deg B+3.
\]

They have opposite parity, so they cannot be equal and their leading terms
cannot cancel. Hence

\[
\operatorname{ord}^{-}_O(f)
=\max(2\deg A,\,2\deg B+3).
\]

A single function carrying \(N\) rational zeros leaves room for at most
\(N-1\) independent classes. Therefore \(N=31\) is minimal for a possible
rank-30 span.

Exact pole order \(31\) forces

\[
\deg B=14,\qquad \deg A\le15.
\]

The norm is

\[
\operatorname{Norm}_{k(E)/k(x)}(A+yB)
=A(x)^2-R(x)B(x)^2,
\]

whose generic degree is

\[
\max(2\deg A,\,3+2\deg B)=31.
\]

This yields the minimal one-function construction target

\[
\boxed{
A(x)^2-R(x)B(x)^2
=c\prod_{i=1}^{31}(x-r_i),
\qquad
\deg A\le15,\quad \deg B=14.
}
\]

## Exact point-recovery conditions

To turn such an identity into thirty-one explicit rational points, the
following conditions are sufficient:

1. \(R\) has nonzero discriminant;
2. \(A^2-RB^2\) has thirty-one distinct rational roots \(r_i\);
3. \(B(r_i)\ne0\) for every \(i\).

Then

\[
P_i=\left(r_i,-\frac{A(r_i)}{B(r_i)}\right)
\]

lies on \(E\), because the norm identity implies

\[
\left(\frac{A(r_i)}{B(r_i)}\right)^2=R(r_i),
\]

and \(A(r_i)+y_iB(r_i)=0\). These points satisfy the forced relation

\[
P_1+\cdots+P_{31}=O.
\]

The identity leaves room for rank 30, but it does **not** prove it. Further
relations can exist, and an exact independence certificate remains necessary.

## Riemann--Roch audit

For \(N\ge1\), Riemann--Roch on a genus-one curve gives

\[
\ell(N(O))=N.
\]

A standard basis is formed by \(1\), powers \(x^i\) with \(2i\le N\), and
terms \(x^j y\) with \(2j+3\le N\). Its pole orders are

\[
0,2,3,4,5,\ldots,N.
\]

The executable certificate checks this basis dimension for every
\(1\le N\le32\), enumerates every degree pattern of exact pole order 30, 31,
and 32, and confirms:

- order 30: \(\deg A=15\), \(\deg B\le13\);
- order 31: \(\deg B=14\), \(\deg A\le15\);
- degree-15 near-square: visible span at most 29;
- minimal one-function rank-30-compatible order: 31.
- minimal one-function rank-31-compatible order: 32.

Run

```bash
python3 research/riemann_roch_visibility_certificate.py \
  --compare certificates/riemann_roch_visibility_barrier.json
```

The proof certificate uses only exact integer arithmetic and the Python
standard library.

## Research consequence

The useful reformulation is not “find thirty split values of a near-square.”
That representation is blocked. The smallest viable single-function problem
is instead:

> Construct a nonsingular cubic \(R\), polynomials
> \(\deg A\le15\), \(\deg B=14\), and thirty-one distinct rational numbers
> \(r_i\) such that \(A^2-RB^2\) splits exactly at the \(r_i\), while the
> resulting point classes have no relation beyond the unavoidable principal
> divisor relation.

This is a sharply defined auxiliary variety and a legitimate inverse-
construction thread. Its difficulty must still be measured against the
record-K3 specialization route.

## The rank-31 target and an exact equivalence

To leave room for 31 independent classes, a single function must have at
least 32 rational zeros.  Exact pole order 32 forces

\[
\deg A=16,\qquad \deg B\le14,
\]

and therefore

\[
\boxed{
A(x)^2-R(x)B(x)^2
=c\prod_{i=1}^{32}(x-r_i).
}
\]

If the roots are distinct and \(B(r_i)\ne0\), the recovered points satisfy
the unavoidable relation

\[
P_1+\cdots+P_{32}=O,
\]

so their span can have rank at most 31.  Certifying any 31 of them as
independent would prove rank at least 31.

There is also a converse that is important for judging the search space.
Suppose \(P_1,\ldots,P_{31}\) are independent rational points and put

\[
P_{32}=-(P_1+\cdots+P_{31}).
\]

Then the 32 points are nonzero and have distinct \(x\)-coordinates: any
collision would give a nontrivial integral relation among the first 31.
Their sum is zero, so

\[
P_1+\cdots+P_{32}-32(O)
\]

is principal.  Its defining function belongs to \(L(32O)\), has exact pole
order 32, and hence has precisely the degree pattern above.  Thus the split
norm-32 formulation is an exact representation of the rank-31 construction
problem—not, by itself, a reduction of its dimension.  It becomes a useful
search method only after imposing additional structure (symmetry, a
low-dimensional family, local conditions, or overlapping point packets).

## Exact reverse experiment on ICARM curve #273

The script `research/icarm273_rr_reverse_certificate.py` applies the converse
construction to the new rank-30 record curve.  Starting from its 30 published
independent points, it computes

\[
P_{31}=-(P_1+\cdots+P_{30})
\]

with the exact group law.  The 31 resulting points have distinct
\(x\)-coordinates and sum to \(O\).  On the model

\[
Y=2y+x,\qquad
Y^2=4x^3+x^2+4a_4x+4a_6,
\]

exact rational linear algebra reconstructs the unique kernel vector for the
evaluation basis

\[
1,x,\ldots,x^{15},Y,Yx,\ldots,Yx^{14}.
\]

The evaluation matrix has rank 30 and nullity one.  The recovered function
has \(\deg A=15\), \(\deg B=14\), and exact coefficient comparison verifies

\[
A^2-(4x^3+x^2+4a_4x+4a_6)B^2
=c\prod_{i=1}^{31}(x-x(P_i)).
\]

This is a successful falsifiable test of the reverse-construction machinery,
not a new independent point.  The committed certificate stores the completed
point, degree/rank data, and cryptographic digests of the large exact
coefficient vectors.
