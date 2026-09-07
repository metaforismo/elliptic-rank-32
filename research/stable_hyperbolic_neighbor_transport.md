# Explicit integral hyperbolic transport of the seven-edge rootless bridge

## Verified result and boundary

The exact seven-edge bridge found by
[run 34164507132](https://github.com/metaforismo/elliptic-rank-31/actions/runs/34164507132)
connects the E8+A2^3 essential lattice to the frozen determinant-948 rootless
lattice. Its path matrices, meeting isometry and rational 17-dimensional
composite have been independently replayed.

A separate construction now gives an **integral unimodular 19-dimensional
isometry** after adjoining a hyperbolic plane U. All seven edge extensions
have determinant -1; the complete matrix has determinant +1. Standard-library
Fraction/Bareiss arithmetic and independent SymPy matrix checks agree.

This is an abstract lattice transport, not a geometric change of K3
fibration. No effective or nef representatives, compatible periods or ample
cones, Riemann-Roch pencils, birational maps, rational P3, or rank-32 curve
have been certified. The fiber and section vectors below are **classes**,
not equations of divisors or coordinates of rational points.

## Construction for one edge

Let G be the positive parent Gram matrix, q(y)=y^t G y/2, and let T be the
recorded child basis, so T^t G T=G_child. Work in

\[
S=U\oplus N(-1),\qquad U=\langle e,f\rangle,\quad e^2=f^2=0,\ e.f=1.
\]

Lift the selected isotropic vector modulo p to an integer y with
q(y) divisible by p squared, using the even lift in the pinned
[Sage neighbor construction](https://github.com/sagemath/sage/blob/10.9/src/sage/quadratic_forms/quadratic_form__neighbors.py).
Independently check that T inverse applied to y/p is integral. The bad-prime
radical case is rejected, not silently interpreted as an even neighbor.

Set F'=(q(y)/p,p,y). Then F'^2=0 and e.F'=p. Compute a Bezout witness Z1
with F'.Z1=1, rejecting divisibility greater than one. Since S is even,

\[
Z=Z_1-\frac{Z_1^2}{2}F'
\]

is integral, Z^2=0 and F'.Z=1. Thus F',Z span an integral U; Z-F' is an
abstract square-minus-two class meeting F' once.

For each child-basis column w of T, choose an index j with y_j nonzero
modulo p and an integer u inverse to y_j modulo p. Put t=p w_j u. Define

\[
H(w)=\left(\frac{(w,y)}p-t\frac{q(y)}{p^2},\ -t,\ w-\frac tp y\right).
\]

The neighbor congruences make these entries integral; this is also checked
entry by entry. Direct expansion gives F'.H(w)=0 and
H(w).H(w')=-(w,w'). Project once more:

\[
H'(w)=H(w)-(Z.H(w))F'.
\]

Then H'(w) is orthogonal to both F' and Z and has the same pairings as H(w).
The columns F', Z and the 17 vectors H'(w) form the matrix Phi. We explicitly
verify Phi^t S_parent Phi=S_child and det(Phi)=+1 or -1. This proves an
integral basis change, rather than only an equality of rational quadratic
spaces.

## Complete path

Compose the three origin-side extensions and the four rootless-side
extensions. Lift the 17-dimensional meeting isometry Q as diag(I_2,Q).
The end-to-end matrix is C_origin diag(I_2,Q) C_rootless inverse. Its integral
entries, determinant +1, inverse integrality, and full Gram identity are
verified independently.

The resulting basis is not optimized for geometry: its largest absolute
entry is 15,171,339,103,561,459,844. The old fiber intersects the transported
rootless fiber class in 5,715,299,020,407,871. These values warn against
mistaking an abstract integral marking for a practical nef divisor or a
small Riemann-Roch pencil. The individual edge intersections remain 2 or 5.

## Reproduction and provenance

```sh
python3 research/stable_hyperbolic_neighbor_transport.py \
  --compare certificates/stable_hyperbolic_rootless_bridge_34164507132.json
python3 -m unittest discover -s tests \
  -p 'test_stable_hyperbolic_neighbor_transport.py' -v
```

The committed raw bridge is retained byte-for-byte from the checked artifact,
SHA-256 `0714dd21746e9be94e53e91b675d445487411e6c9b140af17523543daaf60139`.
The new certificate record hash is
`40f65721819c5e4d4725814a7732f1bea77543b52afe778e0e6bdc84debad5a0`;
the final stable-matrix hash is
`ee16d707181177ffab193826dd2ab6612c7eab612c0d65f0befe8bb22a85f1b5`.

All eight dedicated tests passed in 58.592 seconds, including independent
SymPy checks of every edge and the full composite/inverse, tampering
rejection, a radical p=2 failure, and two additional frozen p=2/p=5 examples.
Those two additional examples come from the earlier inconclusive search;
their local matrix identities are proved without turning that stopped run
into a finite-negative certificate. Experiment dates refer to UTC.

The complete local regression suite then ran 196 tests in 444.335 seconds:
195 passed, with only the pinned-Sage API test skipped on ordinary CPython.
That unchanged API test passed separately in Sage run `34164060337`.
Both workflow definitions also passed `actionlint`; the staged diff passed
`git diff --check`.
