# From Kneser neighbors to explicit elliptic-K3 transformations

## 2026-09-07 update

The earlier missing basis-data stage is now resolved for a new seven-edge
E8+A2^3-to-rootless chain: its rational edge maps and exact meeting isometry
are archived and independently replayed. We have also constructed explicit
integral 19-by-19 isometries after adjoining U, for every edge and for the
composite. See [the stable transport certificate](stable_hyperbolic_neighbor_transport.md).

This does **not** resolve the geometric requirements below. The new fiber
and section vectors are abstract lattice classes, with no effective/nef
representatives, compatible K3 periods, local expansions, or Riemann-Roch
pencils certified. The subsequent discussion of absent matrices describes
the original discovery workflow, not the new archived chain.

## Conclusion

The present positive-definite Kneser neighbor chain is **not enough** to
produce a sequence of Weierstrass transformations.  It gives useful exact
information about possible essential lattices, but it does not specify the
two hyperbolic-plane embeddings in the Neron--Severi lattice, the effective
fiber divisors, or the rational functions that define the new pencils.

There are two distinct notions that have been compressed into the word
"neighbor":

1. a Kneser `p`-neighbor of positive-definite lattices, computed from an
   isotropic line modulo `p`; and
2. an `r`-neighbor step between two elliptic fibrations on the *same* K3
   surface, where their fiber classes have intersection `r`.

An elliptic `r`-neighbor step induces an `r`-neighbor relation between the
corresponding frame lattices.  The converse requires extra marked geometric
data and is not supplied by two abstract Gram matrices.  Elkies--Kumar state
the forward relation explicitly: the geometric step is an explicit
isomorphism between two presentations of `NS(X)`, and the common orthogonal
submodule projects to index `r` in both frame lattices
([Section 5, pp. 15--16](https://arxiv.org/pdf/1209.3527#page=15)).

This distinction changes the next milestone.  Replaying the missing Gram
matrices is useful, but the first genuinely geometric certificate must record
a **marked Neron--Severi transport and a Riemann--Roch pencil**, not merely the
prime and the child Gram matrix.

## Evidence boundary

- **Primary-source facts:** Sage's exact Kneser construction and returned
  basis matrix; Elkies--Kumar's definition of a geometric `r`-neighbor and
  their Riemann--Roch construction; Kumar's degree-two cases; the explicit
  `jac2`/`jac3` auxiliary code; and Elkies' published description of the
  `X(6,79)` anchor.
- **Repository observations:** the present workflow retains Gram matrices,
  primes, and vectors but not the Kneser/LLL basis maps; the frozen path begins
  with two prime-five moves.
- **Inferences from those facts:** the positive-lattice path is a routing
  certificate only; a marked geometric lift needs the additional data listed
  below; and a 2-/3-only marked path would avoid implementing a new degree-five
  surface algorithm.
- **Open gaps:** no explicit `X(6,79)` Weierstrass seed or moduli map was found
  in the public primary package surveyed here, no geometric lift of the frozen
  Kneser path is known, and no explicit elliptic-K3 5-neighbor implementation
  was found.  These are negative audit results, not nonexistence theorems.

## 1. What the Kneser computation records

For a lattice `(N,q)` and a suitable vector `y`, Sage constructs

\[
 N_y = \mathbf Z\frac{y}{p}
       + \{x\in N:(x,y)\equiv0\pmod p\}.
\]

The official Sage implementation forms a rational basis matrix `B` for this
lattice and can return it through `return_matrix=True`
([Sage 10.9 source, lines 141--248](https://github.com/sagemath/sage/blob/10.9/src/sage/quadratic_forms/quadratic_form__neighbors.py#L141-L248)).

The current search in
`.github/workflows/probe-rank29-height-lattice.yml` instead calls
`find_p_neighbor_from_vec(p, vector)` without `return_matrix=True`, then LLL
reduces the child Gram matrix and retains only:

- the parent and child Gram matrices/hashes;
- `p` and the input vector in the current reduced basis;
- root and norm-four counts.

Consequently, the workflow discards both:

- the rational Kneser basis matrix `B`; and
- the integral LLL change of basis `U`.

If `G_i` is the parent Gram matrix, a replay should retain a total rational
matrix `M_i` satisfying

\[
 M_i^{\mathsf T}G_iM_i=G_{i+1}.
\]

With Sage's conventions this is obtained by composing the matrix returned by
`return_matrix=True` with the LLL transformation.  Saving `M_i`, its inverse,
and the common index-`p` sublattice would make the *positive-lattice* chain
hermetic.  It would still not make it a geometric K3 transformation.

## 2. The marked lattice data needed for geometry

Fix one abstract Neron--Severi lattice `S` of signature `(1,18)`.  An elliptic
fibration with section is not specified by its positive essential lattice
alone.  It requires a primitive embedding

\[
 U=\langle F,O+F\rangle\hookrightarrow S,
 \qquad F^2=0,\quad O^2=-2,\quad F\cdot O=1,
\]

where `F` is nef.  The positive essential lattice is only the negative of the
orthogonal complement of this marked `U`.

For every node in a geometric chain one therefore needs:

1. a fixed basis of `S` and its full intersection matrix;
2. the vectors `F` and `O` in that basis;
3. an explicit identification of `U^perp` with the displayed essential Gram
   matrix;
4. the roots represented by actual effective fiber components, including
   multiplicities and the identity component;
5. the section classes, their component intersections, and their fields of
   definition;
6. a nef/ample chamber witness that chooses effective roots rather than their
   negatives.

For an edge `i -> i+1`, the lattice certificate should then contain an
integral change of marking

\[
 \Phi_i:U\oplus N_i(-1)\longrightarrow U\oplus N_{i+1}(-1)
\]

and the images of the old and new fiber and zero-section classes.  Before an
algebraic birational map is known, one must also check that the marking has
the correct discriminant action, period orientation, and ample chamber.  An
abstract isometry of indefinite lattices need not preserve the chosen K3
period or effective cone.  Once explicit mutually inverse rational maps of
surfaces have been verified, those maps themselves provide the stronger
geometric identification.

This is why a Kneser matrix `M_i` is necessary but not sufficient: it relates
two complements inside a rational quadratic space, while `Phi_i` must explain
how the hyperbolic plane changes and must select actual divisor classes on the
same surface.

## 3. The geometric neighbor algorithm

Let the current fibration be

\[
 y^2=x^3+a_2(t)x^2+a_4(t)x+a_6(t)
\]

with fiber class `F`, and let `F'` be the proposed new primitive effective nef
isotropic divisor.  Elkies--Kumar's construction is

\[
 H^0(X,\mathcal O_X(F'))
   =\langle s_0,s_1\rangle,
 \qquad u=s_1/s_0.
\]

To calculate this pencil, decompose `F'=F'_hor+F'_ver` relative to the old
fibration.  If `r=F.F'`, then `F'_hor` has degree `r` on the generic elliptic
fiber.  Choose a basis of the degree-`r` Riemann--Roch space there and impose
the zero/pole conditions forced by every component of `F'_ver`.  These are
linear conditions and leave a two-dimensional space
([Elkies--Kumar, Section 5](https://arxiv.org/pdf/1209.3527#page=15)).

A complete geometric edge must therefore record:

1. `F'` as an explicit sum of old sections and resolved fiber components;
2. exact checks that `F'` is primitive, effective, nef, and has square zero;
3. a divisor class `D` with `D.F'=1`, or two classes whose intersections with
   `F'` are coprime, proving that the new genus-one fibration has a section;
4. the horizontal/vertical decomposition of `F'`;
5. a basis of the generic-fiber Riemann--Roch space;
6. local expansions at every relevant resolved fiber component and the
   resulting linear pole-cancellation conditions;
7. the two surviving functions `s_0,s_1` and the new base parameter `u`;
8. the resulting genus-one equation over the field `K(u)`;
9. a chosen rational point/section and an explicit isomorphism, not only the
   Jacobian invariants, to a minimal Weierstrass model;
10. forward and inverse rational maps in `(t,x,y)` and `(u,X,Y)`;
11. Tate-algorithm checks for the new fibers and exact transport of all
    sections and height pairings;
12. the denominator and discriminant bad locus on the moduli base.

The section criterion is genuinely extra information.  A primitive nef
isotropic class gives a genus-one fibration, while a divisor intersecting it
once gives a section
([Elkies--Kumar, Lemmas 1--2 and Corollary 3](https://arxiv.org/pdf/1209.3527#page=4)).
Replacing a genus-one curve by its Jacobian without tracking this point can
lose the torsor and does not by itself give a birational model of the original
K3 surface.

### Necessary versus sufficient data

There are three different certification levels, and they should not be
conflated.

1. **Positive-lattice replay.**  The two Gram matrices, `p`, the isotropic
   line, `B`, `U`, and the common index-`p` sublattice are sufficient to prove
   a Kneser-neighbor relation.  They are neither necessary nor sufficient to
   prove that two elliptic fibrations occur on the same K3 surface.
2. **Neighbor construction on a fixed K3.**  Within the geometric neighbor
   method, one must at least know the source surface and the actual divisor
   `F'`, compute the two-dimensional pencil `H^0(X,O(F'))`, and trivialize the
   generic genus-one torsor by a section (or an equivalent rational point).
   Abstract classes without effective representatives do not determine the
   rational functions in the pencil.  A marked `NS(X)` and resolved component
   data are the practical exact inputs from which these objects are computed.
3. **Geometric sufficiency.**  Equations for the two minimal surfaces, the new
   parameter `u` as a rational function on the old surface, and mutually
   inverse rational maps on dense open subsets are sufficient to prove that
   the step is a birational change of elliptic fibration.  Symbolic
   substitution in both directions, together with the fiber and section
   divisor checks, then certifies the advertised marking.  At this point a
   separate period-orientation or ample-cone test is a useful independent
   check, not an additional logical prerequisite for identifying the two
   explicit surfaces.

For a **moduli map**, the same data must be rational in the moduli parameters
and accompanied by a Zariski-open bad locus on which denominators,
discriminants, or inverse maps fail.  For transport at the single rational
`X(6,79)` point, it is sufficient to supply and verify the specialized data;
one does not first need a universal family over the whole Shimura curve.

## 4. Explicit algorithms for 2- and 3-neighbors

### 4.1 Two-neighbor

Kumar gives all horizontal degree-two cases after translating sections
([Appendix A](https://arxiv.org/pdf/1105.1715#page=44)):

- `2O`, with generic-fiber basis `{1,x}`;
- `O+P` for a non-2-torsion section `P=(x_0,y_0)`, with basis
  `{1,(y+y_0)/(x-x_0)}`;
- `O+T` for a 2-torsion section `T=(0,0)`, with basis `{1,y/x}`.

After the vertical pole conditions select `u`, substitution gives a double
cover

\[
 z^2=g(t,u)
\]

with `g` cubic or quartic in `t`.  A rational point converts this quartic to a
Weierstrass equation.  Kumar also writes out the corresponding birational
conversion rather than only the new invariants
([Appendix A.2](https://arxiv.org/pdf/1105.1715#page=45)).

The primary auxiliary source for Elkies--Kumar's Hilbert-modular-surface
paper contains a compact executable version.  For example, `12/12.txt` uses

```text
u = (x + e/(f+1)*t^2*(t-1))/(t^3*(t-1))
```

and then calls `jac2` on the resulting quartic.  `jacobians` computes the
classical binary-quartic invariants exactly.  These files are distributed in
the paper's [arXiv source archive](https://export.arxiv.org/e-print/1209.3527),
whose `README.txt` explicitly says that the per-discriminant files contain the
2- and 3-neighbor steps.

### 4.2 Three-neighbor

For horizontal divisor `3O`, the generic-fiber basis is `{1,x,y}`.  The new
parameter has the form

\[
 u=c(t)+d(t)x+e(t)y.
\]

Solving for `y` and substituting produces a plane cubic; its Jacobian gives
the next Weierstrass equation
([Elkies--Kumar, Section 5](https://arxiv.org/pdf/1209.3527#page=15)).

The auxiliary file `17/17.txt` makes the missing local calculation concrete:
it solves for a polynomial expression

```text
y + quadratic(t)*x + quintic(t)
```

whose poles cancel on the prescribed components, substitutes it as the new
parameter, and calls `jac3`.  The file `41/41.txt` similarly states which
component the candidate function must be regular on before solving its
coefficients.  Thus the component marking and local expansions are not
optional implementation details; they determine the pencil.

### 4.3 Five-neighbor

No explicit elliptic-K3 **5-neighbor** implementation was found in the
surveyed primary Elkies sources or in the complete auxiliary archive for
arXiv:1209.3527.  Elkies describes chains of 2-neighbors and occasional
3-neighbors for the high-rank K3 computation
([2007 lectures, Lecture III](https://arxiv.org/pdf/0709.2908#page=21)), and
Elkies--Kumar's executable helper contains only `jac2` and `jac3`.

There is nevertheless a general route.  If `F.F'=5`, then the horizontal
divisor gives a degree-five complete linear system on the old generic fiber.
After the vertical conditions, the new generic fiber can be represented as a
genus-one normal quintic in `P^4`, defined by the five `4 x 4` Pfaffians of a
`5 x 5` alternating matrix of linear forms.  Fisher gives a practical
algorithm for its invariants and hence its Jacobian
([The invariants of a genus one curve](https://arxiv.org/abs/math/0610318);
[elliptic normal quintic model](https://arxiv.org/abs/1110.3520)).

That supplies only the genus-one/Jacobian algebra.  A K3 5-neighbor routine
would still have to compute `H^0(X,O(F'))`, prove the existence and field of
definition of a section, construct the isomorphism to the Jacobian, and track
the surface maps and sections.  It is therefore substantially more work than
calling the degree-five invariant algorithm.

Most importantly, the two prime-5 moves in the current
`5,5,2,2,2,2,2` Kneser chain were introduced as a search heuristic to escape
local minima.  They are not currently certified geometric 5-neighbor steps.
A preferable route, if a suitable marked chain exists, is to search for a
2-/3-neighbor geometric path between the desired `U`-embeddings.  This is a
computational proposal, not a theorem that every present 5-edge can be
replaced.

## 5. Data missing specifically for `X(6,79)`

The public primary sources establish the moduli anchor but do not publish the
record-specific transport.  Elkies says that the rational non-CM point on

\[
 X(6,79)/\langle w_{474}\rangle
\]

produced the rank-17 elliptic K3, while the detailed examples in the 2008
paper are `N=6,14,57,206`, not `(6,79)`
([introduction](https://arxiv.org/pdf/0802.1301#page=2)).  The 2007 lectures
give the quotient curve and rational point but only describe the construction
architecture, not the equations and transformations
([arXiv:0709.2908](https://arxiv.org/abs/0709.2908)).

To transport the known point

\[
 (q,v)=\left(\frac{14}{13},\frac{16064}{2197}\right)
\]

to the `E8+A2^3` chart, one still needs the following exact package.

### Moduli-to-seed data

- A Weierstrass family over the function field of the genus-two Shimura
  quotient, or at minimum the specialized K3 over `Q` at the displayed point.
- A marked realization of the transparent seed with root lattice
  `A11+A3+A2`, including the rank-one section of height `79/12`.
- The actual Kodaira symbol and local component data for the ambiguous `A2`
  factor (`I3` versus `IV` cannot be read from the lattice alone).
- The rational functions from the Shimura coordinates to all Weierstrass
  coefficients, together with the Atkin--Lehner/polarization and quadratic
  twist choice.
- Exact nonvanishing checks for every denominator at the non-CM point.

The 2008 paper explains how such a missing moduli map was recovered in other
levels: parametrize a smaller lattice-polarized family, find a finite-field
seed, lift p-adically while fixing simple moduli functions, reconstruct
algebraic relations, and verify the full lattice embedding symbolically
([computational method](https://arxiv.org/pdf/0802.1301#page=9)).

### Seed-to-target data

- The complete positive-lattice chain with `B`, LLL, and common-sublattice
  matrices for every edge.
- More importantly, compatible primitive `U`-embeddings in one marked
  Neron--Severi lattice and effective representatives of each new fiber.
- A 2-/3-only geometric route, or a new rigorously implemented degree-five
  route for the first two edges.
- At every hop, the Riemann--Roch pencil and mutually inverse birational maps.
- Transported section formulas, not just Mordell--Weil ranks and regulators.
- A final coordinate comparison with the two-split `E8+A2^3` chart and exact
  recovery of the three target sections.

For the rank-31 objective, reconstructing only the single K3 at the non-CM
point would be enough to test the lattice bridge and obtain a fixed
rank-17-over-`Q(t)` surface.  Reconstructing the full Shimura family is needed
for a genuine moduli map and for systematic deformation to other Shimura
points, but not logically necessary for the first transported surface.

## 6. Recommended next certificate

The smallest useful new artifact is a **marked neighbor certificate**, with
one record per edge containing:

```text
parent_gram, child_gram, p, input_vector
kneser_basis_matrix_B, lll_matrix_U, total_matrix_M
common_index_p_sublattice_bases
stable_NS_isometry_Phi
old_F, old_O, new_F, new_O
new_F_component_decomposition
section_index_witness_D
riemann_roch_basis_and_local_conditions
new_parameter_u
genus_one_model, weierstrass_model
forward_map, inverse_map
transported_sections
bad_locus
```

The exact validation gates should include:

\[
 M^{\mathsf T}G_iM=G_{i+1},\qquad
 \Phi^{\mathsf T}(U\oplus-G_{i+1})\Phi=U\oplus-G_i,
\]

the fiber/section intersection identities, symbolic substitution in both
directions, minimal-discriminant and Kodaira checks, and exact height-pairing
preservation.

The immediate implementation order suggested by the sources is:

1. rerun the Kneser chain while retaining `return_matrix=True` and the LLL
   matrices;
2. search in the marked rank-19 Neron--Severi lattice for a 2-/3-neighbor
   path connecting the accessible seed fibration to the target `U`-embedding;
3. reconstruct the specialized `X(6,79)` seed K3 and its marked curves;
4. execute the Riemann--Roch pencils and verify the birational maps;
5. only then evaluate/compare the two-split chart and certify the third
   section.

This keeps the current lattice result as a strong routing signal while
preventing an abstract Kneser edge from being mistaken for an elliptic-surface
transformation.

## Source snapshot used for the auxiliary-code audit

The exact Sage 10.9 source file cited above had SHA-256
`2ee38ce4f398e816ed1634c7737bb97ae42e22337b7aeeafc1eae472cb3e24c8`.
At lines 233--247 it constructs row generators for the neighbor, returns
their transpose when `return_matrix=True`, and otherwise computes
`B * G * B.transpose()`.  Thus, for the returned matrix `T`, the checked
identity is `T.transpose() * G * T = G_new`.

The arXiv v3 source archive for Elkies--Kumar 1209.3527 was downloaded from
`https://export.arxiv.org/e-print/1209.3527` and had SHA-256
`82f3d657100f5d6f6abd79a001af0579439e8370511861d06e40199659e205bc`.
Relevant files and hashes were:

```text
README.txt   7d2583b06ecec2d7537c542835f8072bbf955768af9583e489ad9e2703600bda
12/12.txt    c02e7f3b3ba25153556e1be94cd148143eafb107328557cdd6ba9fad08690848
17/17.txt    1663d7012d02870ada73f7b6a14e9c64e4824dc743e6a4f822b14800df173b25
41/41.txt    466ecc56aea15c4ac0cac5bbe4e1e8eab0ea320d75fb1dc6467ff4b797f698f6
jacobians    b76e650e0f06f4ed668fc01190e21fa533bf04fd50583af324b9d0f432060593
```

This snapshot supports the statements about the actual `jac2`/`jac3`
implementation and local pole-cancellation calculations.  The deductions
about the missing `X(6,79)` package and the proposed marked certificate are
inferences from those primary methods and the current repository artifacts,
not claims made verbatim in the papers.
