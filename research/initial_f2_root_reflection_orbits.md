# Exact initial F2 line orbits under coordinate-root reflections

## Result

Each of the three initial Gram matrices has exactly 65,279 nonzero vectors
with `q(v)=v^t G v/2 = 0` modulo two. Since F2 has only one nonzero scalar,
these vectors are precisely the isotropic projective lines.

For each coordinate vector `e_i` of norm two we construct the integral
reflection `s_i(v)=v-(v,e_i)e_i`. The certificate verifies its matrix identity
`S_i^t G S_i=G` and `S_i^2=I` over the integers. It then exhausts all 131,071
nonzero binary vectors, and computes the complete orbits of the isotropic
subset under these generators, with closure, disjointness, and count checks.

| Initial lattice | Coordinate-root generators | Isotropic lines | Subgroup orbits | Orbits touched by the 32-line hash window |
| --- | ---: | ---: | ---: | ---: |
| Origin E8+A2^3 | 14 | 65,279 | 95 | 23 |
| Transparent A11+K6 | 16 | 65,279 | 34 | 10 |
| Frozen rootless | 0 | 65,279 | 65,279 | 32 |

This is the orbit partition under the **listed reflection subgroup**, not a
computation of the full integral automorphism group, nor a count of lattice
isometry classes. The rootless row has a trivial reflection subgroup; it
does not claim that its full automorphism group is trivial.

## Why this matters to the next experiment

The origin has many very small orbits: five of size one, six of size three,
and eighteen of size nine. A full-coordinate sample removes the earlier
coordinate concentration but remains weighted toward large orbits; its first
32 lines touch only 23 of the 95 origin orbits and 10 of the 34 transparent
orbits. Testing one representative per listed orbit is a concrete next
selection rule, with 129 representatives across the two rootful endpoints.

For the origin these reflections preserve coordinates 14, 15, and 16
(zero-based). The orbit counts in their eight parity classes are 11 for
`000` and 12 in every other class. The class `001` contains 8,128 lines;
the earlier initial lexicographic windows, whose last coordinate was one
and whose other two section coordinates were zero, stayed within this
class. This statement concerns the **initial fixed Gram basis** only. LLL
changes the basis at later search states, so it is not a global confinement
theorem for an entire neighbor path.

The orbit partition alone does not certify a complete list of even
2-neighbors. Two divides the determinant, 948; radical lines and the
non-maximal/even construction boundary must still be treated explicitly.
No representative neighbors, geometric K3 maps, or new rank-32 points have
been certified by this orbit audit.

## Reproduction

```sh
python3 research/projective_root_reflection_orbits.py \
  --compare certificates/initial_f2_root_reflection_orbits.json
python3 -m unittest discover -s tests \
  -p 'test_projective_root_reflection_orbits.py' -v
```

The certificate stores every generator matrix, all 95 and 34 representatives
and orbit sizes for the rootful lattices, a singleton rule for the rootless
case, hashes of the full enumerated sets and partitions, source hashes, and
the exact sampler-window comparison. Its record SHA-256 is
`86add90a1d2f7c80b3f786cbc03a2e265e82e748a7441ec7403f208aca885233`.
All six dedicated tests passed in 2.371 seconds, including a separate
integer-quadratic evaluation and matrix-action check on toy examples.
