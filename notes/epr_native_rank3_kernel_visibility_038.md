# Native rank-three kernel visibility 038

## Result

Audit pass:

    True

After Audit 036,

    E_ii = -1

and

    E_ij = -v S_ij.

Define

    K(v) = -E(v) = I + v C.

Audit 037 gives

    C^2 = 5 I

with two three-dimensional eigenspaces.

Therefore

    spec K(v)
      =
    1 + v sqrt(5)  multiplicity 3
    1 - v sqrt(5)  multiplicity 3.

Positivity requires

    |v| <= 1/sqrt(5).

Inside that interval K has rank six except at the endpoints.

A positive rank-three realization therefore forces

    v = +/- 1/sqrt(5).

The native B_plus3 analyzer sector has Gram kernel

    I + C/sqrt(5),

so it selects

    v = 1/sqrt(5).

For the frozen quartet,

    CHSH = 1 + 3/sqrt(5).

## Boundary

No Born frequency rule is used.

No numerical target fitting is used.

This audit does not prove that the finite receipt correlation is a
positive rank-three kernel.

That is now the remaining native weighting gate.
