# Unique binary weight closure 041

## Result

Audit pass:

    True

The native joint source supplies

    Gamma(A,B)
      =
    Tr[P_odd (A tensor B)].

The native one-wing signed means vanish.

For binary outcomes the four unknown weights are fixed by

    sum w_ab = 1

    sum a w_ab = 0

    sum b w_ab = 0

    sum ab w_ab = Gamma.

The coefficient matrix has rank

    4.

Therefore the solution is unique:

    w_ab
      =
    (1 + a b Gamma) / 4.

For the frozen quartet the resulting algebraic weights reproduce

    CHSH = 1 + 3/5*sqrt(5).

They are positive and normalized and have exact one-half local
marginals.

Audit035 is used only as an independent comparison after the weights
are derived.

## Boundary

This is an algebraic weight theorem.

It does not yet prove that repeated registered trial frequencies equal
these weights.

No Born frequency postulate is used as an input.

The remaining operational gate is:

    registered frequency
        ?=
    canonical algebraic source weight.
