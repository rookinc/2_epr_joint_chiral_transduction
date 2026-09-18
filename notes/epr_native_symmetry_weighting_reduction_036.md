# Native symmetry weighting reduction 036

## Result

Audit pass:

    True

The native within-sector projective analyzer symmetry has order

    60.

The exact covariant one-wing marginal space has dimension

    0.

Thus native signed switching covariance forces

    m_i = 0

for all six settings.

For a factor-exchange-symmetric two-wing correlation table, the exact
covariant space has dimension

    2.

It is spanned by

    D = diagonal identity

and

    G = signed off-diagonal native Gram pattern.

Therefore

    E = alpha D + beta G.

The antisymmetric same-setting condition gives

    alpha = -1.

Writing

    v = -beta

gives the complete remaining family

    E_ii = -1

    E_ij = -v S_ij
        for i != j.

No bilinear correlation ansatz is used.

With zero local marginals, binary moment inversion gives

    P(a,b|x,y)
      =
    (1 + a b E_xy) / 4.

Thus there is no additional joint-weight freedom after v is fixed.

For the frozen Audit-035 quartet

    Alice = (0,1)
    Bob   = (0,2)

the CHSH family is

    S(v) = 1 + 3v

for v >= 0.

Bell violation occurs iff

    v > 1/3.

The conditional Hilbert control corresponds to

    v = 1/sqrt(5)

and reproduces

    S = 1 + 3/sqrt(5).

## Boundary

This audit does not derive v.

The common adjoint line system remains a conditional interface.

The antisymmetric same-axis closure is used.

No Born target value is used to obtain the one-parameter family.

The remaining native weighting problem is one-dimensional.

## Next gate

Derive v from finite native receipt mechanics on the h=1 heralded
source without using 1/sqrt(5) as selector input.
