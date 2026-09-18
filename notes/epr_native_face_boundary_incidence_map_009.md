# Native face boundary incidence map 009

## Result

Audit pass:

    True

The correct Program-01 local boundary domain is not bare G60.

It is the incidence set

    I_F = {(u,B): u belongs to signed face block B}.

There are

    120

such incidences.

Each of the sixty G60 states occurs exactly twice:

    once in a positive signed block;
    once in a negative signed block.

Therefore a G60 state alone does not specify one local face vector.

## Local face module

Each of the six face carriers has exactly four signed blocks:

    two positive;
    two negative.

These four blocks are the canonical basis objects of the four-real-
dimensional Program-01 face permutation module.

For carrier c define

    V_c = R[B_c].

The native boundary presentation is

    beta_c(u,B) = e_B.

No arbitrary within-face coordinate numbering is required.

## Chiral covariance

The native deck involution

    a = Aut(G60)[326]

closes exactly on the signed incidence domain.

On each face it acts as two transpositions:

    one on the positive pair;
    one on the negative pair.

This is the basis-free signed-block action identified in Program 01 as

    S1 = sigma_x K.

Therefore all 120 incidences satisfy

    beta(a u,a B)
      =
    S1 beta(u,B).

## State-only ambiguity

States whose two signed incidences use the same carrier:

    0

States whose two signed incidences use different carriers:

    60

Either way, the signed boundary incidence is required to specify the
basis object.

## Next gate

Construct the two-wing incidence presentation on

    G1800 = (G60 x G60)/diag(a).

The diagonal quotient must descend through the local S1 actions.

If the construction is correct, the chiral compatibility structure found
abstractly in Audit 003 should reappear from this native incidence map
rather than being imposed afterward.
