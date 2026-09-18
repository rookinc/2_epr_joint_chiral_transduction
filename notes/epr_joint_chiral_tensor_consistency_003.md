# Joint chiral tensor consistency 003

## Result

Audit pass:

    True

Program 02 already established that the native G1800 relative-lift
involution built from deck element a has local Program-01 face action

    S1 = sigma_x K.

The present audit asks whether the two point-level presentations

    [a u, v]
    [u, a v]

automatically become the same operator on the full two-face Hilbert
tensor product.

They do not.

## Two chiral tensor sheets

Let V be the four-real-dimensional Program-01 face doublet with complex
structure J.

The same-orientation tensor sheet is

    Q_plus =
      (V_A tensor_R V_B)
      /
      <J_A u tensor v - u tensor J_B v>.

The opposite-orientation sheet is

    Q_minus =
      (V_A tensor_R V_B)
      /
      <J_A u tensor v + u tensor J_B v>.

Both have real dimension 8, hence complex dimension 4.

A one-sided S1 exchanges these sheets.

## Local presentation comparison

The Alice-side map

    S1_A tensor I

and Bob-side map

    I tensor S1_B

both map Q_plus to Q_minus.

However, they are not identical on the full quotient.

Their induced difference has real rank

    4.

The states on which the two local presentations agree form a kernel of
real dimension

    4,

It is not a complex subspace.

The two-sided anti-complex involution

    C_AB = S1_A tensor S1_B

has this compatibility space as its exact fixed real form.

Because C_AB anticommutes with the joint complex structure,
multiplication by i sends the fixed real form to the anti-fixed
real form.

Projectively, the compatible rays form

    RP^3 inside CP^3.

## Meaning

Audit 002 established a point-level joint chiral torsor.

Audit 003 shows that Hilbert coherence is an additional condition.

The native equality

    [a u,v] = [u,a v]

does not by itself license the operator identity

    S1_A tensor I = I tensor S1_B

on every joint Hilbert state.

Instead it exposes a smaller compatibility space on which the two local
presentations can agree.

That space is now the next native preparation target.

## Boundary

No unique joint projective line is selected here.

No singlet has been inserted.

No probability or frequency law is assumed.

No Bell violation or no-signaling claim is made.
