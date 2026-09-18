# Literal phase-decorated incidence table 016

## Result

Audit pass:

    True

Program-01 201FK gives two local relative-complex-orientation states on
each of six native face objects.

Serialize them as

    (face, epsilon)

with

    epsilon = +1 or -1

and

    J_F(face,epsilon)
      =
    J_S direct sum epsilon J_E.

This gives exactly

    6 x 2 = 12

local phase states.

The global exchange of the names +1 and -1 is only a presentation
convention.

## Native event modes

Audit 015 already closed four native signed event-mode objects on each
face:

    two positive
    two negative.

Each mode is identified by its exact signed face block and unique
registered-event anchor label.

The integer local mode slot written into the JSON is only a deterministic
serialization label.

## Decorated domain

The literal native preparation domain is

    I_F_tilde
      =
    {(s,u,B):
      s in S_12,
      face(s)=carrier(B),
      u in B}.

It has exactly

    240

rows.

Profiles:

    20 rows per phase state
    40 rows per face
    48 phase-mode pairs
    5 G60 states per phase-mode pair.

Forgetting the phase maps I_F_tilde exactly two-to-one onto the 120-row
bare boundary-incidence domain.

## Closed interface

The corrected eta map is now fully serialized:

    eta_tilde(s,u,B) = (s,m_B).

The remaining local preparation problem is ell:

    ell:
      S_12 x M_4
        ->
      CP1.

Once ell is native, the two-wing preparation census can classify
projective-line equality and therefore zero versus nonzero
antisymmetric wedge before analyzer settings exist.
