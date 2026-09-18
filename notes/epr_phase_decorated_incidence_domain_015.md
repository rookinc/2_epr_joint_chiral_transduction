# Phase-decorated incidence domain 015

## Result

Audit pass:

    True

The 24 registered seam rows contain 48 event sides.

Every event side matches one of the 24 signed face blocks exactly.

Every signed block occurs exactly twice:

    once in registration 0
    once in registration 1.

Each signed block also has one unique event-anchor label.

Thus each face carries exactly four native mode objects:

    two positive
    two negative.

Registration 0 pairs positive and negative modes with the same anchor
label.

Registration 1 swaps the two anchor labels inside each face.

## Domain correction

The bare boundary-incidence domain has

    120

rows.

Every one of those 120 incidences is compatible with both registered
seam sectors.

Program 01 independently supplies two local relative-complex-orientation
states per face:

    6 faces x 2 phase states = 12 local states.

Therefore a bare incidence cannot be used to manufacture the phase bit.

The correct preparation domain is

    I_F_tilde
      =
    {(s,u,B):
      s in S_12,
      face(s)=carrier(B),
      u in B}.

Each face has

    4 modes x 5 states = 20

bare incidences and two phase states, giving

    6 x 20 x 2 = 240

phase-decorated incidences.

## Next target

Recover the explicit twelve local phase-state rows from the current
201FK face-field checkpoint.

Then construct I_F_tilde literally and derive the CP1 projective line
map ell.

No registration bit or signed-block sign is being silently identified
with the independent complex-orientation phase.
