# Analyzer native character identification 026

## Result

Audit pass:

    True

All 480 native Aut(G60) elements were evaluated on one common domain.

Analyzer character:

    parity of the six-face permutation

Deck character:

    preserve / exchange b and ab

Face character:

    preserve / reverse native signed-face sign

Match counts:

    chi_deck                 480
    chi_face                 240
    chi_deck * chi_face      240

Exact native match:

    chi_deck

## Local phase consequence

Program-01 FR already proves

    delta_f
      =
    chi_deck restricted to Stab(f)

for every native face.

Therefore the local interpretation depends on the exact character
identified above:

    analyzer sector exchange equals chi_deck globally; by FR its restriction to every face stabilizer equals the local relative-complex-orientation character delta_f

## Boundary

This is a finite native character theorem.

It does not identify the binary character with physical spatial parity,
weak-interaction chirality, or a physical detector response law.

The next gate is apparatus construction on the complex2 face carrier.
