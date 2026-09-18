# Native six-axis singlet quartet census 022

## Result

Audit pass:

    True

Audit 021 tested one restricted quartet: the same analyzer pair on both
wings.

This audit keeps the canonical antisymmetric EPR line fixed and exhausts
the complete six-axis native analyzer geometry.

Each native analyzer character sector has six projective lines with

    |n_i dot n_j|^2 = 1/5

for distinct axes.

Under a common adjoint lift into the face su(2), the singlet correlation
is

    E(i,j) = - n_i dot n_j.

A common adjoint-frame rotation changes no Gram entries and therefore no
value in this census.

## Exhaustive native setting census

Alice chooses two distinct axes.

Bob chooses two distinct axes.

Thus there are

    6*5*6*5 = 900

ordered setting quartets per native character sector.

For each quartet all

    16

local binary outcome-sign relabelings are exhausted.

Both character sectors give the same exact profile.

Bell-capable quartets per sector:

    240 / 900.

Exact maximum:

    1 + 3/sqrt(5)

which is greater than

    2.

Every violating quartet has exactly one native analyzer axis shared
between Alice's two-setting set and Bob's two-setting set.

No arbitrary continuous Alice/Bob frame rotation is optimized.

## Meaning

The Audit-021 negative result was not a state/analyzer incompatibility.

It was a restricted-setting result.

The full six-axis native analyzer geometry contains setting quartets
compatible with Bell violation by the canonical antisymmetric EPR line,
provided Alice and Bob inherit the required common adjoint lift.

That relative-frame admission is now the next theorem gate.

The coherent local instrument and frequency law remain downstream.
