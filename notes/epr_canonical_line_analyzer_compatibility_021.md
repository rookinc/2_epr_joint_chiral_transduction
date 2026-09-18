# Canonical EPR line versus analyzer capacity 021

## Result

Audit pass:

    True

A distinct native analyzer pair has the exact qubit normal form

    Z = [[1,0],[0,-1]]

    R = [[-3/5,4/5],
         [ 4/5,3/5]].

The corresponding rank-one positive projectors have overlap

    1/5.

The CHSH operator

    B = Z tensor (Z+R)
        + R tensor (Z-R)

has squared spectral values

    36/25
    164/25

and therefore operator norm

    2 sqrt(41) / 5.

## Canonical preparation test

The Program-02 nonzero preparation branch is the antisymmetric line

    |01> - |10>.

In the common pair normal form it satisfies exactly

    B (|01>-|10>)
      =
    (6/5) (|01>-|10>).

So the canonical prepared line does not attain the Bell-capable
operator norm.

Its CHSH value in this frame is

    6/5.

## Outcome-label stress

All sixteen independent local binary outcome-sign relabelings were
tested.

The absolute CHSH values have profile

    |S| = 6/5 : 8 cases
    |S| = 2   : 8 cases.

Thus local outcome relabeling alone never gives

    |S| > 2.

## Consequence

The inherited analyzer geometry is Bell-capable as an operator family,
but that fact does not yet imply Bell violation for the canonical
source state we actually derived.

The missing object is now sharper:

    a native relative Alice/Bob analyzer-frame law
    or a different native admissible setting quartet.

An arbitrary relative frame may not be chosen to optimize CHSH.

The coherent local instrument remains available, but instrument
integration is paused until this state/operator compatibility question
is resolved.
