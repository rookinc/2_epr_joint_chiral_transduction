# 2 - EPR Joint Chiral Transduction

## Purpose

Construct one native EPR experiment in which a single setting-independent joint preparation encounters two separately chosen local transducers and produces a complete nonsignalling Bell-violating receipt distribution.

The source should correlate relational or chiral structure, not preassign classical local answers.

## Starting point

Already established upstream:

- Native six-axis analyzer geometry exists.
- Its signed conference matrix satisfies

    C^2 = 5 I.

- The selected three-dimensional analyzer sector has Gram kernel

    G_plus = I + C/sqrt(5).

- Native projective switching symmetry forces zero one-wing signed means.
- A pre-setting source receipt isolates the heralded h=1 EPR branch.
- The canonical joint preparation is the exchange-odd line with projector

    P_odd = (I-X)/2.

- The same projector has the exact quaternionic form

    P_odd =
      (I + K1 tensor K1
         + K2 tensor K2
         + K3 tensor K3) / 4.

- The canonical joint contraction is

    Gamma(A,B) = Tr[P_odd (A tensor B)].

- On traceless local observables,

    Gamma(A,B) = -1/2 Tr(AB).

- Native analyzer geometry therefore fixes the distinct-setting visibility to

    v = 1/sqrt(5).

- For binary outcomes the joint algebraic weights are uniquely determined:

    w_ab = (1 + a b Gamma) / 4.

- The frozen native quartet gives

    CHSH = 1 + 3/sqrt(5) > 2.

- No setting-independent positive mixture of preexisting local binary answer
  tables can reproduce this joint table.

## Algebraic EPR assembly lock

Program 02 has closed its algebraic EPR assembly.

The native chain is

    common-domain preparation
      -> pre-setting herald h
      -> h=1 canonical antisymmetric preparation
      -> P_odd=(I-X)/2
      -> Gamma(A,B)=Tr[P_odd(A tensor B)]
      -> Gamma(A,B)=-1/2 Tr(AB)
      -> native six-axis Gram geometry
      -> v=1/sqrt(5)
      -> zero local signed means
      -> unique binary algebraic weights
      -> exact half marginals
      -> CHSH=1+3/sqrt(5).

The final assembly certificate is

    epr_algebraic_assembly_lock_042.

## Native joint preparation

The herald is fixed before analyzer settings are chosen.

The h=1 branch carries the canonical exchange-odd joint line. Its projector is

    P_odd = (I-X)/2.

This is exactly equal to

    P_odd =
      (I + K1 tensor K1
         + K2 tensor K2
         + K3 tensor K3) / 4.

The joint object is therefore intrinsic to the common preparation geometry and
is not assembled from separate local answer tables after settings are chosen.

## Joint contraction law

The canonical projector supplies the exact bilinear joint operator contraction

    Gamma(A,B)
      =
    Tr[P_odd (A tensor B)].

For arbitrary 2x2 local matrices,

    Gamma(A,B)
      =
    1/2 [Tr(A)Tr(B) - Tr(AB)].

On traceless local observables this reduces to

    Gamma(A,B)
      =
    -1/2 Tr(AB).

For normalized Pauli directions this is the negative Euclidean inner product.

## Native analyzer geometry

The six native analyzer lines are governed by a signed conference matrix C with

    C^2 = 5 I.

Its two eigenspaces have dimension three.

For the selected B_plus3 sector,

    G_plus = I + C/sqrt(5).

Hence the native joint contraction gives

    Gamma_ii = -1

and, for distinct settings,

    Gamma_ij = -S_ij/sqrt(5).

Thus the native visibility is

    v = 1/sqrt(5).

No visibility parameter is fitted.

## Binary algebraic measure

Native projective switching symmetry forces the one-wing signed means to zero.

For binary outcomes a,b in {+1,-1}, normalization, the two zero local means,
and the joint product mean Gamma form a full-rank four-equation moment system.

Its unique solution is

    w_ab
      =
    (1 + a b Gamma) / 4.

For the frozen quartet,

    E_00 = -1
    E_01 = -1/sqrt(5)
    E_10 = -1/sqrt(5)
    E_11 = +1/sqrt(5).

The corresponding weights are positive and normalized.

Each local marginal is exactly one half, so the algebraic table is
nonsignalling.

The resulting Bell value is

    CHSH
      =
    1 + 3/sqrt(5)
      =
    2.341640786499874...

which exceeds the Bell-local bound 2.

## Bell-local obstruction

Every deterministic local binary answer table

    (A0,A1,B0,B1)

has

    |CHSH| = 2.

Any setting-independent positive mixture of such tables therefore obeys

    |CHSH| <= 2.

The native algebraic table lies outside that local polytope.

The missing structure is therefore not a hidden census of preexisting local
answers.

## What is closed

Program 02 has closed:

- the common-domain heralded preparation interface,
- the canonical antisymmetric joint object,
- the native joint operator contraction,
- the native 1/sqrt(5) visibility,
- zero one-wing signed means,
- the unique positive normalized binary algebraic measure,
- exact nonsignalling marginals,
- and an exact CHSH value above 2.

No target probability table is inserted as a derivation input.

No Bell-local hidden answer weighting reproduces the result.

## Remaining interfaces

Two interfaces remain explicit.

### Common-adjoint face binding

The native analyzer sector is connected to the Program-01 face adjoint through
the previously declared conditional common-adjoint interface.

Program 02 does not reopen that calibration problem.

### Operational frequency correspondence

The finite construction derives a unique positive normalized algebraic measure

    w_ab.

It does not yet identify that measure with a microscopic repeated-trial
counting law.

This is an operational correspondence problem rather than unresolved freedom
in the EPR table.

## Program boundary

No further internal weighting hunt is required in Program 02.

In particular, do not return to:

- source-state counting as a replacement for the joint measure,
- preexisting local binary answer tables,
- coherent-pulse optimization,
- visibility fitting,
- or further parity labels.

Future work belongs in a separate operational-correspondence or
instrument-realization program.

## Final status

    PROGRAM 02 ALGEBRAIC EPR ASSEMBLY: CLOSED

Canonical result:

    P_odd
      -> Gamma
      -> w_ab
      -> CHSH = 1 + 3/sqrt(5) > 2.

The remaining work is correspondence, not reconstruction of the algebraic EPR
table.
