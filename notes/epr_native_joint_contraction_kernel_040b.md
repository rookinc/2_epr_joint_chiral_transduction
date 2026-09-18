# Native joint contraction kernel 040B

## Result

Audit pass:

    True

The canonical native source projector is

    P_odd = (I-X)/2

and exactly

    P_odd =
      (I + K1 tensor K1
         + K2 tensor K2
         + K3 tensor K3) / 4.

Define the joint operator contraction

    Gamma(A,B)
      =
    Tr[P_odd (A tensor B)].

Exact evaluation on the full 2x2 matrix-unit basis gives

    Gamma(A,B)
      =
    1/2 [Tr(A)Tr(B) - Tr(AB)].

For traceless local observables,

    Gamma(A,B)
      =
    -1/2 Tr(AB).

For normalized Pauli directions this is

    Gamma(n.sigma,m.sigma)
      =
    - n dot m.

The one-wing contractions vanish:

    Gamma(A,I) = 0
    Gamma(I,B) = 0

for traceless A and B.

Audit 037 supplies the native analyzer Gram kernel

    G_plus = I + C/sqrt(5).

Therefore, under the already-declared common-adjoint analyzer interface,

    Gamma_ii = -1

and for distinct settings

    Gamma_ij = -S_ij/sqrt(5).

Thus the canonical source operator itself supplies

    v = 1/sqrt(5)

at the joint operator-kernel level.

## Boundary

This closes native joint operator bilinearity.

It does not identify the operator contraction with empirical trial
frequency.

It does not derive a Born frequency rule.

The common face-adjoint binding remains conditional.

The remaining gate is operational:

    registered binary correlation
        ?=
    canonical joint contraction Gamma.
