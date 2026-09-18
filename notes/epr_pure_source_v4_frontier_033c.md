# Pure-source V4 frontier 033C

The complete Audit-019 branch anatomy was reconstructed directly from
the sealed G1800 and decorated-face interfaces.

Each G1800 state has exactly 16 quotient-level decorated presentations.

The exact branch profile is:

    (zero, nonzero) -> G1800 state count

    (0,16)  -> 16
    (4,12)  -> 112
    (6,10)  -> 196
    (8,8)   -> 1152
    (10,6)  -> 196
    (12,4)  -> 112
    (16,0)  -> 16

The reconstructed profile matches Audit 019 exactly.

Define:

    Omega_NZ
        = the 16 G1800 states with profile (0,16)

    Omega_0
        = the 16 G1800 states with profile (16,0).

Omega_NZ is setting-independent and defined by an intrinsic
pre-setting wedge property.

However, Omega_NZ is not invariant under the full native source V4.

For the native source actions:

    identity:
        Omega_NZ -> Omega_NZ
        Omega_0  -> Omega_0

    X:
        Omega_NZ -> Omega_NZ
        Omega_0  -> Omega_0

    tau:
        Omega_NZ <-> Omega_0

    Xtau:
        Omega_NZ <-> Omega_0.

Thus the native object is the 32-state extreme doublet

    Omega_ext = Omega_NZ union Omega_0,

with a C2 grading epsilon satisfying

    epsilon(X s)   =  epsilon(s)
    epsilon(tau s) = -epsilon(s).

The pure-nonzero 16-state subset is therefore not promoted as a native
source by itself.

Selecting Omega_NZ would require an additional setting-independent
source orientation or registration law.

No probability weighting is introduced to favor Omega_NZ.

Jackalope 2 supplies a structurally related negative-a carrier and an
orientation-sensitive receiver program, but the prior crosswalk does
not derive the missing EPR preparation selector.

Next target:

    test whether an already-native registered G1800 orientation datum
    acts as the epsilon grading on Omega_ext.

If no such same-domain native registration exists, leave source
preparation open rather than selecting Omega_NZ by hand.
