#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import contextlib
import io
import runpy

HERE = Path(__file__).resolve().parents[2]

SOURCE = (
    HERE
    / "scripts/probes"
    / "probe_epr_pure_source_v4_gate_033c.py"
)

print("== 034A EXTREME SOURCE SHEET FACTORIZATION ==")

# Reuse the already-tested 033C reconstruction without copying its
# thousand-line construction into another probe.
sink = io.StringIO()

with contextlib.redirect_stdout(sink):
    ns = runpy.run_path(
        str(SOURCE)
    )

a = ns["a"]
g1800_reps = ns["g1800_reps"]
state_class = ns["state_class"]
pure_nonzero = set(
    ns["pure_nonzero"]
)
pure_zero = set(
    ns["pure_zero"]
)
X = ns["X"]
T = ns["T"]
XT = ns["XT"]

print()
print("PURE_NONZERO_COUNT:", len(pure_nonzero))
print("PURE_ZERO_COUNT:", len(pure_zero))

nz_reps = [
    g1800_reps[sid]
    for sid in sorted(pure_nonzero)
]

z_reps = [
    g1800_reps[sid]
    for sid in sorted(pure_zero)
]

nz_left = sorted({
    u
    for u, v in nz_reps
})

nz_right = sorted({
    v
    for u, v in nz_reps
})

z_left = sorted({
    u
    for u, v in z_reps
})

z_right = sorted({
    v
    for u, v in z_reps
})

print()
print("NZ_LEFT:", nz_left)
print("NZ_RIGHT:", nz_right)
print("ZERO_LEFT:", z_left)
print("ZERO_RIGHT:", z_right)

A = set(z_left)
B = set(nz_right)

checks = {}

checks[
    "zero_representatives_are_A_cross_A"
] = (
    set(z_left) == A
    and set(z_right) == A
    and {
        (u, v)
        for u in A
        for v in A
    }
    == set(z_reps)
)

checks[
    "nonzero_representatives_are_A_cross_B"
] = (
    set(nz_left) == A
    and set(nz_right) == B
    and {
        (u, v)
        for u in A
        for v in B
    }
    == set(nz_reps)
)

aA = {
    a[u]
    for u in A
}

aB = {
    a[u]
    for u in B
}

checks[
    "a_maps_A_to_B"
] = (
    aA == B
)

checks[
    "a_maps_B_to_A"
] = (
    aB == A
)

checks[
    "A_and_B_disjoint"
] = not (
    A & B
)

H = A | B

checks[
    "local_extreme_support_size_8"
] = (
    len(H) == 8
)

print()
print("A:", sorted(A))
print("B:", sorted(B))
print("a(A):", sorted(aA))
print("a(B):", sorted(aB))
print("H:", sorted(H))

# The diagonal-a quotient of H x H should contain exactly 32 states.
extreme_from_H = {
    state_class(u, v)
    for u in H
    for v in H
}

extreme_actual = (
    pure_nonzero
    | pure_zero
)

checks[
    "H_cross_H_diag_a_quotient_has_32_states"
] = (
    len(extreme_from_H) == 32
)

checks[
    "H_cross_H_quotient_equals_extreme_domain"
] = (
    extreme_from_H
    == extreme_actual
)

# Define the local a-sheet coordinate only on H.
#
# A = +1 sheet
# B = -1 sheet
#
# Global reversal of both labels changes no relative parity.
sheet = {
    u: (
        1
        if u in A
        else -1
    )
    for u in H
}

checks[
    "a_flips_local_sheet"
] = all(
    sheet[
        a[u]
    ]
    == -sheet[u]
    for u in H
)

def relative_sheet(sid):
    u, v = g1800_reps[sid]

    if (
        u not in H
        or v not in H
    ):
        # A quotient representative can be on either diagonal sheet.
        # Apply diagonal a once if required.
        u2 = a[u]
        v2 = a[v]

        if (
            u2 not in H
            or v2 not in H
        ):
            raise RuntimeError(
                "extreme state left H support"
            )

        u = u2
        v = v2

    return (
        sheet[u]
        * sheet[v]
    )

zero_parities = Counter(
    relative_sheet(sid)
    for sid in pure_zero
)

nonzero_parities = Counter(
    relative_sheet(sid)
    for sid in pure_nonzero
)

print()
print(
    "ZERO_RELATIVE_SHEET_PROFILE:",
    dict(sorted(zero_parities.items())),
)

print(
    "NONZERO_RELATIVE_SHEET_PROFILE:",
    dict(sorted(nonzero_parities.items())),
)

checks[
    "zero_is_same_sheet"
] = (
    zero_parities
    == Counter({
        1: 16,
    })
)

checks[
    "nonzero_is_opposite_sheet"
] = (
    nonzero_parities
    == Counter({
        -1: 16,
    })
)

# Verify the relative sheet coordinate is well defined on the diagonal-a
# quotient. Applying a to both factors must preserve it.
diag_a_failures = 0

for sid in extreme_actual:
    u, v = g1800_reps[sid]

    p1 = (
        sheet[u]
        * sheet[v]
    )

    p2 = (
        sheet[a[u]]
        * sheet[a[v]]
    )

    if p1 != p2:
        diag_a_failures += 1

checks[
    "relative_sheet_descends_through_diag_a"
] = (
    diag_a_failures == 0
)

# Test the source V4 action directly on this relative coordinate.
X_failures = 0
tau_failures = 0
Xtau_failures = 0

for sid in extreme_actual:
    e = relative_sheet(sid)

    if (
        relative_sheet(
            X[sid]
        )
        != e
    ):
        X_failures += 1

    if (
        relative_sheet(
            T[sid]
        )
        != -e
    ):
        tau_failures += 1

    if (
        relative_sheet(
            XT[sid]
        )
        != -e
    ):
        Xtau_failures += 1

checks[
    "X_preserves_relative_sheet"
] = (
    X_failures == 0
)

checks[
    "tau_flips_relative_sheet"
] = (
    tau_failures == 0
)

checks[
    "Xtau_flips_relative_sheet"
] = (
    Xtau_failures == 0
)

print()
print("X_RELATIVE_SHEET_FAILURES:", X_failures)
print("TAU_RELATIVE_SHEET_FAILURES:", tau_failures)
print(
    "XTAU_RELATIVE_SHEET_FAILURES:",
    Xtau_failures,
)

# Count quotient states by relative sheet directly from H x H.
quotient_sheet = {}

for u in H:
    for v in H:
        sid = state_class(
            u,
            v,
        )

        e = (
            sheet[u]
            * sheet[v]
        )

        if sid in quotient_sheet:
            assert (
                quotient_sheet[sid]
                == e
            )

        quotient_sheet[sid] = e

sheet_profile = Counter(
    quotient_sheet.values()
)

print()
print(
    "EXTREME_QUOTIENT_SHEET_PROFILE:",
    dict(
        sorted(
            sheet_profile.items()
        )
    ),
)

checks[
    "extreme_sheet_profile_16_16"
] = (
    sheet_profile
    == Counter({
        -1: 16,
        1: 16,
    })
)

failed = [
    key
    for key, value
    in checks.items()
    if not value
]

audit_pass = not failed

print()
print("== CHECKS ==")

for key, value in checks.items():
    print(
        "CHECK",
        key,
        "=",
        value,
    )

print()
print("AUDIT_PASS:", audit_pass)
print(
    "FAILED_CHECK_COUNT:",
    len(failed),
)

print(
    "FAILED_CHECKS:",
    failed,
)

if audit_pass:
    print(
        "VERDICT:",
        "extreme_G1800_source_domain_is_exactly_the_diag_a_quotient_"
        "of_an_8_state_two_sheet_local_support_and_zero_vs_nonzero_"
        "wedge_is_exactly_same_vs_opposite_relative_a_sheet_parity",
    )

    print()
    print(
        "KEEPER:",
        "The extreme source does not hide an extra bit. "
        "It retains exactly the relative a-sheet relation."
    )

    print()
    print(
        "NEXT_GATE:",
        "Test whether the already-earned antisymmetric preparation "
        "operator selects the opposite-sheet relation natively, "
        "without postselection or setting dependence."
    )

print()
print("BOUNDARY:")
print(
    "  This does not yet derive why a source prepares the "
    "opposite-sheet rather than same-sheet relation."
)
print(
    "  It identifies the missing preparation choice exactly."
)
