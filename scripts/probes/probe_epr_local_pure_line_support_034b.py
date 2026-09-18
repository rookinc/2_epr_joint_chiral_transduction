#!/usr/bin/env python3

from collections import Counter
from pathlib import Path
import contextlib
import io
import runpy

HERE = Path(__file__).resolve().parents[2]

SOURCE = (
    HERE
    / "scripts/probes"
    / "probe_epr_pure_source_v4_gate_033c.py"
)

print("== 034B LOCAL PURE-LINE SUPPORT GATE ==")

sink = io.StringIO()

with contextlib.redirect_stdout(sink):
    ns = runpy.run_path(
        str(SOURCE)
    )

a = ns["a"]
g1800_reps = ns["g1800_reps"]
state_class = ns["state_class"]
decorated_rows = ns["decorated_rows"]
row_line = ns["row_line"]
rows_by_g60 = ns["rows_by_g60"]
pure_nonzero = set(
    ns["pure_nonzero"]
)
pure_zero = set(
    ns["pure_zero"]
)

checks = {}

# ------------------------------------------------------------
# Compute the four-row projective-line profile of every G60 state.
# ------------------------------------------------------------

line_profile = {}

for u in range(60):
    counts = Counter(
        row_line(row)
        for row in rows_by_g60[u]
    )

    line_profile[u] = (
        int(counts[0]),
        int(counts[1]),
    )

profile_hist = Counter(
    line_profile.values()
)

print()
print("LOCAL_LINE_PROFILE_HISTOGRAM:")

for profile, count in sorted(
    profile_hist.items()
):
    print(
        " ",
        profile,
        "->",
        count,
    )

checks[
    "four_decorated_rows_per_G60_state"
] = all(
    sum(profile) == 4
    for profile in line_profile.values()
)


# ------------------------------------------------------------
# Define pure-line support independently of the EPR branch.
# ------------------------------------------------------------

pure_line_0 = {
    u
    for u, profile
    in line_profile.items()
    if profile == (
        4,
        0,
    )
}

pure_line_1 = {
    u
    for u, profile
    in line_profile.items()
    if profile == (
        0,
        4,
    )
}

H_pure = (
    pure_line_0
    | pure_line_1
)

print()
print(
    "PURE_LINE_0_COUNT:",
    len(pure_line_0),
)

print(
    "PURE_LINE_1_COUNT:",
    len(pure_line_1),
)

print(
    "PURE_LINE_SUPPORT_COUNT:",
    len(H_pure),
)

print(
    "PURE_LINE_0_STATES:",
    sorted(pure_line_0),
)

print(
    "PURE_LINE_1_STATES:",
    sorted(pure_line_1),
)

print(
    "H_PURE:",
    sorted(H_pure),
)

checks[
    "pure_line_support_has_8_states"
] = (
    len(H_pure) == 8
)

checks[
    "pure_line_halves_have_4_states_each"
] = (
    len(pure_line_0) == 4
    and len(pure_line_1) == 4
)


# ------------------------------------------------------------
# Recover H independently from the already-observed extreme G1800
# domain, then compare.
# ------------------------------------------------------------

extreme = (
    pure_nonzero
    | pure_zero
)

H_extreme = set()

for sid in extreme:
    u, v = g1800_reps[sid]

    H_extreme.add(u)
    H_extreme.add(v)

    H_extreme.add(a[u])
    H_extreme.add(a[v])

print()
print(
    "H_FROM_EXTREME_DOMAIN:",
    sorted(H_extreme),
)

checks[
    "independent_pure_line_support_equals_extreme_support"
] = (
    H_pure == H_extreme
)


# ------------------------------------------------------------
# Test native a action on the independently-defined local support.
# ------------------------------------------------------------

checks[
    "a_preserves_pure_line_support"
] = all(
    a[u] in H_pure
    for u in H_pure
)

checks[
    "a_exchanges_the_two_pure_line_halves"
] = (
    {
        a[u]
        for u in pure_line_0
    }
    == pure_line_1
    and {
        a[u]
        for u in pure_line_1
    }
    == pure_line_0
)

a_line_flip_failures = 0

for u in H_pure:
    p = line_profile[u]
    q = line_profile[
        a[u]
    ]

    if p == (4, 0):
        expected = (
            0,
            4,
        )
    elif p == (0, 4):
        expected = (
            4,
            0,
        )
    else:
        raise RuntimeError(
            "non-pure state entered H_pure"
        )

    if q != expected:
        a_line_flip_failures += 1

print()
print(
    "A_PURE_LINE_FLIP_FAILURES:",
    a_line_flip_failures,
)

checks[
    "a_flips_pure_projective_line"
] = (
    a_line_flip_failures == 0
)


# ------------------------------------------------------------
# Reconstruct the G1800 quotient from H_pure x H_pure only.
# ------------------------------------------------------------

pure_product_states = {
    state_class(
        u,
        v,
    )
    for u in H_pure
    for v in H_pure
}

print()
print(
    "PURE_LINE_PRODUCT_QUOTIENT_COUNT:",
    len(pure_product_states),
)

checks[
    "pure_line_product_quotient_has_32_states"
] = (
    len(pure_product_states)
    == 32
)

checks[
    "pure_line_product_quotient_equals_extreme_domain"
] = (
    pure_product_states
    == extreme
)


# ------------------------------------------------------------
# Define local projective-line sign. Its overall sign is gauge.
# The relative product is not.
# ------------------------------------------------------------

line_sign = {}

for u in H_pure:
    if u in pure_line_0:
        line_sign[u] = 1
    else:
        line_sign[u] = -1


def relative_line(sid):
    u, v = g1800_reps[sid]

    # Representatives of states in the pure product remain in H_pure,
    # but retain the diagonal-a fallback explicitly.
    if (
        u not in H_pure
        or v not in H_pure
    ):
        u = a[u]
        v = a[v]

    return (
        line_sign[u]
        * line_sign[v]
    )


zero_relative_profile = Counter(
    relative_line(sid)
    for sid in pure_zero
)

nonzero_relative_profile = Counter(
    relative_line(sid)
    for sid in pure_nonzero
)

print()
print(
    "ZERO_RELATIVE_LINE_PROFILE:",
    dict(
        sorted(
            zero_relative_profile.items()
        )
    ),
)

print(
    "NONZERO_RELATIVE_LINE_PROFILE:",
    dict(
        sorted(
            nonzero_relative_profile.items()
        )
    ),
)

checks[
    "zero_branch_is_same_pure_line"
] = (
    zero_relative_profile
    == Counter({
        1: 16,
    })
)

checks[
    "nonzero_branch_is_opposite_pure_line"
] = (
    nonzero_relative_profile
    == Counter({
        -1: 16,
    })
)


# ------------------------------------------------------------
# Global S2 relabeling of line 0 <-> line 1 changes every local
# sign but leaves the relative product unchanged.
# ------------------------------------------------------------

global_line_swap_failures = 0

for sid in extreme:
    u, v = g1800_reps[sid]

    if (
        u not in H_pure
        or v not in H_pure
    ):
        u = a[u]
        v = a[v]

    original = (
        line_sign[u]
        * line_sign[v]
    )

    swapped = (
        (-line_sign[u])
        * (-line_sign[v])
    )

    if original != swapped:
        global_line_swap_failures += 1

print()
print(
    "GLOBAL_LINE_SWAP_FAILURES:",
    global_line_swap_failures,
)

checks[
    "relative_line_relation_is_global_S2_gauge_invariant"
] = (
    global_line_swap_failures
    == 0
)


# ------------------------------------------------------------
# Classification.
# ------------------------------------------------------------

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
print(
    "AUDIT_PASS:",
    audit_pass,
)

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
        "the_8_state_local_support_underlying_the_extreme_G1800_"
        "domain_is_independently_characterized_as_the_complete_"
        "pure_projective_line_support_and_native_a_exchanges_its_"
        "two_four_state_line_halves"
    )

    print()
    print(
        "KEEPER:",
        "The extreme EPR domain is not defined by its favorable "
        "joint outcome. It is the diagonal-a product of the "
        "independently recognizable local pure-line support."
    )

    print()
    print(
        "NEXT_GATE:",
        "Determine whether the source construction natively "
        "prepares opposite pure-line relation, or instead yields "
        "a setting-independent null/EPR preparation receipt."
    )

print()
print("BOUNDARY:")
print(
    "  Pure-line support does not by itself choose the opposite-line "
    "joint relation."
)
print(
    "  No probability law, Born rule, or analyzer setting enters "
    "this support definition."
)
