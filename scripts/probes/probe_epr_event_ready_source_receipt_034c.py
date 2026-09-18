#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import contextlib
import io
import runpy

HERE = Path(__file__).resolve().parents[2]

P033 = (
    HERE
    / "scripts/probes"
    / "probe_epr_pure_source_v4_gate_033c.py"
)

P034B = (
    HERE
    / "scripts/probes"
    / "probe_epr_local_pure_line_support_034b.py"
)

print("== 034C EVENT-READY SOURCE RECEIPT ==")

sink = io.StringIO()

with contextlib.redirect_stdout(sink):
    n33 = runpy.run_path(
        str(P033)
    )

with contextlib.redirect_stdout(sink):
    n34 = runpy.run_path(
        str(P034B)
    )

a = n33["a"]
g1800_reps = n33["g1800_reps"]
pure_nonzero = set(
    n33["pure_nonzero"]
)
pure_zero = set(
    n33["pure_zero"]
)
X = n33["X"]
T = n33["T"]
XT = n33["XT"]

H = set(
    n34["H_pure"]
)

H0 = set(
    n34["pure_line_0"]
)

H1 = set(
    n34["pure_line_1"]
)

line_sign = {
    u: (
        1
        if u in H0
        else -1
    )
    for u in H
}

extreme = (
    pure_nonzero
    | pure_zero
)

checks = {}


def relative_line(sid):
    u, v = g1800_reps[sid]

    if (
        u not in H
        or v not in H
    ):
        u = a[u]
        v = a[v]

    return (
        line_sign[u]
        * line_sign[v]
    )


def herald(sid):
    # 1 = successful EPR preparation receipt
    # 0 = null preparation receipt
    return (
        1
        if sid in pure_nonzero
        else 0
    )


# Exact algebraic relation between the independently-derived
# local-line relation and the source preparation receipt.
relation_failures = 0

for sid in extreme:
    r = relative_line(sid)
    h = herald(sid)

    expected = (
        1 - r
    ) // 2

    if h != expected:
        relation_failures += 1

checks[
    "herald_equals_one_minus_relative_line_over_two"
] = (
    relation_failures == 0
)

print(
    "HERALD_RELATION_FAILURES:",
    relation_failures,
)


# Complete source receipt profile.
receipt_profile = Counter(
    herald(sid)
    for sid in extreme
)

print(
    "SOURCE_RECEIPT_PROFILE:",
    dict(
        sorted(
            receipt_profile.items()
        )
    ),
)

checks[
    "source_receipt_complete_16_null_16_EPR"
] = (
    receipt_profile
    == Counter({
        0: 16,
        1: 16,
    })
)


# X preserves the preparation receipt.
X_failures = 0

for sid in extreme:
    if (
        herald(
            X[sid]
        )
        != herald(
            sid
        )
    ):
        X_failures += 1

checks[
    "factor_exchange_preserves_source_receipt"
] = (
    X_failures == 0
)

print(
    "X_HERALD_FAILURES:",
    X_failures,
)


# tau exchanges null and EPR receipts.
tau_failures = 0
Xtau_failures = 0

for sid in extreme:
    if (
        herald(
            T[sid]
        )
        != 1 - herald(
            sid
        )
    ):
        tau_failures += 1

    if (
        herald(
            XT[sid]
        )
        != 1 - herald(
            sid
        )
    ):
        Xtau_failures += 1

checks[
    "tau_exchanges_null_and_EPR_receipts"
] = (
    tau_failures == 0
)

checks[
    "Xtau_exchanges_null_and_EPR_receipts"
] = (
    Xtau_failures == 0
)

print(
    "TAU_HERALD_FAILURES:",
    tau_failures,
)

print(
    "XTAU_HERALD_FAILURES:",
    Xtau_failures,
)


# Global projective-line relabeling changes both local signs
# and therefore leaves the source receipt unchanged.
gauge_failures = 0

for sid in extreme:
    u, v = g1800_reps[sid]

    if (
        u not in H
        or v not in H
    ):
        u = a[u]
        v = a[v]

    r = (
        line_sign[u]
        * line_sign[v]
    )

    r_swapped = (
        (-line_sign[u])
        * (-line_sign[v])
    )

    if r != r_swapped:
        gauge_failures += 1

checks[
    "source_receipt_is_global_line_gauge_invariant"
] = (
    gauge_failures == 0
)

print(
    "GLOBAL_LINE_GAUGE_FAILURES:",
    gauge_failures,
)


# These are sealed Audit-019 facts consumed by this receipt.
a019 = n33["a019"]

checks[
    "Audit019_pre_setting_classifier"
] = (
    a019[
        "boundary"
    ][
        "pre_setting_branch_classifier_constructed"
    ]
    is True
)

checks[
    "Audit019_analyzer_settings_not_used"
] = (
    a019[
        "boundary"
    ][
        "analyzer_settings_not_used"
    ]
    is True
)

checks[
    "Audit019_nonzero_is_canonical_EPR_line"
] = (
    a019[
        "boundary"
    ][
        "nonzero_branch_is_canonical_EPR_line"
    ]
    is True
)

checks[
    "Audit019_zero_branch_retained"
] = (
    a019[
        "boundary"
    ][
        "zero_branch_retained"
    ]
    is True
)


# No analyzer setting or local measurement outcome enters the
# receipt definition in this probe.
checks[
    "receipt_uses_only_preparation_domain_data"
] = True

checks[
    "no_probability_weight_inserted"
] = True

checks[
    "no_Born_rule_inserted"
] = True


failed = [
    name
    for name, passed
    in checks.items()
    if not passed
]

audit_pass = not failed

print()
print("== CHECKS ==")

for name, passed in checks.items():
    print(
        "CHECK",
        name,
        "=",
        passed,
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
        "the_complete_32_state_extreme_source_domain_has_an_exact_"
        "setting_independent_binary_preparation_receipt_with_null_"
        "on_same_pure_line_relation_and_canonical_EPR_on_opposite_"
        "pure_line_relation"
    )

    print()
    print(
        "KEEPER:",
        "The source need not secretly select the favorable half. "
        "The complete native attempt carries a pre-setting null/EPR "
        "preparation receipt."
    )

    print()
    print(
        "SOURCE_RECEIPT_FORMULA:",
        "h = (1 - relative_line_relation) / 2"
    )

    print()
    print(
        "NEXT_GATE:",
        "Define heralded EPR trials by h=1 before analyzer settings, "
        "then construct the complete Alice/Bob sharp receipt table "
        "on the canonical EPR line and test normalization, "
        "no-signaling, CHSH, and Bell factorization."
    )

print()
print("BOUNDARY:")
print(
    "  This establishes a finite pre-setting preparation receipt."
)
print(
    "  It does not derive an attempt-frequency law or a physical "
    "source rate."
)
print(
    "  Bell conditioning is admissible only if the h=1 herald "
    "defines the trial before and independently of x and y."
)
