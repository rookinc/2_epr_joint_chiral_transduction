#!/usr/bin/env python3

from itertools import product
from hashlib import sha256
from pathlib import Path
import json

HERE = Path(__file__).resolve().parents[2]

A035 = (
    HERE
    / "artifacts/json"
    / "epr_heralded_conditional_hilbert_table_035.v1.json"
)

A038 = (
    HERE
    / "artifacts/json"
    / "epr_native_rank3_kernel_visibility_038.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_heralded_local_table_obstruction_039.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_heralded_local_table_obstruction_039.md"
)


def load(path):
    return json.loads(
        path.read_text(
            encoding="utf-8",
        )
    )


def digest_json(obj):
    raw = json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")

    return sha256(raw).hexdigest()


print("== 039 HERALDED LOCAL TABLE OBSTRUCTION ==")

a035 = load(A035)
a038 = load(A038)

checks = {}

checks["Audit035_passes"] = (
    a035["audit_pass"] is True
)

checks["Audit038_passes"] = (
    a038["audit_pass"] is True
)

checks[
    "herald_is_pre_setting"
] = (
    a035[
        "source"
    ][
        "analyzer_settings_used_to_define_herald"
    ]
    is False
)

target_chsh = float(
    a035[
        "CHSH"
    ][
        "numeric"
    ]
)

checks[
    "target_CHSH_exceeds_2"
] = (
    target_chsh > 2.0
)


# ------------------------------------------------------------
# Exhaust all deterministic binary local answer tables
#
#     (A0,A1,B0,B1) in {+1,-1}^4.
#
# For each table:
#
#     S
#       =
#     A0 B0
#       +
#     A0 B1
#       +
#     A1 B0
#       -
#     A1 B1.
#
# Every deterministic vertex must satisfy |S|=2.
# ------------------------------------------------------------

rows = []

for (
    A0,
    A1,
    B0,
    B1,
) in product(
    (+1, -1),
    repeat=4,
):
    chsh = (
        A0 * B0
        + A0 * B1
        + A1 * B0
        - A1 * B1
    )

    rows.append({
        "A0": A0,
        "A1": A1,
        "B0": B0,
        "B1": B1,
        "CHSH": chsh,
        "abs_CHSH": abs(chsh),
    })

profile = {}

for row in rows:
    key = str(
        row[
            "CHSH"
        ]
    )

    profile[key] = (
        profile.get(
            key,
            0,
        )
        + 1
    )

max_abs = max(
    row[
        "abs_CHSH"
    ]
    for row in rows
)

checks[
    "deterministic_table_count_16"
] = (
    len(rows) == 16
)

checks[
    "every_deterministic_local_table_has_abs_CHSH_2"
] = all(
    row[
        "abs_CHSH"
    ]
    == 2
    for row in rows
)

checks[
    "deterministic_local_max_CHSH_2"
] = (
    max_abs == 2
)

print(
    "DETERMINISTIC_TABLE_COUNT:",
    len(rows),
)

print(
    "DETERMINISTIC_CHSH_PROFILE:",
    profile,
)

print(
    "DETERMINISTIC_MAX_ABS_CHSH:",
    max_abs,
)


# ------------------------------------------------------------
# Any stochastic local model with a setting-independent source
# is a convex combination of deterministic local tables.
#
# Since CHSH is affine in the joint probabilities, convexity gives
#
#     |S| <= 2.
#
# This includes arbitrary nonnegative source weights mu(lambda),
# not merely uniform counting over the 16 heralded source states.
# ------------------------------------------------------------

checks[
    "convex_local_models_bound_2"
] = (
    checks[
        "every_deterministic_local_table_has_abs_CHSH_2"
    ]
)

checks[
    "target_outside_local_polytope"
] = (
    target_chsh > max_abs
)

print()
print(
    "TARGET_CHSH:",
    target_chsh,
)

print(
    "TARGET_OUTSIDE_LOCAL_POLYTOPE:",
    checks[
        "target_outside_local_polytope"
    ],
)


# ------------------------------------------------------------
# The target distinct-setting probabilities contain sqrt(5):
#
#     1/4 +/- sqrt(5)/20.
#
# Therefore exact equal-count frequencies over any finite number
# of equiprobable source states cannot reproduce them.
#
# This is recorded as a mathematical boundary rather than used as
# the Bell obstruction: arbitrary algebraic native weights remain
# possible.
# ------------------------------------------------------------

target_probability_form = (
    "1/4 +/- sqrt(5)/20"
)

print()
print(
    "TARGET_DISTINCT_SETTING_PROBABILITY_FORM:",
    target_probability_form,
)

print(
    "FINITE_EQUAL_COUNT_EXACT_TARGET:",
    False,
)


# ------------------------------------------------------------
# Classification.
# ------------------------------------------------------------

failed = [
    name
    for name, passed
    in checks.items()
    if not passed
]

audit_pass = not failed

verdict = (
    "the_pre_setting_h1_herald_cannot_reproduce_the_target_EPR_"
    "table_by_any_setting_independent_positive_mixture_of_"
    "preexisting_local_binary_answer_tables"
    if audit_pass
    else
    "heralded_local_table_obstruction_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_heralded_local_table_obstruction_039",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "heralded_source": {
        "receipt":
            "h=1 from Audit034C",

        "settings_used_to_define_herald":
            False,
    },

    "local_answer_table_test": {
        "deterministic_table_count":
            len(rows),

        "CHSH_profile":
            profile,

        "maximum_absolute_CHSH":
            max_abs,

        "stochastic_extension":
            (
                "all setting-independent positive mixtures are "
                "convex combinations of these deterministic vertices"
            ),
    },

    "target": {
        "CHSH":
            a035[
                "CHSH"
            ][
                "exact"
            ],

        "CHSH_numeric":
            target_chsh,

        "distinct_setting_probability_form":
            target_probability_form,
    },

    "checks":
        checks,

    "boundary": {
        "rules_out_uniform_local_counting":
            True,

        "rules_out_arbitrary_positive_setting_independent_local_hidden_variable_weights":
            True,

        "rules_out_preexisting_local_answer_tables":
            True,

        "rules_out_joint_contextual_receipt":
            False,

        "rules_out_algebraic_amplitude_weights":
            False,

        "rules_out_native_positive_rank3_kernel":
            False,

        "Born_rule_assumed":
            False,

        "remote_setting_dependence_allowed":
            False,
    },

    "earned_statement": (
        "After the setting-independent h=1 preparation herald, any "
        "model assigning preexisting binary answers A0,A1,B0,B1 to "
        "each source state has pointwise CHSH magnitude two. Any "
        "setting-independent positive mixture of such tables therefore "
        "has CHSH at most two by convexity. The exact Audit-035 target "
        "has CHSH 1+3/sqrt(5)>2 and cannot arise from such a local "
        "answer-table census. Thus the remaining native receipt law "
        "must not reduce to preassigned local responses with positive "
        "setting-independent source weights."
    ),

    "next_gate": (
        "Construct a genuinely joint relational receipt on the h=1 "
        "preparation domain whose local marginals remain setting-local "
        "but whose joint kernel is not a convex mixture of preexisting "
        "local answer tables. The Audit038 positive rank-three kernel "
        "is the exact target structure."
    ),
}

artifact[
    "artifact_sha256"
] = digest_json(
    artifact
)

JSON_OUT.parent.mkdir(
    parents=True,
    exist_ok=True,
)

NOTE_OUT.parent.mkdir(
    parents=True,
    exist_ok=True,
)

JSON_OUT.write_text(
    json.dumps(
        artifact,
        indent=2,
        sort_keys=True,
        ensure_ascii=True,
    )
    + "\n",
    encoding="ascii",
)

note = f"""# Heralded local-table obstruction 039

## Result

Audit pass:

    {audit_pass}

The h=1 source herald is fixed before analyzer settings.

Every deterministic local binary answer table

    (A0,A1,B0,B1)

has

    |CHSH| = 2.

There are exactly

    {len(rows)}

such deterministic vertices.

Any setting-independent positive stochastic mixture of these vertices
therefore also obeys

    |CHSH| <= 2.

The frozen target instead has

    CHSH = {a035["CHSH"]["exact"]}

which is greater than two.

Therefore the native h=1 source cannot reach the target by assigning
preexisting local binary answers and counting or positively weighting
those answer tables.

The distinct-setting target probabilities also have the exact form

    1/4 +/- sqrt(5)/20,

so finite equal-count frequencies cannot reproduce them exactly.

## Consequence

The missing native receipt law must be genuinely relational at the
joint level.

It may not be replaced by a hidden local answer census.

Local marginals must remain independent of the remote setting.

The Audit038 positive rank-three kernel remains an admissible target.

## Boundary

This does not rule out:

- a joint contextual receipt,
- algebraic native weights,
- a positive rank-three kernel,
- or another finite mechanism not factorizing into preexisting local
  answer tables.

No Born rule is assumed.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print()
print(
    "AUDIT_PASS:",
    audit_pass,
)

print(
    "VERDICT:",
    verdict,
)

print(
    "FAILED_CHECK_COUNT:",
    len(failed),
)

print(
    "FAILED_CHECKS:",
    failed,
)

print(
    "NEXT_GATE:",
    "native_joint_relational_receipt_kernel",
)

print(
    "JSON_OUT:",
    JSON_OUT,
)

print(
    "NOTE_OUT:",
    NOTE_OUT,
)

print(
    "JSON_SHA256:",
    sha256(
        JSON_OUT.read_bytes()
    ).hexdigest(),
)
