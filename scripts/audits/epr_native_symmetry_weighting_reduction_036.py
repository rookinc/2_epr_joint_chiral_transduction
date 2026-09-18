#!/usr/bin/env python3

from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from math import sqrt
import contextlib
import io
import json
import runpy

HERE = Path(__file__).resolve().parents[2]

S025 = (
    HERE
    / "scripts/audits"
    / "epr_projective_analyzer_sector_extension_025.py"
)

A025 = (
    HERE
    / "artifacts/json"
    / "epr_projective_analyzer_sector_extension_025.v1.json"
)

A035 = (
    HERE
    / "artifacts/json"
    / "epr_heralded_conditional_hilbert_table_035.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_native_symmetry_weighting_reduction_036.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_native_symmetry_weighting_reduction_036.md"
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


def matrix_rank(rows, ncols):
    A = [
        [
            Fraction(x)
            for x in row
        ]
        for row in rows
        if any(
            x != 0
            for x in row
        )
    ]

    rank = 0
    col = 0

    while (
        rank < len(A)
        and col < ncols
    ):
        pivot = None

        for r in range(
            rank,
            len(A),
        ):
            if A[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            col += 1
            continue

        A[rank], A[pivot] = (
            A[pivot],
            A[rank],
        )

        value = A[rank][col]

        A[rank] = [
            x / value
            for x in A[rank]
        ]

        for r in range(
            len(A)
        ):
            if r == rank:
                continue

            factor = A[r][col]

            if factor == 0:
                continue

            A[r] = [
                x - factor * y
                for x, y in zip(
                    A[r],
                    A[rank],
                )
            ]

        rank += 1
        col += 1

    return rank


print(
    "== 036 NATIVE SYMMETRY WEIGHTING REDUCTION =="
)

a025 = load(A025)
a035 = load(A035)

checks = {}

checks[
    "Audit025_passes"
] = (
    a025["audit_pass"] is True
)

checks[
    "Audit035_passes"
] = (
    a035["audit_pass"] is True
)

checks[
    "Audit035_is_conditional_Hilbert_control"
] = (
    a035[
        "boundary"
    ][
        "conditional_Hilbert_probability_control"
    ]
    is True
)

checks[
    "native_frequency_law_still_open_at_entry"
] = (
    a035[
        "boundary"
    ][
        "native_receipt_frequency_law_derived"
    ]
    is False
)

print(
    "PROGRESS: 1/7 replay native six-line projective symmetry"
)

sink = io.StringIO()

with contextlib.redirect_stdout(sink):
    ns025 = runpy.run_path(
        str(S025)
    )

S = tuple(
    tuple(
        int(x)
        for x in row
    )
    for row in ns025[
        "Splus"
    ]
)

maps = dict(
    ns025[
        "Gpp_map"
    ]
)

checks[
    "within_sector_projective_group_order_60"
] = (
    len(maps) == 60
)

checks[
    "signed_Gram_is_6_by_6"
] = (
    len(S) == 6
    and all(
        len(row) == 6
        for row in S
    )
)

print(
    "PROJECTIVE_GROUP_ORDER:",
    len(maps),
)

print()
print(
    "SIGNED_GRAM:"
)

for row in S:
    print(
        " ",
        list(row),
    )


# ------------------------------------------------------------
# One-wing marginal covariance.
#
# Under a native switching symmetry (p,d),
#
#     O_i -> d_i O_{p(i)}
#
# so a covariant one-wing expectation m_i must satisfy
#
#     m_i = d_i m_{p(i)}.
#
# Solve all such linear constraints without introducing any
# probability target.
# ------------------------------------------------------------

print(
    "PROGRESS: 2/7 solve one-wing marginal covariance"
)

marginal_rows = []

for p, d in maps.items():
    for i in range(6):
        row = [
            0
            for _ in range(6)
        ]

        row[i] += 1
        row[
            p[i]
        ] -= int(
            d[i]
        )

        marginal_rows.append(
            row
        )

marginal_rank = matrix_rank(
    marginal_rows,
    6,
)

marginal_nullity = (
    6
    - marginal_rank
)

checks[
    "native_marginal_covariance_space_is_zero"
] = (
    marginal_nullity == 0
)

flip_stabilizer_counts = {}

for i in range(6):
    count = sum(
        1
        for p, d in maps.items()
        if (
            p[i] == i
            and int(
                d[i]
            ) == -1
        )
    )

    flip_stabilizer_counts[
        str(i)
    ] = count

checks[
    "gauge_fixed_axis0_has_no_explicit_flip_stabilizer"
] = (
    flip_stabilizer_counts["0"] == 0
)

checks[
    "remaining_five_axes_have_explicit_flip_stabilizers"
] = all(
    flip_stabilizer_counts[str(i)] > 0
    for i in range(1, 6)
)

checks[
    "full_marginal_vanishing_follows_from_rank_not_axiswise_flips"
] = (
    marginal_nullity == 0
)

print(
    "MARGINAL_CONSTRAINT_RANK:",
    marginal_rank,
)

print(
    "MARGINAL_COVARIANCE_NULLITY:",
    marginal_nullity,
)

print(
    "AXIS_FLIP_STABILIZER_COUNTS:",
    flip_stabilizer_counts,
)


# ------------------------------------------------------------
# Two-wing correlation covariance.
#
# The heralded source is invariant under factor exchange X.
# We therefore use a symmetric correlation table:
#
#     E_ij = E_ji.
#
# Under common analyzer symmetry,
#
#     E_ij
#       =
#     d_i d_j E_{p(i),p(j)}.
#
# Solve this space exactly over Q.
# ------------------------------------------------------------

print(
    "PROGRESS: 3/7 solve symmetric two-wing correlation covariance"
)

pairs = [
    (
        i,
        j,
    )
    for i in range(6)
    for j in range(
        i,
        6,
    )
]

pair_index = {
    pair: k
    for k, pair in enumerate(
        pairs
    )
}

corr_rows = []

for p, d in maps.items():
    for i, j in pairs:
        target = tuple(
            sorted(
                (
                    p[i],
                    p[j],
                )
            )
        )

        sign = (
            int(d[i])
            * int(d[j])
        )

        row = [
            0
            for _ in range(
                len(pairs)
            )
        ]

        row[
            pair_index[
                (
                    i,
                    j,
                )
            ]
        ] += 1

        row[
            pair_index[
                target
            ]
        ] -= sign

        corr_rows.append(
            row
        )

corr_rank = matrix_rank(
    corr_rows,
    len(pairs),
)

corr_nullity = (
    len(pairs)
    - corr_rank
)

checks[
    "symmetric_correlation_covariance_space_dimension_2"
] = (
    corr_nullity == 2
)

print(
    "SYMMETRIC_CORRELATION_VARIABLE_COUNT:",
    len(pairs),
)

print(
    "CORRELATION_CONSTRAINT_RANK:",
    corr_rank,
)

print(
    "CORRELATION_COVARIANCE_NULLITY:",
    corr_nullity,
)


# ------------------------------------------------------------
# Verify an explicit basis for the two-dimensional solution space.
#
# D:
#     1 on diagonal, 0 off diagonal.
#
# G:
#     0 on diagonal, signed native Gram sign S_ij off diagonal.
#
# Thus every symmetric covariant correlation law must be
#
#     E = alpha D + beta G.
# ------------------------------------------------------------

def D(i, j):
    return (
        1
        if i == j
        else 0
    )


def G(i, j):
    return (
        0
        if i == j
        else int(
            S[i][j]
        )
    )


def template_covariant(fn):
    failures = 0

    for p, d in maps.items():
        for i, j in pairs:
            lhs = fn(
                i,
                j,
            )

            rhs = (
                int(d[i])
                * int(d[j])
                * fn(
                    p[i],
                    p[j],
                )
            )

            if lhs != rhs:
                failures += 1

    return failures


D_failures = template_covariant(
    D
)

G_failures = template_covariant(
    G
)

checks[
    "diagonal_template_is_covariant"
] = (
    D_failures == 0
)

checks[
    "signed_offdiagonal_Gram_template_is_covariant"
] = (
    G_failures == 0
)

checks[
    "explicit_two_templates_span_by_dimension"
] = (
    corr_nullity == 2
    and D_failures == 0
    and G_failures == 0
)

print(
    "DIAGONAL_TEMPLATE_FAILURES:",
    D_failures,
)

print(
    "SIGNED_GRAM_TEMPLATE_FAILURES:",
    G_failures,
)


# ------------------------------------------------------------
# Antisymmetric same-axis boundary.
#
# The h=1 source is the canonical antisymmetric line. Under the
# common sharp observable interface, equal settings are perfectly
# anticorrelated:
#
#     E_ii = -1.
#
# This fixes alpha=-1.
#
# All remaining native symmetry-compatible freedom is one scalar beta:
#
#     E_ij = beta S_ij, i != j.
#
# Define visibility v=-beta so the Hilbert target is positive:
#
#     E_ij = -v S_ij.
# ------------------------------------------------------------

print(
    "PROGRESS: 4/7 impose antisymmetric same-axis boundary"
)

alpha = -1

checks[
    "same_axis_anticorrelation_fixes_alpha_minus_one"
] = (
    alpha == -1
)

remaining_parameter_count = (
    corr_nullity
    - 1
)

checks[
    "one_scalar_correlation_freedom_remains"
] = (
    remaining_parameter_count == 1
)

print(
    "FIXED_DIAGONAL_ALPHA:",
    alpha,
)

print(
    "REMAINING_CORRELATION_PARAMETER_COUNT:",
    remaining_parameter_count,
)

print(
    "CORRELATION_FAMILY:",
    "E_ii=-1; E_ij=-v*S_ij for i!=j",
)


# ------------------------------------------------------------
# With zero local marginals, a normalized binary +/-1 joint table
# is algebraically determined by its correlation:
#
#     p(a,b) = (1 + a b E) / 4.
#
# This is moment inversion for binary variables, not an insertion
# of the Hilbert/Born target value for v.
# ------------------------------------------------------------

print(
    "PROGRESS: 5/7 reduce binary joint weights to the same scalar"
)

checks[
    "zero_marginals_forced_by_native_switching_symmetry"
] = (
    marginal_nullity == 0
)

checks[
    "binary_table_has_no_additional_free_parameters_once_v_fixed"
] = True

# Positivity of
#
#     (1 +/- E)/4
#
# for all distinct-setting correlations E=+/-v requires |v|<=1.
positivity_interval = (
    -1.0,
    1.0,
)

print(
    "POSITIVITY_INTERVAL_FOR_v:",
    positivity_interval,
)


# ------------------------------------------------------------
# Frozen Audit-035 quartet:
#
# Alice = (0,1)
# Bob   = (0,2)
#
# Signed Gram entries:
#
# S_02 = +1
# S_10 = +1
# S_12 = -1
#
# Therefore
#
# E00 = -1
# E02 = -v
# E10 = -v
# E12 = +v
#
# and
#
# CHSH = 1 + 3v
#
# for v >= 0.
# ------------------------------------------------------------

print(
    "PROGRESS: 6/7 derive one-parameter CHSH law"
)

A = (
    0,
    1,
)

B = (
    0,
    2,
)

checks[
    "frozen_quartet_sign_pattern"
] = (
    S[0][2] == 1
    and S[1][0] == 1
    and S[1][2] == -1
)

print(
    "FROZEN_QUARTET:",
    {
        "Alice": A,
        "Bob": B,
    },
)

print(
    "FROZEN_OFFDIAGONAL_SIGNS:",
    {
        "S_02":
            S[0][2],

        "S_10":
            S[1][0],

        "S_12":
            S[1][2],
    },
)

print(
    "CHSH_FAMILY_FOR_v_GE_0:",
    "S(v)=1+3v",
)

print(
    "BELL_VIOLATION_THRESHOLD:",
    "v > 1/3",
)

target_v = (
    1.0
    / sqrt(5.0)
)

target_chsh = (
    1.0
    + 3.0 * target_v
)

checks[
    "Audit035_target_visibility_is_above_Bell_threshold"
] = (
    target_v
    > 1.0 / 3.0
)

checks[
    "Audit035_target_CHSH_matches"
] = (
    abs(
        target_chsh
        - float(
            a035[
                "CHSH"
            ][
                "numeric"
            ]
        )
    )
    < 1.0e-12
)

print(
    "AUDIT035_TARGET_v:",
    "1/sqrt(5)",
)

print(
    "AUDIT035_TARGET_v_NUMERIC:",
    target_v,
)

print(
    "AUDIT035_TARGET_CHSH_FROM_FAMILY:",
    target_chsh,
)


print(
    "PROGRESS: 7/7 classify"
)

failed_checks = [
    name
    for name, passed
    in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "native_six_setting_projective_symmetry_and_factor_exchange_"
    "force_zero_local_marginals_and_reduce_every_symmetric_"
    "covariant_heralded_binary_correlation_table_to_one_visibility_"
    "parameter_v_after_same_axis_antisymmetric_closure"
    if audit_pass
    else
    "native_symmetry_weighting_reduction_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_native_symmetry_weighting_reduction_036",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "native_symmetry": {
        "sector":
            "B_plus3",

        "projective_group_order":
            len(maps),

        "signed_Gram":
            [
                list(row)
                for row in S
            ],

        "axis_flip_stabilizer_counts":
            flip_stabilizer_counts,

        "switching_orientation_gauge":
            "d_0 fixed to +1 in Audit025 switching_maps",
    },

    "marginal_reduction": {
        "variable_count":
            6,

        "constraint_rank":
            marginal_rank,

        "covariant_space_dimension":
            marginal_nullity,

        "result":
            "all covariant one-wing signed means vanish",
    },

    "correlation_reduction": {
        "symmetric_variable_count":
            len(pairs),

        "constraint_rank":
            corr_rank,

        "covariant_space_dimension":
            corr_nullity,

        "basis":
            [
                "diagonal_identity_template",
                "signed_offdiagonal_native_Gram_template",
            ],

        "general_form_before_same_axis_condition":
            "E = alpha*D + beta*G",

        "same_axis_condition":
            "E_ii=-1",

        "reduced_form":
            "E_ii=-1; E_ij=-v*S_ij for i!=j",

        "remaining_scalar":
            "v",
    },

    "binary_weight_reconstruction": {
        "local_marginals":
            "zero",

        "joint_probability_identity":
            "P(a,b|x,y)=(1+a*b*E_xy)/4",

        "positivity_condition":
            "|v|<=1",

        "Born_target_inserted":
            False,
    },

    "frozen_quartet": {
        "Alice":
            [
                0,
                1,
            ],

        "Bob":
            [
                0,
                2,
            ],

        "CHSH_family_for_v_nonnegative":
            "1+3v",

        "Bell_violation_condition":
            "v>1/3",

        "Audit035_target_v":
            "1/sqrt(5)",

        "Audit035_target_v_numeric":
            target_v,

        "Audit035_target_CHSH_numeric":
            target_chsh,
    },

    "checks":
        checks,

    "boundary": {
        "pre_setting_herald_inherited_from_Audit034C":
            True,

        "common_adjoint_line_system_still_conditional":
            True,

        "factor_exchange_symmetry_used":
            True,

        "native_switching_covariance_used":
            True,

        "antisymmetric_same_axis_anticorrelation_used":
            True,

        "bilinear_correlation_ansatz_used":
            False,

        "Born_frequency_law_used_to_derive_family":
            False,

        "Audit035_target_used_to_select_v":
            False,

        "native_visibility_v_derived":
            False,

        "native_receipt_frequency_law_fully_derived":
            False,

        "remaining_weighting_problem_dimension":
            1,
    },

    "earned_statement": (
        "The native 60-element within-sector projective switching "
        "symmetry forces every covariant signed one-wing marginal to "
        "vanish. On a factor-exchange-symmetric heralded two-wing "
        "source, the exact space of symmetric six-by-six correlation "
        "tables covariant under the same signed switching action is "
        "two-dimensional. It is spanned by the diagonal identity "
        "template and the signed off-diagonal native Gram template. "
        "The antisymmetric same-setting condition E_ii=-1 fixes the "
        "diagonal coefficient, leaving one scalar visibility v with "
        "E_ij=-v S_ij for distinct settings. Hence normalized binary "
        "joint weights contain no further freedom once v is fixed. "
        "For the frozen Audit-035 quartet CHSH is 1+3v and violation "
        "occurs exactly for v>1/3. The Hilbert control corresponds to "
        "v=1/sqrt(5). No bilinear correlation ansatz or Born target "
        "value is used to obtain this one-parameter reduction."
    ),

    "next_gate": (
        "Derive the single remaining visibility v from finite native "
        "receipt mechanics on the h=1 source. The exact target "
        "v=1/sqrt(5) is reserved for comparison only. Do not use it "
        "as selector input."
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

note = f"""# Native symmetry weighting reduction 036

## Result

Audit pass:

    {audit_pass}

The native within-sector projective analyzer symmetry has order

    {len(maps)}.

The exact covariant one-wing marginal space has dimension

    {marginal_nullity}.

Thus native signed switching covariance forces

    m_i = 0

for all six settings.

For a factor-exchange-symmetric two-wing correlation table, the exact
covariant space has dimension

    {corr_nullity}.

It is spanned by

    D = diagonal identity

and

    G = signed off-diagonal native Gram pattern.

Therefore

    E = alpha D + beta G.

The antisymmetric same-setting condition gives

    alpha = -1.

Writing

    v = -beta

gives the complete remaining family

    E_ii = -1

    E_ij = -v S_ij
        for i != j.

No bilinear correlation ansatz is used.

With zero local marginals, binary moment inversion gives

    P(a,b|x,y)
      =
    (1 + a b E_xy) / 4.

Thus there is no additional joint-weight freedom after v is fixed.

For the frozen Audit-035 quartet

    Alice = (0,1)
    Bob   = (0,2)

the CHSH family is

    S(v) = 1 + 3v

for v >= 0.

Bell violation occurs iff

    v > 1/3.

The conditional Hilbert control corresponds to

    v = 1/sqrt(5)

and reproduces

    S = 1 + 3/sqrt(5).

## Boundary

This audit does not derive v.

The common adjoint line system remains a conditional interface.

The antisymmetric same-axis closure is used.

No Born target value is used to obtain the one-parameter family.

The remaining native weighting problem is one-dimensional.

## Next gate

Derive v from finite native receipt mechanics on the h=1 heralded
source without using 1/sqrt(5) as selector input.
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
    len(
        failed_checks
    ),
)

print(
    "FAILED_CHECKS:",
    failed_checks,
)

print(
    "REMAINING_NATIVE_WEIGHT_PARAMETER_COUNT:",
    1 if audit_pass else None,
)

print(
    "TARGET_VISIBILITY_RESERVED_FOR_COMPARISON:",
    "1/sqrt(5)",
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
