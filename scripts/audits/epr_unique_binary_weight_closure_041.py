#!/usr/bin/env python3

from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from math import sqrt
import json

HERE = Path(__file__).resolve().parents[2]

A035 = (
    HERE
    / "artifacts/json"
    / "epr_heralded_conditional_hilbert_table_035.v1.json"
)

A036 = (
    HERE
    / "artifacts/json"
    / "epr_native_symmetry_weighting_reduction_036.v1.json"
)

A040 = (
    HERE
    / "artifacts/json"
    / "epr_native_joint_contraction_kernel_040b.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_unique_binary_weight_closure_041.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_unique_binary_weight_closure_041.md"
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


# Exact quadratic field:
#
#     x = a + b sqrt(5)
#
def q(a=0, b=0):
    return (
        Fraction(a),
        Fraction(b),
    )


def qadd(x, y):
    return (
        x[0] + y[0],
        x[1] + y[1],
    )


def qscale(c, x):
    c = Fraction(c)

    return (
        c * x[0],
        c * x[1],
    )


def qnum(x):
    return (
        float(x[0])
        + float(x[1]) * sqrt(5.0)
    )


def fstr(x):
    if x.denominator == 1:
        return str(
            x.numerator
        )

    return (
        str(x.numerator)
        + "/"
        + str(x.denominator)
    )


def qstr(x):
    a, b = x

    if b == 0:
        return fstr(a)

    if a == 0:
        if b == 1:
            return "sqrt(5)"
        if b == -1:
            return "-sqrt(5)"

        return (
            fstr(b)
            + "*sqrt(5)"
        )

    sign = "+" if b > 0 else "-"
    bb = abs(b)

    if bb == 1:
        tail = "sqrt(5)"
    else:
        tail = (
            fstr(bb)
            + "*sqrt(5)"
        )

    return (
        fstr(a)
        + " "
        + sign
        + " "
        + tail
    )


def rank_rational(A):
    A = [
        [
            Fraction(x)
            for x in row
        ]
        for row in A
    ]

    rows = len(A)
    cols = len(A[0])

    rank = 0
    col = 0

    while (
        rank < rows
        and col < cols
    ):
        pivot = None

        for r in range(
            rank,
            rows,
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

        p = A[rank][col]

        A[rank] = [
            x / p
            for x in A[rank]
        ]

        for r in range(rows):
            if r == rank:
                continue

            c = A[r][col]

            if c == 0:
                continue

            A[r] = [
                x - c * y
                for x, y in zip(
                    A[r],
                    A[rank],
                )
            ]

        rank += 1
        col += 1

    return rank


print("== 041 UNIQUE BINARY WEIGHT CLOSURE ==")

a035 = load(A035)
a036 = load(A036)
a040 = load(A040)

checks = {}

checks["Audit035_passes"] = (
    a035["audit_pass"] is True
)

checks["Audit036_passes"] = (
    a036["audit_pass"] is True
)

checks["Audit040B_passes"] = (
    a040["audit_pass"] is True
)

checks[
    "native_operator_bilinearity_closed"
] = (
    a040[
        "boundary"
    ][
        "native_joint_operator_bilinearity_closed"
    ]
    is True
)

checks[
    "native_frequency_interpretation_open_at_entry"
] = (
    a040[
        "boundary"
    ][
        "native_trial_frequency_law_derived"
    ]
    is False
)


# ------------------------------------------------------------
# Binary weight system.
#
# Unknown vector:
#
#     p = [p++, p+-, p-+, p--].
#
# Constraints:
#
#     normalization      = 1
#     Alice signed mean  = 0
#     Bob signed mean    = 0
#     product mean       = E
#
# Coefficient matrix:
#
#     [ 1  1  1  1 ]
#     [ 1  1 -1 -1 ]
#     [ 1 -1  1 -1 ]
#     [ 1 -1 -1  1 ]
#
# This is a Hadamard matrix with rank four.
# ------------------------------------------------------------

H = [
    [1, 1, 1, 1],
    [1, 1, -1, -1],
    [1, -1, 1, -1],
    [1, -1, -1, 1],
]

H_rank = rank_rational(
    H
)

checks[
    "binary_moment_system_rank_four"
] = (
    H_rank == 4
)

print(
    "BINARY_MOMENT_MATRIX_RANK:",
    H_rank,
)

print(
    "BINARY_WEIGHT_SOLUTION_UNIQUE:",
    H_rank == 4,
)


def weights_from_E(E):
    one = q(
        1,
        0,
    )

    same = qscale(
        Fraction(
            1,
            4,
        ),
        qadd(
            one,
            E,
        ),
    )

    opposite = qscale(
        Fraction(
            1,
            4,
        ),
        qadd(
            one,
            qscale(
                -1,
                E,
            ),
        ),
    )

    return {
        "++": same,
        "+-": opposite,
        "-+": opposite,
        "--": same,
    }


def verify_weights(E, w):
    total = q()

    for value in w.values():
        total = qadd(
            total,
            value,
        )

    alice = qadd(
        qadd(
            w["++"],
            w["+-"],
        ),
        qscale(
            -1,
            qadd(
                w["-+"],
                w["--"],
            ),
        ),
    )

    bob = qadd(
        qadd(
            w["++"],
            w["-+"],
        ),
        qscale(
            -1,
            qadd(
                w["+-"],
                w["--"],
            ),
        ),
    )

    corr = qadd(
        qadd(
            w["++"],
            w["--"],
        ),
        qscale(
            -1,
            qadd(
                w["+-"],
                w["-+"],
            ),
        ),
    )

    return (
        total == q(1),
        alice == q(0),
        bob == q(0),
        corr == E,
    )


# ------------------------------------------------------------
# Exact native operator correlations for the frozen quartet.
#
# From Audit040B:
#
#     Gamma_ii = -1
#
# and for distinct axes:
#
#     Gamma_ij = -S_ij / sqrt(5)
#              = -(S_ij/5) sqrt(5).
#
# Audit036 gives:
#
#     S_02 = +1
#     S_10 = +1
#     S_12 = -1.
# ------------------------------------------------------------

S = a036[
    "native_symmetry"
][
    "signed_Gram"
]

checks[
    "frozen_native_sign_pattern"
] = (
    S[0][2] == 1
    and S[1][0] == 1
    and S[1][2] == -1
)

correlations = {
    "A0B0":
        q(
            -1,
            0,
        ),

    "A0B1":
        q(
            0,
            Fraction(
                -S[0][2],
                5,
            ),
        ),

    "A1B0":
        q(
            0,
            Fraction(
                -S[1][0],
                5,
            ),
        ),

    "A1B1":
        q(
            0,
            Fraction(
                -S[1][2],
                5,
            ),
        ),
}

tables = {}

all_unique_constraints_pass = True
all_nonnegative = True

for label, E in correlations.items():
    w = weights_from_E(
        E
    )

    (
        norm_ok,
        alice_ok,
        bob_ok,
        corr_ok,
    ) = verify_weights(
        E,
        w,
    )

    if not all(
        (
            norm_ok,
            alice_ok,
            bob_ok,
            corr_ok,
        )
    ):
        all_unique_constraints_pass = False

    numeric = {
        key:
            qnum(value)
        for key, value in w.items()
    }

    if any(
        value < -1.0e-12
        for value in numeric.values()
    ):
        all_nonnegative = False

    tables[label] = {
        "correlation_exact":
            qstr(E),

        "correlation_numeric":
            qnum(E),

        "weights_exact": {
            key:
                qstr(value)
            for key, value in w.items()
        },

        "weights_numeric":
            numeric,

        "normalization":
            norm_ok,

        "Alice_mean_zero":
            alice_ok,

        "Bob_mean_zero":
            bob_ok,

        "correlation_recovered":
            corr_ok,
    }


checks[
    "all_native_operator_correlations_have_unique_binary_weight_table"
] = (
    all_unique_constraints_pass
)

checks[
    "all_derived_algebraic_weights_nonnegative"
] = (
    all_nonnegative
)


print()
print("== UNIQUE ALGEBRAIC WEIGHT TABLES ==")

for label, row in tables.items():
    print()
    print(
        label,
        "E=",
        row[
            "correlation_exact"
        ],
    )

    print(
        "  weights=",
        row[
            "weights_exact"
        ],
    )


# ------------------------------------------------------------
# Recover exact CHSH from the source contraction.
# ------------------------------------------------------------

E00 = correlations[
    "A0B0"
]

E01 = correlations[
    "A0B1"
]

E10 = correlations[
    "A1B0"
]

E11 = correlations[
    "A1B1"
]

chsh = qadd(
    qadd(
        qadd(
            qscale(
                -1,
                E00,
            ),
            qscale(
                -1,
                E01,
            ),
        ),
        qscale(
            -1,
            E10,
        ),
    ),
    E11,
)

expected_chsh = q(
    1,
    Fraction(
        3,
        5,
    ),
)

checks[
    "source_contraction_CHSH_exact"
] = (
    chsh
    == expected_chsh
)

checks[
    "source_contraction_CHSH_exceeds_2"
] = (
    qnum(chsh) > 2.0
)

print()
print(
    "CHSH_EXACT:",
    qstr(chsh),
)

print(
    "CHSH_NUMERIC:",
    qnum(chsh),
)


# ------------------------------------------------------------
# Cross-check against Audit035 conditional Hilbert table.
# This is comparison only. Audit035 values are not inputs to the
# weight derivation above.
# ------------------------------------------------------------

target_numeric = float(
    a035[
        "CHSH"
    ][
        "numeric"
    ]
)

checks[
    "independent_algebraic_CHSH_matches_Audit035_control"
] = (
    abs(
        qnum(chsh)
        - target_numeric
    )
    < 1.0e-12
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
    "the_native_joint_operator_contraction_together_with_native_"
    "zero_one_wing_means_uniquely_determines_a_positive_normalized_"
    "four_outcome_algebraic_weight_table_for_each_binary_setting_pair"
    if audit_pass
    else
    "unique_binary_weight_closure_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_unique_binary_weight_closure_041",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "binary_moment_system": {
        "matrix":
            H,

        "rank":
            H_rank,

        "unknowns":
            [
                "w++",
                "w+-",
                "w-+",
                "w--",
            ],

        "constraints":
            [
                "normalization=1",
                "Alice_signed_mean=0",
                "Bob_signed_mean=0",
                "joint_product_mean=Gamma",
            ],

        "unique_solution":
            "w_ab=(1+a*b*Gamma)/4",
    },

    "source_operator": {
        "joint_contraction":
            "Gamma(A,B)=Tr[Podd(A tensor B)]",

        "local_means":
            0,

        "native_visibility":
            "1/sqrt(5)",
    },

    "frozen_quartet_tables":
        tables,

    "CHSH": {
        "exact":
            qstr(chsh),

        "numeric":
            qnum(chsh),

        "exceeds_2":
            qnum(chsh) > 2.0,
    },

    "checks":
        checks,

    "boundary": {
        "native_joint_operator_kernel_closed":
            True,

        "binary_algebraic_weight_table_unique":
            audit_pass,

        "weights_positive_and_normalized":
            all_nonnegative
            and all_unique_constraints_pass,

        "Audit035_target_used_as_derivation_input":
            False,

        "Born_frequency_postulate_used":
            False,

        "registered_trial_frequency_equals_algebraic_weight_derived":
            False,

        "physical_sampling_rule_derived":
            False,

        "common_adjoint_face_binding_still_conditional":
            True,

        "remaining_gate":
            (
                "derive that registered repeated trials realize "
                "the unique algebraic source weights"
            ),
    },

    "earned_statement": (
        "The canonical native joint source supplies the exact "
        "operator contraction Gamma through Audit040B, while Audit036 "
        "forces zero signed one-wing means. For binary outcomes, "
        "normalization, the two zero local means, and the joint "
        "product mean Gamma form a full-rank four-by-four moment "
        "system. The algebraic joint weights are therefore unique: "
        "w_ab=(1+a b Gamma)/4. For the frozen native quartet these "
        "weights are positive, normalized, have exact one-half "
        "marginals, and give CHSH 1+3/sqrt(5). The Audit035 table is "
        "recovered independently rather than used as input. What "
        "remains open is the operational statement that registered "
        "trial frequencies realize these canonical algebraic weights."
    ),

    "next_gate": (
        "Search only for a native repeated-trial or receipt-additivity "
        "law that identifies registered frequencies with the unique "
        "algebraic source weights. Do not revisit source selection, "
        "local hidden answer tables, pulse optimization, or visibility "
        "fitting."
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

note = f"""# Unique binary weight closure 041

## Result

Audit pass:

    {audit_pass}

The native joint source supplies

    Gamma(A,B)
      =
    Tr[P_odd (A tensor B)].

The native one-wing signed means vanish.

For binary outcomes the four unknown weights are fixed by

    sum w_ab = 1

    sum a w_ab = 0

    sum b w_ab = 0

    sum ab w_ab = Gamma.

The coefficient matrix has rank

    {H_rank}.

Therefore the solution is unique:

    w_ab
      =
    (1 + a b Gamma) / 4.

For the frozen quartet the resulting algebraic weights reproduce

    CHSH = {qstr(chsh)}.

They are positive and normalized and have exact one-half local
marginals.

Audit035 is used only as an independent comparison after the weights
are derived.

## Boundary

This is an algebraic weight theorem.

It does not yet prove that repeated registered trial frequencies equal
these weights.

No Born frequency postulate is used as an input.

The remaining operational gate is:

    registered frequency
        ?=
    canonical algebraic source weight.
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
    "BINARY_WEIGHT_SYSTEM_RANK:",
    H_rank,
)

print(
    "UNIQUE_WEIGHT_FORMULA:",
    "w_ab=(1+a*b*Gamma)/4",
)

print(
    "CHSH_EXACT:",
    qstr(chsh),
)

print(
    "REMAINING_GATE:",
    "registered_frequency_equals_algebraic_weight",
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
