#!/usr/bin/env python3

from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json

HERE = Path(__file__).resolve().parents[2]

A011 = (
    HERE
    / "artifacts/json"
    / "epr_native_complex_balanced_boundary_RP3_011.v1.json"
)

A012 = (
    HERE
    / "artifacts/json"
    / "epr_native_boundary_antisymmetrizer_closure_012.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_complete_exchange_parity_closure_013.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_complete_exchange_parity_closure_013.md"
)


def load(path):
    return json.loads(path.read_text())


def eye(n):
    return [
        [Fraction(1) if i == j else Fraction(0) for j in range(n)]
        for i in range(n)
    ]


def matmul(a, b):
    return [
        [
            sum(
                a[i][k] * b[k][j]
                for k in range(len(b))
            )
            for j in range(len(b[0]))
        ]
        for i in range(len(a))
    ]


def matadd(a, b):
    return [
        [
            a[i][j] + b[i][j]
            for j in range(len(a[0]))
        ]
        for i in range(len(a))
    ]


def matsub(a, b):
    return [
        [
            a[i][j] - b[i][j]
            for j in range(len(a[0]))
        ]
        for i in range(len(a))
    ]


def matscale(s, a):
    s = Fraction(s)

    return [
        [
            s * x
            for x in row
        ]
        for row in a
    ]


def rank(matrix):
    a = [
        list(row)
        for row in matrix
    ]

    rows = len(a)
    cols = len(a[0])

    r = 0

    for c in range(cols):
        pivot = None

        for i in range(r, rows):
            if a[i][c] != 0:
                pivot = i
                break

        if pivot is None:
            continue

        a[r], a[pivot] = a[pivot], a[r]

        q = a[r][c]

        a[r] = [
            x / q
            for x in a[r]
        ]

        for i in range(rows):
            if i == r:
                continue

            q = a[i][c]

            if q == 0:
                continue

            a[i] = [
                a[i][j] - q * a[r][j]
                for j in range(cols)
            ]

        r += 1

        if r == rows:
            break

    return r


def equal(a, b):
    return a == b


def digest_json(obj):
    raw = json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("ascii")

    return sha256(raw).hexdigest()


print("PROGRESS: 1/6 load sealed boundary audits")

a011 = load(A011)
a012 = load(A012)

checks = {}

checks["audit_011_passed"] = (
    a011["audit_pass"] is True
)

checks["audit_012_passed"] = (
    a012["audit_pass"] is True
)

checks["native_boundary_real_dimension_4"] = (
    a012[
        "native_boundary_space"
    ]["real_dimension"]
    == 4
)

checks["audit_012_odd_rank_1"] = (
    a012[
        "closure"
    ]["induced_real_rank"]
    == 1
)

checks["audit_012_kernel_dim_3"] = (
    a012[
        "closure"
    ]["kernel_real_dimension"]
    == 3
)

print("PROGRESS: 2/6 construct canonical RP3 coordinates")

# Audit 003/011 Real-form coordinates can be written as
#
# M =
#   [[a+ib, c+id],
#    [c-id, a-ib]]
#
# with real coordinates
#
#   (a,b,c,d).
#
# Factor exchange M -> M^T acts as
#
#   (a,b,c,d) -> (a,b,c,-d).

I = eye(4)

X = [
    [1, 0, 0,  0],
    [0, 1, 0,  0],
    [0, 0, 1,  0],
    [0, 0, 0, -1],
]

checks["X_square_identity"] = equal(
    matmul(X, X),
    I,
)

print("PROGRESS: 3/6 construct complete parity closure")

P_even = matscale(
    Fraction(1, 2),
    matadd(I, X),
)

P_odd = matscale(
    Fraction(1, 2),
    matsub(I, X),
)

checks["P_even_idempotent"] = equal(
    matmul(P_even, P_even),
    P_even,
)

checks["P_odd_idempotent"] = equal(
    matmul(P_odd, P_odd),
    P_odd,
)

checks["P_even_P_odd_zero"] = equal(
    matmul(P_even, P_odd),
    [
        [Fraction(0) for _ in range(4)]
        for _ in range(4)
    ],
)

checks["P_odd_P_even_zero"] = equal(
    matmul(P_odd, P_even),
    [
        [Fraction(0) for _ in range(4)]
        for _ in range(4)
    ],
)

checks["P_even_plus_P_odd_identity"] = equal(
    matadd(P_even, P_odd),
    I,
)

print("PROGRESS: 4/6 measure branch dimensions")

even_rank = rank(P_even)
odd_rank = rank(P_odd)

checks["even_branch_real_rank_3"] = (
    even_rank == 3
)

checks["odd_branch_real_rank_1"] = (
    odd_rank == 1
)

checks["branch_ranks_sum_to_domain_dimension"] = (
    even_rank + odd_rank == 4
)

print("PROGRESS: 5/6 identify branch roles")

# Coordinate d spans the odd ray:
#
# M_d =
#   [[0, i d],
#    [-i d, 0]]
#
# projectively [|01>-|10>].
#
# Coordinates a,b,c span the even branch.

odd_generator = [
    Fraction(0),
    Fraction(0),
    Fraction(0),
    Fraction(1),
]

even_product_example = [
    Fraction(1),
    Fraction(0),
    Fraction(1),
    Fraction(0),
]

checks["odd_generator_survives_P_odd"] = (
    [
        row[0]
        for row in matmul(
            P_odd,
            [[x] for x in odd_generator],
        )
    ]
    == odd_generator
)

checks["odd_generator_killed_by_P_even"] = (
    [
        row[0]
        for row in matmul(
            P_even,
            [[x] for x in odd_generator],
        )
    ]
    == [
        Fraction(0)
        for _ in range(4)
    ]
)

checks["even_example_survives_P_even"] = (
    [
        row[0]
        for row in matmul(
            P_even,
            [[x] for x in even_product_example],
        )
    ]
    == even_product_example
)

checks["even_example_killed_by_P_odd"] = (
    [
        row[0]
        for row in matmul(
            P_odd,
            [[x] for x in even_product_example],
        )
    ]
    == [
        Fraction(0)
        for _ in range(4)
    ]
)

print("PROGRESS: 6/6 classify")

failed_checks = [
    name
    for name, passed in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "native_RP3_boundary_has_complete_setting_independent_"
    "exchange_parity_closure_with_rank1_odd_and_rank3_even_branches"
    if audit_pass
    else
    "complete_exchange_parity_closure_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_complete_exchange_parity_closure_013",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "native_boundary": {
        "real_dimension":
            4,

        "projective_space":
            "RP^3",

        "coordinates":
            "(a,b,c,d)",
    },

    "factor_exchange": {
        "action":
            "(a,b,c,d) -> (a,b,c,-d)",

        "square":
            "I",
    },

    "complete_closure": {
        "even":
            "P_even = (I+X)/2",

        "odd":
            "P_odd = (I-X)/2",

        "sum":
            "P_even + P_odd = I",

        "cross_products":
            "P_even P_odd = P_odd P_even = 0",

        "even_real_rank":
            even_rank,

        "odd_real_rank":
            odd_rank,
    },

    "branch_interpretation": {
        "odd":
            (
                "canonical rank-one EPR line "
                "[|01>-|10>]"
            ),

        "even":
            (
                "three-real-dimensional complementary "
                "exchange-even boundary branch"
            ),
    },

    "earned_statement": (
        "The three-dimensional kernel of the Audit-012 "
        "antisymmetrizer is not an undefined discarded null "
        "space. It is exactly the image of the complementary "
        "exchange-even projector P_even=(I+X)/2. Together "
        "P_even and P_odd form a complete setting-independent "
        "two-branch closure of the native RP3 boundary space. "
        "The odd branch is the canonical rank-one EPR line and "
        "the even branch has real rank three."
    ),

    "checks":
        checks,

    "boundary": {
        "complete_two_branch_closure_constructed":
            True,

        "odd_branch_not_permitted_to_be_postselected_after_settings":
            True,

        "source_may_prepare_odd_before_settings_if_natively_derived":
            True,

        "otherwise_even_branch_must_be_retained":
            True,

        "branch_rank_not_interpreted_as_probability":
            True,

        "no_Born_rule":
            True,

        "no_trial_frequency_law":
            True,

        "no_CHSH_claim":
            True,
    },

    "next_gate": (
        "Determine whether native source preparation before "
        "analyzer selection is constrained to the odd branch. "
        "If not, carry the exchange-even branch forward as a "
        "real preparation outcome. Do not condition the EPR "
        "sample on observing the odd branch after settings are "
        "chosen."
    ),
}

artifact["artifact_sha256"] = digest_json(
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
    )
    + "\n",
    encoding="ascii",
)

note = f"""# Complete exchange-parity closure 013

## Result

Audit pass:

    {audit_pass}

The native balanced boundary is the four-real-dimensional Real form
underlying RP3.

Factor exchange acts in exact Real-form coordinates as

    (a,b,c,d)
      ->
    (a,b,c,-d).

Therefore it has two canonical complementary projectors:

    P_even = (I+X)/2
    P_odd  = (I-X)/2.

They satisfy

    P_even^2 = P_even
    P_odd^2  = P_odd

    P_even P_odd = 0
    P_odd P_even = 0

    P_even + P_odd = I.

## Branch dimensions

The exchange-even branch has real rank

    {even_rank}.

The exchange-odd branch has real rank

    {odd_rank}.

The odd branch is exactly the canonical EPR line

    [|01>-|10>].

The even branch is the complete three-real-dimensional complement.

## No-postselection consequence

The three-dimensional kernel of the odd closure is therefore not an
undefined null region.

It is a genuine complementary preparation branch.

A later EPR protocol has only two honest options:

1. derive a setting-independent native source rule that prepares the odd
   branch before analyzer settings are chosen;

2. retain the even branch in the complete experimental outcome space.

The odd branch may not be selected after seeing settings or outcomes.

No branch rank is interpreted as a probability here.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print(
    "EVEN_BRANCH_REAL_RANK:",
    even_rank,
)
print(
    "ODD_BRANCH_REAL_RANK:",
    odd_rank,
)
print(
    "COMPLETE_CLOSURE:",
    checks["P_even_plus_P_odd_identity"],
)
print(
    "ODD_BRANCH_IS_CANONICAL_EPR_LINE:",
    checks["odd_generator_survives_P_odd"],
)
print(
    "POSTSELECTION_ALLOWED:",
    False,
)
print("FAILED_CHECK_COUNT:", len(failed_checks))
print("FAILED_CHECKS:", failed_checks)
print("JSON_OUT:", JSON_OUT)
print("NOTE_OUT:", NOTE_OUT)
print(
    "JSON_SHA256:",
    sha256(JSON_OUT.read_bytes()).hexdigest(),
)
