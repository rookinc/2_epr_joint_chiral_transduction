#!/usr/bin/env python3

from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json

HERE = Path(__file__).resolve().parents[2]

INPUT = (
    HERE
    / "artifacts/json"
    / "epr_native_deck_face_chiral_transduction_002.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_joint_chiral_tensor_consistency_003.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_joint_chiral_tensor_consistency_003.md"
)


def load(path):
    return json.loads(path.read_text())


def eye(n):
    return [
        [1 if i == j else 0 for j in range(n)]
        for i in range(n)
    ]


def zeros(r, c):
    return [[0 for _ in range(c)] for _ in range(r)]


def matmul(a, b):
    out = zeros(len(a), len(b[0]))

    for i in range(len(a)):
        for k in range(len(b)):
            if a[i][k] == 0:
                continue

            for j in range(len(b[0])):
                out[i][j] += a[i][k] * b[k][j]

    return out


def matscale(s, a):
    return [
        [s * x for x in row]
        for row in a
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
    return matadd(a, matscale(-1, b))


def kron(a, b):
    out = []

    for ar in a:
        rows = [
            []
            for _ in range(len(b))
        ]

        for av in ar:
            for i, br in enumerate(b):
                rows[i].extend(
                    av * bv
                    for bv in br
                )

        out.extend(rows)

    return out


def hcat(*mats):
    return [
        sum(
            (list(m[i]) for m in mats),
            [],
        )
        for i in range(len(mats[0]))
    ]


def rank(matrix):
    a = [
        [Fraction(x) for x in row]
        for row in matrix
    ]

    if not a:
        return 0

    rows = len(a)
    cols = len(a[0])

    r = 0
    c = 0

    while r < rows and c < cols:
        pivot = None

        for i in range(r, rows):
            if a[i][c] != 0:
                pivot = i
                break

        if pivot is None:
            c += 1
            continue

        a[r], a[pivot] = a[pivot], a[r]

        pv = a[r][c]

        a[r] = [
            x / pv
            for x in a[r]
        ]

        for i in range(rows):
            if i == r:
                continue

            factor = a[i][c]

            if factor == 0:
                continue

            a[i] = [
                a[i][j] - factor * a[r][j]
                for j in range(cols)
            ]

        r += 1
        c += 1

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


print("PROGRESS: 1/7 load Audit 002")

source = load(INPUT)

checks = {}

checks["audit_002_passed"] = (
    source["audit_pass"] is True
)

checks["audit_002_selected_native_a"] = (
    source[
        "selected_native_S1_kernel"
    ]["native_name"]
    == "a"
)

checks["audit_002_selected_S1"] = (
    source[
        "selected_native_S1_kernel"
    ]["face_operator"]
    == "S1"
)

print("PROGRESS: 2/7 construct local Program-01 doublet")

# Real coordinates:
#
#   (Re z0, Im z0, Re z1, Im z1)
#
# J is multiplication by i.
#
# S = sigma_x K:
#
#   (z0,z1) -> (conj z1, conj z0)

J = [
    [0, -1, 0, 0],
    [1,  0, 0, 0],
    [0,  0, 0, -1],
    [0,  0, 1,  0],
]

S = [
    [0,  0, 1,  0],
    [0,  0, 0, -1],
    [1,  0, 0,  0],
    [0, -1, 0,  0],
]

I4 = eye(4)

checks["J_square_minus_identity"] = equal(
    matmul(J, J),
    matscale(-1, I4),
)

checks["S_square_identity"] = equal(
    matmul(S, S),
    I4,
)

checks["S_anticommutes_with_J"] = equal(
    matmul(S, J),
    matscale(
        -1,
        matmul(J, S),
    ),
)

print("PROGRESS: 3/7 construct same/opposite chiral tensor sheets")

JI = kron(J, I4)
IJ = kron(I4, J)

# Same complex orientation:
#
#   Ju tensor v = u tensor Jv

R_PLUS = matsub(JI, IJ)

# Opposite complex orientation:
#
#   Ju tensor v = -u tensor Jv

R_MINUS = matadd(JI, IJ)

rank_plus = rank(R_PLUS)
rank_minus = rank(R_MINUS)

q_plus_dim = 16 - rank_plus
q_minus_dim = 16 - rank_minus

checks["same_chiral_relation_rank_8"] = (
    rank_plus == 8
)

checks["opposite_chiral_relation_rank_8"] = (
    rank_minus == 8
)

checks["same_chiral_tensor_real_dim_8"] = (
    q_plus_dim == 8
)

checks["opposite_chiral_tensor_real_dim_8"] = (
    q_minus_dim == 8
)

print("PROGRESS: 4/7 apply one-sided S1")

SA = kron(S, I4)
SB = kron(I4, S)

SA_R_PLUS = matmul(SA, R_PLUS)
SB_R_PLUS = matmul(SB, R_PLUS)

checks["Alice_S_maps_plus_relations_to_minus"] = (
    rank(
        hcat(
            R_MINUS,
            SA_R_PLUS,
        )
    )
    == rank_minus
)

checks["Bob_S_maps_plus_relations_to_minus"] = (
    rank(
        hcat(
            R_MINUS,
            SB_R_PLUS,
        )
    )
    == rank_minus
)

checks["Alice_S_does_not_preserve_plus_sheet"] = (
    rank(
        hcat(
            R_PLUS,
            SA_R_PLUS,
        )
    )
    > rank_plus
)

checks["Bob_S_does_not_preserve_plus_sheet"] = (
    rank(
        hcat(
            R_PLUS,
            SB_R_PLUS,
        )
    )
    > rank_plus
)

print("PROGRESS: 5/7 compare the two local presentations")

# Audit 002 says the same native G1800 torsor swap can be
# represented at point level by acting with a on Alice or Bob.
#
# At Hilbert level this would require SA and SB to agree after
# passing from Q_plus to Q_minus.
#
# Test whether their difference vanishes modulo the target
# relation space.

DIFF = matsub(SA, SB)

difference_augmented_rank = rank(
    hcat(
        R_MINUS,
        DIFF,
    )
)

difference_induced_rank = (
    difference_augmented_rank
    - rank_minus
)

compatibility_real_dimension = (
    q_plus_dim
    - difference_induced_rank
)

checks["full_tensor_local_presentations_are_not_identical"] = (
    difference_induced_rank > 0
)

checks["compatibility_kernel_is_nonzero"] = (
    compatibility_real_dimension > 0
)

checks["compatibility_kernel_real_dimension_is_4"] = (
    compatibility_real_dimension == 4
)

print("PROGRESS: 6/7 construct two-sided comparison structure")

C = matmul(SA, SB)

I16 = eye(16)

checks["SA_and_SB_commute_as_real_maps"] = equal(
    matmul(SA, SB),
    matmul(SB, SA),
)

checks["two_sided_structure_square_identity"] = equal(
    matmul(C, C),
    I16,
)

checks["two_sided_structure_preserves_plus_sheet"] = (
    rank(
        hcat(
            R_PLUS,
            matmul(C, R_PLUS),
        )
    )
    == rank_plus
)

checks["two_sided_structure_anticommutes_with_joint_J"] = (
    equal(
        matmul(C, JI),
        matscale(
            -1,
            matmul(JI, C),
        ),
    )
)

# Fixed and anti-fixed dimensions on Q_plus.

C_MINUS_I = matsub(C, I16)
C_PLUS_I = matadd(C, I16)

fixed_induced_rank = (
    rank(
        hcat(
            R_PLUS,
            C_MINUS_I,
        )
    )
    - rank_plus
)

antifixed_induced_rank = (
    rank(
        hcat(
            R_PLUS,
            C_PLUS_I,
        )
    )
    - rank_plus
)

fixed_real_dimension = (
    q_plus_dim
    - fixed_induced_rank
)

antifixed_real_dimension = (
    q_plus_dim
    - antifixed_induced_rank
)

checks["two_sided_fixed_real_dimension_is_4"] = (
    fixed_real_dimension == 4
)

checks["two_sided_antifixed_real_dimension_is_4"] = (
    antifixed_real_dimension == 4
)

# Exact relation between the disagreement of the two
# one-sided presentations and the two-sided involution:
#
#   SA - SB = SA (I - C)
#
# where
#
#   C = SA SB.
#
# Since SA maps Q_plus isomorphically onto Q_minus,
# the compatibility kernel is exactly the C-fixed locus
# on Q_plus.

checks[
    "presentation_difference_factorizes_through_C"
] = equal(
    DIFF,
    matmul(
        SA,
        matsub(I16, C),
    ),
)

checks[
    "compatibility_kernel_equals_C_fixed_locus_by_dimension"
] = (
    compatibility_real_dimension
    == fixed_real_dimension
)

# C anticommutes with the joint complex structure.
#
# Therefore its fixed locus is not closed under J.
# It is a four-real-dimensional real form of the
# four-complex-dimensional joint tensor fiber.
#
# Its real projective rays form RP^3 inside CP^3.

checks[
    "compatibility_fixed_locus_is_not_complex_linear"
] = (
    fixed_real_dimension > 0
    and checks[
        "two_sided_structure_anticommutes_with_joint_J"
    ]
)

checks[
    "compatibility_fixed_real_dimension_is_4"
] = (
    compatibility_real_dimension == 4
    and fixed_real_dimension == 4
)

print("PROGRESS: 7/7 classify")

failed_checks = [
    name
    for name, passed in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "g1800_chiral_torsor_induces_two_distinct_local_"
    "Hilbert_presentations_with_four_real_dimensional_"
    "consistency_kernel"
    if audit_pass
    else
    "joint_chiral_tensor_consistency_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_joint_chiral_tensor_consistency_003",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "source_audit":
        "epr_native_deck_face_chiral_transduction_002",

    "local_structure": {
        "real_dimension":
            4,

        "complex_dimension":
            2,

        "complex_structure":
            "J_F",

        "anti_complex_reflection":
            "S1 = sigma_x K",
    },

    "tensor_sheets": {
        "same_chirality": {
            "relation":
                "J_A tensor I - I tensor J_B",

            "relation_rank":
                rank_plus,

            "real_dimension":
                q_plus_dim,

            "complex_dimension":
                q_plus_dim // 2,
        },

        "opposite_chirality": {
            "relation":
                "J_A tensor I + I tensor J_B",

            "relation_rank":
                rank_minus,

            "real_dimension":
                q_minus_dim,

            "complex_dimension":
                q_minus_dim // 2,
        },
    },

    "one_sided_transduction": {
        "Alice":
            "S1_A tensor I maps Q_plus to Q_minus",

        "Bob":
            "I tensor S1_B maps Q_plus to Q_minus",

        "full_tensor_presentations_identical":
            False,

        "difference_induced_rank":
            difference_induced_rank,

        "consistency_kernel_real_dimension":
            compatibility_real_dimension,

        "compatibility_space_type":
            "real fixed locus of C_AB",

        "complex_dimension":
            None,

        "J_invariant":
            False,

        "projective_ray_space":
            "RP^3 inside CP^3",
    },

    "two_sided_structure": {
        "operator":
            "S1_A tensor S1_B",

        "square":
            "+I",

        "complex_character":
            "anti-complex on Q_plus",

        "fixed_real_dimension":
            fixed_real_dimension,

        "antifixed_real_dimension":
            antifixed_real_dimension,
    },

    "earned_statement": (
        "The native G1800 chiral torsor does not automatically "
        "promote to equality of the two one-sided S1 actions on "
        "the entire two-face Hilbert tensor product. Each local "
        "S1 exchanges the same-chirality and opposite-chirality "
        "tensor sheets. Their disagreement factors exactly as "
        "SA(I-C_AB), where C_AB=SA SB. Therefore the states on "
        "which the Alice-side and Bob-side presentations agree "
        "are exactly the fixed locus of C_AB on Q_plus. This "
        "compatibility locus has real dimension four. Because "
        "C_AB anticommutes with the joint complex structure, "
        "the fixed locus is a real form rather than a complex "
        "two-dimensional subspace. Its projective rays form "
        "RP^3 inside the ambient CP^3. A further native "
        "preparation condition is required to select a "
        "distinguished projective ray."
    ),

    "checks":
        checks,

    "boundary": {
        "point_level_tau_equality_does_not_imply_full_Hilbert_operator_equality":
            True,

        "compatibility_fixed_locus_is_not_complex_linear":
            True,

        "compatible_projective_ray_space_is_RP3":
            True,

        "compatibility_kernel_is_not_yet_a_unique_joint_line":
            True,

        "compatibility_kernel_is_not_yet_a_Bell_resource":
            True,

        "no_singlet_inserted":
            True,

        "no_Born_rule_assumed":
            True,

        "no_probability_law":
            True,

        "no_CHSH_claim":
            True,
    },

    "next_gate": (
        "Search the exact RP3 compatibility locus for a native "
        "preparation selector. Test whether registered joint "
        "orientation, deck history, or local closure selects a "
        "distinguished projective ray without using Bell score, "
        "singlet resemblance, or Born weighting as the selector."
    ),
}

artifact["artifact_sha256"] = digest_json(artifact)

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

note = f"""# Joint chiral tensor consistency 003

## Result

Audit pass:

    {audit_pass}

Program 02 already established that the native G1800 relative-lift
involution built from deck element a has local Program-01 face action

    S1 = sigma_x K.

The present audit asks whether the two point-level presentations

    [a u, v]
    [u, a v]

automatically become the same operator on the full two-face Hilbert
tensor product.

They do not.

## Two chiral tensor sheets

Let V be the four-real-dimensional Program-01 face doublet with complex
structure J.

The same-orientation tensor sheet is

    Q_plus =
      (V_A tensor_R V_B)
      /
      <J_A u tensor v - u tensor J_B v>.

The opposite-orientation sheet is

    Q_minus =
      (V_A tensor_R V_B)
      /
      <J_A u tensor v + u tensor J_B v>.

Both have real dimension 8, hence complex dimension 4.

A one-sided S1 exchanges these sheets.

## Local presentation comparison

The Alice-side map

    S1_A tensor I

and Bob-side map

    I tensor S1_B

both map Q_plus to Q_minus.

However, they are not identical on the full quotient.

Their induced difference has real rank

    {difference_induced_rank}.

The states on which the two local presentations agree form a kernel of
real dimension

    {compatibility_real_dimension},

It is not a complex subspace.

The two-sided anti-complex involution

    C_AB = S1_A tensor S1_B

has this compatibility space as its exact fixed real form.

Because C_AB anticommutes with the joint complex structure,
multiplication by i sends the fixed real form to the anti-fixed
real form.

Projectively, the compatible rays form

    RP^3 inside CP^3.

## Meaning

Audit 002 established a point-level joint chiral torsor.

Audit 003 shows that Hilbert coherence is an additional condition.

The native equality

    [a u,v] = [u,a v]

does not by itself license the operator identity

    S1_A tensor I = I tensor S1_B

on every joint Hilbert state.

Instead it exposes a smaller compatibility space on which the two local
presentations can agree.

That space is now the next native preparation target.

## Boundary

No unique joint projective line is selected here.

No singlet has been inserted.

No probability or frequency law is assumed.

No Bell violation or no-signaling claim is made.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print("Q_PLUS_REAL_DIM:", q_plus_dim)
print("Q_MINUS_REAL_DIM:", q_minus_dim)
print(
    "LOCAL_PRESENTATION_DIFFERENCE_INDUCED_RANK:",
    difference_induced_rank,
)
print(
    "COMPATIBILITY_KERNEL_REAL_DIM:",
    compatibility_real_dimension,
)
print(
    "COMPATIBILITY_SPACE_TYPE:",
    "real_fixed_locus",
)
print(
    "COMPATIBILITY_SPACE_J_INVARIANT:",
    False,
)
print(
    "PROJECTIVE_COMPATIBILITY_SPACE:",
    "RP^3",
)
print(
    "TWO_SIDED_FIXED_REAL_DIM:",
    fixed_real_dimension,
)
print(
    "TWO_SIDED_ANTIFIXED_REAL_DIM:",
    antifixed_real_dimension,
)
print("FAILED_CHECK_COUNT:", len(failed_checks))
print("FAILED_CHECKS:", failed_checks)
print("JSON_OUT:", JSON_OUT)
print("NOTE_OUT:", NOTE_OUT)
print(
    "JSON_SHA256:",
    sha256(JSON_OUT.read_bytes()).hexdigest(),
)
