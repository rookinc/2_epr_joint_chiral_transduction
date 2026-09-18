#!/usr/bin/env python3

from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json

HERE = Path(__file__).resolve().parents[2]

A004 = (
    HERE
    / "artifacts/json"
    / "epr_factor_exchange_odd_ray_frontier_004.v1.json"
)

A006 = (
    HERE
    / "artifacts/json"
    / "epr_canonical_joint_projector_identity_006.v1.json"
)

A010 = (
    HERE
    / "artifacts/json"
    / "epr_g1800_joint_boundary_incidence_quotient_010.v1.json"
)

A011 = (
    HERE
    / "artifacts/json"
    / "epr_native_complex_balanced_boundary_RP3_011.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_native_boundary_antisymmetrizer_closure_012.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_native_boundary_antisymmetrizer_closure_012.md"
)


def load(path):
    return json.loads(path.read_text())


def eye(n):
    return [
        [1 if i == j else 0 for j in range(n)]
        for i in range(n)
    ]


def zeros(r, c):
    return [
        [0 for _ in range(c)]
        for _ in range(r)
    ]


def matmul(a, b):
    out = zeros(len(a), len(b[0]))

    for i in range(len(a)):
        for k in range(len(b)):
            if a[i][k] == 0:
                continue

            for j in range(len(b[0])):
                out[i][j] += a[i][k] * b[k][j]

    return out


def matadd(a, b):
    return [
        [
            a[i][j] + b[i][j]
            for j in range(len(a[0]))
        ]
        for i in range(len(a))
    ]


def matscale(s, a):
    return [
        [s * x for x in row]
        for row in a
    ]


def matsub(a, b):
    return matadd(
        a,
        matscale(-1, b),
    )


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


def nullspace(matrix):
    a = [
        [Fraction(x) for x in row]
        for row in matrix
    ]

    rows = len(a)
    cols = len(a[0])

    pivot_cols = []
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

        pivot_cols.append(c)
        r += 1

        if r == rows:
            break

    free_cols = [
        c
        for c in range(cols)
        if c not in pivot_cols
    ]

    basis = []

    for free in free_cols:
        v = [
            Fraction(0)
            for _ in range(cols)
        ]

        v[free] = Fraction(1)

        for row_index, pivot_col in reversed(
            list(enumerate(pivot_cols))
        ):
            total = sum(
                a[row_index][j] * v[j]
                for j in range(pivot_col + 1, cols)
            )

            v[pivot_col] = -total

        basis.append(v)

    return basis


def columns(vectors):
    if not vectors:
        return []

    return [
        [
            vectors[j][i]
            for j in range(len(vectors))
        ]
        for i in range(len(vectors[0]))
    ]


def get_col(matrix, j):
    return [
        [matrix[i][j]]
        for i in range(len(matrix))
    ]


def independent_column_basis(matrix):
    if not matrix:
        return []

    selected = []

    for j in range(len(matrix[0])):
        col = get_col(matrix, j)

        if not selected:
            if rank(col) > 0:
                selected.append(j)
            continue

        current = [
            [
                matrix[i][k]
                for k in selected
            ]
            for i in range(len(matrix))
        ]

        if rank(
            hcat(current, col)
        ) > rank(current):
            selected.append(j)

    return [
        [
            matrix[i][j]
            for j in selected
        ]
        for i in range(len(matrix))
    ]


def matvec(a, v):
    return [
        sum(
            a[i][j] * v[j]
            for j in range(len(v))
        )
        for i in range(len(a))
    ]


def in_span(matrix, vector):
    col = [
        [x]
        for x in vector
    ]

    return (
        rank(
            hcat(matrix, col)
        )
        == rank(matrix)
    )


def equal(a, b):
    return a == b


def digest_json(obj):
    raw = json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("ascii")

    return sha256(raw).hexdigest()


print("PROGRESS: 1/8 load sealed interfaces")

a004 = load(A004)
a006 = load(A006)
a010 = load(A010)
a011 = load(A011)

checks = {}

checks["audit_004_passed"] = (
    a004["audit_pass"] is True
)

checks["audit_006_passed"] = (
    a006["audit_pass"] is True
)

checks["audit_010_passed"] = (
    a010["audit_pass"] is True
)

checks["audit_011_passed"] = (
    a011["audit_pass"] is True
)

checks["audit_011_native_RP3"] = (
    a011[
        "result"
    ]["projective_space"]
    == "RP^3"
)

print("PROGRESS: 2/8 construct local Program-01 operators")

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
I16 = eye(16)

C = kron(S, S)

JI = kron(J, I4)
IJ = kron(I4, J)

R = matsub(
    JI,
    IJ,
)

checks["C_square_identity"] = equal(
    matmul(C, C),
    I16,
)

checks["balancing_relation_rank_8"] = (
    rank(R) == 8
)

print("PROGRESS: 3/8 reconstruct native fixed and relation spaces")

C_MINUS_I = matsub(
    C,
    I16,
)

F_basis = nullspace(
    C_MINUS_I
)

F = columns(
    F_basis
)

F_dim = len(
    F_basis
)

checks["native_simultaneous_S1_fixed_dim_8"] = (
    F_dim == 8
)

constraint = matmul(
    C_MINUS_I,
    R,
)

preimage_basis = nullspace(
    constraint
)

preimage_matrix = columns(
    preimage_basis
)

relation_fixed_images = matmul(
    R,
    preimage_matrix,
)

L = independent_column_basis(
    relation_fixed_images
)

L_dim = rank(L)

checks["fixed_relation_dim_4"] = (
    L_dim == 4
)

N_dim = (
    F_dim
    - L_dim
)

checks["native_balanced_boundary_dim_4"] = (
    N_dim == 4
)

print("PROGRESS: 4/8 construct native factor exchange")

# Swap the two local four-real-dimensional factors.

X = zeros(16, 16)

for i in range(4):
    for j in range(4):
        source = i * 4 + j
        target = j * 4 + i
        X[target][source] = 1

checks["factor_exchange_square_identity"] = equal(
    matmul(X, X),
    I16,
)

checks["factor_exchange_commutes_with_C"] = equal(
    matmul(X, C),
    matmul(C, X),
)

XR = matmul(
    X,
    R,
)

checks["factor_exchange_preserves_balancing_relation_space"] = (
    rank(
        hcat(
            R,
            XR,
        )
    )
    == rank(R)
)

XF = matmul(
    X,
    F,
)

checks["factor_exchange_preserves_native_fixed_space"] = (
    rank(
        hcat(
            F,
            XF,
        )
    )
    == rank(F)
)

XL = matmul(
    X,
    L,
)

checks["factor_exchange_preserves_fixed_relation_space"] = (
    rank(
        hcat(
            L,
            XL,
        )
    )
    == rank(L)
)

print("PROGRESS: 5/8 construct descended antisymmetrizer")

# Use A = I-X to avoid fractions.
#
# P_boundary = A/2.
#
# A^2 = 2A is equivalent to P_boundary^2=P_boundary.

A = matsub(
    I16,
    X,
)

checks["scaled_antisymmetrizer_identity"] = equal(
    matmul(A, A),
    matscale(2, A),
)

AF = matmul(
    A,
    F,
)

induced_image_rank = (
    rank(
        hcat(
            L,
            AF,
        )
    )
    - L_dim
)

induced_kernel_dim = (
    N_dim
    - induced_image_rank
)

checks["descended_antisymmetrizer_rank_1"] = (
    induced_image_rank == 1
)

checks["descended_antisymmetrizer_kernel_dim_3"] = (
    induced_kernel_dim == 3
)

print("PROGRESS: 6/8 identify canonical image line")

# Raw representative of the Audit-003 fixed-form ray
#
#     i(|01> - |10>)
#
# using the section
#
#     Im(c_ij) represented by J e_i tensor e_j.
#
# Local real basis:
#
#     e0, J e0, e1, J e1.
#
# Therefore:
#
#     + J e0 tensor e1
#     - J e1 tensor e0

phi = [
    0
    for _ in range(16)
]

phi[1 * 4 + 2] = 1
phi[3 * 4 + 0] = -1

checks["canonical_raw_representative_is_C_fixed"] = (
    matvec(C, phi)
    == phi
)

checks["canonical_raw_representative_not_killed_by_relation"] = (
    not in_span(
        L,
        phi,
    )
)

Xphi = matvec(
    X,
    phi,
)

checks["canonical_ray_is_exchange_odd_mod_relation"] = (
    in_span(
        L,
        [
            Xphi[i] + phi[i]
            for i in range(16)
        ],
    )
)

Aphi = matvec(
    A,
    phi,
)

checks["antisymmetrizer_selects_canonical_ray_mod_relation"] = (
    in_span(
        L,
        [
            Aphi[i] - 2 * phi[i]
            for i in range(16)
        ],
    )
)

checks["induced_image_is_exactly_canonical_line"] = (
    rank(
        hcat(
            L,
            AF,
            [
                [x]
                for x in phi
            ],
        )
    )
    - L_dim
    == 1
)

print("PROGRESS: 7/8 classify null branch")

# Any compatible boundary state in the 3-real-dimensional kernel
# maps to zero under P_boundary.
#
# This null branch is part of the complete closure map and may
# not be postselected away if the operation is later promoted to
# an experimental preparation/instrument.

checks["null_branch_real_dimension_3"] = (
    induced_kernel_dim == 3
)

null_fraction_dimension_only = {
    "kernel_real_dimension":
        induced_kernel_dim,

    "domain_real_dimension":
        N_dim,

    "warning":
        (
            "Dimension ratio is not a probability. "
            "No measure or Born weighting is assumed."
        ),
}

print("PROGRESS: 8/8 classify")

failed_checks = [
    name
    for name, passed in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "native_balanced_boundary_admits_setting_independent_"
    "rank_one_factor_exchange_antisymmetrizer_closure_"
    "with_explicit_three_dimensional_null_branch"
    if audit_pass
    else
    "native_boundary_antisymmetrizer_closure_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_native_boundary_antisymmetrizer_closure_012",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "native_boundary_space": {
        "source":
            "Audit 011 native balanced boundary",

        "real_dimension":
            N_dim,

        "projective_space":
            "RP^3",
    },

    "closure": {
        "factor_exchange":
            "X",

        "operator":
            "P_boundary = (I-X)/2",

        "scaled_operator":
            "A = I-X",

        "idempotence":
            "A^2 = 2A",

        "setting_independent":
            True,

        "induced_real_rank":
            induced_image_rank,

        "kernel_real_dimension":
            induced_kernel_dim,

        "image":
            "canonical Audit-006 joint line",

        "image_representative":
            "i(|01>-|10>)",
    },

    "null_branch":
        null_fraction_dimension_only,

    "earned_statement": (
        "The native G1800 balanced boundary space carries a "
        "canonical setting-independent factor-exchange action. "
        "Its antisymmetrizer P_boundary=(I-X)/2 descends through "
        "both the simultaneous-S1 boundary quotient and the "
        "Program-01 complex-balancing relation. On the resulting "
        "four-real-dimensional RP3 compatibility space it is an "
        "idempotent real-linear map of rank one. Its nonzero "
        "image is exactly the canonical EPR line identified in "
        "Audits 006-007. Its kernel has real dimension three, "
        "giving an explicit null branch that must remain visible "
        "in any later operational preparation protocol."
    ),

    "checks":
        checks,

    "boundary": {
        "finite_boundary_closure_operator_constructed":
            True,

        "closure_is_setting_independent":
            True,

        "closure_image_is_canonical_line":
            True,

        "source_dynamics_do_not_yet_prove_closure_is_executed":
            True,

        "null_branch_must_not_be_postselected_away":
            True,

        "dimension_ratio_not_interpreted_as_probability":
            True,

        "no_Born_rule":
            True,

        "no_trial_frequency_law":
            True,

        "no_no_signaling_claim":
            True,

        "no_CHSH_claim":
            True,
    },

    "next_gate": (
        "Determine whether the native source/common preparation "
        "actually executes or guarantees this boundary "
        "antisymmetrizer closure before analyzer settings are "
        "chosen. If not, retain the three-dimensional null/even "
        "branch explicitly. A preparation theorem must specify "
        "how trials enter the nonzero image without using "
        "setting-dependent or outcome-dependent postselection."
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

note = f"""# Native boundary antisymmetrizer closure 012

## Result

Audit pass:

    {audit_pass}

Audit 011 produced the native complex-balanced boundary compatibility
space

    RP^3

with underlying real dimension

    {N_dim}.

Factor exchange X descends exactly to this space.

Define

    P_boundary = (I-X)/2.

Equivalently use the integer-scaled operator

    A = I-X,

which satisfies

    A^2 = 2A.

Therefore P_boundary is idempotent.

## Rank

On the native balanced boundary space,

    rank_R(P_boundary) = {induced_image_rank}

and

    dim_R ker(P_boundary) = {induced_kernel_dim}.

The nonzero image is exactly the canonical projective line

    [|01>-|10>]

with Audit-003 Real-form representative

    i(|01>-|10>).

Thus the canonical line is now the image of an explicit finite
setting-independent boundary closure operator.

## Null branch

The closure has a three-real-dimensional kernel.

That kernel is an explicit null/even branch.

If this closure is later promoted to an operational source-preparation
rule, null trials may not be discarded after analyzer settings or
outcomes are known.

No dimension ratio is interpreted as a probability.

## Remaining gap

The operator exists canonically on the native boundary apparatus.

What is not yet proved is that the native source dynamics actually
execute or guarantee this closure before Alice and Bob choose their
settings.

That is now the preparation theorem frontier.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print(
    "NATIVE_BALANCED_BOUNDARY_REAL_DIM:",
    N_dim,
)
print(
    "ANTISYMMETRIZER_INDUCED_RANK:",
    induced_image_rank,
)
print(
    "ANTISYMMETRIZER_KERNEL_REAL_DIM:",
    induced_kernel_dim,
)
print(
    "IMAGE_IS_CANONICAL_LINE:",
    checks[
        "induced_image_is_exactly_canonical_line"
    ],
)
print(
    "NULL_BRANCH_EXPLICIT:",
    True,
)
print("FAILED_CHECK_COUNT:", len(failed_checks))
print("FAILED_CHECKS:", failed_checks)
print("JSON_OUT:", JSON_OUT)
print("NOTE_OUT:", NOTE_OUT)
print(
    "JSON_SHA256:",
    sha256(JSON_OUT.read_bytes()).hexdigest(),
)
