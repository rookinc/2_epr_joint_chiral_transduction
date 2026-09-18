#!/usr/bin/env python3

from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json

HERE = Path(__file__).resolve().parents[2]

A003 = (
    HERE
    / "artifacts/json"
    / "epr_joint_chiral_tensor_consistency_003.v1.json"
)

A004 = (
    HERE
    / "artifacts/json"
    / "epr_factor_exchange_odd_ray_frontier_004.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_diagonal_quaternionic_ray_selector_005.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_diagonal_quaternionic_ray_selector_005.md"
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


def matvec(a, v):
    return [
        sum(
            a[i][j] * v[j]
            for j in range(len(v))
        )
        for i in range(len(a))
    ]


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


def realify_complex_matrix(m):
    n = len(m)

    out = zeros(2 * n, 2 * n)

    for i in range(n):
        for j in range(n):
            z = m[i][j]

            a = int(z.real)
            b = int(z.imag)

            out[2*i][2*j] = a
            out[2*i][2*j + 1] = -b
            out[2*i + 1][2*j] = b
            out[2*i + 1][2*j + 1] = a

    return out


def vcat(*mats):
    out = []

    for m in mats:
        out.extend(m)

    return out


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


def nullity(matrix):
    return len(matrix[0]) - rank(matrix)


def equal(a, b):
    return a == b


def digest_json(obj):
    raw = json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("ascii")

    return sha256(raw).hexdigest()


print("PROGRESS: 1/7 load sealed interfaces")

a003 = load(A003)
a004 = load(A004)

checks = {}

checks["audit_003_passed"] = (
    a003["audit_pass"] is True
)

checks["audit_004_passed"] = (
    a004["audit_pass"] is True
)

checks["audit_003_compatibility_is_RP3"] = (
    a003[
        "one_sided_transduction"
    ]["projective_ray_space"]
    == "RP^3 inside CP^3"
)

checks["audit_004_unique_odd_ray"] = (
    a004[
        "joint_Hilbert_geometry"
    ]["compatible_odd_real_dimension"]
    == 1
)

print("PROGRESS: 2/7 construct local quaternionic doublet")

# Fundamental SU(2) quaternionic units:
#
# K1 = i sigma_x
# K2 = -i sigma_y
# K3 = i sigma_z
#
# Each squares to -I and
#
# K1 K2 = K3,
# K2 K3 = K1,
# K3 K1 = K2.
#
# These are one concrete event-adapted realization of the
# Program-01 three-axis quaternionic algebra. The result
# below is invariant under permutation of the three axes.

K1 = [
    [0 + 0j, 0 + 1j],
    [0 + 1j, 0 + 0j],
]

K2 = [
    [0 + 0j, -1 + 0j],
    [1 + 0j, 0 + 0j],
]

K3 = [
    [0 + 1j, 0 + 0j],
    [0 + 0j, 0 - 1j],
]


def cmul(a, b):
    out = [
        [0j for _ in range(len(b[0]))]
        for _ in range(len(a))
    ]

    for i in range(len(a)):
        for k in range(len(b)):
            for j in range(len(b[0])):
                out[i][j] += a[i][k] * b[k][j]

    return out


I2c = [
    [1 + 0j, 0 + 0j],
    [0 + 0j, 1 + 0j],
]

minus_I2c = [
    [-1 + 0j, 0 + 0j],
    [0 + 0j, -1 + 0j],
]

checks["K1_square_minus_I"] = (
    cmul(K1, K1) == minus_I2c
)

checks["K2_square_minus_I"] = (
    cmul(K2, K2) == minus_I2c
)

checks["K3_square_minus_I"] = (
    cmul(K3, K3) == minus_I2c
)

checks["K1_K2_equals_K3"] = (
    cmul(K1, K2) == K3
)

checks["K2_K3_equals_K1"] = (
    cmul(K2, K3) == K1
)

checks["K3_K1_equals_K2"] = (
    cmul(K3, K1) == K2
)

print("PROGRESS: 3/7 construct diagonal quaternionic actions")

D1c = kron(K1, K1)
D2c = kron(K2, K2)
D3c = kron(K3, K3)

D1 = realify_complex_matrix(D1c)
D2 = realify_complex_matrix(D2c)
D3 = realify_complex_matrix(D3c)

I8 = eye(8)

checks["D1_square_identity"] = equal(
    matmul(D1, D1),
    I8,
)

checks["D2_square_identity"] = equal(
    matmul(D2, D2),
    I8,
)

checks["D3_square_identity"] = equal(
    matmul(D3, D3),
    I8,
)

checks["diagonal_actions_commute"] = (
    equal(
        matmul(D1, D2),
        matmul(D2, D1),
    )
    and equal(
        matmul(D2, D3),
        matmul(D3, D2),
    )
    and equal(
        matmul(D3, D1),
        matmul(D1, D3),
    )
)

print("PROGRESS: 4/7 reconstruct Audit-003 Real form")

# C(M) = sigma_x conjugate(M) sigma_x
#
# coefficient order:
#
# Re c00, Im c00,
# Re c01, Im c01,
# Re c10, Im c10,
# Re c11, Im c11.

C = [
    [0, 0, 0, 0, 0, 0, 1,  0],
    [0, 0, 0, 0, 0, 0, 0, -1],

    [0, 0, 0, 0, 1,  0, 0,  0],
    [0, 0, 0, 0, 0, -1, 0,  0],

    [0, 0, 1,  0, 0, 0, 0,  0],
    [0, 0, 0, -1, 0, 0, 0,  0],

    [1,  0, 0, 0, 0, 0, 0,  0],
    [0, -1, 0, 0, 0, 0, 0,  0],
]

C_MINUS_I = matsub(C, I8)

compat_dim = nullity(
    C_MINUS_I
)

checks["compatibility_real_dimension_4"] = (
    compat_dim == 4
)

checks["each_diagonal_axis_preserves_real_structure"] = (
    equal(
        matmul(C, D1),
        matmul(D1, C),
    )
    and equal(
        matmul(C, D2),
        matmul(D2, C),
    )
    and equal(
        matmul(C, D3),
        matmul(D3, C),
    )
)

print("PROGRESS: 5/7 measure common invariant locus")

D1_MINUS_I = matsub(D1, I8)
D2_MINUS_I = matsub(D2, I8)
D3_MINUS_I = matsub(D3, I8)

common_diag_fixed_dim = nullity(
    vcat(
        D1_MINUS_I,
        D2_MINUS_I,
        D3_MINUS_I,
    )
)

compatible_common_diag_fixed_dim = nullity(
    vcat(
        C_MINUS_I,
        D1_MINUS_I,
        D2_MINUS_I,
        D3_MINUS_I,
    )
)

checks["common_diagonal_quaternionic_fixed_real_dim_2"] = (
    common_diag_fixed_dim == 2
)

checks["compatible_common_diagonal_fixed_real_dim_1"] = (
    compatible_common_diag_fixed_dim == 1
)

print("PROGRESS: 6/7 compare with exchange-odd ray")

# C-fixed representative of the alternating projective ray:
#
#     i(|01> - |10>)

odd = [
    0, 0,
    0, 1,
    0, -1,
    0, 0,
]

checks["alternating_ray_is_C_fixed"] = (
    matvec(C, odd) == odd
)

checks["alternating_ray_fixed_by_D1"] = (
    matvec(D1, odd) == odd
)

checks["alternating_ray_fixed_by_D2"] = (
    matvec(D2, odd) == odd
)

checks["alternating_ray_fixed_by_D3"] = (
    matvec(D3, odd) == odd
)

checks["common_frame_neutral_ray_matches_004_odd_ray"] = all([
    checks["alternating_ray_is_C_fixed"],
    checks["alternating_ray_fixed_by_D1"],
    checks["alternating_ray_fixed_by_D2"],
    checks["alternating_ray_fixed_by_D3"],
    compatible_common_diag_fixed_dim == 1,
])

print("PROGRESS: 7/7 classify")

failed_checks = [
    name
    for name, passed in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "within_the_earned_RP3_the_unique_ray_neutral_under_"
    "all_three_diagonal_quaternionic_face_axes_is_exactly_"
    "the_exchange_odd_alternating_ray"
    if audit_pass
    else
    "diagonal_quaternionic_ray_selector_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_diagonal_quaternionic_ray_selector_005",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "input_geometry": {
        "ambient_joint_projective_space":
            "CP^3",

        "Audit003_compatible_locus":
            "RP^3",

        "Audit004_exchange_odd_frontier":
            "one compatible projective ray",
    },

    "finite_face_frame": {
        "local_generators":
            [
                "K1 = i sigma_x",
                "K2 = -i sigma_y",
                "K3 = i sigma_z",
            ],

        "algebra":
            "quaternionic fundamental SU(2) axes",

        "axis_slot_S3_invariant":
            True,

        "joint_actions":
            [
                "K1 tensor K1",
                "K2 tensor K2",
                "K3 tensor K3",
            ],
    },

    "dimensions": {
        "Audit003_compatibility_real_dimension":
            compat_dim,

        "common_diagonal_quaternionic_fixed_real_dimension":
            common_diag_fixed_dim,

        "compatible_common_diagonal_fixed_real_dimension":
            compatible_common_diag_fixed_dim,

        "compatible_common_diagonal_projective_ray_count":
            (
                1
                if compatible_common_diag_fixed_dim == 1
                else None
            ),
    },

    "ray": {
        "representative":
            "i(|01>-|10>)",

        "projective_form":
            "[|01>-|10>]",

        "exchange_parity":
            "odd",

        "product_status":
            "non-product",
    },

    "earned_statement": (
        "The three quaternionic Program-01 face axes provide a "
        "finite common-frame neutrality test independent of the "
        "unresolved face/Weyl S3 axis-slot crosswalk. On the "
        "full complex joint fiber their simultaneous diagonal "
        "fixed space is one complex line. Intersecting with the "
        "Audit-003 Real compatibility locus leaves one real "
        "line, hence exactly one compatible projective ray. "
        "That ray is the same alternating exchange-odd ray "
        "isolated independently by Audit 004."
    ),

    "checks":
        checks,

    "boundary": {
        "finite_selector_geometry_only":
            True,

        "source_invariance_under_diagonal_face_axes_not_yet_derived":
            True,

        "ray_not_yet_declared_prepared":
            True,

        "no_singlet_probability_law_inserted":
            True,

        "no_Born_rule":
            True,

        "no_trial_weighting":
            True,

        "no_no_signaling_claim":
            True,

        "no_CHSH_claim":
            True,
    },

    "next_gate": (
        "Test whether the actual setting-independent native "
        "source preparation is neutral under the common "
        "quaternionic face-frame action, or independently "
        "selects odd factor-exchange parity. Either native "
        "selector would promote the unique ray from a symmetry "
        "candidate to a prepared joint state. Do not infer the "
        "selector from the ray's Bell properties."
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

note = f"""# Diagonal quaternionic ray selector 005

## Result

Audit pass:

    {audit_pass}

Program 01 supplies three quaternionic derived face axes.

Using one standard fundamental-doublet realization,

    K1 = i sigma_x
    K2 = -i sigma_y
    K3 = i sigma_z,

their diagonal joint actions are

    K1 tensor K1
    K2 tensor K2
    K3 tensor K3.

The result is invariant under permutation of the three axis names.

## Common-frame neutrality

On the full complex two-doublet tensor fiber, the common fixed space of
the three diagonal quaternionic actions has real dimension

    {common_diag_fixed_dim},

equivalently one complex line.

Audit 003 supplied a four-real-dimensional compatible Real form.

Intersecting the two gives real dimension

    {compatible_common_diag_fixed_dim}.

Therefore exactly one compatible projective ray is neutral under all
three common face-axis actions.

It is represented by

    i (|01> - |10>),

and projectively by

    [|01> - |10>].

This is exactly the exchange-odd ray independently isolated by Audit 004.

## Meaning

Two independent finite characterizations now converge on the same ray:

1. unique compatible exchange-odd ray;
2. unique compatible common-quaternionic-frame-neutral ray.

Neither characterization is yet a source-preparation theorem.

The remaining question is whether the native source actually carries
either selector before Alice and Bob choose their analyzer settings.

## Boundary

No prepared singlet is declared.

No Born rule is assumed.

No probability table, no-signaling theorem, Bell factorization result,
or CHSH value is claimed.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print(
    "COMPATIBILITY_REAL_DIM:",
    compat_dim,
)
print(
    "COMMON_DIAGONAL_FIXED_REAL_DIM:",
    common_diag_fixed_dim,
)
print(
    "COMPATIBLE_COMMON_DIAGONAL_FIXED_REAL_DIM:",
    compatible_common_diag_fixed_dim,
)
print(
    "UNIQUE_COMPATIBLE_PROJECTIVE_RAY:",
    compatible_common_diag_fixed_dim == 1,
)
print(
    "RAY:",
    "[|01>-|10>]",
)
print(
    "MATCHES_AUDIT004_EXCHANGE_ODD_RAY:",
    checks[
        "common_frame_neutral_ray_matches_004_odd_ray"
    ],
)
print("FAILED_CHECK_COUNT:", len(failed_checks))
print("FAILED_CHECKS:", failed_checks)
print("JSON_OUT:", JSON_OUT)
print("NOTE_OUT:", NOTE_OUT)
print(
    "JSON_SHA256:",
    sha256(JSON_OUT.read_bytes()).hexdigest(),
)
