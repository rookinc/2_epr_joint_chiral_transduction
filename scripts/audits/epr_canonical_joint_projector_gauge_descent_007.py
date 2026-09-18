#!/usr/bin/env python3

from fractions import Fraction
from hashlib import sha256
from itertools import permutations, product
from pathlib import Path
import json

HERE = Path(__file__).resolve().parents[2]

A006 = (
    HERE
    / "artifacts/json"
    / "epr_canonical_joint_projector_identity_006.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_canonical_joint_projector_gauge_descent_007.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_canonical_joint_projector_gauge_descent_007.md"
)


def load(path):
    return json.loads(path.read_text())


def eye(n):
    return [
        [1 + 0j if i == j else 0 + 0j for j in range(n)]
        for i in range(n)
    ]


def zeros(r, c):
    return [
        [0 + 0j for _ in range(c)]
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


def dagger(a):
    return [
        [
            a[j][i].conjugate()
            for j in range(len(a))
        ]
        for i in range(len(a[0]))
    ]


def equal(a, b):
    return a == b


def digest_json(obj):
    raw = json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("ascii")

    return sha256(raw).hexdigest()


print("PROGRESS: 1/7 load Audit 006")

a006 = load(A006)

checks = {}

checks["audit_006_passed"] = (
    a006["audit_pass"] is True
)

checks["audit_006_exact_projector_identity"] = (
    a006["checks"]["exact_projector_identity"]
    is True
)

print("PROGRESS: 2/7 reconstruct canonical projector")

I2 = eye(2)
I4 = eye(4)

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

KS = [K1, K2, K3]

X = [
    [1 + 0j, 0 + 0j, 0 + 0j, 0 + 0j],
    [0 + 0j, 0 + 0j, 1 + 0j, 0 + 0j],
    [0 + 0j, 1 + 0j, 0 + 0j, 0 + 0j],
    [0 + 0j, 0 + 0j, 0 + 0j, 1 + 0j],
]

P = matscale(
    Fraction(1, 2),
    matsub(I4, X),
)

Qsum = I4

for K in KS:
    Qsum = matadd(
        Qsum,
        kron(K, K),
    )

Q = matscale(
    Fraction(1, 4),
    Qsum,
)

checks["reconstructed_projectors_equal"] = (
    P == Q
)

print("PROGRESS: 3/7 test all S3 axis permutations")

s3_failure_count = 0

for perm in permutations(range(3)):
    qsum = I4

    for i in perm:
        qsum = matadd(
            qsum,
            kron(KS[i], KS[i]),
        )

    qp = matscale(
        Fraction(1, 4),
        qsum,
    )

    if qp != P:
        s3_failure_count += 1

checks["all_six_axis_permutations_preserve_projector"] = (
    s3_failure_count == 0
)

print("PROGRESS: 4/7 test independent axis sign conventions")

# Ki -> +/- Ki does not alter Ki tensor Ki.
#
# This covers all eight independent sign choices, including
# orientation reversals of the displayed quaternionic frame.

sign_failure_count = 0

for signs in product((-1, 1), repeat=3):
    qsum = I4

    for sign, K in zip(signs, KS):
        KP = matscale(sign, K)

        qsum = matadd(
            qsum,
            kron(KP, KP),
        )

    qp = matscale(
        Fraction(1, 4),
        qsum,
    )

    if qp != P:
        sign_failure_count += 1

checks["all_eight_axis_sign_choices_preserve_projector"] = (
    sign_failure_count == 0
)

print("PROGRESS: 5/7 test representative common U2 frame changes")

# Exact representative common U(2) transformations.
#
# Since the exchange projector commutes with U tensor U for
# every U, these finite exact examples are regression witnesses.
# The general theorem is algebraic and recorded below.

U_LIST = {
    "identity": [
        [1 + 0j, 0 + 0j],
        [0 + 0j, 1 + 0j],
    ],

    "sigma_x": [
        [0 + 0j, 1 + 0j],
        [1 + 0j, 0 + 0j],
    ],

    "phase_z": [
        [0 + 1j, 0 + 0j],
        [0 + 0j, 0 - 1j],
    ],

    "quarter_turn": [
        [0 + 0j, 1 + 0j],
        [-1 + 0j, 0 + 0j],
    ],
}

u2_failure_count = 0

for name, U in U_LIST.items():
    UU = kron(U, U)

    transformed = matmul(
        matmul(
            UU,
            P,
        ),
        dagger(UU),
    )

    if transformed != P:
        u2_failure_count += 1

checks["representative_common_U2_changes_preserve_projector"] = (
    u2_failure_count == 0
)

print("PROGRESS: 6/7 verify general common-frame commutator identity")

# The swap satisfies
#
#     X (A tensor B) = (B tensor A) X.
#
# Hence for A=B=U,
#
#     X (U tensor U) = (U tensor U) X.
#
# This is independent of matrix entries.
#
# We encode the theorem as an exact symbolic contract and verify
# the finite basis action of X used in the derivation.

basis_pairs = [
    (0, 0),
    (0, 1),
    (1, 0),
    (1, 1),
]

swap_basis_exact = True

for col, (i, j) in enumerate(basis_pairs):
    expected = basis_pairs.index((j, i))

    for row in range(4):
        want = (
            1 + 0j
            if row == expected
            else 0 + 0j
        )

        if X[row][col] != want:
            swap_basis_exact = False

checks["factor_exchange_is_exact_tensor_swap"] = (
    swap_basis_exact
)

checks["general_common_U2_invariance_theorem_applies"] = (
    checks["factor_exchange_is_exact_tensor_swap"]
    and checks["reconstructed_projectors_equal"]
)

print("PROGRESS: 7/7 classify")

failed_checks = [
    name
    for name, passed in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "canonical_joint_projector_descends_through_S3_axis_"
    "ambiguity_sign_conventions_and_common_U2_frame_gauge"
    if audit_pass
    else
    "canonical_joint_projector_gauge_descent_failed"
)

artifact = {
    "artifact_id":
        "epr_canonical_joint_projector_gauge_descent_007",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "projector": {
        "exchange_form":
            "P_joint = (I-X)/2",

        "quaternionic_form":
            (
                "P_joint = "
                "(I + sum_i Ki tensor Ki)/4"
            ),

        "complex_rank":
            1,
    },

    "gauge_descent": {
        "S3_axis_permutation_failure_count":
            s3_failure_count,

        "axis_sign_failure_count":
            sign_failure_count,

        "representative_common_U2_failure_count":
            u2_failure_count,

        "common_U2_general_reason":
            (
                "X commutes with U tensor U for every common "
                "local frame change U, so any polynomial in X "
                "including P_joint=(I-X)/2 is common-U2 "
                "invariant."
            ),
    },

    "earned_statement": (
        "The canonical rank-one joint projector from Audit 006 "
        "does not depend on the unresolved Program-01 S3 "
        "axis-slot crosswalk, on quaternionic axis sign "
        "conventions, or on a common local U(2) frame choice. "
        "Its exchange form P=(I-X)/2 makes common-U2 "
        "invariance exact, while its quaternionic form makes "
        "axis-frame independence explicit. Therefore the "
        "canonical joint line descends through the remaining "
        "local presentation freedoms of Program 01."
    ),

    "checks":
        checks,

    "boundary": {
        "projector_is_interface_canonical":
            True,

        "face_Weyl_axis_crosswalk_not_needed":
            True,

        "source_preparation_execution_not_yet_derived":
            True,

        "projector_not_yet_declared_physical_operation":
            True,

        "no_Born_rule":
            True,

        "no_probability_law":
            True,

        "no_CHSH_claim":
            True,
    },

    "next_gate": (
        "Construct a setting-independent native preparation "
        "operation or source-sector rule whose image is the "
        "gauge-descended canonical line. The remaining problem "
        "is execution, not state selection."
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

note = f"""# Canonical joint projector gauge descent 007

## Result

Audit pass:

    {audit_pass}

Audit 006 established the exact identity

    P_joint
      =
    (I-X)/2

and

    P_joint
      =
    (I + K1 tensor K1
       + K2 tensor K2
       + K3 tensor K3)/4.

The present audit tests whether this apparently canonical projector still
depends on local presentation choices.

It does not.

## Axis-slot ambiguity

All six permutations of

    K1, K2, K3

leave the quaternionic expression unchanged.

Therefore the unresolved Program-01 S3 face/Weyl axis-slot crosswalk
cannot move P_joint.

## Axis signs

All eight independent sign choices

    Ki -> +/- Ki

also leave the projector unchanged because every term occurs as

    Ki tensor Ki.

## Common local U(2) frame

The exchange operator satisfies

    X (U tensor U)
      =
    (U tensor U) X

for every common local frame change U.

Therefore

    (U tensor U)
    P_joint
    (U tensor U)^dagger
      =
    P_joint.

The canonical line is consequently independent of the common Schur phase,
common SU(2) frame, and the unresolved naming of the three quaternionic
axes.

## Meaning

Program 01 may be consumed exactly at its earned boundary.

Program 02 does not need to solve the historical face/Weyl S3 axis
crosswalk in order to define the joint preparation line.

The remaining preparation problem is now execution:

    what native setting-independent operation or source rule
    places the preparation on this line?

No probability law is introduced here.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print(
    "S3_AXIS_PERMUTATION_FAILURE_COUNT:",
    s3_failure_count,
)
print(
    "AXIS_SIGN_FAILURE_COUNT:",
    sign_failure_count,
)
print(
    "COMMON_U2_REPRESENTATIVE_FAILURE_COUNT:",
    u2_failure_count,
)
print(
    "GENERAL_COMMON_U2_INVARIANCE:",
    checks[
        "general_common_U2_invariance_theorem_applies"
    ],
)
print(
    "FACE_WEYL_S3_CROSSWALK_REQUIRED:",
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
