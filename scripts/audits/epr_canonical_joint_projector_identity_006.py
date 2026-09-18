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

A005 = (
    HERE
    / "artifacts/json"
    / "epr_diagonal_quaternionic_ray_selector_005.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_canonical_joint_projector_identity_006.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_canonical_joint_projector_identity_006.md"
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


def equal(a, b):
    return a == b


def real_integer_matrix(a):
    out = []

    for row in a:
        rr = []

        for z in row:
            if z.imag != 0:
                raise RuntimeError(
                    "matrix expected to be real"
                )

            value = z.real

            if int(value) != value:
                raise RuntimeError(
                    "matrix expected to be integral"
                )

            rr.append(int(value))

        out.append(rr)

    return out


def rank_rational(matrix):
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


def c_real_structure(v):
    # C(M) = sigma_x conjugate(M) sigma_x
    #
    # M coefficient order:
    #
    # c00, c01,
    # c10, c11

    return [
        v[3].conjugate(),
        v[2].conjugate(),
        v[1].conjugate(),
        v[0].conjugate(),
    ]


def digest_json(obj):
    raw = json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("ascii")

    return sha256(raw).hexdigest()


print("PROGRESS: 1/7 load sealed selector audits")

a004 = load(A004)
a005 = load(A005)

checks = {}

checks["audit_004_passed"] = (
    a004["audit_pass"] is True
)

checks["audit_005_passed"] = (
    a005["audit_pass"] is True
)

checks["audit_004_unique_odd_ray"] = (
    a004[
        "joint_Hilbert_geometry"
    ]["compatible_odd_real_dimension"]
    == 1
)

checks["audit_005_unique_neutral_ray"] = (
    a005[
        "dimensions"
    ][
        "compatible_common_diagonal_fixed_real_dimension"
    ]
    == 1
)

print("PROGRESS: 2/7 construct quaternionic doublet")

# Standard oriented quaternionic basis:
#
# K1 = i sigma_x
# K2 = -i sigma_y
# K3 = i sigma_z
#
# K1 K2 = K3
# K2 K3 = K1
# K3 K1 = K2

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

I2 = eye(2)

checks["K1_square_minus_I"] = equal(
    matmul(K1, K1),
    matscale(-1, I2),
)

checks["K2_square_minus_I"] = equal(
    matmul(K2, K2),
    matscale(-1, I2),
)

checks["K3_square_minus_I"] = equal(
    matmul(K3, K3),
    matscale(-1, I2),
)

checks["K1_K2_equals_K3"] = equal(
    matmul(K1, K2),
    K3,
)

checks["K2_K3_equals_K1"] = equal(
    matmul(K2, K3),
    K1,
)

checks["K3_K1_equals_K2"] = equal(
    matmul(K3, K1),
    K2,
)

print("PROGRESS: 3/7 construct factor exchange")

# Tensor basis:
#
# |00>, |01>, |10>, |11>

X = [
    [1 + 0j, 0 + 0j, 0 + 0j, 0 + 0j],
    [0 + 0j, 0 + 0j, 1 + 0j, 0 + 0j],
    [0 + 0j, 1 + 0j, 0 + 0j, 0 + 0j],
    [0 + 0j, 0 + 0j, 0 + 0j, 1 + 0j],
]

I4 = eye(4)

checks["factor_exchange_square_identity"] = equal(
    matmul(X, X),
    I4,
)

A = matsub(
    I4,
    X,
)

checks["scaled_odd_projector_identity"] = equal(
    matmul(A, A),
    matscale(2, A),
)

print("PROGRESS: 4/7 construct quaternionic neutral projector")

D1 = kron(K1, K1)
D2 = kron(K2, K2)
D3 = kron(K3, K3)

checks["D1_square_identity"] = equal(
    matmul(D1, D1),
    I4,
)

checks["D2_square_identity"] = equal(
    matmul(D2, D2),
    I4,
)

checks["D3_square_identity"] = equal(
    matmul(D3, D3),
    I4,
)

checks["D1_D2_equals_D3"] = equal(
    matmul(D1, D2),
    D3,
)

checks["D2_D3_equals_D1"] = equal(
    matmul(D2, D3),
    D1,
)

checks["D3_D1_equals_D2"] = equal(
    matmul(D3, D1),
    D2,
)

Q = matadd(
    matadd(
        I4,
        D1,
    ),
    matadd(
        D2,
        D3,
    ),
)

checks["scaled_quaternionic_projector_identity"] = equal(
    matmul(Q, Q),
    matscale(4, Q),
)

print("PROGRESS: 5/7 compare the two selectors exactly")

# P_odd = (I-X)/2
#
# P_Q = (I+D1+D2+D3)/4
#
# Exact equality is equivalent to
#
# 2(I-X) = I+D1+D2+D3.

checks["exact_projector_identity"] = equal(
    matscale(2, A),
    Q,
)

A_real = real_integer_matrix(A)
Q_real = real_integer_matrix(Q)

odd_rank = rank_rational(
    A_real
)

quaternionic_rank = rank_rational(
    Q_real
)

checks["odd_projector_complex_rank_one"] = (
    odd_rank == 1
)

checks["quaternionic_projector_complex_rank_one"] = (
    quaternionic_rank == 1
)

print("PROGRESS: 6/7 verify canonical ray and Real form")

psi = [
    0 + 0j,
    1 + 0j,
    -1 + 0j,
    0 + 0j,
]

phi = [
    1j * x
    for x in psi
]

checks["psi_exchange_odd"] = (
    matvec(X, psi)
    == [
        -x
        for x in psi
    ]
)

checks["psi_D1_fixed"] = (
    matvec(D1, psi)
    == psi
)

checks["psi_D2_fixed"] = (
    matvec(D2, psi)
    == psi
)

checks["psi_D3_fixed"] = (
    matvec(D3, psi)
    == psi
)

checks["odd_scaled_projector_selects_psi"] = (
    matvec(A, psi)
    == [
        2 * x
        for x in psi
    ]
)

checks["quaternionic_scaled_projector_selects_psi"] = (
    matvec(Q, psi)
    == [
        4 * x
        for x in psi
    ]
)

checks["psi_is_C_antifixed"] = (
    c_real_structure(psi)
    == [
        -x
        for x in psi
    ]
)

checks["i_psi_is_C_fixed"] = (
    c_real_structure(phi)
    == phi
)

print("PROGRESS: 7/7 classify")

failed_checks = [
    name
    for name, passed in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "factor_exchange_odd_projector_equals_common_"
    "quaternionic_neutral_projector_exactly"
    if audit_pass
    else
    "canonical_joint_projector_identity_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_canonical_joint_projector_identity_006",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "projectors": {
        "factor_exchange_odd":
            "P_odd = (I-X)/2",

        "common_quaternionic_neutral":
            (
                "P_Q = "
                "(I + K1xK1 + K2xK2 + K3xK3)/4"
            ),

        "exact_identity":
            "P_odd = P_Q",

        "integer_scaled_identity":
            (
                "2(I-X) = "
                "I + K1xK1 + K2xK2 + K3xK3"
            ),

        "complex_rank":
            1,
    },

    "canonical_line": {
        "projective_representative":
            "[|01>-|10>]",

        "exchange_character":
            "-1",

        "diagonal_quaternionic_character":
            "+1,+1,+1",

        "Audit003_real_form_representative":
            "i(|01>-|10>)",
    },

    "earned_statement": (
        "The factor-exchange antisymmetrizer and the common "
        "quaternionic-frame neutral projector are not merely "
        "two selectors with the same observed ray. They are "
        "exactly the same rank-one complex projector on the "
        "two-doublet joint fiber. Thus exchange oddness and "
        "common quaternionic neutrality define one canonical "
        "setting-independent projective line before any Bell "
        "score or probability law is introduced."
    ),

    "checks":
        checks,

    "boundary": {
        "canonical_line_is_mathematically_selected":
            True,

        "native_preparation_dynamics_not_yet_constructed":
            True,

        "projector_not_yet_admitted_as_physical_operation":
            True,

        "no_Born_rule":
            True,

        "no_trial_weighting":
            True,

        "no_probability_table":
            True,

        "no_no_signaling_claim":
            True,

        "no_CHSH_claim":
            True,
    },

    "next_gate": (
        "Construct or identify a native setting-independent "
        "preparation operation whose image is this canonical "
        "rank-one line. The projector is now intrinsic; the "
        "remaining problem is whether the finite mechanics "
        "actually prepares it rather than merely allowing it."
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

note = f"""# Canonical joint projector identity 006

## Result

Audit pass:

    {audit_pass}

Two independently derived selectors from Audits 004 and 005 are
exactly the same operator.

The factor-exchange odd projector is

    P_odd = (I-X)/2.

The common quaternionic-frame neutral projector is

    P_Q =
      (I + K1 tensor K1
         + K2 tensor K2
         + K3 tensor K3) / 4.

Exact finite arithmetic gives

    P_odd = P_Q.

Equivalently,

    2(I-X)
      =
    I
      + K1 tensor K1
      + K2 tensor K2
      + K3 tensor K3.

Both have complex rank one.

## Canonical line

Their common image is the projective line

    [|01> - |10>].

It is:

    exchange odd;

    fixed by each diagonal quaternionic action;

    represented inside the Audit-003 Real form by

        i(|01> - |10>).

Thus the line was not selected by a Bell score.

It is the exact common image of two native-interface symmetry
constructions established before any probability law.

## What is now closed

The joint-state selection geometry is no longer ambiguous once both
earned structural conditions are imposed.

There is one canonical rank-one joint line.

## What remains open

A canonical line is not yet a preparation process.

Program 02 must still construct or identify a setting-independent
native preparation operation whose image is this line.

The present audit does not declare the projector to be an admitted
physical operation.

No Born rule, trial weighting, probability table, no-signaling result,
or CHSH value is asserted.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print("ODD_PROJECTOR_COMPLEX_RANK:", odd_rank)
print(
    "QUATERNIONIC_PROJECTOR_COMPLEX_RANK:",
    quaternionic_rank,
)
print(
    "EXACT_PROJECTOR_IDENTITY:",
    checks["exact_projector_identity"],
)
print(
    "CANONICAL_PROJECTIVE_LINE:",
    "[|01>-|10>]",
)
print(
    "REAL_FORM_REPRESENTATIVE:",
    "i(|01>-|10>)",
)
print("FAILED_CHECK_COUNT:", len(failed_checks))
print("FAILED_CHECKS:", failed_checks)
print("JSON_OUT:", JSON_OUT)
print("NOTE_OUT:", NOTE_OUT)
print(
    "JSON_SHA256:",
    sha256(JSON_OUT.read_bytes()).hexdigest(),
)
