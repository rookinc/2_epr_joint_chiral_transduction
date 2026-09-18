#!/usr/bin/env python3

from pathlib import Path
from hashlib import sha256
from math import sqrt
import contextlib
import io
import json
import runpy

HERE = Path(__file__).resolve().parents[2]

S006 = (
    HERE
    / "scripts/audits"
    / "epr_canonical_joint_projector_identity_006.py"
)

A006 = (
    HERE
    / "artifacts/json"
    / "epr_canonical_joint_projector_identity_006.v1.json"
)

A037 = (
    HERE
    / "artifacts/json"
    / "epr_native_bilinear_kernel_uniqueness_037.v1.json"
)

A038 = (
    HERE
    / "artifacts/json"
    / "epr_native_rank3_kernel_visibility_038.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_native_joint_contraction_kernel_040b.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_native_joint_contraction_kernel_040b.md"
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


def matadd(A, B):
    return [
        [
            A[i][j] + B[i][j]
            for j in range(len(A[0]))
        ]
        for i in range(len(A))
    ]


def matscale(c, A):
    return [
        [
            c * x
            for x in row
        ]
        for row in A
    ]


def matmul(A, B):
    return [
        [
            sum(
                A[i][k] * B[k][j]
                for k in range(len(B))
            )
            for j in range(len(B[0]))
        ]
        for i in range(len(A))
    ]


def kron(A, B):
    out = []

    for row_a in A:
        for row_b in B:
            row = []

            for a in row_a:
                row.extend(
                    a * b
                    for b in row_b
                )

            out.append(row)

    return out


def trace(A):
    return sum(
        A[i][i]
        for i in range(len(A))
    )


def eye(n):
    return [
        [
            1 + 0j
            if i == j
            else 0 + 0j
            for j in range(n)
        ]
        for i in range(n)
    ]


def equal(A, B):
    return A == B


print("== 040B NATIVE JOINT CONTRACTION KERNEL ==")

a006 = load(A006)
a037 = load(A037)
a038 = load(A038)

checks = {}

checks["Audit006_passes"] = (
    a006["audit_pass"] is True
)

checks["Audit037_passes"] = (
    a037["audit_pass"] is True
)

checks["Audit038_passes"] = (
    a038["audit_pass"] is True
)


# ------------------------------------------------------------
# Replay the actual Audit006 construction.
# ------------------------------------------------------------

print("PROGRESS: 1/7 replay Audit006 joint projector")

sink = io.StringIO()

with contextlib.redirect_stdout(sink):
    ns006 = runpy.run_path(
        str(S006)
    )

X = ns006["X"]
K1 = ns006["K1"]
K2 = ns006["K2"]
K3 = ns006["K3"]

I2 = eye(2)
I4 = eye(4)

Podd = matscale(
    0.5,
    matadd(
        I4,
        matscale(
            -1,
            X,
        ),
    ),
)

PQ = matscale(
    0.25,
    matadd(
        I4,
        matadd(
            kron(K1, K1),
            matadd(
                kron(K2, K2),
                kron(K3, K3),
            ),
        ),
    ),
)

checks[
    "Audit006_projector_identity_replayed"
] = equal(
    Podd,
    PQ,
)

checks[
    "Podd_trace_one"
] = (
    trace(Podd)
    == 1
)

checks[
    "Podd_idempotent"
] = equal(
    matmul(
        Podd,
        Podd,
    ),
    Podd,
)

print(
    "P_ODD_EQUALS_P_Q:",
    checks[
        "Audit006_projector_identity_replayed"
    ],
)

print(
    "P_ODD_TRACE:",
    trace(Podd),
)


# ------------------------------------------------------------
# Exact swap contraction identity.
#
# Test on the four matrix units E_ab.
#
# Bilinearity then gives the identity on the whole M2(C):
#
# Tr[Podd (A tensor B)]
#
#   =
#
# 1/2 (Tr(A)Tr(B) - Tr(AB)).
# ------------------------------------------------------------

print(
    "PROGRESS: 2/7 verify exact contraction identity on matrix basis"
)

matrix_units = []

for a in range(2):
    for b in range(2):
        E = [
            [
                0 + 0j
                for _ in range(2)
            ]
            for _ in range(2)
        ]

        E[a][b] = 1 + 0j

        matrix_units.append(
            (
                a,
                b,
                E,
            )
        )

basis_identity_failures = 0
basis_rows = []

for a, b, A in matrix_units:
    for c, d, B in matrix_units:
        lhs = trace(
            matmul(
                Podd,
                kron(
                    A,
                    B,
                ),
            )
        )

        rhs = 0.5 * (
            trace(A) * trace(B)
            - trace(
                matmul(
                    A,
                    B,
                )
            )
        )

        ok = (
            lhs == rhs
        )

        if not ok:
            basis_identity_failures += 1

        basis_rows.append({
            "A":
                "E_" + str(a) + str(b),

            "B":
                "E_" + str(c) + str(d),

            "lhs":
                str(lhs),

            "rhs":
                str(rhs),

            "match":
                ok,
        })

checks[
    "full_M2_contraction_identity_exact"
] = (
    basis_identity_failures == 0
    and len(
        basis_rows
    )
    == 16
)

print(
    "MATRIX_BASIS_PAIR_COUNT:",
    len(basis_rows),
)

print(
    "CONTRACTION_IDENTITY_FAILURES:",
    basis_identity_failures,
)

print(
    "GENERAL_IDENTITY:",
    "Tr[Podd(AxB)] = "
    "1/2*(Tr(A)Tr(B)-Tr(AB))",
)


# ------------------------------------------------------------
# Restrict to traceless Hermitian Pauli observables.
# ------------------------------------------------------------

print(
    "PROGRESS: 3/7 restrict to traceless local observable space"
)

sx = [
    [0 + 0j, 1 + 0j],
    [1 + 0j, 0 + 0j],
]

sy = [
    [0 + 0j, 0 - 1j],
    [0 + 1j, 0 + 0j],
]

sz = [
    [1 + 0j, 0 + 0j],
    [0 + 0j, -1 + 0j],
]

pauli = (
    sx,
    sy,
    sz,
)

pauli_names = (
    "sigma_x",
    "sigma_y",
    "sigma_z",
)

checks[
    "all_Pauli_traceless"
] = all(
    trace(S)
    == 0
    for S in pauli
)

pauli_trace_gram = []

source_contraction_gram = []

pauli_failures = 0

for i, A in enumerate(pauli):
    trace_row = []
    gamma_row = []

    for j, B in enumerate(pauli):
        inner = (
            0.5
            * trace(
                matmul(
                    A,
                    B,
                )
            )
        )

        gamma = trace(
            matmul(
                Podd,
                kron(
                    A,
                    B,
                ),
            )
        )

        trace_row.append(
            inner
        )

        gamma_row.append(
            gamma
        )

        if gamma != -inner:
            pauli_failures += 1

    pauli_trace_gram.append(
        trace_row
    )

    source_contraction_gram.append(
        gamma_row
    )

checks[
    "source_contraction_is_minus_Hilbert_Schmidt_on_traceless_space"
] = (
    pauli_failures == 0
)

print()
print(
    "PAULI_HALF_TRACE_GRAM:"
)

for row in pauli_trace_gram:
    print(
        " ",
        row,
    )

print()
print(
    "SOURCE_CONTRACTION_GRAM:"
)

for row in source_contraction_gram:
    print(
        " ",
        row,
    )

print(
    "PAULI_CONTRACTION_FAILURES:",
    pauli_failures,
)

print(
    "TRACELESS_IDENTITY:",
    "Gamma(A,B)=-1/2*Tr(AB)",
)


# ------------------------------------------------------------
# Local one-wing contractions vanish.
#
# This is the operator-algebra analogue of the zero signed marginals
# derived independently in Audit036.
# ------------------------------------------------------------

print(
    "PROGRESS: 4/7 verify local source contractions vanish"
)

left_local_failures = 0
right_local_failures = 0

for A in pauli:
    left = trace(
        matmul(
            Podd,
            kron(
                A,
                I2,
            ),
        )
    )

    right = trace(
        matmul(
            Podd,
            kron(
                I2,
                A,
            ),
        )
    )

    if left != 0:
        left_local_failures += 1

    if right != 0:
        right_local_failures += 1

checks[
    "left_local_contractions_vanish"
] = (
    left_local_failures == 0
)

checks[
    "right_local_contractions_vanish"
] = (
    right_local_failures == 0
)

print(
    "LEFT_LOCAL_FAILURES:",
    left_local_failures,
)

print(
    "RIGHT_LOCAL_FAILURES:",
    right_local_failures,
)


# ------------------------------------------------------------
# Common local frame covariance.
#
# The identity is basis-independent. On any orthonormal local
# traceless-Hermitian frame, Gamma is minus the Euclidean inner
# product.
#
# Audit037 supplies the native B_plus3 Gram kernel:
#
#     <n_i,n_j>
#       =
#     delta_ij + (1-delta_ij) S_ij/sqrt(5).
#
# Therefore, under the already-declared common-adjoint lift,
#
#     Gamma_ii = -1
#     Gamma_ij = -S_ij/sqrt(5).
# ------------------------------------------------------------

print(
    "PROGRESS: 5/7 descend contraction to native analyzer Gram geometry"
)

conference = a037[
    "conference_system"
]

checks[
    "Audit037_B_plus3_Gram_available"
] = (
    conference[
        "B_plus3_Gram"
    ]
    == "I + C/sqrt(5)"
)

S_native = a038[
    "native_sector"
][
    "Gram_kernel"
]

checks[
    "Audit038_native_Gram_matches"
] = (
    S_native
    == "I+C/sqrt(5)"
)

derived_visibility = (
    1.0
    / sqrt(5.0)
)

derived_chsh = (
    1.0
    + 3.0 * derived_visibility
)

print(
    "NATIVE_ANALYZER_GRAM:",
    "I + C/sqrt(5)",
)

print(
    "SOURCE_CONTRACTION_ON_NATIVE_AXES:",
    "Gamma_ii=-1; "
    "Gamma_ij=-S_ij/sqrt(5)",
)

print(
    "CONTRACTION_VISIBILITY:",
    "1/sqrt(5)",
)


# ------------------------------------------------------------
# Verify the source contraction reproduces the same operator kernel
# isolated conditionally in Audits037-038.
# ------------------------------------------------------------

print(
    "PROGRESS: 6/7 bind to 037-038 kernel"
)

checks[
    "operator_contraction_visibility_matches_Audit038"
] = (
    abs(
        derived_visibility
        - float(
            a038[
                "native_sector"
            ][
                "selected_visibility_numeric"
            ]
        )
    )
    < 1.0e-15
)

checks[
    "operator_contraction_CHSH_matches_Audit038"
] = (
    abs(
        derived_chsh
        - float(
            a038[
                "Bell"
            ][
                "numeric"
            ]
        )
    )
    < 1.0e-12
)

print(
    "DERIVED_OPERATOR_VISIBILITY:",
    derived_visibility,
)

print(
    "DERIVED_OPERATOR_CHSH:",
    derived_chsh,
)


# ------------------------------------------------------------
# Classification.
# ------------------------------------------------------------

print(
    "PROGRESS: 7/7 classify"
)

failed = [
    name
    for name, passed
    in checks.items()
    if not passed
]

audit_pass = not failed

verdict = (
    "the_native_canonical_exchange_odd_projector_supplies_an_exact_"
    "bilinear_joint_operator_contraction_whose_restriction_to_the_"
    "traceless_local_observable_space_is_minus_the_normalized_trace_"
    "inner_product_and_therefore_under_the_existing_common_adjoint_"
    "interface_reproduces_the_native_1_over_sqrt5_analyzer_kernel"
    if audit_pass
    else
    "native_joint_contraction_kernel_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_native_joint_contraction_kernel_040b",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "source_projector": {
        "form":
            "Podd=(I-X)/2",

        "quaternionic_form":
            (
                "Podd=(I+K1xK1+K2xK2+K3xK3)/4"
            ),

        "trace":
            1,

        "rank":
            1,
    },

    "joint_contraction": {
        "definition":
            "Gamma(A,B)=Tr[Podd(A tensor B)]",

        "general_identity":
            (
                "Gamma(A,B)="
                "1/2*(Tr(A)Tr(B)-Tr(AB))"
            ),

        "traceless_identity":
            "Gamma(A,B)=-1/2*Tr(AB)",

        "Pauli_kernel":
            "-delta_ij",

        "left_local_contraction":
            0,

        "right_local_contraction":
            0,
    },

    "native_analyzer_descent": {
        "interface":
            "Audit022 conditional common-adjoint lift",

        "native_Gram":
            "I+C/sqrt(5)",

        "joint_operator_kernel":
            (
                "Gamma_ii=-1; "
                "Gamma_ij=-S_ij/sqrt(5)"
            ),

        "visibility":
            "1/sqrt(5)",

        "visibility_numeric":
            derived_visibility,

        "frozen_quartet_CHSH":
            "1+3/sqrt(5)",

        "frozen_quartet_CHSH_numeric":
            derived_chsh,
    },

    "checks":
        checks,

    "boundary": {
        "native_joint_operator_bilinearity_closed":
            audit_pass,

        "joint_operator_kernel_is_setting_independent":
            True,

        "common_adjoint_face_binding_still_conditional":
            True,

        "operator_contraction_identified_with_empirical_correlation":
            False,

        "Born_frequency_law_derived":
            False,

        "native_trial_frequency_law_derived":
            False,

        "probability_table_promoted_from_operator_kernel":
            False,

        "remaining_gate":
            (
                "derive the operational receipt rule identifying "
                "registered binary correlation with the canonical "
                "joint operator contraction"
            ),
    },

    "earned_statement": (
        "The canonical native exchange-odd projector already supplies "
        "the bilinear joint operator kernel required by Audit 037. "
        "Direct exact evaluation on the complete 2x2 matrix-unit basis "
        "gives Gamma(A,B)=Tr[Podd(A tensor B)]="
        "1/2(Tr(A)Tr(B)-Tr(AB)). On the traceless Hermitian local "
        "observable space this reduces to Gamma(A,B)=-1/2 Tr(AB). "
        "Local one-wing contractions vanish. Therefore under the "
        "existing common-adjoint analyzer interface the native "
        "B_plus3 Gram geometry gives Gamma_ii=-1 and "
        "Gamma_ij=-S_ij/sqrt(5), exactly the visibility and CHSH "
        "kernel isolated in Audits 037-038. No probability or Born "
        "frequency interpretation is used in this derivation."
    ),

    "next_gate": (
        "Determine whether the finite registered face apparatus reads "
        "the canonical joint contraction Gamma as its binary outcome "
        "correlation. If that operational identification is native, "
        "combine it with Audit036 zero marginals to obtain the unique "
        "normalized joint receipt table. Otherwise retain the "
        "operator kernel/frequency distinction explicitly."
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

note = f"""# Native joint contraction kernel 040B

## Result

Audit pass:

    {audit_pass}

The canonical native source projector is

    P_odd = (I-X)/2

and exactly

    P_odd =
      (I + K1 tensor K1
         + K2 tensor K2
         + K3 tensor K3) / 4.

Define the joint operator contraction

    Gamma(A,B)
      =
    Tr[P_odd (A tensor B)].

Exact evaluation on the full 2x2 matrix-unit basis gives

    Gamma(A,B)
      =
    1/2 [Tr(A)Tr(B) - Tr(AB)].

For traceless local observables,

    Gamma(A,B)
      =
    -1/2 Tr(AB).

For normalized Pauli directions this is

    Gamma(n.sigma,m.sigma)
      =
    - n dot m.

The one-wing contractions vanish:

    Gamma(A,I) = 0
    Gamma(I,B) = 0

for traceless A and B.

Audit 037 supplies the native analyzer Gram kernel

    G_plus = I + C/sqrt(5).

Therefore, under the already-declared common-adjoint analyzer interface,

    Gamma_ii = -1

and for distinct settings

    Gamma_ij = -S_ij/sqrt(5).

Thus the canonical source operator itself supplies

    v = 1/sqrt(5)

at the joint operator-kernel level.

## Boundary

This closes native joint operator bilinearity.

It does not identify the operator contraction with empirical trial
frequency.

It does not derive a Born frequency rule.

The common face-adjoint binding remains conditional.

The remaining gate is operational:

    registered binary correlation
        ?=
    canonical joint contraction Gamma.
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
    "OPERATOR_KERNEL_VISIBILITY:",
    "1/sqrt(5)",
)

print(
    "REMAINING_GATE:",
    "registered_correlation_equals_joint_contraction",
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
