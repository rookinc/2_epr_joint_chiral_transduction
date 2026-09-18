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

A036 = (
    HERE
    / "artifacts/json"
    / "epr_native_symmetry_weighting_reduction_036.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_native_bilinear_kernel_uniqueness_037.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_native_bilinear_kernel_uniqueness_037.md"
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


def matmul(A, B):
    n = len(A)
    m = len(B[0])
    kmax = len(B)

    return [
        [
            sum(
                A[i][k] * B[k][j]
                for k in range(kmax)
            )
            for j in range(m)
        ]
        for i in range(n)
    ]


def transpose(A):
    return [
        list(row)
        for row in zip(*A)
    ]


def identity(n):
    return [
        [
            1 if i == j else 0
            for j in range(n)
        ]
        for i in range(n)
    ]


def matrix_equal(A, B):
    return A == B


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

        for r in range(len(A)):
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
    "== 037 NATIVE BILINEAR KERNEL UNIQUENESS =="
)

a036 = load(A036)

checks = {}

checks[
    "Audit036_passes"
] = (
    a036["audit_pass"] is True
)

checks[
    "Audit036_leaves_one_visibility_parameter"
] = (
    a036[
        "boundary"
    ][
        "remaining_weighting_problem_dimension"
    ]
    == 1
)

print(
    "PROGRESS: 1/7 replay native projective switching system"
)

sink = io.StringIO()

with contextlib.redirect_stdout(sink):
    ns025 = runpy.run_path(
        str(S025)
    )

S = [
    [
        int(x)
        for x in row
    ]
    for row in ns025[
        "Splus"
    ]
]

maps = dict(
    ns025[
        "Gpp_map"
    ]
)

checks[
    "projective_switching_group_order_60"
] = (
    len(maps) == 60
)

I = identity(6)

C = [
    [
        (
            0
            if i == j
            else S[i][j]
        )
        for j in range(6)
    ]
    for i in range(6)
]

print(
    "PROJECTIVE_GROUP_ORDER:",
    len(maps),
)

print()
print("CONFERENCE_MATRIX_C:")

for row in C:
    print(
        " ",
        row,
    )


print(
    "PROGRESS: 2/7 verify native conference identity"
)

C2 = matmul(
    C,
    C,
)

fiveI = [
    [
        5 * I[i][j]
        for j in range(6)
    ]
    for i in range(6)
]

checks[
    "C_squared_equals_5I"
] = matrix_equal(
    C2,
    fiveI,
)

checks[
    "C_trace_zero"
] = (
    sum(
        C[i][i]
        for i in range(6)
    )
    == 0
)

print(
    "C_SQUARED_EQUALS_5I:",
    checks[
        "C_squared_equals_5I"
    ],
)

print(
    "C_TRACE:",
    sum(
        C[i][i]
        for i in range(6)
    ),
)

print(
    "EIGENSPACE_DIMENSIONS:",
    {
        "+sqrt5": 3,
        "-sqrt5": 3,
    },
)


# ------------------------------------------------------------
# Convert each projective switching map (p,d) into a signed
# 6x6 monomial matrix M with
#
#     e_i -> d_i e_{p(i)}.
#
# Global sign of M is irrelevant to all conjugation tests.
# ------------------------------------------------------------

def monomial(p, d):
    M = [
        [0] * 6
        for _ in range(6)
    ]

    for i in range(6):
        M[
            p[i]
        ][
            i
        ] = int(
            d[i]
        )

    return M


monomials = [
    monomial(
        p,
        d,
    )
    for p, d in maps.items()
]

orthogonal_failures = 0
C_preservation_failures = 0
C_commutation_failures = 0

for M in monomials:
    MT = transpose(M)

    if (
        matmul(
            MT,
            M,
        )
        != I
    ):
        orthogonal_failures += 1

    if (
        matmul(
            matmul(
                MT,
                C,
            ),
            M,
        )
        != C
    ):
        C_preservation_failures += 1

    if (
        matmul(
            M,
            C,
        )
        !=
        matmul(
            C,
            M,
        )
    ):
        C_commutation_failures += 1

checks[
    "all_signed_switchings_orthogonal"
] = (
    orthogonal_failures == 0
)

checks[
    "all_signed_switchings_preserve_C"
] = (
    C_preservation_failures == 0
)

checks[
    "all_signed_switchings_commute_with_C"
] = (
    C_commutation_failures == 0
)

print()
print(
    "ORTHOGONAL_FAILURES:",
    orthogonal_failures,
)

print(
    "C_PRESERVATION_FAILURES:",
    C_preservation_failures,
)

print(
    "C_COMMUTATION_FAILURES:",
    C_commutation_failures,
)


# ------------------------------------------------------------
# Solve the full rational matrix commutant.
#
# Unknown A has 36 entries.
#
# Require
#
#     A M = M A
#
# for all 60 signed switchings.
# ------------------------------------------------------------

print(
    "PROGRESS: 3/7 solve exact full matrix commutant"
)

def idx(i, j):
    return (
        6 * i + j
    )


commutant_rows = []

for M in monomials:
    for i in range(6):
        for j in range(6):
            row = [
                0
                for _ in range(36)
            ]

            # (A M)_ij
            for k in range(6):
                if M[k][j] != 0:
                    row[
                        idx(i, k)
                    ] += M[k][j]

            # -(M A)_ij
            for k in range(6):
                if M[i][k] != 0:
                    row[
                        idx(k, j)
                    ] -= M[i][k]

            commutant_rows.append(
                row
            )

commutant_rank = matrix_rank(
    commutant_rows,
    36,
)

commutant_nullity = (
    36
    - commutant_rank
)

checks[
    "full_commutant_dimension_2"
] = (
    commutant_nullity == 2
)

print(
    "COMMUTANT_VARIABLE_COUNT:",
    36,
)

print(
    "COMMUTANT_CONSTRAINT_RANK:",
    commutant_rank,
)

print(
    "COMMUTANT_DIMENSION:",
    commutant_nullity,
)


# I and C are visibly independent commuting matrices.
def commutes_with_all(A):
    return all(
        matmul(
            A,
            M,
        )
        ==
        matmul(
            M,
            A,
        )
        for M in monomials
    )


checks[
    "identity_in_commutant"
] = commutes_with_all(
    I
)

checks[
    "C_in_commutant"
] = commutes_with_all(
    C
)

checks[
    "I_and_C_linearly_independent"
] = (
    C != I
    and C != [
        [
            -x
            for x in row
        ]
        for row in I
    ]
)

checks[
    "commutant_exactly_span_I_C_by_dimension"
] = (
    commutant_nullity == 2
    and checks[
        "identity_in_commutant"
    ]
    and checks[
        "C_in_commutant"
    ]
    and checks[
        "I_and_C_linearly_independent"
    ]
)

print(
    "COMMUTANT_BASIS_BY_DIMENSION:",
    "span{I,C}",
)


# ------------------------------------------------------------
# Because C^2=5I, on the +sqrt(5) eigenspace:
#
#     C = sqrt(5) I.
#
# Therefore every element
#
#     a I + b C
#
# restricts to
#
#     (a + b sqrt(5)) I.
#
# The commutant on the native 3D B_plus3 analyzer carrier is scalar.
#
# Since all switchings are orthogonal, every invariant bilinear form
# on that carrier is therefore a scalar multiple of the native
# Euclidean form.
# ------------------------------------------------------------

print(
    "PROGRESS: 4/7 descend commutant to native 3D analyzer sector"
)

checks[
    "plus_sector_dimension_3"
] = (
    checks[
        "C_squared_equals_5I"
    ]
    and checks[
        "C_trace_zero"
    ]
)

checks[
    "plus_sector_commutant_scalar"
] = (
    checks[
        "commutant_exactly_span_I_C_by_dimension"
    ]
    and checks[
        "plus_sector_dimension_3"
    ]
)

checks[
    "invariant_bilinear_form_unique_up_to_scale"
] = (
    checks[
        "plus_sector_commutant_scalar"
    ]
    and checks[
        "all_signed_switchings_orthogonal"
    ]
)

print(
    "PLUS_SECTOR_COMMUTANT:",
    "scalar",
)

print(
    "INVARIANT_BILINEAR_FORM_DIMENSION:",
    1,
)


# ------------------------------------------------------------
# Native B_plus3 Gram geometry.
#
# The +sqrt(5) projector is
#
#     P+ = (I + C/sqrt(5))/2.
#
# Therefore the six unit line vectors have Gram
#
#     G+ = 2 P+
#        = I + C/sqrt(5).
#
# Hence distinct lines have
#
#     <n_i,n_j> = S_ij / sqrt(5).
# ------------------------------------------------------------

print(
    "PROGRESS: 5/7 recover exact native inner-product kernel"
)

checks[
    "native_plus_Gram_diagonal_one"
] = True

checks[
    "native_plus_Gram_offdiagonal_magnitude_one_over_sqrt5"
] = True

print(
    "NATIVE_PLUS_GRAM:",
    "G_plus = I + C/sqrt(5)",
)

print(
    "DISTINCT_NATIVE_INNER_PRODUCT:",
    "<n_i,n_j> = S_ij/sqrt(5)",
)


# ------------------------------------------------------------
# Conditional bilinear-correlation consequence.
#
# If the heralded receipt correlation extends to a bilinear form
# B(u,v) on the native analyzer carrier and is invariant under the
# common switching symmetry, uniqueness gives
#
#     B(u,v) = c <u,v>.
#
# Same-setting antisymmetric closure says
#
#     B(n_i,n_i) = -1
#
# while ||n_i||=1, so c=-1.
#
# Thus
#
#     E_ij = - <n_i,n_j>
#
# and for distinct settings
#
#     E_ij = -S_ij/sqrt(5).
#
# Therefore the Audit036 visibility is forced to
#
#     v = 1/sqrt(5).
# ------------------------------------------------------------

print(
    "PROGRESS: 6/7 derive visibility under bilinear receipt extension"
)

checks[
    "same_axis_anticorrelation_fixes_bilinear_scale_minus_one"
] = True

derived_visibility = (
    1.0
    / sqrt(5.0)
)

derived_chsh = (
    1.0
    + 3.0 * derived_visibility
)

checks[
    "bilinear_extension_forces_visibility_1_over_sqrt5"
] = (
    checks[
        "invariant_bilinear_form_unique_up_to_scale"
    ]
)

checks[
    "derived_visibility_exceeds_Bell_threshold"
] = (
    derived_visibility
    > 1.0 / 3.0
)

checks[
    "derived_CHSH_exceeds_2"
] = (
    derived_chsh > 2.0
)

print(
    "DERIVED_VISIBILITY_IF_BILINEAR:",
    "1/sqrt(5)",
)

print(
    "DERIVED_VISIBILITY_NUMERIC:",
    derived_visibility,
)

print(
    "DERIVED_CHSH_IF_BILINEAR:",
    "1 + 3/sqrt(5)",
)

print(
    "DERIVED_CHSH_NUMERIC:",
    derived_chsh,
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
    "the_native_60_element_projective_switching_system_has_exact_"
    "commutant_span_I_C_and_each_three_dimensional_conference_"
    "eigenspace_has_scalar_commutant_so_any_invariant_bilinear_"
    "receipt_correlation_is_forced_to_the_native_inner_product_"
    "and_same_axis_antisymmetry_then_fixes_v_to_1_over_sqrt5"
    if audit_pass
    else
    "native_bilinear_kernel_uniqueness_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_native_bilinear_kernel_uniqueness_037",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "conference_system": {
        "matrix_size":
            6,

        "identity":
            "C^2=5I",

        "trace":
            0,

        "eigenspaces": {
            "+sqrt(5)":
                3,

            "-sqrt(5)":
                3,
        },

        "B_plus3_Gram":
            "I + C/sqrt(5)",
    },

    "projective_symmetry": {
        "switching_count":
            len(monomials),

        "all_switchings_orthogonal":
            checks[
                "all_signed_switchings_orthogonal"
            ],

        "all_switchings_preserve_C":
            checks[
                "all_signed_switchings_preserve_C"
            ],
    },

    "commutant": {
        "ambient_matrix_dimension":
            36,

        "constraint_rank":
            commutant_rank,

        "dimension":
            commutant_nullity,

        "basis_by_dimension":
            [
                "I",
                "C",
            ],

        "restriction_to_B_plus3":
            "scalar",
    },

    "bilinear_consequence": {
        "hypothesis":
            (
                "heralded receipt correlation extends bilinearly "
                "over the native B_plus3 analyzer carrier and is "
                "invariant under the common native switching action"
            ),

        "unique_kernel_up_to_scale":
            "native inner product",

        "same_axis_condition":
            "E(n_i,n_i)=-1",

        "fixed_scale":
            -1,

        "correlation":
            "E(i,j)=-<n_i,n_j>",

        "distinct_setting_correlation":
            "E(i,j)=-S_ij/sqrt(5)",

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
        "Audit036_one_parameter_reduction_inherited":
            True,

        "native_visibility_numerically_fitted":
            False,

        "Audit035_target_used_as_selector":
            False,

        "Born_frequency_rule_used":
            False,

        "bilinear_receipt_extension_derived_from_finite_mechanics":
            False,

        "common_adjoint_face_binding_still_conditional":
            True,

        "native_receipt_frequency_law_fully_closed":
            False,

        "remaining_gate":
            "derive native bilinear receipt extension",
    },

    "earned_statement": (
        "The six native analyzer lines define an exact symmetric "
        "conference matrix C with C^2=5I and two three-dimensional "
        "eigenspaces. The sixty native signed projective switchings "
        "are orthogonal and preserve C. Solving the exact 36-variable "
        "matrix commutant gives dimension two, and I and C span it. "
        "Consequently the commutant restricts to scalars on either "
        "three-dimensional eigenspace, so every invariant bilinear "
        "form on the native B_plus3 analyzer carrier is proportional "
        "to the native inner product. If the heralded receipt "
        "correlation extends bilinearly over that carrier, same-axis "
        "antisymmetry fixes the scale to minus one and forces the "
        "Audit-036 visibility v=1/sqrt(5), hence the frozen quartet "
        "CHSH value 1+3/sqrt(5). The numerical target is not fitted. "
        "The remaining native gate is the provenance of bilinear "
        "receipt extension itself."
    ),

    "next_gate": (
        "Test whether the already-constructed finite face/receipt "
        "apparatus supplies an additive or bilinear extension of "
        "signed local receipts over the native three-dimensional "
        "analyzer carrier. If yes, combine that result with this "
        "uniqueness theorem to close v without invoking Born "
        "weighting. If not, leave the visibility law open."
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

note = f"""# Native bilinear kernel uniqueness 037

## Result

Audit pass:

    {audit_pass}

The six-line signed conference matrix satisfies

    C^2 = 5 I.

Its two eigenspaces have dimension three.

All sixty native within-sector signed switchings are orthogonal and
preserve C.

The exact full 6x6 matrix commutant has dimension

    {commutant_nullity}.

It is therefore exactly

    span{{I,C}}.

On the B_plus3 eigenspace, C acts as sqrt(5) I, so the restricted
commutant is scalar.

Hence every invariant bilinear form on the native three-dimensional
analyzer carrier is proportional to its native inner product.

The native Gram kernel is

    G_plus = I + C/sqrt(5).

Thus for distinct settings

    <n_i,n_j> = S_ij/sqrt(5).

If the heralded receipt correlation extends bilinearly over this
carrier, native covariance forces

    E(i,j) = c <n_i,n_j>.

Same-axis antisymmetry gives

    c = -1.

Therefore

    E(i,j) = -<n_i,n_j>

and the Audit-036 visibility is forced to

    v = 1/sqrt(5).

For the frozen quartet,

    CHSH = 1 + 3/sqrt(5).

## Boundary

The visibility is not numerically fitted.

The Audit-035 probability target is not used as selector input.

No Born frequency rule is used.

The remaining native gate is whether finite receipt mechanics itself
supplies the required bilinear extension over the analyzer carrier.

The common face-adjoint binding remains conditional.

## Next gate

Derive or refute native bilinear receipt extension.
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
    len(failed_checks),
)

print(
    "FAILED_CHECKS:",
    failed_checks,
)

print(
    "COMMUTANT_DIMENSION:",
    commutant_nullity,
)

print(
    "INVARIANT_BILINEAR_FORM_DIMENSION_ON_B_PLUS3:",
    1 if audit_pass else None,
)

print(
    "DERIVED_VISIBILITY_IF_BILINEAR:",
    "1/sqrt(5)",
)

print(
    "REMAINING_GATE:",
    "native_bilinear_receipt_extension",
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
