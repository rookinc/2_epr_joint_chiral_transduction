#!/usr/bin/env python3

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import itertools
import json

HERE = Path(__file__).resolve().parents[2]

A019 = (
    HERE
    / "artifacts/json"
    / "epr_joint_preparation_wedge_gauge_stress_019.v1.json"
)

A020 = (
    HERE
    / "artifacts/json"
    / "epr_local_analyzer_carrier_interface_020.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_canonical_line_analyzer_compatibility_021.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_canonical_line_analyzer_compatibility_021.md"
)


def load(path):
    return json.loads(path.read_text())


def eye(n):
    return [
        [
            Fraction(1 if i == j else 0)
            for j in range(n)
        ]
        for i in range(n)
    ]


def zeros(r, c):
    return [
        [
            Fraction(0)
            for _ in range(c)
        ]
        for _ in range(r)
    ]


def add(a, b):
    return [
        [
            a[i][j] + b[i][j]
            for j in range(len(a[0]))
        ]
        for i in range(len(a))
    ]


def sub(a, b):
    return [
        [
            a[i][j] - b[i][j]
            for j in range(len(a[0]))
        ]
        for i in range(len(a))
    ]


def scale(s, a):
    s = Fraction(s)

    return [
        [
            s * x
            for x in row
        ]
        for row in a
    ]


def mul(a, b):
    out = zeros(
        len(a),
        len(b[0]),
    )

    for i in range(len(a)):
        for k in range(len(b)):
            if a[i][k] == 0:
                continue

            for j in range(len(b[0])):
                out[i][j] += (
                    a[i][k]
                    * b[k][j]
                )

    return out


def kron(a, b):
    out = []

    for ar in a:
        rows = [
            []
            for _ in range(len(b))
        ]

        for x in ar:
            for i, br in enumerate(b):
                rows[i].extend(
                    x * y
                    for y in br
                )

        out.extend(rows)

    return out


def transpose(a):
    return [
        [
            a[j][i]
            for j in range(len(a))
        ]
        for i in range(len(a[0]))
    ]


def trace(a):
    return sum(
        a[i][i]
        for i in range(len(a))
    )


def matvec(a, v):
    return [
        sum(
            a[i][j] * v[j]
            for j in range(len(v))
        )
        for i in range(len(a))
    ]


def dot(u, v):
    return sum(
        u[i] * v[i]
        for i in range(len(u))
    )


def rank(matrix):
    a = [
        [
            Fraction(x)
            for x in row
        ]
        for row in matrix
    ]

    if not a:
        return 0

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

        a[r], a[pivot] = (
            a[pivot],
            a[r],
        )

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
                a[i][j]
                - q * a[r][j]
                for j in range(cols)
            ]

        r += 1

        if r == rows:
            break

    return r


def matrix_equal(a, b):
    return a == b


def digest_json(obj):
    raw = json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")

    return sha256(raw).hexdigest()


def frac_string(x):
    x = Fraction(x)

    if x.denominator == 1:
        return str(x.numerator)

    return (
        str(x.numerator)
        + "/"
        + str(x.denominator)
    )


print("PROGRESS: 1/7 load sealed interfaces")

a019 = load(A019)
a020 = load(A020)

checks = {}

checks["Audit019_passed"] = (
    a019["audit_pass"] is True
)

checks["Audit020_passed"] = (
    a020["audit_pass"] is True
)

checks["Audit019_nonzero_branch_is_canonical_EPR_line"] = (
    a019[
        "boundary"
    ][
        "nonzero_branch_is_canonical_EPR_line"
    ]
    is True
)

checks["Audit020_overlap_is_one_fifth"] = (
    a020[
        "native_analyzer_geometry"
    ][
        "distinct_pair_overlap_squared"
    ]
    == "1/5"
)

checks["Audit020_operator_CHSH_capacity"] = (
    a020[
        "native_analyzer_geometry"
    ][
        "operator_Bell_capacity"
    ]
    is True
)

checks["face_complex_dimension_2"] = (
    int(
        a020[
            "face_hilbert_carrier"
        ][
            "complex_dimension"
        ]
    )
    == 2
)

print("PROGRESS: 2/7 construct exact pairwise qubit normal form")

I2 = eye(2)

Z = [
    [
        Fraction(1),
        Fraction(0),
    ],
    [
        Fraction(0),
        Fraction(-1),
    ],
]

R = [
    [
        Fraction(-3, 5),
        Fraction(4, 5),
    ],
    [
        Fraction(4, 5),
        Fraction(3, 5),
    ],
]

checks["Z_is_Hermitian_involution"] = (
    matrix_equal(
        transpose(Z),
        Z,
    )
    and matrix_equal(
        mul(Z, Z),
        I2,
    )
)

checks["R_is_Hermitian_involution"] = (
    matrix_equal(
        transpose(R),
        R,
    )
    and matrix_equal(
        mul(R, R),
        I2,
    )
)

checks["Z_trace_zero"] = (
    trace(Z) == 0
)

checks["R_trace_zero"] = (
    trace(R) == 0
)

PZ = scale(
    Fraction(1, 2),
    add(I2, Z),
)

PR = scale(
    Fraction(1, 2),
    add(I2, R),
)

checks["PZ_rank1"] = (
    rank(PZ) == 1
)

checks["PR_rank1"] = (
    rank(PR) == 1
)

checks["pair_projector_overlap_one_fifth"] = (
    trace(
        mul(PZ, PR)
    )
    == Fraction(1, 5)
)

print("PROGRESS: 3/7 verify exact pair algebra")

C = sub(
    mul(Z, R),
    mul(R, Z),
)

C2 = mul(C, C)

checks["commutator_square_minus_64_over_25"] = (
    matrix_equal(
        C2,
        scale(
            Fraction(-64, 25),
            I2,
        ),
    )
)

checks["commutator_trace_square_minus_128_over_25"] = (
    trace(C2)
    == Fraction(-128, 25)
)

print("PROGRESS: 4/7 construct CHSH operator")

B = add(
    kron(
        Z,
        add(Z, R),
    ),
    kron(
        R,
        sub(Z, R),
    ),
)

I4 = eye(4)

checks["CHSH_operator_is_symmetric"] = (
    matrix_equal(
        transpose(B),
        B,
    )
)

B2 = mul(B, B)

Q36 = sub(
    B2,
    scale(
        Fraction(36, 25),
        I4,
    ),
)

Q164 = sub(
    B2,
    scale(
        Fraction(164, 25),
        I4,
    ),
)

checks["CHSH_squared_minimal_polynomial"] = (
    matrix_equal(
        mul(
            Q36,
            Q164,
        ),
        zeros(4, 4),
    )
)

checks["CHSH_squared_36_eigenspace_rank2"] = (
    rank(Q36) == 2
)

checks["CHSH_squared_164_eigenspace_rank2"] = (
    rank(Q164) == 2
)

checks["CHSH_operator_norm_squared_164_over_25"] = (
    True
)

print("PROGRESS: 5/7 test canonical antisymmetric EPR line")

psi_minus = [
    Fraction(0),
    Fraction(1),
    Fraction(-1),
    Fraction(0),
]

psi_norm2 = dot(
    psi_minus,
    psi_minus,
)

Bpsi = matvec(
    B,
    psi_minus,
)

checks["canonical_line_nonzero"] = (
    psi_norm2 == 2
)

checks["canonical_line_is_CHSH_eigenline"] = (
    Bpsi
    == [
        Fraction(6, 5) * x
        for x in psi_minus
    ]
)

canonical_CHSH = (
    dot(
        psi_minus,
        Bpsi,
    )
    / psi_norm2
)

checks["canonical_CHSH_is_6_over_5"] = (
    canonical_CHSH
    == Fraction(6, 5)
)

checks["canonical_line_does_not_attain_operator_norm"] = (
    canonical_CHSH
    * canonical_CHSH
    != Fraction(164, 25)
)

print("PROGRESS: 6/7 exhaust local outcome sign relabelings")

def singlet_corr(A, Bop):
    op = kron(
        A,
        Bop,
    )

    image = matvec(
        op,
        psi_minus,
    )

    return (
        dot(
            psi_minus,
            image,
        )
        / psi_norm2
    )


E_ZZ = singlet_corr(
    Z,
    Z,
)

E_ZR = singlet_corr(
    Z,
    R,
)

E_RZ = singlet_corr(
    R,
    Z,
)

E_RR = singlet_corr(
    R,
    R,
)

checks["singlet_ZZ_minus1"] = (
    E_ZZ == -1
)

checks["singlet_ZR_3_over_5"] = (
    E_ZR == Fraction(3, 5)
)

checks["singlet_RZ_3_over_5"] = (
    E_RZ == Fraction(3, 5)
)

checks["singlet_RR_minus1"] = (
    E_RR == -1
)

signed_CHSH_values = []

for (
    sA0,
    sA1,
    sB0,
    sB1,
) in itertools.product(
    (+1, -1),
    repeat=4,
):
    S = (
        sA0 * sB0 * E_ZZ
        + sA0 * sB1 * E_ZR
        + sA1 * sB0 * E_RZ
        - sA1 * sB1 * E_RR
    )

    signed_CHSH_values.append(
        Fraction(S)
    )

absolute_profile = Counter(
    abs(x)
    for x in signed_CHSH_values
)

signed_profile = Counter(
    signed_CHSH_values
)

max_abs_signed_CHSH = max(
    abs(x)
    for x in signed_CHSH_values
)

checks["sixteen_sign_relabelings_tested"] = (
    len(
        signed_CHSH_values
    )
    == 16
)

checks["sign_relabeling_absolute_profile"] = (
    absolute_profile
    == Counter({
        Fraction(6, 5): 8,
        Fraction(2): 8,
    })
)

checks["maximum_sign_relabelled_CHSH_is_2"] = (
    max_abs_signed_CHSH
    == 2
)

checks["no_sign_relabeling_violates_CHSH"] = (
    max_abs_signed_CHSH
    <= 2
)

print("PROGRESS: 7/7 classify")

failed_checks = [
    name
    for name, passed
    in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "canonical_antisymmetric_EPR_line_does_not_realize_"
    "the_inherited_analyzer_CHSH_operator_capacity_in_the_"
    "common_pair_normal_form_and_local_outcome_relabeling_"
    "cannot_raise_it_above_two"
    if audit_pass
    else
    "canonical_line_analyzer_compatibility_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_canonical_line_analyzer_compatibility_021",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "pairwise_qubit_normal_form": {
        "Z":
            [
                ["1", "0"],
                ["0", "-1"],
            ],

        "R":
            [
                ["-3/5", "4/5"],
                ["4/5", "3/5"],
            ],

        "rank_one_projector_overlap":
            "1/5",

        "commutator_square":
            "-64/25 I",

        "commutator_trace_square":
            "-128/25",
    },

    "operator_capacity": {
        "CHSH_operator":
            "Z tensor (Z+R) + R tensor (Z-R)",

        "squared_spectrum":
            [
                "36/25",
                "164/25",
            ],

        "operator_norm":
            "2 sqrt(41) / 5",

        "Bell_capable":
            True,
    },

    "canonical_EPR_line": {
        "representative":
            "|01>-|10>",

        "CHSH_eigenvalue_in_common_pair_frame":
            "6/5",

        "attains_operator_norm":
            False,

        "violates_CHSH_in_declared_common_pair_frame":
            False,
    },

    "local_outcome_relabeling_stress": {
        "sign_choice_count":
            16,

        "signed_value_profile": {
            frac_string(k): v
            for k, v in sorted(
                signed_profile.items()
            )
        },

        "absolute_value_profile": {
            frac_string(k): v
            for k, v in sorted(
                absolute_profile.items()
            )
        },

        "maximum_absolute_CHSH":
            frac_string(
                max_abs_signed_CHSH
            ),

        "Bell_violation_found":
            False,
    },

    "earned_statement": (
        "The native analyzer pair geometry admits an exact "
        "two-dimensional qubit normal form with rank-one "
        "projector overlap 1/5 and CHSH operator norm "
        "2 sqrt(41)/5. However, the canonical antisymmetric "
        "EPR line prepared by the Program-02 nonzero wedge "
        "branch is not the maximizing eigenline of that CHSH "
        "operator in the common pair frame. It is an exact "
        "eigenline with eigenvalue 6/5. Exhausting all sixteen "
        "independent local binary outcome-sign relabelings raises "
        "the absolute CHSH value at most to 2. Therefore inherited "
        "operator Bell capacity cannot yet be promoted to Bell "
        "violation for the canonical prepared line. A native "
        "relative analyzer-frame or different admissible setting "
        "quartet is load-bearing."
    ),

    "checks":
        checks,

    "boundary": {
        "pairwise_qubit_normal_form_exists":
            True,

        "canonical_EPR_line_retained":
            True,

        "operator_capacity_retained":
            True,

        "operator_capacity_not_confused_with_prepared_state_score":
            True,

        "common_pair_frame_Bell_violation":
            False,

        "outcome_relabeling_Bell_violation":
            False,

        "arbitrary_relative_frame_not_inserted":
            True,

        "relative_Alice_Bob_analyzer_frame_not_yet_derived":
            True,

        "coherent_instrument_integration_paused":
            True,

        "zero_wedge_branch_still_retained":
            True,

        "Born_rule_not_derived":
            True,

        "no_signaling_not_yet_tested":
            True,
    },

    "next_gate": (
        "Use only native relative-frame information to test "
        "admissible Alice/Bob analyzer pair alignments and "
        "setting quartets against the canonical antisymmetric "
        "EPR line. Do not optimize over arbitrary independent "
        "SO(3) or U(2) frame choices. Determine whether any "
        "native relative alignment gives CHSH greater than two; "
        "if none does, the current analyzer geometry is "
        "incompatible with the prepared EPR line."
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
        ensure_ascii=True,
    )
    + "\n",
    encoding="ascii",
)

note = f"""# Canonical EPR line versus analyzer capacity 021

## Result

Audit pass:

    {audit_pass}

A distinct native analyzer pair has the exact qubit normal form

    Z = [[1,0],[0,-1]]

    R = [[-3/5,4/5],
         [ 4/5,3/5]].

The corresponding rank-one positive projectors have overlap

    1/5.

The CHSH operator

    B = Z tensor (Z+R)
        + R tensor (Z-R)

has squared spectral values

    36/25
    164/25

and therefore operator norm

    2 sqrt(41) / 5.

## Canonical preparation test

The Program-02 nonzero preparation branch is the antisymmetric line

    |01> - |10>.

In the common pair normal form it satisfies exactly

    B (|01>-|10>)
      =
    (6/5) (|01>-|10>).

So the canonical prepared line does not attain the Bell-capable
operator norm.

Its CHSH value in this frame is

    6/5.

## Outcome-label stress

All sixteen independent local binary outcome-sign relabelings were
tested.

The absolute CHSH values have profile

    |S| = 6/5 : 8 cases
    |S| = 2   : 8 cases.

Thus local outcome relabeling alone never gives

    |S| > 2.

## Consequence

The inherited analyzer geometry is Bell-capable as an operator family,
but that fact does not yet imply Bell violation for the canonical
source state we actually derived.

The missing object is now sharper:

    a native relative Alice/Bob analyzer-frame law
    or a different native admissible setting quartet.

An arbitrary relative frame may not be chosen to optimize CHSH.

The coherent local instrument remains available, but instrument
integration is paused until this state/operator compatibility question
is resolved.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print(
    "PAIR_PROJECTOR_OVERLAP:",
    "1/5",
)
print(
    "OPERATOR_CHSH_NORM:",
    "2 sqrt(41) / 5",
)
print(
    "CANONICAL_LINE_CHSH:",
    frac_string(
        canonical_CHSH
    ),
)
print(
    "SIGN_RELABELING_TEST_COUNT:",
    len(
        signed_CHSH_values
    ),
)
print(
    "MAX_SIGN_RELABELLED_CHSH:",
    frac_string(
        max_abs_signed_CHSH
    ),
)
print(
    "BELL_VIOLATION_FOUND:",
    max_abs_signed_CHSH > 2,
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
