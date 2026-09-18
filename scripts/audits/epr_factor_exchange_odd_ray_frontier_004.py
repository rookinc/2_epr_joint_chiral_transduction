#!/usr/bin/env python3

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json

HERE = Path(__file__).resolve().parents[2]

P41 = (
    Path.home()
    / "dev/cori/research/mathematics"
    / "41-order-4-dodecahedral-residue"
)

A002 = (
    HERE
    / "artifacts/json"
    / "epr_native_deck_face_chiral_transduction_002.v1.json"
)

A003 = (
    HERE
    / "artifacts/json"
    / "epr_joint_chiral_tensor_consistency_003.v1.json"
)

AUT = (
    P41
    / "artifacts/json"
    / "native_g60_automorphism_deck_cache.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_factor_exchange_odd_ray_frontier_004.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_factor_exchange_odd_ray_frontier_004.md"
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
    if not matrix:
        return 0

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


def complex_matrix_from_real_vector(v):
    return [
        [
            complex(v[0], v[1]),
            complex(v[2], v[3]),
        ],
        [
            complex(v[4], v[5]),
            complex(v[6], v[7]),
        ],
    ]


def det2(m):
    return (
        m[0][0] * m[1][1]
        - m[0][1] * m[1][0]
    )


print("PROGRESS: 1/7 load sealed interfaces")

a002 = load(A002)
a003 = load(A003)
aut = load(AUT)

checks = {}

checks["audit_002_passed"] = (
    a002["audit_pass"] is True
)

checks["audit_003_passed"] = (
    a003["audit_pass"] is True
)

checks["audit_002_selected_native_a"] = (
    a002[
        "selected_native_S1_kernel"
    ]["native_name"]
    == "a"
)

checks["audit_003_compatibility_is_real_fixed_locus"] = (
    a003[
        "one_sided_transduction"
    ]["compatibility_space_type"]
    == "real fixed locus of C_AB"
)

print("PROGRESS: 2/7 reconstruct native G1800 factor exchange")

a_index = int(
    a002[
        "selected_native_S1_kernel"
    ]["automorphism_index"]
)

aut_rows = {
    int(row["automorphism_index"]): row
    for row in aut["measurements"]["automorphism_rows"]
}

a_perm = tuple(
    int(v)
    for v in aut_rows[a_index]["permutation"]
)

checks["native_a_has_60_points"] = (
    len(a_perm) == 60
)

checks["native_a_is_fixed_point_free"] = all(
    a_perm[i] != i
    for i in range(60)
)

checks["native_a_is_involution"] = all(
    a_perm[a_perm[i]] == i
    for i in range(60)
)


def diag_rep(u, v):
    return min(
        (u, v),
        (a_perm[u], a_perm[v]),
    )


joint_reps = sorted({
    diag_rep(u, v)
    for u in range(60)
    for v in range(60)
})

joint_index = {
    rep: i
    for i, rep in enumerate(joint_reps)
}


def joint_class(u, v):
    return joint_index[diag_rep(u, v)]


swap = {}
tau = {}

swap_well_defined = True
tau_well_defined = True

for cid, (u, v) in enumerate(joint_reps):
    s0 = joint_class(v, u)
    s1 = joint_class(
        a_perm[v],
        a_perm[u],
    )

    if s0 != s1:
        swap_well_defined = False

    t0 = joint_class(
        a_perm[u],
        v,
    )

    t1 = joint_class(
        u,
        a_perm[v],
    )

    if t0 != t1:
        tau_well_defined = False

    swap[cid] = s0
    tau[cid] = t0

checks["g1800_class_count_1800"] = (
    len(joint_reps) == 1800
)

checks["factor_swap_well_defined"] = (
    swap_well_defined
)

checks["tau_well_defined"] = (
    tau_well_defined
)

checks["factor_swap_is_involution"] = all(
    swap[swap[i]] == i
    for i in swap
)

checks["factor_swap_commutes_with_tau"] = all(
    swap[tau[i]]
    == tau[swap[i]]
    for i in swap
)

swap_fixed_count = sum(
    swap[i] == i
    for i in swap
)

swap_orbit_sizes = []
seen = set()

for i in range(len(joint_reps)):
    if i in seen:
        continue

    orbit = {
        i,
        swap[i],
    }

    seen.update(orbit)
    swap_orbit_sizes.append(len(orbit))

swap_orbit_profile = Counter(
    swap_orbit_sizes
)

checks["factor_swap_fixed_count_60"] = (
    swap_fixed_count == 60
)

checks["factor_swap_orbit_profile_exact"] = (
    swap_orbit_profile
    == Counter({
        1: 60,
        2: 870,
    })
)

print("PROGRESS: 3/7 construct joint CP3 presentation")

# Real coordinates for the complex coefficient matrix
#
#     [[c00, c01],
#      [c10, c11]]
#
# ordered as
#
# Re c00, Im c00,
# Re c01, Im c01,
# Re c10, Im c10,
# Re c11, Im c11.

I8 = eye(8)

J = [
    [0, -1, 0,  0, 0,  0, 0,  0],
    [1,  0, 0,  0, 0,  0, 0,  0],

    [0,  0, 0, -1, 0,  0, 0,  0],
    [0,  0, 1,  0, 0,  0, 0,  0],

    [0,  0, 0,  0, 0, -1, 0,  0],
    [0,  0, 0,  0, 1,  0, 0,  0],

    [0,  0, 0,  0, 0,  0, 0, -1],
    [0,  0, 0,  0, 0,  0, 1,  0],
]

# Canonical factor exchange:
#
#     M -> M^T.

P = [
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],

    [0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0],

    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0],

    [0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 1],
]

# Audit 003 Real structure:
#
#     C(M) = sigma_x conjugate(M) sigma_x.

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

checks["J_square_minus_identity"] = equal(
    matmul(J, J),
    matscale(-1, I8),
)

checks["factor_exchange_square_identity"] = equal(
    matmul(P, P),
    I8,
)

checks["real_structure_square_identity"] = equal(
    matmul(C, C),
    I8,
)

checks["factor_exchange_commutes_with_J"] = equal(
    matmul(P, J),
    matmul(J, P),
)

checks["real_structure_anticommutes_with_J"] = equal(
    matmul(C, J),
    matscale(
        -1,
        matmul(J, C),
    ),
)

checks["factor_exchange_commutes_with_real_structure"] = equal(
    matmul(P, C),
    matmul(C, P),
)

print("PROGRESS: 4/7 measure exchange sectors")

P_MINUS_I = matsub(P, I8)
P_PLUS_I = matadd(P, I8)

C_MINUS_I = matsub(C, I8)

exchange_even_real_dim = nullity(
    P_MINUS_I
)

exchange_odd_real_dim = nullity(
    P_PLUS_I
)

compat_real_dim = nullity(
    C_MINUS_I
)

compat_even_real_dim = nullity(
    vcat(
        C_MINUS_I,
        P_MINUS_I,
    )
)

compat_odd_real_dim = nullity(
    vcat(
        C_MINUS_I,
        P_PLUS_I,
    )
)

checks["exchange_even_complex_dimension_3"] = (
    exchange_even_real_dim == 6
)

checks["exchange_odd_complex_dimension_1"] = (
    exchange_odd_real_dim == 2
)

checks["compatibility_real_dimension_4"] = (
    compat_real_dim == 4
)

checks["compatible_even_real_dimension_3"] = (
    compat_even_real_dim == 3
)

checks["compatible_odd_real_dimension_1"] = (
    compat_odd_real_dim == 1
)

print("PROGRESS: 5/7 verify unique odd compatible ray")

# Generator of the C-fixed, P-odd real line:
#
#     M =
#       [[ 0,  i],
#        [-i,  0]]
#
# This is projectively the usual alternating tensor
#
#     |01> - |10>
#
# but no quantum preparation rule is assumed here.

odd_generator = [
    0, 0,
    0, 1,
    0, -1,
    0, 0,
]

checks["odd_generator_is_C_fixed"] = (
    matvec(C, odd_generator)
    == odd_generator
)

checks["odd_generator_is_exchange_odd"] = (
    matvec(P, odd_generator)
    == [
        -x
        for x in odd_generator
    ]
)

odd_matrix = complex_matrix_from_real_vector(
    odd_generator
)

odd_det = det2(
    odd_matrix
)

checks["odd_generator_is_nonproduct"] = (
    odd_det != 0
)

# Exhibit a compatible exchange-even product ray so that
# compatibility alone is visibly insufficient.

even_product = [
    1, 0,
    1, 0,
    1, 0,
    1, 0,
]

checks["even_product_is_C_fixed"] = (
    matvec(C, even_product)
    == even_product
)

checks["even_product_is_exchange_even"] = (
    matvec(P, even_product)
    == even_product
)

even_product_matrix = (
    complex_matrix_from_real_vector(
        even_product
    )
)

checks["even_compatible_example_is_product"] = (
    det2(even_product_matrix) == 0
)

print("PROGRESS: 6/7 classify selector frontier")

checks["odd_parity_would_select_unique_compatible_projective_ray"] = (
    compat_odd_real_dim == 1
)

checks["compatibility_alone_does_not_select_entanglement"] = (
    checks[
        "even_compatible_example_is_product"
    ]
    and checks[
        "odd_generator_is_nonproduct"
    ]
)

print("PROGRESS: 7/7 seal result")

failed_checks = [
    name
    for name, passed in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "canonical_g1800_factor_exchange_reduces_RP3_to_"
    "one_unique_odd_compatible_projective_ray_but_"
    "native_odd_parity_selection_remains_open"
    if audit_pass
    else
    "factor_exchange_odd_ray_frontier_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_factor_exchange_odd_ray_frontier_004",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "native_carrier": {
        "carrier":
            "(G60 x G60)/diag(a)",

        "deck_kernel":
            "a",

        "deck_automorphism_index":
            a_index,

        "factor_exchange":
            "X([u,v]) = [v,u]",

        "factor_exchange_well_defined":
            swap_well_defined,

        "factor_exchange_fixed_joint_state_count":
            swap_fixed_count,

        "factor_exchange_orbit_profile": {
            str(k): v
            for k, v in sorted(
                swap_orbit_profile.items()
            )
        },

        "factor_exchange_commutes_with_tau_a":
            checks[
                "factor_swap_commutes_with_tau"
            ],
    },

    "joint_Hilbert_geometry": {
        "ambient":
            "CP^3",

        "compatibility_locus":
            "RP^3",

        "factor_exchange_even_sector":
            "Sym^2(C^2), complex dimension 3",

        "factor_exchange_odd_sector":
            "Lambda^2(C^2), complex dimension 1",

        "compatible_even_real_dimension":
            compat_even_real_dim,

        "compatible_odd_real_dimension":
            compat_odd_real_dim,

        "compatible_odd_projective_dimension":
            0,

        "compatible_odd_ray":
            (
                "projective line represented by "
                "i(|01>-|10>)"
            ),
    },

    "earned_statement": (
        "The canonical exchange of the two identical G60 "
        "factors descends to the a-kernel G1800 carrier, is an "
        "involution, and commutes with the native relative-lift "
        "exchange tau_a. On the joint Program-01 Hilbert fiber "
        "it commutes with both the joint complex structure and "
        "the Audit-003 Real structure. The Audit-003 compatible "
        "RP3 therefore splits into an exchange-even RP2-like "
        "real locus of dimension three and an exchange-odd real "
        "line of dimension one. The latter contains exactly one "
        "compatible projective ray, represented by the "
        "alternating tensor. Compatibility alone does not "
        "select this ray because the even compatible locus "
        "contains product states."
    ),

    "checks":
        checks,

    "boundary": {
        "factor_exchange_symmetry_is_native_but_exchange_parity_is_not_yet_selected":
            True,

        "odd_ray_is_not_declared_prepared":
            True,

        "alternating_ray_not_inserted_as_target_probability_law":
            True,

        "no_Born_rule_assumed":
            True,

        "no_trial_weighting_rule":
            True,

        "no_no_signaling_claim":
            True,

        "no_CHSH_claim":
            True,

        "no_Bell_nonfactorizability_claim":
            True,
    },

    "next_gate": (
        "Determine whether the native source preparation, "
        "registered orientation, or closure law selects odd "
        "factor-exchange parity before analyzer settings are "
        "chosen. If odd parity is natively selected, the joint "
        "preparation ray is then unique. If parity is not "
        "selected, do not promote the alternating ray by "
        "recognition or Bell score."
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

note = f"""# Factor exchange odd-ray frontier 004

## Result

Audit pass:

    {audit_pass}

The a-kernel common carrier

    G1800 = (G60 x G60)/diag(a)

has a canonical factor exchange

    X([u,v]) = [v,u].

It is well defined, involutive, and commutes with the native
relative-lift exchange tau_a.

At the Program-01 joint Hilbert level, factor exchange commutes with
the joint complex structure and with the Audit-003 Real structure.

## Compatibility locus

Audit 003 established the compatible projective locus

    RP^3 inside CP^3.

Factor exchange splits the ambient joint fiber into:

    exchange-even:
        Sym^2(C^2)
        complex dimension 3

    exchange-odd:
        Lambda^2(C^2)
        complex dimension 1.

Inside the Audit-003 Real form, the exchange-even fixed locus has real
dimension

    {compat_even_real_dim}

while the exchange-odd fixed locus has real dimension

    {compat_odd_real_dim}.

Therefore the exchange-odd compatible sector contains exactly one real
projective ray.

A representative is

    i (|01> - |10>).

The phase i places the alternating projective ray inside the selected
Real form.

## Important negative result

The compatible RP^3 also contains exchange-even product rays.

Therefore joint chiral compatibility by itself does not select
entanglement.

Factor exchange gives a clean one-bit preparation frontier:

    even
    or
    odd.

If a native setting-independent source rule selects odd exchange parity,
the compatible preparation ray becomes unique.

No such parity-selection theorem is claimed here.

## Boundary

The alternating ray has not been inserted as a prepared state.

No Born rule, frequency law, probability table, no-signaling theorem,
Bell factorization result, or CHSH value is claimed.

The next question is purely native:

    Does the source preparation select odd factor-exchange parity
    before Alice and Bob choose analyzer settings?
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print("G1800_CLASS_COUNT:", len(joint_reps))
print("FACTOR_SWAP_FIXED_COUNT:", swap_fixed_count)
print(
    "FACTOR_SWAP_ORBIT_PROFILE:",
    dict(sorted(swap_orbit_profile.items())),
)
print(
    "FACTOR_SWAP_COMMUTES_WITH_TAU:",
    checks[
        "factor_swap_commutes_with_tau"
    ],
)
print(
    "COMPATIBILITY_REAL_DIM:",
    compat_real_dim,
)
print(
    "COMPATIBLE_EVEN_REAL_DIM:",
    compat_even_real_dim,
)
print(
    "COMPATIBLE_ODD_REAL_DIM:",
    compat_odd_real_dim,
)
print(
    "ODD_COMPATIBLE_PROJECTIVE_RAY_COUNT:",
    1 if compat_odd_real_dim == 1 else None,
)
print(
    "ODD_GENERATOR_NONPRODUCT:",
    checks[
        "odd_generator_is_nonproduct"
    ],
)
print(
    "EVEN_COMPATIBLE_PRODUCT_EXISTS:",
    checks[
        "even_compatible_example_is_product"
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
