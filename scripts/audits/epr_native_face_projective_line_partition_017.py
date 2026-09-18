#!/usr/bin/env python3

from collections import Counter, defaultdict
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

AUT = (
    P41
    / "artifacts/json"
    / "native_g60_automorphism_deck_cache.v1.json"
)

SIGNED = (
    P41
    / "artifacts/json"
    / "signed_24_face_block_closure_audit_025.json"
)

FV = (
    P41
    / "artifacts/provenance"
    / "project41_positive_pair_hermitian_geometry_201fv.v1.json"
)

FK = (
    P41
    / "artifacts/provenance"
    / "project41_face_field_higgs_bridge_checkpoint_receipt_201fk.v1.json"
)

A016 = (
    HERE
    / "artifacts/json"
    / "epr_literal_phase_decorated_incidence_table_016.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_native_face_projective_line_partition_017.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_native_face_projective_line_partition_017.md"
)


def load(path):
    return json.loads(path.read_text())


def digest_json(obj):
    raw = json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")

    return sha256(raw).hexdigest()


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


def transpose(a):
    return [
        [
            a[j][i]
            for j in range(len(a))
        ]
        for i in range(len(a[0]))
    ]


def matmul(a, b):
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


def outer(u, v):
    return [
        [
            u[i] * v[j]
            for j in range(len(v))
        ]
        for i in range(len(u))
    ]


def matrix_equal(a, b):
    return a == b


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


def perm_order(p):
    q = tuple(range(len(p)))
    current = tuple(range(len(p)))

    for n in range(1, 33):
        current = tuple(
            p[current[i]]
            for i in range(len(p))
        )

        if current == q:
            return n

    raise RuntimeError(
        "permutation order too large"
    )


def permutation_matrix(p):
    n = len(p)
    out = zeros(n, n)

    for i, j in enumerate(p):
        out[j][i] = Fraction(1)

    return out


def frac_string(x):
    x = Fraction(x)

    if x.denominator == 1:
        return str(x.numerator)

    return (
        str(x.numerator)
        + "/"
        + str(x.denominator)
    )


def matrix_strings(a):
    return [
        [
            frac_string(x)
            for x in row
        ]
        for row in a
    ]


def basis_vector(n, i):
    return [
        Fraction(1 if j == i else 0)
        for j in range(n)
    ]


def matvec(a, v):
    return [
        sum(
            a[i][j] * v[j]
            for j in range(len(v))
        )
        for i in range(len(a))
    ]


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


def same_complex_line(J, i, j):
    ei = basis_vector(4, i)
    ej = basis_vector(4, j)

    Jei = matvec(J, ei)
    Jej = matvec(J, ej)

    M = columns([
        ei,
        Jei,
        ej,
        Jej,
    ])

    return rank(M) == 2


def complex_line_partition(J):
    unused = set(range(4))
    parts = []

    while unused:
        i = min(unused)

        line = sorted(
            j
            for j in range(4)
            if same_complex_line(
                J,
                i,
                j,
            )
        )

        parts.append(tuple(line))

        for j in line:
            unused.discard(j)

    return tuple(
        sorted(parts)
    )


print("PROGRESS: 1/8 load canonical objects")

aut = load(AUT)
signed = load(SIGNED)
fv = load(FV)
fk = load(FK)
a016 = load(A016)

checks = {}

aut_rows = aut[
    "measurements"
]["automorphism_rows"]

checks["AutG60_cache_has_480_automorphisms"] = (
    len(aut_rows) == 480
)

checks["AutG60_rows_are_60_point_permutations"] = all(
    len(row["permutation"]) == 60
    and sorted(
        int(x)
        for x in row["permutation"]
    ) == list(range(60))
    for row in aut_rows
)

checks["signed_block_audit_passes"] = (
    signed["audit_pass"] is True
)

checks["201FV_passes"] = (
    fv["audit_pass"] is True
)

checks["201FK_passes"] = (
    fk["audit_pass"] is True
)

checks["Audit016_passes"] = (
    a016["audit_pass"] is True
)

print("PROGRESS: 2/8 reconstruct signed block action")

blocks = signed[
    "measurements"
]["signed_blocks"]

block_key_to_index = {}

for i, block in enumerate(blocks):
    key = tuple(
        sorted(
            int(x)
            for x in block[
                "five_state_block"
            ]
        )
    )

    block_key_to_index[key] = i

carrier_to_blocks = defaultdict(list)

for i, block in enumerate(blocks):
    carrier_to_blocks[
        int(block["carrier_index"])
    ].append(i)

canonical_carrier = int(
    fv[
        "canonical_face"
    ]["carrier_index"]
)

local_blocks = sorted(
    carrier_to_blocks[
        canonical_carrier
    ]
)

checks["canonical_carrier_is_zero"] = (
    canonical_carrier == 0
)

checks["canonical_face_has_four_blocks"] = (
    len(local_blocks) == 4
)

local_position = {
    block_index: pos
    for pos, block_index in enumerate(
        local_blocks
    )
}

closure_failure_count = 0
carrier_map_failure_count = 0

stabilizer_permutations = []

for row in aut_rows:
    g = tuple(
        int(x)
        for x in row["permutation"]
    )

    image_blocks = []

    for block_index in local_blocks:
        key = tuple(
            sorted(
                g[int(u)]
                for u in blocks[
                    block_index
                ]["five_state_block"]
            )
        )

        image = block_key_to_index.get(
            key
        )

        if image is None:
            closure_failure_count += 1
            image_blocks = None
            break

        image_blocks.append(image)

    if image_blocks is None:
        continue

    target_carriers = {
        int(
            blocks[i][
                "carrier_index"
            ]
        )
        for i in image_blocks
    }

    if len(target_carriers) != 1:
        carrier_map_failure_count += 1
        continue

    target_carrier = next(
        iter(target_carriers)
    )

    if target_carrier != canonical_carrier:
        continue

    if any(
        i not in local_position
        for i in image_blocks
    ):
        carrier_map_failure_count += 1
        continue

    p = tuple(
        local_position[i]
        for i in image_blocks
    )

    stabilizer_permutations.append(p)

checks["signed_block_action_closes"] = (
    closure_failure_count == 0
)

checks["carrier_action_is_well_defined"] = (
    carrier_map_failure_count == 0
)

print("PROGRESS: 3/8 certify local D8 image")

face_stabilizer_count = len(
    stabilizer_permutations
)

perm_multiplicity = Counter(
    stabilizer_permutations
)

local_image = sorted(
    perm_multiplicity.keys()
)

order_profile = Counter(
    perm_order(p)
    for p in local_image
)

checks["face_stabilizer_order_80"] = (
    face_stabilizer_count == 80
)

checks["local_D8_image_order_8"] = (
    len(local_image) == 8
)

checks["local_action_kernel_order_10"] = all(
    count == 10
    for count in perm_multiplicity.values()
)

checks["local_D8_order_profile"] = (
    order_profile
    == Counter({
        1: 1,
        2: 5,
        4: 2,
    })
)

checks["matches_201FV_local_D8_order"] = (
    len(local_image)
    ==
    int(
        fv[
            "canonical_face"
        ]["local_D8_image_order"]
    )
)

checks["matches_201FV_kernel_order"] = (
    10
    ==
    int(
        fv[
            "canonical_face"
        ]["local_action_kernel_order"]
    )
)

print("PROGRESS: 4/8 select order-four square rotation")

order4_perms = sorted(
    p
    for p in local_image
    if perm_order(p) == 4
)

checks["exactly_two_order4_local_rotations"] = (
    len(order4_perms) == 2
)

R_perm = order4_perms[0]
R = permutation_matrix(
    R_perm
)

R2 = matmul(
    R,
    R,
)

I4 = eye(4)

checks["R_fourth_identity"] = matrix_equal(
    matmul(R2, R2),
    I4,
)

# Follow the chosen four-cycle.
cycle = [0]

for _ in range(3):
    cycle.append(
        R_perm[
            cycle[-1]
        ]
    )

checks["chosen_rotation_is_literal_four_cycle"] = (
    len(set(cycle)) == 4
    and R_perm[
        cycle[-1]
    ] == cycle[0]
)

print("PROGRESS: 5/8 derive 1 + chi + E2 decomposition")

# The square permutation representation of D8 decomposes as
#
#     1 + chi + E2.
#
# For an order-four square rotation R:
#
#     R = +1 on the trivial axis,
#     R = -1 on the chi axis,
#     R^2 = -I on E2.
#
# Choose:
#
#     t = (1,1,1,1)
#
# and an alternating R=-1 vector s along the four-cycle.
#
# The sign of s is an orientation convention only.

t = [
    Fraction(1)
    for _ in range(4)
]

s = [
    Fraction(0)
    for _ in range(4)
]

for k, i in enumerate(cycle):
    s[i] = Fraction(
        1 if k % 2 == 0 else -1
    )

Rt = matvec(R, t)
Rs = matvec(R, s)

checks["trivial_axis_R_eigenvalue_plus1"] = (
    Rt == t
)

checks["chi_axis_R_eigenvalue_minus1"] = (
    Rs
    == [
        -x
        for x in s
    ]
)

checks["trivial_and_chi_axes_orthogonal"] = (
    sum(
        t[i] * s[i]
        for i in range(4)
    )
    == 0
)

P_t = matscale(
    Fraction(1, 4),
    outer(t, t),
)

P_s = matscale(
    Fraction(1, 4),
    outer(s, s),
)

P_scalar = matadd(
    P_t,
    P_s,
)

P_E = matsub(
    I4,
    P_scalar,
)

checks["scalar_projector_rank_2"] = (
    rank(P_scalar) == 2
)

checks["E2_projector_rank_2"] = (
    rank(P_E) == 2
)

checks["R2_plus_on_scalar_plane"] = matrix_equal(
    matmul(R2, P_scalar),
    P_scalar,
)

checks["R2_minus_on_E2_plane"] = matrix_equal(
    matmul(R2, P_E),
    matscale(-1, P_E),
)

print("PROGRESS: 6/8 construct four native complex structures")

# J_S rotates the scalar plane t -> s -> -t.
#
# J_E is the selected square rotation restricted to E2.
#
# Independent orientation signs on J_S and J_E produce
# the four exact candidates already recorded by 201FV.

J_S = matscale(
    Fraction(1, 4),
    matsub(
        outer(s, t),
        outer(t, s),
    ),
)

J_E = matmul(
    R,
    P_E,
)

checks["J_S_square_minus_scalar_projector"] = matrix_equal(
    matmul(J_S, J_S),
    matscale(-1, P_scalar),
)

checks["J_E_square_minus_E2_projector"] = matrix_equal(
    matmul(J_E, J_E),
    matscale(-1, P_E),
)

checks["J_S_J_E_zero"] = matrix_equal(
    matmul(J_S, J_E),
    zeros(4, 4),
)

checks["J_E_J_S_zero"] = matrix_equal(
    matmul(J_E, J_S),
    zeros(4, 4),
)

candidate_rows = []

for scalar_orientation in (-1, +1):
    for E2_orientation in (-1, +1):
        J = matadd(
            matscale(
                scalar_orientation,
                J_S,
            ),
            matscale(
                E2_orientation,
                J_E,
            ),
        )

        J2 = matmul(
            J,
            J,
        )

        JT = transpose(J)

        orthogonal = matrix_equal(
            matmul(JT, J),
            I4,
        )

        skew = matrix_equal(
            JT,
            matscale(-1, J),
        )

        partition = complex_line_partition(
            J
        )

        candidate_rows.append({
            "scalar_orientation":
                scalar_orientation,

            "E2_orientation":
                E2_orientation,

            "J_squared_minus_identity":
                matrix_equal(
                    J2,
                    matscale(-1, I4),
                ),

            "J_orthogonal":
                orthogonal,

            "J_skew":
                skew,

            "line_partition_local_positions":
                [
                    list(pair)
                    for pair in partition
                ],

            "matrix":
                matrix_strings(J),
        })

checks["four_complex_structure_candidates"] = (
    len(candidate_rows) == 4
)

checks["all_candidates_square_to_minus_identity"] = all(
    row[
        "J_squared_minus_identity"
    ]
    for row in candidate_rows
)

checks["all_candidates_orthogonal"] = all(
    row["J_orthogonal"]
    for row in candidate_rows
)

checks["all_candidates_skew"] = all(
    row["J_skew"]
    for row in candidate_rows
)

print("PROGRESS: 7/8 classify native projective line partitions")

partition_counts = Counter(
    tuple(
        tuple(pair)
        for pair in row[
            "line_partition_local_positions"
        ]
    )
    for row in candidate_rows
)

checks["exactly_two_distinct_projective_line_partitions"] = (
    len(partition_counts) == 2
)

checks["each_partition_occurs_for_two_overall_orientations"] = all(
    count == 2
    for count in partition_counts.values()
)

block_rows_016 = {
    int(row["signed_block_index"]): row
    for row in a016[
        "native_mode_space"
    ]["block_rows"]
}

local_mode_rows = []

for pos, block_index in enumerate(
    local_blocks
):
    source = block_rows_016[
        block_index
    ]

    local_mode_rows.append({
        "local_position":
            pos,

        "signed_block_index":
            block_index,

        "sign":
            source["sign"],

        "event_anchor_label":
            source[
                "event_anchor_label"
            ],

        "five_state_block":
            source[
                "five_state_block"
            ],
    })

position_sign = {
    row["local_position"]: row["sign"]
    for row in local_mode_rows
}

position_block = {
    row["local_position"]:
        row["signed_block_index"]
    for row in local_mode_rows
}

position_anchor = {
    row["local_position"]:
        row["event_anchor_label"]
    for row in local_mode_rows
}

partition_rows = []

cross_sign_failure_count = 0
positive_separation_failure_count = 0
negative_separation_failure_count = 0

for partition, multiplicity in sorted(
    partition_counts.items()
):
    lines = []

    for pair in partition:
        signs = [
            position_sign[i]
            for i in pair
        ]

        if set(signs) != {
            "positive",
            "negative",
        }:
            cross_sign_failure_count += 1

        lines.append({
            "local_positions":
                list(pair),

            "signed_block_indices":
                [
                    position_block[i]
                    for i in pair
                ],

            "signs":
                signs,

            "event_anchor_labels":
                [
                    position_anchor[i]
                    for i in pair
                ],
        })

    positive_positions = [
        i
        for i in range(4)
        if position_sign[i]
        == "positive"
    ]

    negative_positions = [
        i
        for i in range(4)
        if position_sign[i]
        == "negative"
    ]

    def same_partition_line(i, j):
        return any(
            i in pair
            and j in pair
            for pair in partition
        )

    if same_partition_line(
        positive_positions[0],
        positive_positions[1],
    ):
        positive_separation_failure_count += 1

    if same_partition_line(
        negative_positions[0],
        negative_positions[1],
    ):
        negative_separation_failure_count += 1

    partition_rows.append({
        "multiplicity_among_four_J_candidates":
            multiplicity,

        "lines":
            lines,
    })

checks["every_complex_line_pairs_one_positive_one_negative_mode"] = (
    cross_sign_failure_count == 0
)

checks["positive_event_modes_are_distinct_complex_lines_for_both_phases"] = (
    positive_separation_failure_count == 0
)

checks["negative_event_modes_are_distinct_complex_lines_for_both_phases"] = (
    negative_separation_failure_count == 0
)

checks["agrees_with_201FV_positive_orthogonality"] = (
    fv[
        "checks"
    ][
        "positive_pair_hermitian_overlap_zero"
    ]
    is True
    and positive_separation_failure_count == 0
)

print("PROGRESS: 8/8 classify")

failed_checks = [
    name
    for name, passed in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "canonical_native_face_D8_permutation_module_derives_"
    "exactly_two_relative_complex_line_partitions_each_"
    "pairing_one_positive_and_one_negative_event_mode"
    if audit_pass
    else
    "native_face_projective_line_partition_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_native_face_projective_line_partition_017",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "canonical_face": {
        "carrier_index":
            canonical_carrier,

        "signed_block_indices":
            local_blocks,

        "mode_rows":
            local_mode_rows,
    },

    "native_local_D8": {
        "face_stabilizer_order":
            face_stabilizer_count,

        "image_order":
            len(local_image),

        "kernel_order":
            10,

        "order_profile": {
            str(k): v
            for k, v in sorted(
                order_profile.items()
            )
        },

        "image_permutations":
            [
                {
                    "permutation":
                        list(p),

                    "order":
                        perm_order(p),

                    "native_realization_count":
                        perm_multiplicity[p],
                }
                for p in local_image
            ],

        "chosen_order4_rotation":
            list(R_perm),

        "chosen_cycle":
            cycle,

        "orientation_choice_is_presentation_only":
            True,
    },

    "representation_decomposition": {
        "real_module":
            "1 + chi + E2",

        "trivial_axis":
            [frac_string(x) for x in t],

        "chi_axis":
            [frac_string(x) for x in s],

        "scalar_plane_dimension":
            2,

        "E2_plane_dimension":
            2,

        "J_S":
            matrix_strings(J_S),

        "J_E":
            matrix_strings(J_E),
    },

    "complex_structure_candidates":
        candidate_rows,

    "projective_line_result": {
        "distinct_partition_count":
            len(partition_rows),

        "partitions":
            partition_rows,

        "all_lines_pair_opposite_event_signs":
            cross_sign_failure_count == 0,

        "positive_modes_distinct_in_both_partitions":
            positive_separation_failure_count == 0,

        "phase_assignment_status":
            (
                "two native relative-orientation partitions "
                "derived on canonical face; assignment of "
                "epsilon=+1 versus epsilon=-1 is defined only "
                "up to the already-admitted global phase-label swap"
            ),
    },

    "earned_statement": (
        "The four native signed event-mode basis objects on the "
        "canonical face carry the exact four-point permutation "
        "representation of a D8 local image of the order-80 "
        "face stabilizer with kernel order ten. Its real "
        "permutation module decomposes intrinsically as "
        "1 + chi + E2. The order-four square rotation determines "
        "the E2 complex orientation while the trivial and chi "
        "axes determine the scalar-plane complex orientation. "
        "The four sign choices produce the four Program-01 "
        "orthogonal complex-structure candidates but only two "
        "distinct projective complex-line partitions. Every "
        "line in both partitions contains exactly one positive "
        "and one negative native event mode; consequently the "
        "two positive event modes remain on distinct complex "
        "lines for either relative complex orientation, in "
        "agreement with 201FV. This derives the canonical-face "
        "projective-line content without selecting a transverse "
        "Pauli frame or a raw Hilbert transport matrix."
    ),

    "checks":
        checks,

    "boundary": {
        "canonical_face_projective_line_partitions_derived":
            audit_pass,

        "raw_Pauli_frame_not_selected":
            True,

        "raw_Hilbert_transport_matrix_not_selected":
            True,

        "common_g_not_identified_with_S1":
            True,

        "r1_not_identified_with_S1":
            True,

        "epsilon_partition_assignment_only_up_to_global_label_swap":
            True,

        "six_face_global_transport_not_yet_closed":
            True,

        "full_48_entry_ell_not_yet_serialized":
            True,

        "joint_wedge_census_not_yet_run":
            True,

        "no_probability_law":
            True,

        "no_Born_rule":
            True,

        "no_CHSH_claim":
            True,
    },

    "next_gate": (
        "Transport the two canonical-face projective-line "
        "partitions over all six native face carriers using "
        "the current global/local phase bridge. Test that the "
        "kernel of the phase character makes the transported "
        "partition independent of transporter choice, leaving "
        "only one global epsilon-label swap. Then serialize "
        "ell on all forty-eight phase-mode pairs."
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

note = f"""# Native face projective-line partition 017

## Result

Audit pass:

    {audit_pass}

The canonical face stabilizer contributes

    {face_stabilizer_count}

native automorphisms.

Its action on the four signed event-mode basis objects has image size

    {len(local_image)}

and kernel size

    10.

The local image has order profile

    {dict(sorted(order_profile.items()))}

which is the D8 square action.

## Native real-module decomposition

Choose either native order-four square rotation.

The four-real-dimensional permutation module splits intrinsically as

    1 + chi + E2.

The trivial and chi axes form the scalar plane.

The remaining two-plane is E2.

A scalar-plane complex structure J_S and the square rotation restricted
to E2 give J_E.

Independent orientation signs produce the four Program-01 complex
structures

    +/- J_S direct_sum +/- J_E.

All four square to -I, are orthogonal, and are skew.

## Projective lines

Those four complex structures produce exactly

    {len(partition_rows)}

distinct partitions of the four native real event-mode basis objects
into complex lines.

Each partition occurs twice because simultaneous reversal of both
complex orientations changes J to -J without changing CP1 lines.

For every derived line partition:

    each complex line contains one positive event mode
    and one negative event mode;

    the two positive event modes lie on distinct complex lines;

    the two negative event modes lie on distinct complex lines.

This agrees with Program-01 201FV, where the positive event pair is
Hermitian-orthogonal for every native relative-complex-orientation
choice.

## Boundary

No transverse Pauli frame has been selected.

No raw Hilbert transport matrix has been selected.

The two canonical-face phase partitions are derived only up to the
already-admitted global naming swap of epsilon=+1 and epsilon=-1.

The next step is to transport this two-partition object over all six
faces and serialize the full 48-entry ell map.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print(
    "FACE_STABILIZER_COUNT:",
    face_stabilizer_count,
)
print(
    "LOCAL_D8_IMAGE_ORDER:",
    len(local_image),
)
print(
    "LOCAL_D8_ORDER_PROFILE:",
    dict(sorted(
        order_profile.items()
    )),
)
print(
    "ORDER4_ROTATION_COUNT:",
    len(order4_perms),
)
print(
    "DISTINCT_PROJECTIVE_LINE_PARTITIONS:",
    len(partition_rows),
)
print(
    "CROSS_SIGN_LINE_FAILURE_COUNT:",
    cross_sign_failure_count,
)
print(
    "POSITIVE_SEPARATION_FAILURE_COUNT:",
    positive_separation_failure_count,
)
print(
    "NEGATIVE_SEPARATION_FAILURE_COUNT:",
    negative_separation_failure_count,
)
print(
    "FULL_ELL_SERIALIZED:",
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
