#!/usr/bin/env python3

from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path
import json

HERE = Path(__file__).resolve().parents[2]

P41 = (
    Path.home()
    / "dev/cori/research/mathematics"
    / "41-order-4-dodecahedral-residue"
)

SIGNED = (
    P41
    / "artifacts/json"
    / "signed_24_face_block_closure_audit_025.json"
)

AUT = (
    P41
    / "artifacts/json"
    / "native_g60_automorphism_deck_cache.v1.json"
)

A009 = (
    HERE
    / "artifacts/json"
    / "epr_native_face_boundary_incidence_map_009.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_g1800_joint_boundary_incidence_quotient_010.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_g1800_joint_boundary_incidence_quotient_010.md"
)


def load(path):
    return json.loads(path.read_text())


def digest_json(obj):
    raw = json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("ascii")

    return sha256(raw).hexdigest()


print("PROGRESS: 1/8 load canonical inputs")

signed = load(SIGNED)
aut = load(AUT)
a009 = load(A009)

checks = {}

checks["audit_009_passed"] = (
    a009["audit_pass"] is True
)

checks["audit_009_induced_operator_S1"] = (
    a009[
        "native_a_action"
    ]["induced_face_operator"]
    == "S1"
)

blocks = signed["measurements"]["signed_blocks"]

checks["signed_block_count_24"] = (
    len(blocks) == 24
)

print("PROGRESS: 2/8 reconstruct one-wing incidence domain")

incidences = []
state_to_blocks = defaultdict(list)
carrier_to_blocks = defaultdict(list)

for block_index, block in enumerate(blocks):
    carrier = int(block["carrier_index"])
    sign = block["sign"]

    carrier_to_blocks[carrier].append(
        block_index
    )

    for state in block["five_state_block"]:
        state = int(state)

        row = (
            state,
            block_index,
        )

        incidences.append(row)
        state_to_blocks[state].append(
            block_index
        )

checks["one_wing_incidence_count_120"] = (
    len(incidences) == 120
)

checks["every_state_has_two_boundary_incidences"] = all(
    len(state_to_blocks[u]) == 2
    for u in range(60)
)

checks["four_blocks_per_carrier"] = all(
    len(carrier_to_blocks[c]) == 4
    for c in range(6)
)

print("PROGRESS: 3/8 reconstruct native a action")

a_index = int(
    a009[
        "native_a_action"
    ]["automorphism_index"]
)

aut_rows = {
    int(row["automorphism_index"]): row
    for row in aut["measurements"]["automorphism_rows"]
}

a = tuple(
    int(v)
    for v in aut_rows[a_index]["permutation"]
)

checks["native_a_has_60_points"] = (
    len(a) == 60
)

checks["native_a_is_fixed_point_free"] = all(
    a[u] != u
    for u in range(60)
)

checks["native_a_is_involution"] = all(
    a[a[u]] == u
    for u in range(60)
)

block_by_states = {}

for i, block in enumerate(blocks):
    key = tuple(
        sorted(
            int(v)
            for v in block["five_state_block"]
        )
    )

    block_by_states[key] = i


def block_image(block_index):
    image_key = tuple(
        sorted(
            a[int(v)]
            for v in blocks[
                block_index
            ]["five_state_block"]
        )
    )

    return block_by_states[
        image_key
    ]


a_block = {
    i: block_image(i)
    for i in range(24)
}

checks["a_block_action_is_involution"] = all(
    a_block[a_block[i]] == i
    for i in range(24)
)

checks["a_block_action_has_no_fixed_blocks"] = all(
    a_block[i] != i
    for i in range(24)
)

checks["a_preserves_block_carrier"] = all(
    int(
        blocks[
            a_block[i]
        ]["carrier_index"]
    )
    ==
    int(
        blocks[i]["carrier_index"]
    )
    for i in range(24)
)

checks["a_preserves_block_sign"] = all(
    blocks[
        a_block[i]
    ]["sign"]
    ==
    blocks[i]["sign"]
    for i in range(24)
)

incidence_set = set(
    incidences
)

checks["a_closes_on_incidence_domain"] = all(
    (
        a[u],
        a_block[B],
    )
    in incidence_set
    for u, B in incidences
)

print("PROGRESS: 4/8 construct native G1800")

def diag_state_rep(u, v):
    return min(
        (u, v),
        (a[u], a[v]),
    )


joint_state_reps = sorted({
    diag_state_rep(u, v)
    for u in range(60)
    for v in range(60)
})

joint_state_index = {
    rep: i
    for i, rep in enumerate(
        joint_state_reps
    )
}


def joint_state_class(u, v):
    return joint_state_index[
        diag_state_rep(u, v)
    ]


checks["G1800_state_count_1800"] = (
    len(joint_state_reps) == 1800
)

print("PROGRESS: 5/8 construct joint boundary quotient")

raw_joint_boundary_count = (
    len(incidences)
    * len(incidences)
)

checks["raw_joint_boundary_count_14400"] = (
    raw_joint_boundary_count == 14400
)


def diag_incidence_image(row):
    u, B, v, C = row

    return (
        a[u],
        a_block[B],
        a[v],
        a_block[C],
    )


def diag_incidence_rep(row):
    image = diag_incidence_image(
        row
    )

    return min(
        row,
        image,
    )


joint_boundary_classes = set()

for u, B in incidences:
    for v, C in incidences:
        row = (
            u,
            B,
            v,
            C,
        )

        joint_boundary_classes.add(
            diag_incidence_rep(row)
        )

joint_boundary_reps = sorted(
    joint_boundary_classes
)

checks["joint_boundary_class_count_7200"] = (
    len(joint_boundary_reps) == 7200
)

checks["diagonal_incidence_action_is_free"] = all(
    diag_incidence_image(row)
    != row
    for row in joint_boundary_reps
)

print("PROGRESS: 6/8 project joint boundary to G1800")

state_to_boundary_classes = defaultdict(
    list
)

for boundary_id, row in enumerate(
    joint_boundary_reps
):
    u, B, v, C = row

    sid = joint_state_class(
        u,
        v,
    )

    state_to_boundary_classes[
        sid
    ].append(
        boundary_id
    )

boundary_fiber_profile = Counter(
    len(rows)
    for rows in state_to_boundary_classes.values()
)

checks["all_1800_joint_states_have_boundary_fibers"] = (
    len(state_to_boundary_classes)
    == 1800
)

checks["every_G1800_state_has_four_boundary_presentations"] = (
    boundary_fiber_profile
    == Counter({
        4: 1800,
    })
)

sign_pair_failure_count = 0
carrier_pair_failure_count = 0

expected_sign_pairs = {
    ("positive", "positive"),
    ("positive", "negative"),
    ("negative", "positive"),
    ("negative", "negative"),
}

for sid, boundary_ids in (
    state_to_boundary_classes.items()
):
    sign_pairs = []
    carrier_pairs = []

    for bid in boundary_ids:
        u, B, v, C = (
            joint_boundary_reps[bid]
        )

        sign_pairs.append(
            (
                blocks[B]["sign"],
                blocks[C]["sign"],
            )
        )

        carrier_pairs.append(
            (
                int(
                    blocks[B][
                        "carrier_index"
                    ]
                ),
                int(
                    blocks[C][
                        "carrier_index"
                    ]
                ),
            )
        )

    if set(sign_pairs) != expected_sign_pairs:
        sign_pair_failure_count += 1

    if len(set(carrier_pairs)) != 4:
        carrier_pair_failure_count += 1

checks["every_G1800_state_has_all_four_sign_pairs"] = (
    sign_pair_failure_count == 0
)

checks["every_G1800_state_has_four_distinct_carrier_pairs"] = (
    carrier_pair_failure_count == 0
)

print("PROGRESS: 7/8 verify induced S1 tensor quotient")

# On a fixed ordered pair of face carriers,
#
# V_A tensor_R V_B
#
# has 4 x 4 = 16 canonical basis tensors.
#
# Diagonal a acts by
#
#     (B,C) -> (aB,aC),
#
# which is precisely the basis permutation of
#
#     S1_A tensor S1_B.
#
# Since a has no fixed block on either local face, every
# local tensor basis orbit has size 2, giving 8 quotient
# basis classes.

carrier_pair_profiles = {}

for ca in range(6):
    for cb in range(6):
        BA = sorted(
            carrier_to_blocks[ca]
        )

        BB = sorted(
            carrier_to_blocks[cb]
        )

        basis_pairs = [
            (B, C)
            for B in BA
            for C in BB
        ]

        basis_set = set(
            basis_pairs
        )

        seen = set()
        orbit_sizes = []

        closure_ok = True

        for pair in basis_pairs:
            image = (
                a_block[pair[0]],
                a_block[pair[1]],
            )

            if image not in basis_set:
                closure_ok = False

        for pair in basis_pairs:
            if pair in seen:
                continue

            image = (
                a_block[pair[0]],
                a_block[pair[1]],
            )

            orbit = {
                pair,
                image,
            }

            seen.update(
                orbit
            )

            orbit_sizes.append(
                len(orbit)
            )

        profile = Counter(
            orbit_sizes
        )

        carrier_pair_profiles[
            f"{ca},{cb}"
        ] = {
            "basis_count":
                len(basis_pairs),

            "closure_ok":
                closure_ok,

            "orbit_profile": {
                str(k): v
                for k, v in sorted(
                    profile.items()
                )
            },

            "coinvariant_real_dimension":
                len(orbit_sizes),

            "fixed_real_dimension":
                len(orbit_sizes),
        }

checks["all_36_carrier_pairs_close_under_diagonal_a"] = all(
    row["closure_ok"]
    for row in carrier_pair_profiles.values()
)

checks["all_local_tensor_basis_actions_are_eight_two_cycles"] = all(
    row["orbit_profile"]
    == {"2": 8}
    for row in carrier_pair_profiles.values()
)

checks["native_simultaneous_S1_coinvariant_dimension_is_8"] = all(
    row[
        "coinvariant_real_dimension"
    ]
    == 8
    for row in carrier_pair_profiles.values()
)

checks["native_simultaneous_S1_fixed_dimension_is_8"] = all(
    row[
        "fixed_real_dimension"
    ]
    == 8
    for row in carrier_pair_profiles.values()
)

print("PROGRESS: 8/8 classify")

failed_checks = [
    name
    for name, passed in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "native_G1800_joint_boundary_incidence_quotient_has_"
    "7200_classes_four_presentations_per_state_and_exact_"
    "S1_tensor_S1_eight_dimensional_local_quotient"
    if audit_pass
    else
    "g1800_joint_boundary_incidence_quotient_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_g1800_joint_boundary_incidence_quotient_010",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "one_wing": {
        "incidence_count":
            len(incidences),

        "domain":
            "I_F = {(u,B): u in signed face block B}",
    },

    "two_wing_raw_domain": {
        "definition":
            "I_F x I_F",

        "raw_count":
            raw_joint_boundary_count,
    },

    "diagonal_quotient": {
        "kernel":
            "a",

        "relation":
            (
                "(u,B;v,C) ~ "
                "(a u,a B;a v,a C)"
            ),

        "class_count":
            len(
                joint_boundary_reps
            ),

        "G1800_state_count":
            len(
                joint_state_reps
            ),

        "boundary_presentations_per_G1800_state":
            4,

        "sign_pair_set":
            [
                "++",
                "+-",
                "-+",
                "--",
            ],
    },

    "local_tensor_quotient": {
        "raw_real_dimension":
            16,

        "induced_operator":
            "S1_A tensor S1_B",

        "basis_action":
            (
                "(B,C) -> "
                "(aB,aC)"
            ),

        "orbit_profile":
            "eight 2-cycles",

        "coinvariant_real_dimension":
            8,

        "fixed_real_dimension":
            8,

        "carrier_pair_count":
            36,
    },

    "earned_statement": (
        "The native two-wing boundary domain is the diagonal-a "
        "quotient of the product of the two 120-element signed "
        "state-face incidence domains. It contains exactly "
        "7200 boundary classes over 1800 G1800 states, with "
        "exactly four boundary presentations per common state. "
        "Those four presentations realize all four local sign "
        "pairs ++,+-,-+,-- and four distinct ordered face-"
        "carrier pairs. On every ordered pair of local four-"
        "real-dimensional face modules, the same diagonal-a "
        "quotient induces exactly the simultaneous Program-01 "
        "operator S1_A tensor S1_B. Its sixteen canonical "
        "basis tensors form eight two-cycles, yielding an "
        "eight-real-dimensional native coinvariant/fixed "
        "sector before complex balancing."
    ),

    "checks":
        checks,

    "boundary": {
        "full_two_wing_native_boundary_domain_constructed":
            True,

        "simultaneous_S1_structure_is_native_not_imposed":
            True,

        "complex_balancing_not_yet_applied":
            True,

        "Audit003_RP3_not_yet_claimed_recovered":
            True,

        "canonical_rank_one_projector_not_yet_executed":
            True,

        "no_Born_rule":
            True,

        "no_probability_table":
            True,

        "no_CHSH_claim":
            True,
    },

    "next_gate": (
        "Apply the Program-01 complex structures to the "
        "eight-real-dimensional simultaneous-S1 boundary "
        "quotient. Impose the same-chirality balancing law "
        "J_A tensor I = I tensor J_B and test whether the "
        "resulting native boundary space has real dimension "
        "four and is exactly the Audit-003 C_AB fixed real "
        "form, recovering RP3 from native incidence data."
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

note = f"""# G1800 joint boundary incidence quotient 010

## Result

Audit pass:

    {audit_pass}

The one-wing native boundary domain has

    120

signed state-face incidences.

The two-wing raw product therefore contains

    14400

ordered incidence pairs.

Quotienting by the diagonal native deck involution

    (u,B;v,C)
      ~
    (a u,a B;a v,a C)

gives exactly

    {len(joint_boundary_reps)}

joint boundary-incidence classes.

## Projection to G1800

The same diagonal action on the underlying state pair gives

    1800

G1800 states.

Every G1800 state has exactly

    4

joint boundary presentations.

They realize exactly the four local sign combinations

    ++
    +-
    -+
    --.

The four presentations also use four distinct ordered face-carrier
pairs.

## Local tensor action

Each ordered pair of local face carriers supplies

    4 x 4 = 16

canonical real tensor basis objects.

The native diagonal-a action is

    (B,C) -> (aB,aC).

Audit 009 identifies the one-wing induced action as S1.

Therefore the two-wing action is exactly

    S1_A tensor S1_B.

For every one of the 36 ordered carrier pairs, the sixteen tensor basis
objects split into

    eight 2-cycles.

Thus the native simultaneous-S1 coinvariant space has real dimension

    8.

Over R this has the same dimension as the fixed sector obtained by
symmetrizing each two-cycle.

## Meaning

The simultaneous chiral compatibility structure has now emerged directly
from the native G1800 boundary-incidence quotient.

It was not imposed as an abstract Hilbert constraint.

The next step is to introduce the already-earned Program-01 complex
structures and test whether complex balancing reduces this native
eight-real-dimensional space to the four-real-dimensional RP3 found
abstractly in Audit 003.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print(
    "RAW_JOINT_BOUNDARY_COUNT:",
    raw_joint_boundary_count,
)
print(
    "JOINT_BOUNDARY_CLASS_COUNT:",
    len(joint_boundary_reps),
)
print(
    "G1800_STATE_COUNT:",
    len(joint_state_reps),
)
print(
    "BOUNDARY_FIBER_PROFILE:",
    dict(sorted(
        boundary_fiber_profile.items()
    )),
)
print(
    "SIGN_PAIR_FAILURE_COUNT:",
    sign_pair_failure_count,
)
print(
    "CARRIER_PAIR_FAILURE_COUNT:",
    carrier_pair_failure_count,
)
print(
    "LOCAL_TENSOR_QUOTIENT_REAL_DIM:",
    8,
)
print(
    "INDUCED_OPERATOR:",
    "S1_A tensor S1_B",
)
print("FAILED_CHECK_COUNT:", len(failed_checks))
print("FAILED_CHECKS:", failed_checks)
print("JSON_OUT:", JSON_OUT)
print("NOTE_OUT:", NOTE_OUT)
print(
    "JSON_SHA256:",
    sha256(JSON_OUT.read_bytes()).hexdigest(),
)
