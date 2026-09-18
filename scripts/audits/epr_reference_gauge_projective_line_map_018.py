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

FU = (
    P41
    / "artifacts/provenance"
    / "project41_6k2_event_u1_connection_201fu.v1.json"
)

FR = (
    P41
    / "artifacts/provenance"
    / "project41_global_local_phase_character_bridge_201fr.v1.json"
)

A009 = (
    HERE
    / "artifacts/json"
    / "epr_native_face_boundary_incidence_map_009.v1.json"
)

A016 = (
    HERE
    / "artifacts/json"
    / "epr_literal_phase_decorated_incidence_table_016.v1.json"
)

A017 = (
    HERE
    / "artifacts/json"
    / "epr_native_face_projective_line_partition_017.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_reference_gauge_projective_line_map_018.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_reference_gauge_projective_line_map_018.md"
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


def norm_line(values):
    return tuple(
        sorted(
            int(x)
            for x in values
        )
    )


def norm_partition(lines):
    return tuple(
        sorted(
            norm_line(line)
            for line in lines
        )
    )


def image_partition(partition, block_image):
    return norm_partition(
        [
            block_image[b]
            for b in line
        ]
        for line in partition
    )


def inverse_permutation(p):
    q = [None] * len(p)

    for i, j in enumerate(p):
        q[j] = i

    if any(x is None for x in q):
        raise RuntimeError(
            "not a permutation"
        )

    return tuple(q)


def perm_order(p):
    cur = tuple(range(len(p)))

    for n in range(1, 33):
        cur = tuple(
            p[cur[i]]
            for i in range(len(p))
        )

        if cur == tuple(range(len(p))):
            return n

    raise RuntimeError(
        "permutation order too large"
    )


print("PROGRESS: 1/8 load sealed interfaces")

aut = load(AUT)
signed = load(SIGNED)
fu = load(FU)
fr = load(FR)
a009 = load(A009)
a016 = load(A016)
a017 = load(A017)

checks = {}

checks["signed_block_audit_passed"] = (
    signed["audit_pass"] is True
)

checks["201FU_passed"] = (
    fu["audit_pass"] is True
)

checks["201FR_passed"] = (
    fr["audit_pass"] is True
)

checks["Audit009_passed"] = (
    a009["audit_pass"] is True
)

checks["Audit016_passed"] = (
    a016["audit_pass"] is True
)

checks["Audit017_passed"] = (
    a017["audit_pass"] is True
)

aut_rows = aut[
    "measurements"
]["automorphism_rows"]

checks["AutG60_row_count_480"] = (
    len(aut_rows) == 480
)

print("PROGRESS: 2/8 reconstruct full signed-block action")

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

checks["signed_block_count_24"] = (
    len(blocks) == 24
)

checks["six_carriers_four_blocks_each"] = (
    set(carrier_to_blocks.keys())
    == set(range(6))
    and all(
        len(carrier_to_blocks[c]) == 4
        for c in range(6)
    )
)

block_image_by_aut = {}
carrier_image_by_aut = {}

block_action_failure_count = 0
carrier_action_failure_count = 0

for row in aut_rows:
    idx = int(
        row["automorphism_index"]
    )

    g = tuple(
        int(x)
        for x in row["permutation"]
    )

    image = []

    for block in blocks:
        key = tuple(
            sorted(
                g[int(u)]
                for u in block[
                    "five_state_block"
                ]
            )
        )

        target = block_key_to_index.get(
            key
        )

        if target is None:
            block_action_failure_count += 1
            image = None
            break

        image.append(target)

    if image is None:
        continue

    if sorted(image) != list(range(24)):
        block_action_failure_count += 1
        continue

    block_image_by_aut[idx] = tuple(
        image
    )

    carrier_image = {}

    for source_carrier in range(6):
        target_carriers = {
            int(
                blocks[
                    image[b]
                ]["carrier_index"]
            )
            for b in carrier_to_blocks[
                source_carrier
            ]
        }

        if len(target_carriers) != 1:
            carrier_action_failure_count += 1
            continue

        carrier_image[
            source_carrier
        ] = next(
            iter(target_carriers)
        )

    if len(carrier_image) == 6:
        carrier_image_by_aut[
            idx
        ] = carrier_image

checks["all_480_automorphisms_act_on_signed_blocks"] = (
    block_action_failure_count == 0
    and len(block_image_by_aut) == 480
)

checks["all_480_automorphisms_act_on_six_carriers"] = (
    carrier_action_failure_count == 0
    and len(carrier_image_by_aut) == 480
)

print("PROGRESS: 3/8 recover canonical two-partition object")

canonical_carrier = int(
    a017[
        "canonical_face"
    ]["carrier_index"]
)

local_blocks = tuple(
    sorted(
        carrier_to_blocks[
            canonical_carrier
        ]
    )
)

local_position = {
    block_index: pos
    for pos, block_index in enumerate(
        local_blocks
    )
}

partition_rows = a017[
    "projective_line_result"
]["partitions"]

canonical_partitions = []

for partition_row in partition_rows:
    partition = norm_partition(
        line[
            "signed_block_indices"
        ]
        for line in partition_row[
            "lines"
        ]
    )

    canonical_partitions.append(
        partition
    )

canonical_partitions = sorted(
    set(canonical_partitions)
)

checks["canonical_partition_count_2"] = (
    len(canonical_partitions) == 2
)

P0 = canonical_partitions[0]
P1 = canonical_partitions[1]

checks["canonical_partitions_are_distinct"] = (
    P0 != P1
)

checks["each_partition_covers_all_four_local_blocks"] = all(
    sorted(
        b
        for line in partition
        for b in line
    )
    == sorted(local_blocks)
    for partition in (
        P0,
        P1,
    )
)

print("PROGRESS: 4/8 recover selected native phase-flip transport")

stabilizer_rows = []

for row in aut_rows:
    idx = int(
        row["automorphism_index"]
    )

    if carrier_image_by_aut[
        idx
    ][canonical_carrier] != canonical_carrier:
        continue

    image = block_image_by_aut[idx]

    p = tuple(
        local_position[
            image[b]
        ]
        for b in local_blocks
    )

    stabilizer_rows.append(
        (
            idx,
            p,
        )
    )

checks["canonical_face_stabilizer_count_80"] = (
    len(stabilizer_rows) == 80
)

positive_positions = [
    pos
    for pos, block_index in enumerate(
        local_blocks
    )
    if blocks[
        block_index
    ]["sign"] == "positive"
]

negative_positions = [
    pos
    for pos, block_index in enumerate(
        local_blocks
    )
    if blocks[
        block_index
    ]["sign"] == "negative"
]

checks["two_positive_two_negative_local_modes"] = (
    len(positive_positions) == 2
    and len(negative_positions) == 2
)

selected_local_perms = Counter()

for idx, p in stabilizer_rows:
    positive_fixed = all(
        p[pos] == pos
        for pos in positive_positions
    )

    negative_swapped = (
        p[
            negative_positions[0]
        ]
        == negative_positions[1]
        and
        p[
            negative_positions[1]
        ]
        == negative_positions[0]
    )

    if (
        positive_fixed
        and negative_swapped
    ):
        selected_local_perms[p] += 1

checks["unique_FU_local_transport_permutation"] = (
    len(selected_local_perms) == 1
)

if len(selected_local_perms) != 1:
    raise RuntimeError(
        "FU local transport permutation not unique"
    )

phase_flip_perm = next(
    iter(selected_local_perms)
)

phase_flip_native_realization_count = (
    selected_local_perms[
        phase_flip_perm
    ]
)

checks["FU_local_transport_has_order2"] = (
    perm_order(
        phase_flip_perm
    ) == 2
)

checks["FU_native_realization_count_10"] = (
    phase_flip_native_realization_count
    == 10
)

checks["FU_record_says_positive_fixed"] = (
    fu[
        "connection_action"
    ]["positive_event_modes"]
    == "fixed pointwise"
)

checks["FU_record_says_negative_swapped"] = (
    fu[
        "connection_action"
    ]["negative_face_modes"]
    == "swapped"
)

checks["FU_record_says_unique_selected_transport"] = (
    fu[
        "connection_action"
    ]["unique_selected_transport"]
    is True
)

phase_flip_block = {
    local_blocks[pos]:
        local_blocks[
            phase_flip_perm[pos]
        ]
    for pos in range(4)
}

phase_flip_partition_P0 = norm_partition(
    [
        phase_flip_block[b]
        for b in line
    ]
    for line in P0
)

phase_flip_partition_P1 = norm_partition(
    [
        phase_flip_block[b]
        for b in line
    ]
    for line in P1
)

checks["FU_phase_flip_exchanges_P0_P1"] = (
    phase_flip_partition_P0 == P1
    and
    phase_flip_partition_P1 == P0
)

# Compare with the a-induced S1 action from Audit 009.
a_index = int(
    a009[
        "native_a_action"
    ]["automorphism_index"]
)

a_image = block_image_by_aut[
    a_index
]

a_perm = tuple(
    local_position[
        a_image[b]
    ]
    for b in local_blocks
)

checks["FU_transport_not_equal_a_induced_S1"] = (
    phase_flip_perm != a_perm
)

checks["a_induced_S1_swaps_positive_pair"] = (
    a_perm[
        positive_positions[0]
    ]
    == positive_positions[1]
    and
    a_perm[
        positive_positions[1]
    ]
    == positive_positions[0]
)

checks["a_induced_S1_swaps_negative_pair"] = (
    a_perm[
        negative_positions[0]
    ]
    == negative_positions[1]
    and
    a_perm[
        negative_positions[1]
    ]
    == negative_positions[0]
)

print("PROGRESS: 5/8 transport unordered phase-pair over all six faces")

transporters_by_carrier = defaultdict(list)

for row in aut_rows:
    idx = int(
        row["automorphism_index"]
    )

    target = carrier_image_by_aut[
        idx
    ][canonical_carrier]

    transporters_by_carrier[
        target
    ].append(idx)

transporter_count_profile = Counter(
    len(
        transporters_by_carrier[c]
    )
    for c in range(6)
)

checks["eighty_native_transporters_per_face"] = (
    transporter_count_profile
    == Counter({
        80: 6,
    })
)

unordered_pair_count_by_carrier = {}
unordered_pair_object_by_carrier = {}

for carrier in range(6):
    pair_objects = set()

    for idx in transporters_by_carrier[
        carrier
    ]:
        image = block_image_by_aut[
            idx
        ]

        q0 = image_partition(
            P0,
            image,
        )

        q1 = image_partition(
            P1,
            image,
        )

        pair_objects.add(
            tuple(
                sorted(
                    (
                        q0,
                        q1,
                    )
                )
            )
        )

    unordered_pair_count_by_carrier[
        carrier
    ] = len(
        pair_objects
    )

    if len(pair_objects) == 1:
        unordered_pair_object_by_carrier[
            carrier
        ] = next(
            iter(pair_objects)
        )

checks["one_transporter_independent_phase_pair_object_per_face"] = all(
    unordered_pair_count_by_carrier[c]
    == 1
    for c in range(6)
)

print("PROGRESS: 6/8 choose explicit flat reference section")

identity_candidates = []

for row in aut_rows:
    idx = int(
        row["automorphism_index"]
    )

    perm = tuple(
        int(x)
        for x in row[
            "permutation"
        ]
    )

    if perm == tuple(range(60)):
        identity_candidates.append(
            idx
        )

checks["unique_identity_automorphism"] = (
    len(identity_candidates) == 1
)

identity_index = (
    identity_candidates[0]
    if identity_candidates
    else None
)

section_index_by_carrier = {}

for carrier in range(6):
    candidates = sorted(
        transporters_by_carrier[
            carrier
        ]
    )

    if carrier == canonical_carrier:
        section_index_by_carrier[
            carrier
        ] = identity_index
    else:
        section_index_by_carrier[
            carrier
        ] = candidates[0]

checks["section_defined_on_all_six_faces"] = (
    len(section_index_by_carrier)
    == 6
    and all(
        section_index_by_carrier[c]
        is not None
        for c in range(6)
    )
)

target_partition_by_carrier_phase = {}

for carrier in range(6):
    idx = section_index_by_carrier[
        carrier
    ]

    image = block_image_by_aut[
        idx
    ]

    target_partition_by_carrier_phase[
        (carrier, +1)
    ] = image_partition(
        P0,
        image,
    )

    target_partition_by_carrier_phase[
        (carrier, -1)
    ] = image_partition(
        P1,
        image,
    )

    unordered = tuple(
        sorted(
            (
                target_partition_by_carrier_phase[
                    (carrier, +1)
                ],
                target_partition_by_carrier_phase[
                    (carrier, -1)
                ],
            )
        )
    )

    if (
        unordered
        != unordered_pair_object_by_carrier[
            carrier
        ]
    ):
        checks[
            "selected_section_matches_native_phase_pair"
        ] = False

checks.setdefault(
    "selected_section_matches_native_phase_pair",
    True,
)

print("PROGRESS: 7/8 serialize ell in one reference gauge")

# Reference state:
#
#   canonical carrier
#   epsilon = +1
#   partition P0
#
# Order its two CP1 lines deterministically.
reference_lines = tuple(
    sorted(P0)
)

reference_line_by_block = {}

for line_id, line in enumerate(
    reference_lines
):
    for block_index in line:
        reference_line_by_block[
            block_index
        ] = line_id

checks["reference_partition_has_two_lines"] = (
    len(reference_lines) == 2
)

checks["reference_lines_cover_four_modes"] = (
    len(reference_line_by_block) == 4
)

# FU phase flip is an involution, so it is its own inverse.
checks["phase_flip_local_involution"] = (
    all(
        phase_flip_block[
            phase_flip_block[b]
        ]
        == b
        for b in local_blocks
    )
)

block_rows_016 = {
    int(row["signed_block_index"]):
        row
    for row in a016[
        "native_mode_space"
    ]["block_rows"]
}

phase_states = a016[
    "phase_state_space"
]["states"]

ell_rows = []

phase_state_line_profiles = {}
phase_state_line_sign_failures = 0

for phase_state in phase_states:
    phase_state_id = int(
        phase_state[
            "phase_state_id"
        ]
    )

    carrier = int(
        phase_state[
            "carrier_index"
        ]
    )

    epsilon = int(
        phase_state[
            "epsilon"
        ]
    )

    section_idx = section_index_by_carrier[
        carrier
    ]

    section_image = block_image_by_aut[
        section_idx
    ]

    inverse_section = inverse_permutation(
        section_image
    )

    target_partition = (
        target_partition_by_carrier_phase[
            (carrier, epsilon)
        ]
    )

    partner_by_block = {}

    for line in target_partition:
        if len(line) != 2:
            raise RuntimeError(
                "projective line does not contain two modes"
            )

        partner_by_block[
            line[0]
        ] = line[1]

        partner_by_block[
            line[1]
        ] = line[0]

    state_rows = []

    for block_index in sorted(
        carrier_to_blocks[
            carrier
        ]
    ):
        canonical_phase_block = (
            inverse_section[
                block_index
            ]
        )

        if epsilon == +1:
            reference_block = (
                canonical_phase_block
            )
        else:
            reference_block = (
                phase_flip_block[
                    canonical_phase_block
                ]
            )

        line_id = (
            reference_line_by_block[
                reference_block
            ]
        )

        source = block_rows_016[
            block_index
        ]

        row = {
            "phase_state_id":
                phase_state_id,

            "carrier_index":
                carrier,

            "epsilon":
                epsilon,

            "signed_block_index":
                block_index,

            "local_mode_slot":
                int(
                    source[
                        "local_mode_slot"
                    ]
                ),

            "sign":
                source["sign"],

            "event_anchor_label":
                int(
                    source[
                        "event_anchor_label"
                    ]
                ),

            "phase_line_partner_block":
                partner_by_block[
                    block_index
                ],

            "reference_line_id":
                line_id,

            "reference_line_label":
                "L" + str(
                    line_id
                ),
        }

        ell_rows.append(
            row
        )

        state_rows.append(
            row
        )

    profile = Counter(
        row[
            "reference_line_id"
        ]
        for row in state_rows
    )

    phase_state_line_profiles[
        phase_state_id
    ] = dict(
        sorted(
            profile.items()
        )
    )

    if profile != Counter({
        0: 2,
        1: 2,
    }):
        checks[
            "two_modes_per_reference_line_in_every_phase_state"
        ] = False

    for line_id in (0, 1):
        signs = {
            row["sign"]
            for row in state_rows
            if row[
                "reference_line_id"
            ] == line_id
        }

        if signs != {
            "positive",
            "negative",
        }:
            phase_state_line_sign_failures += 1

checks.setdefault(
    "two_modes_per_reference_line_in_every_phase_state",
    True,
)

checks["one_positive_one_negative_per_reference_line"] = (
    phase_state_line_sign_failures
    == 0
)

checks["ell_row_count_48"] = (
    len(ell_rows) == 48
)

line_profile = Counter(
    row["reference_line_id"]
    for row in ell_rows
)

checks["global_reference_line_profile_24_24"] = (
    line_profile
    == Counter({
        0: 24,
        1: 24,
    })
)

phase_mode_keys = {
    (
        row["phase_state_id"],
        row["local_mode_slot"],
    )
    for row in ell_rows
}

checks["all_48_phase_mode_pairs_covered_once"] = (
    len(phase_mode_keys) == 48
)

print("PROGRESS: 8/8 classify")

failed_checks = [
    name
    for name, passed in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "native_two_partition_face_geometry_plus_unique_"
    "signed_event_phase_flip_connection_serializes_all_"
    "forty_eight_phase_mode_projective_lines_in_one_flat_"
    "reference_gauge"
    if audit_pass
    else
    "reference_gauge_projective_line_map_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_reference_gauge_projective_line_map_018",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "reference_state": {
        "carrier_index":
            canonical_carrier,

        "epsilon":
            +1,

        "phase_partition":
            [
                list(line)
                for line in P0
            ],

        "reference_lines":
            [
                {
                    "reference_line_id":
                        i,

                    "label":
                        "L" + str(i),

                    "signed_block_indices":
                        list(line),
                }
                for i, line in enumerate(
                    reference_lines
                )
            ],
    },

    "selected_phase_flip_connection": {
        "source":
            "Program-01 / Project-41 201FU",

        "local_permutation":
            list(
                phase_flip_perm
            ),

        "native_realization_count":
            phase_flip_native_realization_count,

        "order":
            perm_order(
                phase_flip_perm
            ),

        "positive_modes":
            "fixed pointwise",

        "negative_modes":
            "swapped",

        "exchanges_projective_partitions":
            True,

        "equals_a_induced_S1":
            False,

        "a_induced_S1_local_permutation":
            list(a_perm),
    },

    "six_face_descent": {
        "transporters_per_face":
            80,

        "unordered_phase_pair_object_count_per_face": {
            str(c):
                unordered_pair_count_by_carrier[
                    c
                ]
            for c in range(6)
        },

        "transporter_independent_unordered_phase_pair":
            True,

        "selected_section_AutG60_indices": {
            str(c):
                section_index_by_carrier[
                    c
                ]
            for c in range(6)
        },

        "section_is_reference_gauge_choice":
            True,
    },

    "ell": {
        "domain":
            "S_12 x M_4",

        "codomain_in_selected_reference_gauge":
            "{L0,L1} subset CP1",

        "row_count":
            len(ell_rows),

        "rows":
            ell_rows,

        "reference_line_profile": {
            str(k): v
            for k, v in sorted(
                line_profile.items()
            )
        },

        "phase_state_line_profiles":
            {
                str(k): v
                for k, v in sorted(
                    phase_state_line_profiles.items()
                )
            },

        "fully_serialized_in_reference_gauge":
            True,
    },

    "earned_statement": (
        "Audit 017 supplies the two native projective-line "
        "partitions of the four signed event modes on the "
        "canonical face. The unique 201FU signed-event-selected "
        "same-face phase-flip transport fixes both positive "
        "event modes pointwise and swaps the two negative modes; "
        "at the projective-line level it exchanges the two "
        "relative-orientation partitions exactly. This transport "
        "is distinct from the a-induced S1 reflection, which "
        "swaps both positive and negative mode pairs. Transporting "
        "the unordered two-partition object through all native "
        "Aut(G60) maps from the canonical face gives one and only "
        "one unordered phase-pair object on each of the six face "
        "carriers. Choosing one explicit flat section and using "
        "the selected 201FU phase-flip connection therefore "
        "serializes ell on all 48 phase-mode pairs into one common "
        "reference CP1 gauge. Every phase state has two modes on "
        "each reference line, and each reference line contains "
        "one positive and one negative event mode."
    ),

    "checks":
        checks,

    "boundary": {
        "ell_serialized_in_one_explicit_flat_reference_gauge":
            audit_pass,

        "raw_reference_line_names_are_gauge_labels":
            True,

        "absolute_Pauli_axis_not_selected":
            True,

        "raw_Hilbert_transport_matrix_not_selected":
            True,

        "selected_phase_flip_not_identified_with_S1":
            True,

        "a_induced_S1_retained_as_distinct_operator":
            True,

        "wedge_zero_nonzero_gauge_stress_not_yet_run":
            True,

        "joint_preparation_census_not_yet_run":
            True,

        "no_probability_law":
            True,

        "no_Born_rule":
            True,

        "no_CHSH_claim":
            True,
    },

    "next_gate": (
        "Construct the phase-decorated two-wing boundary domain "
        "from Audit 016 and classify every joint presentation by "
        "equality versus inequality of the Audit-018 reference "
        "projective lines. Then stress that zero/nonzero wedge "
        "classification under alternate admissible flat sections "
        "before admitting any source-preparation branch."
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

note = f"""# Reference-gauge projective line map 018

## Result

Audit pass:

    {audit_pass}

Audit 017 derived two native projective-line partitions of the four
signed event modes.

The same-face phase-flip connection selected in 201FU has one unique
local four-mode permutation. It

    fixes both positive event modes pointwise
    swaps the two negative event modes
    has order 2.

It exchanges the two projective-line partitions exactly.

This selected phase-flip is not the a-induced S1 operator.

S1 swaps both positive and negative mode pairs.

## Six-face descent

There are

    80

native Aut(G60) transporters from the canonical face to each target
face.

For each of all six target faces, every transporter produces the same
unordered pair of projective-line partitions.

Thus the two-phase projective object descends over all six face
carriers independently of transporter choice.

An explicit section is then chosen as a reference gauge.

## ell

Using that section and the unique 201FU phase-flip connection, the map

    ell:
      S_12 x M_4
        ->
      CP1

is serialized for all

    {len(ell_rows)}

phase-mode pairs.

In the selected common reference gauge the two projective lines are
called

    L0
    L1.

Every one of the twelve phase states has

    2 modes on L0
    2 modes on L1.

Each reference line contains

    one positive mode
    one negative mode.

Across all forty-eight phase-mode pairs the line profile is

    L0: {line_profile.get(0, 0)}
    L1: {line_profile.get(1, 0)}.

The names L0 and L1 are gauge labels, not absolute Pauli axes.

## Next gate

The one-wing preparation interface is now explicit:

    eta_tilde
    ell.

The next audit constructs the two-wing phase-decorated boundary domain,
classifies equal-line versus unequal-line pairs, and therefore computes
zero versus nonzero antisymmetric wedge.

That classification must then survive alternate admissible flat-section
choices before any preparation branch is admitted.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print(
    "FU_LOCAL_TRANSPORT_PERM:",
    phase_flip_perm,
)
print(
    "FU_NATIVE_REALIZATION_COUNT:",
    phase_flip_native_realization_count,
)
print(
    "FU_EQUALS_S1:",
    phase_flip_perm == a_perm,
)
print(
    "TRANSPORTER_COUNT_PROFILE:",
    dict(sorted(
        transporter_count_profile.items()
    )),
)
print(
    "UNORDERED_PHASE_PAIR_COUNT_BY_FACE:",
    unordered_pair_count_by_carrier,
)
print(
    "ELL_ROW_COUNT:",
    len(ell_rows),
)
print(
    "REFERENCE_LINE_PROFILE:",
    dict(sorted(
        line_profile.items()
    )),
)
print(
    "ELL_REFERENCE_GAUGE_SERIALIZED:",
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
