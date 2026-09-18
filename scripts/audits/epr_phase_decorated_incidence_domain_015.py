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

BD = (
    P41
    / "artifacts/provenance"
    / "project41_registered_seam_world_address_crosswalk_201bd.v1.json"
)

FR = (
    P41
    / "artifacts/provenance"
    / "project41_global_local_phase_character_bridge_201fr.v1.json"
)

FS = (
    P41
    / "artifacts/provenance"
    / "project41_signed_event_doublet_vector_lift_201fs.v1.json"
)

SIGNED = (
    P41
    / "artifacts/json"
    / "signed_24_face_block_closure_audit_025.json"
)

FP7 = (
    P41
    / "artifacts/provenance"
    / "project41_u2_gauge_bridge_connection_frontier_checkpoint_201fp7.v1.json"
)

A014 = (
    HERE
    / "artifacts/json"
    / "epr_flat_groupoid_preparation_interface_014.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_phase_decorated_incidence_domain_015.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_phase_decorated_incidence_domain_015.md"
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


print("PROGRESS: 1/7 load sealed interfaces")

bd = load(BD)
fr = load(FR)
fs = load(FS)
signed = load(SIGNED)
fp7 = load(FP7)
a014 = load(A014)

checks = {}

checks["201BD_passed"] = bd["audit_pass"] is True
checks["201FR_passed"] = fr["audit_pass"] is True
checks["201FS_passed"] = fs["audit_pass"] is True
checks["signed_block_audit_passed"] = signed["audit_pass"] is True
checks["FP7_passed"] = fp7["audit_pass"] is True
checks["Audit014_passed"] = a014["audit_pass"] is True

print("PROGRESS: 2/7 reconstruct signed block lookup")

blocks = signed["measurements"]["signed_blocks"]

block_key_to_index = {}

for i, block in enumerate(blocks):
    key = tuple(
        sorted(
            int(x)
            for x in block["five_state_block"]
        )
    )

    if key in block_key_to_index:
        raise RuntimeError(
            "duplicate signed block state set"
        )

    block_key_to_index[key] = i

checks["signed_block_count_24"] = (
    len(blocks) == 24
)

carrier_to_blocks = defaultdict(list)

for i, block in enumerate(blocks):
    carrier_to_blocks[
        int(block["carrier_index"])
    ].append(i)

checks["six_face_carriers"] = (
    set(carrier_to_blocks.keys())
    == set(range(6))
)

checks["four_signed_modes_per_carrier"] = all(
    len(carrier_to_blocks[c]) == 4
    for c in range(6)
)

print("PROGRESS: 3/7 crosswalk all registered event sides")

rows = bd["anchor_rows"]

side_records = []

match_failures = 0

for row_index, row in enumerate(rows):
    for side in ("inward", "outward"):
        key = tuple(
            sorted(
                int(x)
                for x in row[
                    side + "_states"
                ]
            )
        )

        block_index = block_key_to_index.get(
            key
        )

        if block_index is None:
            match_failures += 1
            continue

        block = blocks[block_index]

        side_records.append({
            "row_index":
                row_index,

            "registration_index":
                int(row["registration_index"]),

            "anchor_index":
                int(row["anchor_index"]),

            "anchor_positive_index":
                int(row["anchor_positive_index"]),

            "anchor_negative_index":
                int(row["anchor_negative_index"]),

            "side":
                side,

            "block_index":
                block_index,

            "carrier_index":
                int(block["carrier_index"]),

            "sign":
                block["sign"],
        })

checks["registered_event_row_count_24"] = (
    len(rows) == 24
)

checks["event_side_count_48"] = (
    len(side_records) == 48
)

checks["all_event_sides_match_exact_signed_blocks"] = (
    match_failures == 0
)

print("PROGRESS: 4/7 measure block occurrence anatomy")

block_records = defaultdict(list)

for rec in side_records:
    block_records[
        rec["block_index"]
    ].append(rec)

occurrence_profile = Counter(
    len(v)
    for v in block_records.values()
)

checks["every_signed_block_occurs_twice"] = (
    occurrence_profile
    == Counter({2: 24})
)

registration_profile_failures = 0
anchor_label_failures = 0
side_sign_failures = 0

block_anchor_label = {}

for block_index in range(24):
    recs = block_records[block_index]

    registrations = {
        rec["registration_index"]
        for rec in recs
    }

    if registrations != {0, 1}:
        registration_profile_failures += 1

    sign = blocks[block_index]["sign"]

    if sign == "positive":
        labels = {
            rec["anchor_positive_index"]
            for rec in recs
        }

        if any(
            rec["side"] != "outward"
            for rec in recs
        ):
            side_sign_failures += 1

    elif sign == "negative":
        labels = {
            rec["anchor_negative_index"]
            for rec in recs
        }

        if any(
            rec["side"] != "inward"
            for rec in recs
        ):
            side_sign_failures += 1

    else:
        raise RuntimeError(
            "unexpected block sign"
        )

    if len(labels) != 1:
        anchor_label_failures += 1
    else:
        block_anchor_label[
            block_index
        ] = next(iter(labels))

checks["every_block_occurs_once_per_registration"] = (
    registration_profile_failures == 0
)

checks["outward_is_exact_positive_block_side"] = (
    side_sign_failures == 0
)

checks["each_signed_block_has_unique_event_anchor_label"] = (
    anchor_label_failures == 0
    and len(block_anchor_label) == 24
)

print("PROGRESS: 5/7 recover four native mode objects per face")

carrier_positive_anchors = defaultdict(set)
carrier_negative_anchors = defaultdict(set)

for block_index, anchor in (
    block_anchor_label.items()
):
    block = blocks[block_index]
    carrier = int(
        block["carrier_index"]
    )

    if block["sign"] == "positive":
        carrier_positive_anchors[
            carrier
        ].add(anchor)
    else:
        carrier_negative_anchors[
            carrier
        ].add(anchor)

checks["two_positive_anchor_modes_per_face"] = all(
    len(
        carrier_positive_anchors[c]
    )
    == 2
    for c in range(6)
)

checks["two_negative_anchor_modes_per_face"] = all(
    len(
        carrier_negative_anchors[c]
    )
    == 2
    for c in range(6)
)

checks["positive_and_negative_anchor_label_sets_match_per_face"] = all(
    carrier_positive_anchors[c]
    ==
    carrier_negative_anchors[c]
    for c in range(6)
)

# Registration 0 and 1 describe two ways to pair the
# positive and negative mode objects on the same face.

pairing_by_registration = {
    0: defaultdict(dict),
    1: defaultdict(dict),
}

carrier_pairing_failure_count = 0

for rec in side_records:
    if rec["side"] != "outward":
        continue

    row = rows[
        rec["row_index"]
    ]

    reg = int(
        row["registration_index"]
    )

    pos_anchor = int(
        row["anchor_positive_index"]
    )

    neg_anchor = int(
        row["anchor_negative_index"]
    )

    carrier = rec["carrier_index"]

    # Verify inward side is on the same carrier.
    inward_key = tuple(
        sorted(
            int(x)
            for x in row["inward_states"]
        )
    )

    inward_block = block_key_to_index[
        inward_key
    ]

    inward_carrier = int(
        blocks[inward_block][
            "carrier_index"
        ]
    )

    if inward_carrier != carrier:
        carrier_pairing_failure_count += 1

    pairing_by_registration[
        reg
    ][carrier][
        pos_anchor
    ] = neg_anchor

checks["registered_pairings_stay_inside_face_carrier"] = (
    carrier_pairing_failure_count == 0
)

reg0_identity = True
reg1_swap = True

for carrier in range(6):
    anchors = (
        carrier_positive_anchors[
            carrier
        ]
    )

    map0 = pairing_by_registration[
        0
    ][carrier]

    map1 = pairing_by_registration[
        1
    ][carrier]

    if set(map0.keys()) != anchors:
        reg0_identity = False

    if set(map1.keys()) != anchors:
        reg1_swap = False

    for a in anchors:
        if map0.get(a) != a:
            reg0_identity = False

        if map1.get(a) == a:
            reg1_swap = False

        if map1.get(
            map1.get(a)
        ) != a:
            reg1_swap = False

checks["registration0_pairs_same_anchor_modes"] = (
    reg0_identity
)

checks["registration1_swaps_two_anchor_modes_within_each_face"] = (
    reg1_swap
)

print("PROGRESS: 6/7 test whether bare incidence contains registration")

incidence_registration_sets = {}

for block_index, block in enumerate(blocks):
    regs = {
        rec["registration_index"]
        for rec in block_records[
            block_index
        ]
    }

    for state in block[
        "five_state_block"
    ]:
        incidence_registration_sets[
            (
                int(state),
                block_index,
            )
        ] = set(regs)

checks["one_wing_incidence_count_120"] = (
    len(
        incidence_registration_sets
    )
    == 120
)

registration_support_profile = Counter(
    len(regs)
    for regs in incidence_registration_sets.values()
)

checks["every_bare_incidence_supports_both_registrations"] = (
    registration_support_profile
    == Counter({2: 120})
    and all(
        regs == {0, 1}
        for regs in incidence_registration_sets.values()
    )
)

print("PROGRESS: 7/7 classify corrected eta domain")

local_state_count = int(
    fp7[
        "results"
    ]["201FLB"][
        "local_order_parameter_state_count"
    ]
)

face_count = int(
    fr[
        "native_data"
    ]["face_count"]
)

checks["local_order_parameter_state_count_12"] = (
    local_state_count == 12
)

checks["local_state_count_is_two_per_face"] = (
    local_state_count
    == 2 * face_count
)

# Each face carrier has 20 bare incidences:
#
# 4 signed mode objects x 5 G60 states.
#
# A phase-decorated domain with two local order states
# per face therefore has 40 rows per face and 240 total.

bare_incidences_per_face = {
    c: sum(
        len(
            blocks[i][
                "five_state_block"
            ]
        )
        for i in carrier_to_blocks[c]
    )
    for c in range(6)
}

checks["twenty_bare_incidences_per_face"] = all(
    bare_incidences_per_face[c] == 20
    for c in range(6)
)

phase_decorated_incidence_count = sum(
    2 * bare_incidences_per_face[c]
    for c in range(6)
)

checks["phase_decorated_domain_expected_count_240"] = (
    phase_decorated_incidence_count
    == 240
)

failed_checks = [
    name
    for name, passed in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "registered_event_crosswalk_recovers_exact_four_mode_"
    "objects_per_face_but_bare_I_F_does_not_carry_the_"
    "independent_local_phase_state_so_eta_requires_a_"
    "phase_decorated_240_row_domain"
    if audit_pass
    else
    "phase_decorated_incidence_domain_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_phase_decorated_incidence_domain_015",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "registered_event_crosswalk": {
        "event_row_count":
            len(rows),

        "event_side_count":
            len(side_records),

        "signed_block_count":
            len(blocks),

        "occurrences_per_signed_block":
            2,

        "registration_support_per_signed_block":
            [0, 1],

        "registration0_mode_pairing":
            "same-anchor",

        "registration1_mode_pairing":
            "within-face anchor swap",
    },

    "native_mode_objects": {
        "face_count":
            6,

        "signed_modes_per_face":
            4,

        "positive_modes_per_face":
            2,

        "negative_modes_per_face":
            2,

        "mode_identifier":
            (
                "exact signed face block, equivalently "
                "sign plus its unique event-anchor label "
                "inside the carrier"
            ),
    },

    "bare_incidence_domain": {
        "definition":
            "I_F = {(u,B): u belongs to signed block B}",

        "count":
            120,

        "registrations_supported_per_incidence":
            2,

        "contains_independent_local_phase_bit":
            False,
    },

    "local_phase_interface": {
        "local_order_parameter_state_count":
            local_state_count,

        "face_count":
            face_count,

        "states_per_face":
            2,

        "source":
            (
                "Program-01 phase character delta_f = "
                "chi_deck restricted to Stab(f)"
            ),
    },

    "corrected_eta_domain": {
        "name":
            "I_F_tilde",

        "definition":
            (
                "{(s,u,B): s in S_12, "
                "face(s)=carrier(B), u in B}"
            ),

        "expected_count":
            phase_decorated_incidence_count,

        "map":
            "eta_tilde(s,u,B) = (s,m_B)",
    },

    "earned_statement": (
        "The registered-event crosswalk closes the native "
        "mode-object half of the preparation interface. Every "
        "one of the 48 inward/outward event sides matches one "
        "of the 24 signed face blocks exactly; every signed "
        "block occurs twice, once in each registered seam "
        "sector, and carries one unique event-anchor label. "
        "Each face therefore has exactly four native signed "
        "mode objects, two positive and two negative. The two "
        "registered seam sectors pair those objects differently: "
        "registration 0 pairs equal anchor labels while "
        "registration 1 swaps the two anchor labels inside each "
        "face. Crucially, every one of the 120 bare boundary "
        "incidences is compatible with both registrations. "
        "Program-01 independently supplies a binary local "
        "relative-complex-orientation state on each of six "
        "faces. Therefore the bare incidence domain I_F does "
        "not by itself carry the independent phase state needed "
        "for the EPR projective-line map. The correct native "
        "preparation domain is the phase-decorated incidence "
        "domain I_F_tilde with expected size 240."
    ),

    "checks":
        checks,

    "boundary": {
        "Audit014_eta_domain_corrected":
            True,

        "bare_incidence_not_used_to_invent_phase":
            True,

        "registration_index_not_identified_with_delta_f":
            True,

        "signed_block_sign_not_identified_with_delta_f":
            True,

        "mode_object_crosswalk_closed":
            True,

        "explicit_S12_phase_state_rows_not_yet_serialized":
            True,

        "projective_line_table_ell_not_yet_derived":
            True,

        "wedge_census_not_yet_run":
            True,

        "no_probability_law":
            True,

        "no_Born_rule":
            True,

        "no_CHSH_claim":
            True,
    },

    "next_gate": (
        "Recover the explicit twelve local phase-state rows "
        "from the current Program-01 face-field checkpoint "
        "201FK, not from historical 201FN matrices. Construct "
        "the 240-row phase-decorated incidence domain and then "
        "derive the projective line map ell for its four native "
        "mode objects per face."
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

note = f"""# Phase-decorated incidence domain 015

## Result

Audit pass:

    {audit_pass}

The 24 registered seam rows contain 48 event sides.

Every event side matches one of the 24 signed face blocks exactly.

Every signed block occurs exactly twice:

    once in registration 0
    once in registration 1.

Each signed block also has one unique event-anchor label.

Thus each face carries exactly four native mode objects:

    two positive
    two negative.

Registration 0 pairs positive and negative modes with the same anchor
label.

Registration 1 swaps the two anchor labels inside each face.

## Domain correction

The bare boundary-incidence domain has

    120

rows.

Every one of those 120 incidences is compatible with both registered
seam sectors.

Program 01 independently supplies two local relative-complex-orientation
states per face:

    6 faces x 2 phase states = 12 local states.

Therefore a bare incidence cannot be used to manufacture the phase bit.

The correct preparation domain is

    I_F_tilde
      =
    {{(s,u,B):
      s in S_12,
      face(s)=carrier(B),
      u in B}}.

Each face has

    4 modes x 5 states = 20

bare incidences and two phase states, giving

    6 x 20 x 2 = {phase_decorated_incidence_count}

phase-decorated incidences.

## Next target

Recover the explicit twelve local phase-state rows from the current
201FK face-field checkpoint.

Then construct I_F_tilde literally and derive the CP1 projective line
map ell.

No registration bit or signed-block sign is being silently identified
with the independent complex-orientation phase.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print(
    "EVENT_SIDE_COUNT:",
    len(side_records),
)
print(
    "SIGNED_BLOCK_OCCURRENCE_PROFILE:",
    dict(sorted(
        occurrence_profile.items()
    )),
)
print(
    "BARE_INCIDENCE_COUNT:",
    len(
        incidence_registration_sets
    ),
)
print(
    "BARE_INCIDENCE_REGISTRATION_SUPPORT_PROFILE:",
    dict(sorted(
        registration_support_profile.items()
    )),
)
print(
    "LOCAL_ORDER_PARAMETER_STATE_COUNT:",
    local_state_count,
)
print(
    "PHASE_DECORATED_INCIDENCE_COUNT:",
    phase_decorated_incidence_count,
)
print(
    "CORRECTED_ETA_DOMAIN:",
    "I_F_tilde",
)
print("FAILED_CHECK_COUNT:", len(failed_checks))
print("FAILED_CHECKS:", failed_checks)
print("JSON_OUT:", JSON_OUT)
print("NOTE_OUT:", NOTE_OUT)
print(
    "JSON_SHA256:",
    sha256(JSON_OUT.read_bytes()).hexdigest(),
)
