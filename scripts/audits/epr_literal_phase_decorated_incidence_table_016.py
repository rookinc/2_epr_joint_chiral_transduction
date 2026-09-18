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

FK = (
    P41
    / "artifacts/provenance"
    / "project41_face_field_higgs_bridge_checkpoint_receipt_201fk.v1.json"
)

FR = (
    P41
    / "artifacts/provenance"
    / "project41_global_local_phase_character_bridge_201fr.v1.json"
)

BD = (
    P41
    / "artifacts/provenance"
    / "project41_registered_seam_world_address_crosswalk_201bd.v1.json"
)

SIGNED = (
    P41
    / "artifacts/json"
    / "signed_24_face_block_closure_audit_025.json"
)

A015 = (
    HERE
    / "artifacts/json"
    / "epr_phase_decorated_incidence_domain_015.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_literal_phase_decorated_incidence_table_016.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_literal_phase_decorated_incidence_table_016.md"
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

fk = load(FK)
fr = load(FR)
bd = load(BD)
signed = load(SIGNED)
a015 = load(A015)

checks = {}

checks["201FK_passed"] = fk["audit_pass"] is True
checks["201FR_passed"] = fr["audit_pass"] is True
checks["201BD_passed"] = bd["audit_pass"] is True
checks["signed_block_audit_passed"] = signed["audit_pass"] is True
checks["Audit015_passed"] = a015["audit_pass"] is True

print("PROGRESS: 2/7 certify local two-state phase object")

local_face = fk["local_face_fiber"]
order_parameter = fk["local_order_parameter"]
doublet = fk["symmetry_reduced_complex_doublet"]

checks["native_face_object_count_6"] = (
    int(local_face["native_face_object_count"]) == 6
)

checks["signed_blocks_per_face_4"] = (
    int(local_face["signed_blocks_per_face_object"]) == 4
)

checks["local_phase_state_count_2"] = (
    int(order_parameter["state_count"]) == 2
)

checks["order_parameter_is_relative_complex_orientation"] = (
    order_parameter["order_parameter_type"]
    == "relative complex orientation"
)

checks["raw_face_block_does_not_select_phase"] = (
    order_parameter["raw_face_block_is_selector"] is False
)

checks["raw_face_vector_does_not_select_phase"] = (
    order_parameter["raw_face_vector_is_selector"] is False
)

checks["state_space_is_two_relative_JE_orientations"] = (
    order_parameter["state_space"]
    == "{[J_S direct_sum J_E], [J_S direct_sum (-J_E)]}"
)

checks["full_complex_structure_recorded"] = (
    doublet["full_complex_structure"]
    == "J_F = J_S direct_sum J_E"
)

checks["phase_character_matches_201FR"] = (
    order_parameter["swap_character"]
    == fr["native_data"]["local_order_parameter_character"]
)

print("PROGRESS: 3/7 reconstruct signed mode objects")

blocks = signed["measurements"]["signed_blocks"]

block_key_to_index = {}

for i, block in enumerate(blocks):
    key = tuple(
        sorted(
            int(x)
            for x in block["five_state_block"]
        )
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

checks["four_blocks_on_each_carrier"] = all(
    len(carrier_to_blocks[c]) == 4
    for c in range(6)
)

print("PROGRESS: 4/7 recover native mode labels from registered events")

block_anchor_labels = defaultdict(set)

for row in bd["anchor_rows"]:
    for side in ("inward", "outward"):
        states = tuple(
            sorted(
                int(x)
                for x in row[side + "_states"]
            )
        )

        block_index = block_key_to_index[states]
        block = blocks[block_index]

        if block["sign"] == "positive":
            label = int(
                row["anchor_positive_index"]
            )
        else:
            label = int(
                row["anchor_negative_index"]
            )

        block_anchor_labels[
            block_index
        ].add(label)

checks["every_block_has_one_native_anchor_label"] = all(
    len(block_anchor_labels[i]) == 1
    for i in range(24)
)

block_anchor = {
    i: next(iter(block_anchor_labels[i]))
    for i in range(24)
}

carrier_anchor_sets = {}

for carrier in range(6):
    anchors = {
        block_anchor[i]
        for i in carrier_to_blocks[carrier]
    }

    carrier_anchor_sets[carrier] = anchors

checks["two_anchor_labels_per_face"] = all(
    len(carrier_anchor_sets[c]) == 2
    for c in range(6)
)

# Native mode identity is the signed block itself.
#
# For readable serialization, also attach:
#
#   sign
#   event-anchor label
#
# The integer local_mode_slot below is only a deterministic
# serialization convention. It is not claimed as additional
# native structure.

block_to_slot = {}

for carrier in range(6):
    anchors = sorted(
        carrier_anchor_sets[carrier]
    )

    mode_keys = []

    for sign in ("negative", "positive"):
        for anchor in anchors:
            candidates = [
                i
                for i in carrier_to_blocks[carrier]
                if blocks[i]["sign"] == sign
                and block_anchor[i] == anchor
            ]

            if len(candidates) != 1:
                raise RuntimeError(
                    "mode-key lookup is not unique"
                )

            mode_keys.append(candidates[0])

    for slot, block_index in enumerate(mode_keys):
        block_to_slot[block_index] = slot

checks["all_24_blocks_receive_local_slots"] = (
    len(block_to_slot) == 24
)

print("PROGRESS: 5/7 construct literal twelve phase states")

phase_states = []

for carrier in range(6):
    for epsilon in (+1, -1):
        phase_state_id = (
            2 * carrier
            + (0 if epsilon == +1 else 1)
        )

        phase_states.append({
            "phase_state_id":
                phase_state_id,

            "carrier_index":
                carrier,

            "epsilon":
                epsilon,

            "complex_structure":
                (
                    "J_S direct sum J_E"
                    if epsilon == +1
                    else
                    "J_S direct sum (-J_E)"
                ),

            "order_parameter":
                "relative complex orientation",
        })

checks["literal_phase_state_count_12"] = (
    len(phase_states) == 12
)

checks["two_phase_states_per_carrier"] = all(
    sum(
        1
        for row in phase_states
        if row["carrier_index"] == c
    ) == 2
    for c in range(6)
)

checks["phase_state_ids_are_0_through_11"] = (
    sorted(
        row["phase_state_id"]
        for row in phase_states
    )
    == list(range(12))
)

print("PROGRESS: 6/7 construct literal 240-row decorated incidence table")

decorated_rows = []

for phase in phase_states:
    carrier = phase["carrier_index"]

    for block_index in sorted(
        carrier_to_blocks[carrier]
    ):
        block = blocks[block_index]

        for state in sorted(
            int(x)
            for x in block["five_state_block"]
        ):
            decorated_rows.append({
                "phase_state_id":
                    phase["phase_state_id"],

                "carrier_index":
                    carrier,

                "epsilon":
                    phase["epsilon"],

                "complex_structure":
                    phase["complex_structure"],

                "g60_state":
                    state,

                "signed_block_index":
                    block_index,

                "block_sign":
                    block["sign"],

                "event_anchor_label":
                    block_anchor[block_index],

                "local_mode_slot":
                    block_to_slot[block_index],

                "local_mode_identity":
                    {
                        "signed_block_index":
                            block_index,

                        "sign":
                            block["sign"],

                        "event_anchor_label":
                            block_anchor[block_index],
                    },
            })

checks["decorated_row_count_240"] = (
    len(decorated_rows) == 240
)

phase_row_profile = Counter(
    row["phase_state_id"]
    for row in decorated_rows
)

checks["twenty_rows_per_phase_state"] = all(
    phase_row_profile[s] == 20
    for s in range(12)
)

carrier_row_profile = Counter(
    row["carrier_index"]
    for row in decorated_rows
)

checks["forty_rows_per_face_carrier"] = all(
    carrier_row_profile[c] == 40
    for c in range(6)
)

phase_mode_profile = Counter(
    (
        row["phase_state_id"],
        row["local_mode_slot"],
    )
    for row in decorated_rows
)

checks["five_G60_states_per_phase_mode"] = all(
    count == 5
    for count in phase_mode_profile.values()
)

checks["phase_mode_pair_count_48"] = (
    len(phase_mode_profile) == 48
)

print("PROGRESS: 7/7 verify forgetting phase gives Audit015 domain")

bare_projection = Counter(
    (
        row["g60_state"],
        row["signed_block_index"],
    )
    for row in decorated_rows
)

checks["bare_projection_has_120_rows"] = (
    len(bare_projection) == 120
)

checks["forgetting_phase_is_exactly_two_to_one"] = all(
    count == 2
    for count in bare_projection.values()
)

checks["matches_Audit015_expected_count"] = (
    len(decorated_rows)
    ==
    int(
        a015[
            "corrected_eta_domain"
        ]["expected_count"]
    )
)

failed_checks = [
    name
    for name, passed in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "literal_twelve_state_relative_complex_orientation_table_"
    "and_exact_240_row_phase_decorated_native_boundary_"
    "incidence_domain_constructed"
    if audit_pass
    else
    "literal_phase_decorated_incidence_table_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_literal_phase_decorated_incidence_table_016",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "phase_state_space": {
        "face_count":
            6,

        "states_per_face":
            2,

        "state_count":
            len(phase_states),

        "state_formula":
            "(face, epsilon), epsilon in {+1,-1}",

        "complex_structure_formula":
            "J_F(face,epsilon) = J_S direct sum epsilon J_E",

        "swap_character":
            order_parameter["swap_character"],

        "global_label_swap_is_presentation_convention":
            True,

        "states":
            phase_states,
    },

    "native_mode_space": {
        "signed_block_count":
            24,

        "modes_per_face":
            4,

        "mode_identity":
            (
                "signed face block plus its unique "
                "registered-event anchor label"
            ),

        "local_mode_slot_is_serialization_only":
            True,

        "block_rows": [
            {
                "signed_block_index":
                    i,

                "carrier_index":
                    int(blocks[i]["carrier_index"]),

                "sign":
                    blocks[i]["sign"],

                "event_anchor_label":
                    block_anchor[i],

                "local_mode_slot":
                    block_to_slot[i],

                "five_state_block":
                    sorted(
                        int(x)
                        for x in blocks[i]["five_state_block"]
                    ),
            }
            for i in range(24)
        ],
    },

    "decorated_incidence_domain": {
        "name":
            "I_F_tilde",

        "definition":
            (
                "{(s,u,B): s in S_12, "
                "face(s)=carrier(B), u in B}"
            ),

        "row_count":
            len(decorated_rows),

        "rows_per_phase_state":
            20,

        "rows_per_face":
            40,

        "phase_mode_pair_count":
            48,

        "states_per_phase_mode":
            5,

        "forget_phase_projection":
            "2-to-1 onto I_F",

        "rows":
            decorated_rows,
    },

    "eta_tilde": {
        "formula":
            "eta_tilde(s,u,B) = (s,m_B)",

        "domain_row_count":
            len(decorated_rows),

        "codomain_phase_mode_pair_count":
            len(phase_mode_profile),

        "fully_serialized":
            True,
    },

    "earned_statement": (
        "Program-01 201FK supplies exactly two relative-complex-"
        "orientation states on each of the six native face "
        "objects, represented as J_S direct sum epsilon J_E "
        "with epsilon = plus or minus one. Combined with the "
        "four exact signed event-mode objects per face closed "
        "in Audit 015, this constructs the literal twelve-state "
        "local phase table and the exact 240-row phase-decorated "
        "native boundary-incidence domain I_F_tilde. Every "
        "phase state carries twenty incidence rows, every face "
        "carries forty, every one of the forty-eight phase-mode "
        "pairs carries five native G60 states, and forgetting "
        "the phase is exactly two-to-one onto the 120-row bare "
        "incidence domain."
    ),

    "checks":
        checks,

    "boundary": {
        "phase_and_event_coordinates_both_explicit":
            True,

        "raw_face_block_does_not_select_phase":
            True,

        "raw_face_vector_does_not_select_phase":
            True,

        "registration_index_not_identified_with_phase":
            True,

        "signed_block_sign_not_identified_with_phase":
            True,

        "epsilon_label_global_swap_is_conventional":
            True,

        "projective_line_map_ell_not_yet_derived":
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
        "Derive ell on the forty-eight explicit phase-mode "
        "pairs. Determine the CP1 projective line occupied by "
        "each signed event mode under each relative complex "
        "orientation, using current Program-01 face Hilbert "
        "interfaces only. Then classify the phase-decorated "
        "two-wing boundary domain by projective-line equality "
        "before analyzer settings are introduced."
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

note = f"""# Literal phase-decorated incidence table 016

## Result

Audit pass:

    {audit_pass}

Program-01 201FK gives two local relative-complex-orientation states on
each of six native face objects.

Serialize them as

    (face, epsilon)

with

    epsilon = +1 or -1

and

    J_F(face,epsilon)
      =
    J_S direct sum epsilon J_E.

This gives exactly

    6 x 2 = {len(phase_states)}

local phase states.

The global exchange of the names +1 and -1 is only a presentation
convention.

## Native event modes

Audit 015 already closed four native signed event-mode objects on each
face:

    two positive
    two negative.

Each mode is identified by its exact signed face block and unique
registered-event anchor label.

The integer local mode slot written into the JSON is only a deterministic
serialization label.

## Decorated domain

The literal native preparation domain is

    I_F_tilde
      =
    {{(s,u,B):
      s in S_12,
      face(s)=carrier(B),
      u in B}}.

It has exactly

    {len(decorated_rows)}

rows.

Profiles:

    20 rows per phase state
    40 rows per face
    48 phase-mode pairs
    5 G60 states per phase-mode pair.

Forgetting the phase maps I_F_tilde exactly two-to-one onto the 120-row
bare boundary-incidence domain.

## Closed interface

The corrected eta map is now fully serialized:

    eta_tilde(s,u,B) = (s,m_B).

The remaining local preparation problem is ell:

    ell:
      S_12 x M_4
        ->
      CP1.

Once ell is native, the two-wing preparation census can classify
projective-line equality and therefore zero versus nonzero
antisymmetric wedge before analyzer settings exist.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print(
    "PHASE_STATE_COUNT:",
    len(phase_states),
)
print(
    "PHASE_MODE_PAIR_COUNT:",
    len(phase_mode_profile),
)
print(
    "DECORATED_INCIDENCE_ROW_COUNT:",
    len(decorated_rows),
)
print(
    "ROWS_PER_PHASE_STATE_PROFILE:",
    dict(sorted(
        Counter(
            phase_row_profile.values()
        ).items()
    )),
)
print(
    "FORGET_PHASE_MULTIPLICITY_PROFILE:",
    dict(sorted(
        Counter(
            bare_projection.values()
        ).items()
    )),
)
print(
    "ETA_TILDE_FULLY_SERIALIZED:",
    True,
)
print(
    "NEXT_MISSING_INTERFACE:",
    "ell",
)
print("FAILED_CHECK_COUNT:", len(failed_checks))
print("FAILED_CHECKS:", failed_checks)
print("JSON_OUT:", JSON_OUT)
print("NOTE_OUT:", NOTE_OUT)
print(
    "JSON_SHA256:",
    sha256(JSON_OUT.read_bytes()).hexdigest(),
)
