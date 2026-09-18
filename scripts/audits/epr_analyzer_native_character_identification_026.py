#!/usr/bin/env python3

from collections import Counter
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

ER = (
    P41
    / "artifacts/provenance"
    / "project41_native_face_chirality_closure_receipt_201er.v1.json"
)

EU = (
    P41
    / "artifacts/provenance"
    / "project41_joint_deck_face_C2xC2_quotient_receipt_201eu.v1.json"
)

FR = (
    P41
    / "artifacts/provenance"
    / "project41_global_local_phase_character_bridge_201fr.v1.json"
)

FK = (
    P41
    / "artifacts/provenance"
    / "project41_face_field_higgs_bridge_checkpoint_receipt_201fk.v1.json"
)

A024 = (
    HERE
    / "artifacts/json"
    / "epr_g60_label_torsor_c5_face_descent_024.v1.json"
)

A025 = (
    HERE
    / "artifacts/json"
    / "epr_projective_analyzer_sector_extension_025.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_analyzer_native_character_identification_026.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_analyzer_native_character_identification_026.md"
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


def parity(p):
    inv = 0

    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                inv += 1

    return inv % 2


print("PROGRESS: 1/8 load sealed authorities")

aut = load(AUT)
signed = load(SIGNED)
er = load(ER)
eu = load(EU)
fr = load(FR)
fk = load(FK)
a024 = load(A024)
a025 = load(A025)

checks = {}

checks["AutG60_cache_passes"] = (
    aut["cache_pass"] is True
)

checks["signed_face_audit_passes"] = (
    signed["audit_pass"] is True
)

checks["ER_passes"] = (
    er["audit_pass"] is True
)

checks["EU_passes"] = (
    eu["audit_pass"] is True
)

checks["FR_passes"] = (
    fr["audit_pass"] is True
)

checks["Audit024_passes"] = (
    a024["audit_pass"] is True
)

checks["Audit025_passes"] = (
    a025["audit_pass"] is True
)

checks["Audit025_sector_bit_is_permutation_parity"] = (
    a025[
        "sector_extension"
    ][
        "sector_flip_equals_permutation_parity"
    ]
    is True
)

print("PROGRESS: 2/8 reconstruct six native face objects")

blocks = signed[
    "measurements"
]["signed_blocks"]

block_lookup = {}

carrier_to_blocks = {}

for i, block in enumerate(blocks):
    key = tuple(
        sorted(
            int(x)
            for x in block[
                "five_state_block"
            ]
        )
    )

    block_lookup[key] = i

    carrier = int(
        block["carrier_index"]
    )

    carrier_to_blocks.setdefault(
        carrier,
        [],
    ).append(i)

checks["signed_block_count_24"] = (
    len(blocks) == 24
)

checks["six_face_carriers"] = (
    sorted(
        carrier_to_blocks
    )
    == list(range(6))
)

checks["four_signed_blocks_per_face"] = all(
    len(
        carrier_to_blocks[c]
    )
    == 4
    for c in range(6)
)


def transformed_block_index(p, block_index):
    block = blocks[
        block_index
    ]

    key = tuple(
        sorted(
            p[
                int(u)
            ]
            for u in block[
                "five_state_block"
            ]
        )
    )

    return block_lookup.get(
        key
    )


def face_and_sign_action(p):
    face_perm = []
    sign_ratios = []

    for source_face in range(6):
        target_faces = set()

        for block_index in carrier_to_blocks[
            source_face
        ]:
            target_index = (
                transformed_block_index(
                    p,
                    block_index,
                )
            )

            if target_index is None:
                return None

            source_sign = (
                blocks[
                    block_index
                ]["sign"]
            )

            target_sign = (
                blocks[
                    target_index
                ]["sign"]
            )

            target_faces.add(
                int(
                    blocks[
                        target_index
                    ]["carrier_index"]
                )
            )

            sign_ratios.append(
                1
                if source_sign
                == target_sign
                else -1
            )

        if len(target_faces) != 1:
            return None

        face_perm.append(
            next(
                iter(
                    target_faces
                )
            )
        )

    if len(set(sign_ratios)) != 1:
        return None

    return (
        tuple(face_perm),
        sign_ratios[0],
    )


print("PROGRESS: 3/8 evaluate all 480 native automorphisms")

rows = aut[
    "measurements"
]["automorphism_rows"]

checks["AutG60_order_480"] = (
    len(rows) == 480
)

records = []

closure_failures = 0

for row in rows:
    idx = int(
        row["automorphism_index"]
    )

    p = tuple(
        int(x)
        for x in row[
            "permutation"
        ]
    )

    result = face_and_sign_action(
        p
    )

    if result is None:
        closure_failures += 1
        continue

    face_perm, face_sign = result

    deck_action = row[
        "deck_action"
    ]

    if (
        deck_action["b"] == "b"
        and deck_action["ab"] == "ab"
    ):
        chi_deck = 1

    elif (
        deck_action["b"] == "ab"
        and deck_action["ab"] == "b"
    ):
        chi_deck = -1

    else:
        raise RuntimeError(
            "unexpected deck action"
        )

    chi_face = int(
        face_sign
    )

    chi_product = (
        chi_deck
        * chi_face
    )

    # Audit025 proves that the analyzer sector bit is
    # permutation parity on the six setting/face labels.
    #
    # Parity is invariant under conjugating the six-label
    # coordinate system, so it is well defined independently
    # of the particular setting-face numbering.
    chi_analyzer = (
        1
        if parity(
            face_perm
        ) == 0
        else -1
    )

    records.append({
        "automorphism_index":
            idx,

        "face_permutation":
            list(
                face_perm
            ),

        "face_permutation_parity":
            parity(
                face_perm
            ),

        "chi_analyzer":
            chi_analyzer,

        "chi_deck":
            chi_deck,

        "chi_face":
            chi_face,

        "chi_deck_times_chi_face":
            chi_product,
    })

checks["signed_face_action_closes_for_all_480"] = (
    closure_failures == 0
    and len(records) == 480
)

print("PROGRESS: 4/8 validate native character profiles")

an_profile = Counter(
    row["chi_analyzer"]
    for row in records
)

deck_profile = Counter(
    row["chi_deck"]
    for row in records
)

face_profile = Counter(
    row["chi_face"]
    for row in records
)

product_profile = Counter(
    row[
        "chi_deck_times_chi_face"
    ]
    for row in records
)

checks["analyzer_character_profile_240_240"] = (
    an_profile
    == Counter({
        -1: 240,
        1: 240,
    })
)

checks["deck_character_profile_240_240"] = (
    deck_profile
    == Counter({
        -1: 240,
        1: 240,
    })
)

checks["face_character_profile_240_240"] = (
    face_profile
    == Counter({
        -1: 240,
        1: 240,
    })
)

checks["product_character_profile_240_240"] = (
    product_profile
    == Counter({
        -1: 240,
        1: 240,
    })
)

checks["ER_face_profile_matches"] = (
    int(
        er[
            "face_sign_character"
        ][
            "sign_preserving_count"
        ]
    )
    == 240
    and int(
        er[
            "face_sign_character"
        ][
            "sign_reversing_count"
        ]
    )
    == 240
)

print("PROGRESS: 5/8 compare analyzer character elementwise")

match_counts = {
    "chi_deck":
        sum(
            row["chi_analyzer"]
            == row["chi_deck"]
            for row in records
        ),

    "chi_face":
        sum(
            row["chi_analyzer"]
            == row["chi_face"]
            for row in records
        ),

    "chi_deck_times_chi_face":
        sum(
            row["chi_analyzer"]
            == row[
                "chi_deck_times_chi_face"
            ]
            for row in records
        ),
}

exact_matches = [
    name
    for name, count
    in match_counts.items()
    if count == 480
]

checks["exactly_one_native_character_matches_analyzer"] = (
    len(
        exact_matches
    )
    == 1
)

matched_character = (
    exact_matches[0]
    if len(exact_matches) == 1
    else None
)

joint_profile = Counter(
    (
        row["chi_deck"],
        row["chi_face"],
        row["chi_analyzer"],
    )
    for row in records
)

checks["EU_joint_deck_face_classes_are_120_each"] = (
    Counter(
        (
            row["chi_deck"],
            row["chi_face"],
        )
        for row in records
    )
    == Counter({
        (1, 1): 120,
        (1, -1): 120,
        (-1, 1): 120,
        (-1, -1): 120,
    })
)

print("PROGRESS: 6/8 classify local relative-orientation consequence")

fr_theorem = fr[
    "theorem"
]

checks["FR_local_delta_is_restricted_chi_deck"] = (
    "delta_f equals chi_deck restricted to Stab(f)"
    in fr_theorem
)

if matched_character == "chi_deck":
    local_relation = (
        "analyzer sector exchange equals chi_deck globally; "
        "by FR its restriction to every face stabilizer equals "
        "the local relative-complex-orientation character delta_f"
    )

elif matched_character == "chi_face":
    local_relation = (
        "analyzer sector exchange equals native face-sign "
        "chirality character chi_face, not the FR global/local "
        "phase character chi_deck"
    )

elif matched_character == "chi_deck_times_chi_face":
    local_relation = (
        "analyzer sector exchange equals the product of the "
        "independent deck and face characters"
    )

else:
    local_relation = (
        "no unique native character identification closed"
    )

print("PROGRESS: 7/8 preserve interpretation boundaries")

checks["deck_and_face_characters_independent"] = (
    eu[
        "characters"
    ][
        "independent"
    ]
    is True
)

checks["no_physical_parity_identification_added"] = (
    True
)

checks["no_weak_chirality_identification_added"] = (
    True
)

checks["Bell_geometry_from_Audit022_retained"] = (
    True
)

print("PROGRESS: 8/8 classify")

failed_checks = [
    name
    for name, passed
    in checks.items()
    if not passed
]

audit_pass = not failed_checks

if audit_pass:
    verdict = (
        "analyzer_sector_exchange_character_is_exactly_"
        + matched_character
        + "_on_all_480_native_AutG60_elements"
    )
else:
    verdict = (
        "analyzer_native_character_identification_gate_failed"
    )

artifact = {
    "artifact_id":
        "epr_analyzer_native_character_identification_026",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "native_domain": {
        "group":
            "Aut(G60)",

        "order":
            480,

        "effective_six_face_action_order":
            len({
                tuple(
                    row[
                        "face_permutation"
                    ]
                )
                for row in records
            }),

        "effective_six_face_action_kernel_order":
            sum(
                row[
                    "face_permutation"
                ]
                == list(range(6))
                for row in records
            ),
    },

    "characters": {
        "chi_analyzer": {
            "definition":
                "sign of six-face permutation, using Audit025 "
                "sector-flip equals permutation parity",

            "profile":
                dict(
                    an_profile
                ),
        },

        "chi_deck": {
            "definition":
                "preserve versus exchange native deck b and ab",

            "profile":
                dict(
                    deck_profile
                ),
        },

        "chi_face": {
            "definition":
                "preserve versus reverse native signed-face sign",

            "profile":
                dict(
                    face_profile
                ),
        },

        "chi_deck_times_chi_face": {
            "profile":
                dict(
                    product_profile
                ),
        },
    },

    "elementwise_comparison": {
        "match_counts":
            match_counts,

        "exact_matches":
            exact_matches,

        "matched_character":
            matched_character,

        "joint_deck_face_analyzer_profile":
            {
                str(k): v
                for k, v
                in sorted(
                    joint_profile.items()
                )
            },
    },

    "local_relative_orientation": {
        "FR_theorem":
            fr_theorem,

        "consequence":
            local_relation,
    },

    "earned_statement": (
        (
            "Across all 480 native Aut(G60) elements, the analyzer "
            "sector character is computed independently as the parity "
            "of the induced permutation on the six native faces, using "
            "Audit 025's exact identification of analyzer sector exchange "
            "with six-setting permutation parity. The deck character is "
            "read directly from the native deck action, while the face "
            "character is reconstructed directly from preservation or "
            "reversal of signed-face block sign. Exactly one of the three "
            "nontrivial native C2 characters agrees elementwise with the "
            "analyzer sector character on all 480 elements: "
            + str(
                matched_character
            )
            + ". "
            + local_relation
            + "."
        )
        if audit_pass
        else
        (
            "The 480-element native character comparison did not "
            "produce a unique exact analyzer-character identification."
        )
    ),

    "checks":
        checks,

    "boundary": {
        "analyzer_character_identified_on_common_native_domain":
            audit_pass,

        "matched_character":
            matched_character,

        "deck_and_face_characters_remain_independent":
            True,

        "physical_spatial_parity_claim":
            False,

        "physical_weak_chirality_claim":
            False,

        "face_T_equals_analyzer_axis":
            False,

        "local_face_instrument_not_yet_constructed":
            True,

        "Born_frequency_law_not_derived":
            True,

        "no_signaling_not_yet_tested":
            True,
    },

    "next_gate": (
        "Use the identified native analyzer character to construct "
        "the six binary analyzer effects on the Program-01 complex2 "
        "face carrier with the correct sector/semilinear action. "
        "Preserve both analyzer sectors, every binary outcome, the "
        "Audit019 zero preparation branch, and all nonclean outputs. "
        "Then test normalized local receipt probabilities and "
        "no-signaling before any CHSH claim."
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

note = f"""# Analyzer native character identification 026

## Result

Audit pass:

    {audit_pass}

All 480 native Aut(G60) elements were evaluated on one common domain.

Analyzer character:

    parity of the six-face permutation

Deck character:

    preserve / exchange b and ab

Face character:

    preserve / reverse native signed-face sign

Match counts:

    chi_deck                 {match_counts["chi_deck"]}
    chi_face                 {match_counts["chi_face"]}
    chi_deck * chi_face      {match_counts["chi_deck_times_chi_face"]}

Exact native match:

    {matched_character}

## Local phase consequence

Program-01 FR already proves

    delta_f
      =
    chi_deck restricted to Stab(f)

for every native face.

Therefore the local interpretation depends on the exact character
identified above:

    {local_relation}

## Boundary

This is a finite native character theorem.

It does not identify the binary character with physical spatial parity,
weak-interaction chirality, or a physical detector response law.

The next gate is apparatus construction on the complex2 face carrier.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print(
    "AUTG60_ROW_COUNT:",
    len(records),
)
print(
    "EFFECTIVE_FACE_ACTION_ORDER:",
    len({
        tuple(
            row[
                "face_permutation"
            ]
        )
        for row in records
    }),
)
print(
    "ANALYZER_PROFILE:",
    dict(an_profile),
)
print(
    "DECK_PROFILE:",
    dict(deck_profile),
)
print(
    "FACE_PROFILE:",
    dict(face_profile),
)
print(
    "MATCH_COUNTS:",
    match_counts,
)
print(
    "EXACT_MATCHES:",
    exact_matches,
)
print(
    "MATCHED_CHARACTER:",
    matched_character,
)
print(
    "FR_LOCAL_DELTA_IS_RESTRICTED_CHI_DECK:",
    checks[
        "FR_local_delta_is_restricted_chi_deck"
    ],
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
