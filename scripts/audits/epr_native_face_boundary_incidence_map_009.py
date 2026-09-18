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

A002 = (
    HERE
    / "artifacts/json"
    / "epr_native_deck_face_chiral_transduction_002.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_native_face_boundary_incidence_map_009.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_native_face_boundary_incidence_map_009.md"
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


print("PROGRESS: 1/7 load canonical inputs")

signed = load(SIGNED)
aut = load(AUT)
a002 = load(A002)

checks = {}

checks["audit_002_passed"] = (
    a002["audit_pass"] is True
)

checks["audit_002_selected_native_a"] = (
    a002[
        "selected_native_S1_kernel"
    ]["native_name"]
    == "a"
)

checks["audit_002_local_operator_is_S1"] = (
    a002[
        "selected_native_S1_kernel"
    ]["face_operator"]
    == "S1"
)

blocks = signed["measurements"]["signed_blocks"]

checks["signed_block_count_24"] = (
    len(blocks) == 24
)

print("PROGRESS: 2/7 construct incidence domain")

incidences = []

state_to_blocks = defaultdict(list)
state_to_positive = defaultdict(list)
state_to_negative = defaultdict(list)
carrier_to_blocks = defaultdict(list)

for block_index, block in enumerate(blocks):
    carrier = int(block["carrier_index"])
    sign = block["sign"]

    carrier_to_blocks[carrier].append(
        block_index
    )

    for state in block["five_state_block"]:
        state = int(state)

        row = {
            "state": state,
            "block_index": block_index,
            "carrier": carrier,
            "sign": sign,
        }

        incidences.append(row)
        state_to_blocks[state].append(
            block_index
        )

        if sign == "positive":
            state_to_positive[state].append(
                block_index
            )
        elif sign == "negative":
            state_to_negative[state].append(
                block_index
            )

checks["incidence_count_120"] = (
    len(incidences) == 120
)

checks["all_60_states_present"] = (
    set(state_to_blocks) == set(range(60))
)

checks["every_state_has_two_signed_face_incidences"] = all(
    len(state_to_blocks[u]) == 2
    for u in range(60)
)

checks["every_state_has_one_positive_incidence"] = all(
    len(state_to_positive[u]) == 1
    for u in range(60)
)

checks["every_state_has_one_negative_incidence"] = all(
    len(state_to_negative[u]) == 1
    for u in range(60)
)

print("PROGRESS: 3/7 construct six local face modules")

carrier_profiles = {}

for carrier in range(6):
    rows = sorted(
        carrier_to_blocks[carrier]
    )

    signs = Counter(
        blocks[i]["sign"]
        for i in rows
    )

    carrier_profiles[str(carrier)] = {
        "block_indices": rows,
        "block_count": len(rows),
        "sign_profile": dict(
            sorted(signs.items())
        ),
        "real_module_dimension": len(rows),
    }

checks["six_carriers_present"] = (
    set(carrier_to_blocks)
    == set(range(6))
)

checks["four_basis_blocks_per_carrier"] = all(
    len(carrier_to_blocks[c]) == 4
    for c in range(6)
)

checks["two_positive_two_negative_per_carrier"] = all(
    Counter(
        blocks[i]["sign"]
        for i in carrier_to_blocks[c]
    )
    == Counter({
        "positive": 2,
        "negative": 2,
    })
    for c in range(6)
)

print("PROGRESS: 4/7 measure state-only ambiguity")

same_carrier_count = 0
different_carrier_count = 0

state_rows = []

for u in range(60):
    bi = sorted(
        state_to_blocks[u]
    )

    carriers = [
        int(blocks[i]["carrier_index"])
        for i in bi
    ]

    signs = [
        blocks[i]["sign"]
        for i in bi
    ]

    same_carrier = (
        carriers[0] == carriers[1]
    )

    if same_carrier:
        same_carrier_count += 1
    else:
        different_carrier_count += 1

    state_rows.append({
        "state": u,
        "block_indices": bi,
        "carriers": carriers,
        "signs": signs,
        "same_carrier": same_carrier,
    })

checks["state_only_map_is_not_single_valued"] = all(
    len(state_to_blocks[u]) == 2
    for u in range(60)
)

print("PROGRESS: 5/7 reconstruct native a action")

a_index = int(
    a002[
        "selected_native_S1_kernel"
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

    return block_by_states.get(
        image_key
    )


a_block = {
    i: block_image(i)
    for i in range(24)
}

checks["a_closes_on_all_signed_blocks"] = all(
    a_block[i] is not None
    for i in range(24)
)

checks["a_preserves_carrier"] = all(
    int(
        blocks[a_block[i]][
            "carrier_index"
        ]
    )
    ==
    int(
        blocks[i]["carrier_index"]
    )
    for i in range(24)
)

checks["a_preserves_sign"] = all(
    blocks[a_block[i]]["sign"]
    == blocks[i]["sign"]
    for i in range(24)
)

checks["a_has_no_fixed_signed_block"] = all(
    a_block[i] != i
    for i in range(24)
)

checks["a_block_action_is_involution"] = all(
    a_block[a_block[i]] == i
    for i in range(24)
)

print("PROGRESS: 6/7 verify vector-level incidence covariance")

# For each carrier c, the local real face module V_c
# is the free real vector space on its four signed blocks.
#
# beta(u,B) = e_B.
#
# The induced action rho_c(a) is the permutation of those
# four basis vectors generated by B -> aB.
#
# Covariance is therefore checked on every actual incidence:
#
# beta(a u, a B) = rho_c(a) beta(u,B).

incidence_lookup = {
    (
        row["state"],
        row["block_index"],
    )
    for row in incidences
}

covariance_failures = []

for row in incidences:
    u = row["state"]
    B = row["block_index"]

    au = a[u]
    aB = a_block[B]

    if (au, aB) not in incidence_lookup:
        covariance_failures.append({
            "state": u,
            "block_index": B,
            "image_state": au,
            "image_block_index": aB,
        })

checks["all_120_incidence_vectors_are_a_covariant"] = (
    len(covariance_failures) == 0
)

# Per carrier, a must act as two transpositions:
# one inside the positive pair and one inside the negative pair.
#
# This is exactly the basis-free content of S1=(1,0,3,2).

carrier_a_cycle_profiles = {}

for carrier in range(6):
    local_blocks = sorted(
        carrier_to_blocks[carrier]
    )

    local_set = set(
        local_blocks
    )

    visited = set()
    cycle_sizes = []

    for B in local_blocks:
        if B in visited:
            continue

        cycle = []
        x = B

        while x not in visited:
            visited.add(x)
            cycle.append(x)
            x = a_block[x]

        cycle_sizes.append(
            len(cycle)
        )

    cycle_sizes.sort()

    sign_pairs_preserved = all(
        a_block[B] in local_set
        and blocks[a_block[B]]["sign"]
        == blocks[B]["sign"]
        for B in local_blocks
    )

    carrier_a_cycle_profiles[
        str(carrier)
    ] = {
        "cycle_sizes": cycle_sizes,
        "sign_pairs_preserved":
            sign_pairs_preserved,
    }

checks["a_is_two_transpositions_on_every_face_module"] = all(
    carrier_a_cycle_profiles[
        str(c)
    ]["cycle_sizes"]
    == [2, 2]
    for c in range(6)
)

checks["a_matches_basis_free_S1_action"] = (
    checks[
        "a_is_two_transpositions_on_every_face_module"
    ]
    and checks["a_preserves_sign"]
    and checks["audit_002_local_operator_is_S1"]
)

print("PROGRESS: 7/7 classify")

failed_checks = [
    name
    for name, passed in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "native_boundary_presentation_is_a_120_incidence_map_"
    "into_six_four_real_face_modules_with_exact_a_to_S1_"
    "vector_covariance"
    if audit_pass
    else
    "native_face_boundary_incidence_map_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_native_face_boundary_incidence_map_009",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "domain": {
        "type":
            "signed state-face incidence domain",

        "definition":
            "I_F = {(u,B): u belongs to signed face block B}",

        "incidence_count":
            len(incidences),

        "G60_state_count":
            60,

        "signed_block_count":
            24,
    },

    "local_face_modules":
        carrier_profiles,

    "state_only_ambiguity": {
        "every_state_incidence_count":
            2,

        "positive_incidence_per_state":
            1,

        "negative_incidence_per_state":
            1,

        "same_carrier_state_count":
            same_carrier_count,

        "different_carrier_state_count":
            different_carrier_count,

        "state_only_boundary_map_single_valued":
            False,
    },

    "boundary_map": {
        "definition":
            "beta_c(u,B) = e_B in V_c = R[B_c]",

        "local_real_dimension":
            4,

        "basis_objects":
            "the four canonical signed blocks on the selected carrier",

        "requires_boundary_incidence_context":
            True,
    },

    "native_a_action": {
        "native_name":
            "a",

        "automorphism_index":
            a_index,

        "induced_face_operator":
            "S1",

        "carrier_cycle_profiles":
            carrier_a_cycle_profiles,

        "incidence_covariance_failure_count":
            len(covariance_failures),

        "vector_covariance":
            "beta(a u,a B) = S1 beta(u,B)",
    },

    "earned_statement": (
        "The correct native local face presentation is not a "
        "state-only map from G60 into one doublet. It is a "
        "state-face incidence map on 120 native incidences. "
        "Each carrier supplies four canonical signed face-block "
        "basis objects, giving a four-real-dimensional local "
        "module. Every G60 state has exactly two boundary "
        "presentations, one positive and one negative. The "
        "native deck involution a acts equivariantly on all "
        "120 incidences and induces exactly the Program-01 "
        "S1 action on each local four-real-dimensional face "
        "module."
    ),

    "checks":
        checks,

    "boundary": {
        "state_only_map_rejected":
            True,

        "incidence_map_is_pre_measurement_boundary_presentation":
            True,

        "no_event_outcome_probability_assigned":
            True,

        "complex_doublet_coordinate_not_yet_needed":
            True,

        "joint_G1800_boundary_map_not_yet_constructed":
            True,

        "no_Born_rule":
            True,

        "no_CHSH_claim":
            True,
    },

    "next_gate": (
        "Lift the one-wing incidence presentation to the "
        "a-kernel G1800 common carrier. Construct the joint "
        "boundary-incidence domain carrying one Alice signed "
        "face incidence and one Bob signed face incidence, "
        "then test descent through (u,v)~(au,av) using the "
        "exact S1 covariance established here. The quotient "
        "should recover the Audit-003 chiral compatibility "
        "condition rather than imposing it by hand."
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

note = f"""# Native face boundary incidence map 009

## Result

Audit pass:

    {audit_pass}

The correct Program-01 local boundary domain is not bare G60.

It is the incidence set

    I_F = {{(u,B): u belongs to signed face block B}}.

There are

    {len(incidences)}

such incidences.

Each of the sixty G60 states occurs exactly twice:

    once in a positive signed block;
    once in a negative signed block.

Therefore a G60 state alone does not specify one local face vector.

## Local face module

Each of the six face carriers has exactly four signed blocks:

    two positive;
    two negative.

These four blocks are the canonical basis objects of the four-real-
dimensional Program-01 face permutation module.

For carrier c define

    V_c = R[B_c].

The native boundary presentation is

    beta_c(u,B) = e_B.

No arbitrary within-face coordinate numbering is required.

## Chiral covariance

The native deck involution

    a = Aut(G60)[{a_index}]

closes exactly on the signed incidence domain.

On each face it acts as two transpositions:

    one on the positive pair;
    one on the negative pair.

This is the basis-free signed-block action identified in Program 01 as

    S1 = sigma_x K.

Therefore all 120 incidences satisfy

    beta(a u,a B)
      =
    S1 beta(u,B).

## State-only ambiguity

States whose two signed incidences use the same carrier:

    {same_carrier_count}

States whose two signed incidences use different carriers:

    {different_carrier_count}

Either way, the signed boundary incidence is required to specify the
basis object.

## Next gate

Construct the two-wing incidence presentation on

    G1800 = (G60 x G60)/diag(a).

The diagonal quotient must descend through the local S1 actions.

If the construction is correct, the chiral compatibility structure found
abstractly in Audit 003 should reappear from this native incidence map
rather than being imposed afterward.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print("INCIDENCE_COUNT:", len(incidences))
print(
    "STATE_INCIDENCE_PROFILE:",
    Counter(
        len(state_to_blocks[u])
        for u in range(60)
    ),
)
print(
    "SAME_CARRIER_STATE_COUNT:",
    same_carrier_count,
)
print(
    "DIFFERENT_CARRIER_STATE_COUNT:",
    different_carrier_count,
)
print(
    "A_INCIDENCE_COVARIANCE_FAILURE_COUNT:",
    len(covariance_failures),
)
print(
    "A_MATCHES_BASIS_FREE_S1:",
    checks["a_matches_basis_free_S1_action"],
)
print("FAILED_CHECK_COUNT:", len(failed_checks))
print("FAILED_CHECKS:", failed_checks)
print("JSON_OUT:", JSON_OUT)
print("NOTE_OUT:", NOTE_OUT)
print(
    "JSON_SHA256:",
    sha256(JSON_OUT.read_bytes()).hexdigest(),
)
