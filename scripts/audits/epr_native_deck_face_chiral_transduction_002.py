#!/usr/bin/env python3

from collections import defaultdict, Counter
from hashlib import sha256
from pathlib import Path
import json

HERE = Path(__file__).resolve().parents[2]
QR = HERE.parent

P01 = QR / "1_full_face_hilbert_gluing"
P41 = (
    Path.home()
    / "dev/cori/research/mathematics"
    / "41-order-4-dodecahedral-residue"
)

AUT_PATH = (
    P41
    / "artifacts/json"
    / "native_g60_automorphism_deck_cache.v1.json"
)

SIGNED_PATH = (
    P41
    / "artifacts/json"
    / "signed_24_face_block_closure_audit_025.json"
)

KERNEL_PATH = (
    P41
    / "artifacts/json"
    / "g1800_silent_involution_deck_kernel_069.v1.json"
)

FACE_LOCK_PATH = (
    P01
    / "artifacts/json"
    / "003_surviving_face_spin_frontier_lock.json"
)

FACE_DESCENT_PATH = (
    P01
    / "artifacts/json"
    / "004_face_axis_realization_gauge_descent.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_native_deck_face_chiral_transduction_002.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_native_deck_face_chiral_transduction_002.md"
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


def perm_order(p):
    n = len(p)
    seen = [False] * n
    orders = []

    for i in range(n):
        if seen[i]:
            continue

        j = i
        length = 0

        while not seen[j]:
            seen[j] = True
            j = p[j]
            length += 1

        orders.append(length)

    if all(x == 1 for x in orders):
        return 1

    if all(x in (1, 2) for x in orders):
        return 2

    raise RuntimeError(
        "unexpected permutation order for selected deck row"
    )


print("PROGRESS: 1/7 load canonical interfaces")

aut = load(AUT_PATH)
signed = load(SIGNED_PATH)
kernel = load(KERNEL_PATH)
face_lock = load(FACE_LOCK_PATH)
face_descent = load(FACE_DESCENT_PATH)

aut_rows = {
    int(row["automorphism_index"]): row
    for row in aut["measurements"]["automorphism_rows"]
}

blocks = signed["measurements"]["signed_blocks"]

checks = {}

checks["signed_block_count_is_24"] = (
    len(blocks) == 24
)

checks["sealed_01_connection_residual_is_S1"] = (
    face_lock["checks"][
        "connection_residual_generator_S1"
    ]
    is True
)

checks["sealed_01_S1_is_sigma_x_K"] = (
    face_lock["checks"]["S1_sigma_x_K"]
    is True
)

checks["sealed_01_S1_fixes_T"] = (
    face_descent["checks"]["S1_fixes_T"]
    is True
)

checks["sealed_01_r1_not_S1"] = (
    face_lock["checks"]["r1_not_vertical_S1"]
    is True
)

print("PROGRESS: 2/7 recover native deck names")

deck_name_to_index = {}

for row in kernel["measurements"]["automorphism_rows"]:
    labels = json.loads(row["silent_kernel_labels"])

    if not labels:
        if row["is_identity"]:
            deck_name_to_index["identity"] = int(
                row["automorphism_index"]
            )
        continue

    if len(labels) != 1:
        raise RuntimeError(
            "expected one native deck label per silent involution"
        )

    name = labels[0]

    deck_name_to_index[name] = int(
        row["automorphism_index"]
    )

expected_names = {"identity", "a", "b", "ab"}

checks["native_deck_names_complete"] = (
    set(deck_name_to_index) == expected_names
)

deck_perms = {}

for name, index in deck_name_to_index.items():
    deck_perms[name] = tuple(
        int(v)
        for v in aut_rows[index]["permutation"]
    )

checks["native_deck_permutations_have_60_points"] = all(
    len(p) == 60
    for p in deck_perms.values()
)

checks["three_nonidentity_deck_elements_are_free"] = all(
    all(p[i] != i for i in range(60))
    for name, p in deck_perms.items()
    if name != "identity"
)

checks["three_nonidentity_deck_elements_have_order_2"] = all(
    perm_order(p) == 2
    for name, p in deck_perms.items()
    if name != "identity"
)

print("PROGRESS: 3/7 reconstruct signed face-block action")

block_by_states = {}

for index, block in enumerate(blocks):
    key = tuple(
        sorted(
            int(v)
            for v in block["five_state_block"]
        )
    )

    if key in block_by_states:
        raise RuntimeError(
            "duplicate signed block state set"
        )

    block_by_states[key] = index


def block_image(block_index, perm):
    block = blocks[block_index]

    image = tuple(
        sorted(
            perm[int(v)]
            for v in block["five_state_block"]
        )
    )

    return block_by_states.get(image)


positive_by_carrier = defaultdict(list)
negative_by_carrier = defaultdict(list)

for i, block in enumerate(blocks):
    carrier = int(block["carrier_index"])
    sign = block["sign"]

    if sign == "positive":
        positive_by_carrier[carrier].append(i)
    elif sign == "negative":
        negative_by_carrier[carrier].append(i)
    else:
        raise RuntimeError(
            f"unexpected signed block sign: {sign}"
        )

for carrier in range(6):
    positive_by_carrier[carrier].sort()
    negative_by_carrier[carrier].sort()

checks["six_face_carriers_present"] = (
    set(positive_by_carrier) == set(range(6))
    and set(negative_by_carrier) == set(range(6))
)

checks["two_positive_and_two_negative_blocks_per_face"] = all(
    len(positive_by_carrier[c]) == 2
    and len(negative_by_carrier[c]) == 2
    for c in range(6)
)


def carrier_of(block_index):
    return int(
        blocks[block_index]["carrier_index"]
    )


def sign_of(block_index):
    return blocks[block_index]["sign"]


action_rows = {}

for name, perm in deck_perms.items():
    image_map = {
        i: block_image(i, perm)
        for i in range(len(blocks))
    }

    missing = [
        i
        for i, image in image_map.items()
        if image is None
    ]

    carrier_preserved = []
    sign_preserved = []
    fixed_blocks = []

    for i, image in image_map.items():
        if image is None:
            continue

        if carrier_of(image) == carrier_of(i):
            carrier_preserved.append(i)

        if sign_of(image) == sign_of(i):
            sign_preserved.append(i)

        if image == i:
            fixed_blocks.append(i)

    within_pair_swap_faces = []

    for carrier in range(6):
        pos = positive_by_carrier[carrier]
        neg = negative_by_carrier[carrier]

        pos_images = {
            image_map[i]
            for i in pos
        }

        neg_images = {
            image_map[i]
            for i in neg
        }

        pos_swapped = (
            len(pos) == 2
            and None not in pos_images
            and pos_images == set(pos)
            and all(image_map[i] != i for i in pos)
        )

        neg_swapped = (
            len(neg) == 2
            and None not in neg_images
            and neg_images == set(neg)
            and all(image_map[i] != i for i in neg)
        )

        if pos_swapped and neg_swapped:
            within_pair_swap_faces.append(carrier)

    action_rows[name] = {
        "automorphism_index":
            deck_name_to_index[name],

        "missing_block_image_count":
            len(missing),

        "carrier_preserved_count":
            len(carrier_preserved),

        "sign_preserved_count":
            len(sign_preserved),

        "fixed_block_count":
            len(fixed_blocks),

        "within_positive_and_negative_pair_swap_face_count":
            len(within_pair_swap_faces),

        "all_faces_simultaneous_within_pair_swap":
            len(within_pair_swap_faces) == 6,

        "block_image_map": {
            str(i): image_map[i]
            for i in range(len(blocks))
        },
    }

checks["all_deck_actions_close_on_signed_blocks"] = all(
    row["missing_block_image_count"] == 0
    for row in action_rows.values()
)

print("PROGRESS: 4/7 identify native S1 deck element")

s1_candidates = []

for name, row in action_rows.items():
    if name == "identity":
        continue

    if (
        row["carrier_preserved_count"] == 24
        and row["sign_preserved_count"] == 24
        and row["fixed_block_count"] == 0
        and row[
            "all_faces_simultaneous_within_pair_swap"
        ]
    ):
        s1_candidates.append(name)

checks["exactly_one_native_nonidentity_matches_S1_action"] = (
    len(s1_candidates) == 1
)

s1_native_name = (
    s1_candidates[0]
    if len(s1_candidates) == 1
    else None
)

s1_native_index = (
    deck_name_to_index[s1_native_name]
    if s1_native_name is not None
    else None
)

print("PROGRESS: 5/7 instantiate native G1800 torsor")

native_g1800 = {}

if s1_native_name is not None:
    B = deck_perms[s1_native_name]

    orbit_rep = {}

    for x in range(60):
        orbit_rep[x] = min(x, B[x])

    local_orbits = sorted(
        set(orbit_rep.values())
    )

    local_orbit_index = {
        rep: i
        for i, rep in enumerate(local_orbits)
    }

    def local_q(x):
        return local_orbit_index[orbit_rep[x]]

    def diag_rep(u, v):
        return min(
            (u, v),
            (B[u], B[v]),
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

    fibers = defaultdict(list)
    tau = {}

    tau_left_right_equal = True

    for cid, (u, v) in enumerate(joint_reps):
        address = (
            local_q(u),
            local_q(v),
        )

        fibers[address].append(cid)

        left = joint_class(B[u], v)
        right = joint_class(u, B[v])

        if left != right:
            tau_left_right_equal = False

        tau[cid] = left

    fiber_profile = Counter(
        len(rows)
        for rows in fibers.values()
    )

    tau_fixed_point_free = all(
        tau[cid] != cid
        for cid in tau
    )

    tau_involution = all(
        tau[tau[cid]] == cid
        for cid in tau
    )

    native_g1800 = {
        "kernel_name":
            s1_native_name,

        "kernel_automorphism_index":
            s1_native_index,

        "local_orbit_count":
            len(local_orbits),

        "joint_class_count":
            len(joint_reps),

        "local_address_pair_count":
            len(fibers),

        "fiber_size_profile": {
            str(k): v
            for k, v in sorted(
                fiber_profile.items()
            )
        },

        "tau_left_equals_right":
            tau_left_right_equal,

        "tau_fixed_point_free":
            tau_fixed_point_free,

        "tau_is_involution":
            tau_involution,
    }

    checks["native_S1_kernel_has_30_local_orbits"] = (
        len(local_orbits) == 30
    )

    checks["native_S1_kernel_gives_1800_joint_classes"] = (
        len(joint_reps) == 1800
    )

    checks["native_S1_kernel_gives_900_local_address_pairs"] = (
        len(fibers) == 900
    )

    checks["native_S1_kernel_gives_exact_two_state_fibers"] = (
        fiber_profile == Counter({2: 900})
    )

    checks["native_tau_can_move_left_or_right"] = (
        tau_left_right_equal
    )

    checks["native_tau_is_fixed_point_free"] = (
        tau_fixed_point_free
    )

    checks["native_tau_is_involution"] = (
        tau_involution
    )

else:
    for name in (
        "native_S1_kernel_has_30_local_orbits",
        "native_S1_kernel_gives_1800_joint_classes",
        "native_S1_kernel_gives_900_local_address_pairs",
        "native_S1_kernel_gives_exact_two_state_fibers",
        "native_tau_can_move_left_or_right",
        "native_tau_is_fixed_point_free",
        "native_tau_is_involution",
    ):
        checks[name] = False

print("PROGRESS: 6/7 classify chiral transduction")

chiral_transduction_pass = all([
    checks[
        "exactly_one_native_nonidentity_matches_S1_action"
    ],
    checks[
        "sealed_01_connection_residual_is_S1"
    ],
    checks[
        "sealed_01_S1_is_sigma_x_K"
    ],
    checks[
        "sealed_01_S1_fixes_T"
    ],
    checks[
        "native_tau_can_move_left_or_right"
    ],
])

checks["joint_chiral_transduction_interface_closes"] = (
    chiral_transduction_pass
)

print("PROGRESS: 7/7 seal result")

failed_checks = [
    name
    for name, passed in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "native_g1800_relative_lift_exchange_has_exact_"
    "local_S1_anti_complex_face_transduction"
    if audit_pass
    else
    "native_deck_face_chiral_transduction_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_native_deck_face_chiral_transduction_002",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "native_deck_name_to_automorphism_index":
        deck_name_to_index,

    "signed_face_action":
        action_rows,

    "selected_native_S1_kernel": {
        "native_name":
            s1_native_name,

        "automorphism_index":
            s1_native_index,

        "face_operator":
            "S1",

        "event_adapted_form":
            "sigma_x K",

        "complex_character":
            "anti-complex",

        "J_action":
            "S1 J_F S1^-1 = -J_F",

        "T_action":
            "S1(T) = T",
    },

    "native_g1800_instantiation":
        native_g1800,

    "theorem": (
        "Among the three nonidentity native G60 deck "
        "involutions, exactly one induces on every signed "
        "Program-01 face the H65B connection-residual action: "
        "simultaneous swap within the positive registered pair "
        "and within the negative registered pair. The sealed "
        "Program-01 interface identifies that induced operator "
        "as S1=sigma_x K, which is anti-complex and fixes the "
        "neutral T axis. Using that same native involution as "
        "the diagonal G1800 kernel produces an exact two-state "
        "relative-lift torsor whose swap may be executed on "
        "either factor, [Bu,v]=[u,Bv]. Thus one upstairs joint "
        "relative-lift exchange has equivalent local "
        "Program-01 presentations as the selected anti-complex "
        "face reflection."
    ),

    "interpretation": {
        "upstairs":
            (
                "A two-state relative G1800 lift torsor over "
                "each ordered pair of local quotient addresses."
            ),

        "downstairs":
            (
                "The torsor exchange is presented on either "
                "local functional face as S1=sigma_x K."
            ),

        "preserved_local_structure":
            (
                "The neutral analyzer axis T survives the "
                "local chiral reflection exactly."
            ),

        "changed_local_structure":
            (
                "The local complex orientation J_F reverses."
            ),
    },

    "checks":
        checks,

    "boundary": {
        "S1_is_induced_face_operator_not_literal_G60_permutation":
            True,

        "tau_AB_is_not_identified_with_r1":
            True,

        "selected_transport_itself_is_not_S1":
            True,

        "relative_torsor_is_not_yet_Bell_resource":
            True,

        "no_joint_coherent_state_selected":
            True,

        "no_probability_law_derived":
            True,

        "no_Born_rule_assumed":
            True,

        "no_no_signaling_claim":
            True,

        "no_CHSH_claim":
            True,

        "face_Weyl_axis_S3_crosswalk_remains_open":
            True,
    },

    "next_gate": (
        "Test whether this native joint chiral torsor induces "
        "a coherent relative pairing of the two Program-01 "
        "complex doublets and whether that pairing selects an "
        "intrinsic joint projective line without inserting a "
        "singlet or Born weighting by hand."
    ),
}

artifact["artifact_sha256"] = digest_json(artifact)

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

note = f"""# Native deck to face chiral transduction 002

## Result

Audit pass:

    {audit_pass}

Native deck element inducing the selected Program-01 reflection:

    {s1_native_name}

Native automorphism-cache index:

    {s1_native_index}

The native deck action was applied directly to the canonical twenty-four
signed face blocks.

Exactly one nonidentity native deck involution acts on every one of the
six faces as a simultaneous swap inside both registered sign pairs:

    positive_0 <-> positive_1
    negative_0 <-> negative_1

Program 01 identifies that exact connection-residual action as

    S1 = (1,0,3,2)

and, in the event-adapted complex doublet,

    S1 = sigma_x K.

Therefore

    S1 J_F S1^-1 = -J_F

while

    S1(T) = T.

## G1800 transduction

Use the identified native deck involution as the diagonal quotient
kernel.

Then

    G1800 = (G60 x G60)/diag(B)

has 1800 joint states over 900 ordered local-address pairs, with two
joint lifts over every local pair.

The residual joint swap is

    tau_AB [u,v] = [Bu,v] = [u,Bv].

The equality means that the same upstairs relative-lift exchange can be
presented on Alice's factor or Bob's factor.

The direct signed-face calculation now supplies the local Program-01
meaning of that factor action:

    B -> S1 = sigma_x K.

Thus the G1800 relative-lift swap has an exact local anti-complex face
presentation on either wing.

## Earned statement

One upstairs native G1800 relative-lift exchange admits two equivalent
downstairs presentations as the selected Program-01 anti-complex
reflection.

The local complex orientation reverses.

The neutral local analyzer axis survives.

This is the first finite Joint Chiral Transduction theorem of Program 02.

## Boundary

B and S1 are not literally the same permutation on the same carrier.

B is the native G60 deck involution.

S1 is its induced operator on the sealed Program-01 face apparatus.

The selected transport itself is not S1.

The joint torsor is not yet a Bell resource.

No coherent joint Hilbert state has yet been selected.

No trial weighting law, probability table, no-signaling theorem, or CHSH
value is claimed.

The unresolved Program-01 face/Weyl S3 axis-slot crosswalk remains open.

## Next gate

Determine whether the native chiral torsor supplies a coherent relative
pairing of the two local complex doublets.

In particular, test whether an intrinsic one-dimensional joint projective
subspace is selected by the native relation itself.

Do not insert the quantum singlet as an answer.

Do not assume squared-norm frequencies.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print("DECK_NAME_TO_INDEX:", deck_name_to_index)

for name in ("identity", "a", "b", "ab"):
    row = action_rows[name]

    print(
        "ACTION:",
        name,
        "index=",
        row["automorphism_index"],
        "carrier_preserved=",
        row["carrier_preserved_count"],
        "sign_preserved=",
        row["sign_preserved_count"],
        "fixed_blocks=",
        row["fixed_block_count"],
        "within_pair_swap_faces=",
        row[
            "within_positive_and_negative_pair_swap_face_count"
        ],
    )

print("S1_NATIVE_NAME:", s1_native_name)
print("S1_NATIVE_INDEX:", s1_native_index)

if native_g1800:
    print(
        "G1800_CLASS_COUNT:",
        native_g1800["joint_class_count"],
    )

    print(
        "G1800_FIBER_SIZE_PROFILE:",
        native_g1800["fiber_size_profile"],
    )

    print(
        "TAU_LEFT_EQUALS_RIGHT:",
        native_g1800["tau_left_equals_right"],
    )

print(
    "JOINT_CHIRAL_TRANSDUCTION_INTERFACE_CLOSES:",
    checks[
        "joint_chiral_transduction_interface_closes"
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


if __name__ == "__main__":
    pass
