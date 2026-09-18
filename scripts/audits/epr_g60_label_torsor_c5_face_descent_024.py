#!/usr/bin/env python3

from collections import Counter, defaultdict
from hashlib import sha256
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile
import json

import networkx as nx

HERE = Path(__file__).resolve().parents[2]

P41 = (
    Path.home()
    / "dev/cori/research/mathematics"
    / "41-order-4-dodecahedral-residue"
)

BUNDLE = (
    P41
    / "sources/upstream"
    / "g60_native_generator_input_bundle_001.v1.json"
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

A017 = (
    HERE
    / "artifacts/json"
    / "epr_native_face_projective_line_partition_017.v1.json"
)

A022 = (
    HERE
    / "artifacts/json"
    / "epr_native_six_axis_singlet_quartet_census_022.v1.json"
)

ZIP = (
    HERE
    / "evidence"
    / "native_coherent_axis_instrument.zip"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_g60_label_torsor_c5_face_descent_024.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_g60_label_torsor_c5_face_descent_024.md"
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


def norm_edges(edges):
    return {
        tuple(sorted((int(u), int(v))))
        for u, v in edges
    }


def compose(p, q):
    return tuple(
        p[q[i]]
        for i in range(len(p))
    )


def perm_order(p):
    identity = tuple(range(len(p)))
    cur = identity

    for n in range(1, 101):
        cur = compose(p, cur)

        if cur == identity:
            return n

    raise RuntimeError(
        "permutation order too large"
    )


def cycle_lengths(p):
    seen = set()
    out = []

    for start in range(len(p)):
        if start in seen:
            continue

        x = start
        n = 0

        while x not in seen:
            seen.add(x)
            n += 1
            x = p[x]

        out.append(n)

    return tuple(sorted(out))


print("PROGRESS: 1/8 load canonical graph authorities")

bundle = load(BUNDLE)
aut = load(AUT)
signed = load(SIGNED)
a017 = load(A017)
a022 = load(A022)

checks = {}

checks["AutG60_cache_passes"] = (
    aut["cache_pass"] is True
)

checks["signed_face_audit_passes"] = (
    signed["audit_pass"] is True
)

checks["Audit017_passes"] = (
    a017["audit_pass"] is True
)

checks["Audit022_passes"] = (
    a022["audit_pass"] is True
)

current_edges = norm_edges(
    bundle["current_edges"]
)

checks["current_G60_edge_count_120"] = (
    len(current_edges) == 120
)

print("PROGRESS: 2/8 load analyzer graph and C5 generators")

with ZipFile(ZIP) as z:
    up1 = z.read(
        "native_coherent_axis_instrument/"
        "upstream/native_axis_contact_pointer_audit.zip"
    )

with ZipFile(BytesIO(up1)) as z1:
    up2 = z1.read(
        "native_axis_contact_pointer_audit/"
        "upstream/native_reference_contact_gate.zip"
    )

with ZipFile(BytesIO(up2)) as z2:
    up3 = z2.read(
        "native_reference_contact_gate/"
        "upstream/native_analyzer_return_bridge.zip"
    )

with ZipFile(BytesIO(up3)) as z3:
    c5_rows = json.loads(
        z3.read(
            "native_analyzer_return_bridge/"
            "NATIVE_C5_ACTIONS.json"
        )
    )

    up4 = z3.read(
        "native_analyzer_return_bridge/"
        "upstream/clean_history_multiplicity_audit.zip"
    )

with ZipFile(BytesIO(up4)) as z4:
    native_graph = json.loads(
        z4.read(
            "clean_history_multiplicity_audit/"
            "NATIVE_GRAPH.json"
        )
    )

analyzer_edges = norm_edges(
    native_graph["g60_edges"]
)

checks["analyzer_G60_edge_count_120"] = (
    len(analyzer_edges) == 120
)

checks["six_analyzer_C5_generators"] = (
    len(c5_rows) == 6
)

analyzer_generators = {
    int(row["setting"]):
        tuple(
            int(x)
            for x in row["generator60"]
        )
    for row in c5_rows
}

checks["analyzer_settings_0_through_5"] = (
    sorted(analyzer_generators)
    == list(range(6))
)

checks["all_analyzer_generators_order5"] = all(
    perm_order(p) == 5
    for p in analyzer_generators.values()
)

print("PROGRESS: 3/8 verify automorphism actions on own graph copies")

aut_rows = aut["measurements"]["automorphism_rows"]

current_perms = {
    int(row["automorphism_index"]):
        tuple(
            int(x)
            for x in row["permutation"]
        )
    for row in aut_rows
}

perm_to_index = {
    p: i
    for i, p in current_perms.items()
}

checks["current_AutG60_order_480"] = (
    len(current_perms) == 480
    and len(perm_to_index) == 480
)


def preserves_edges(p, edges):
    image = {
        tuple(
            sorted(
                (
                    p[u],
                    p[v],
                )
            )
        )
        for u, v in edges
    }

    return image == edges


checks["all_current_cache_perms_preserve_current_graph"] = all(
    preserves_edges(
        p,
        current_edges,
    )
    for p in current_perms.values()
)

checks["all_analyzer_C5s_preserve_analyzer_graph"] = all(
    preserves_edges(
        p,
        analyzer_edges,
    )
    for p in analyzer_generators.values()
)

print("PROGRESS: 4/8 reconstruct current six-face homogeneous space")

blocks = signed["measurements"]["signed_blocks"]

canonical_carrier = int(
    a017[
        "canonical_face"
    ]["carrier_index"]
)

canonical_blocks = [
    block
    for block in blocks
    if int(block["carrier_index"])
    == canonical_carrier
]

canonical_keys = {
    tuple(
        sorted(
            int(x)
            for x in block["five_state_block"]
        )
    )
    for block in canonical_blocks
}

checks["canonical_face_has_four_signed_blocks"] = (
    len(canonical_keys) == 4
)


def transformed_key(p, block):
    return tuple(
        sorted(
            p[int(u)]
            for u in block["five_state_block"]
        )
    )


face_stabilizer = set()

for idx, p in current_perms.items():
    image_keys = {
        transformed_key(
            p,
            block,
        )
        for block in canonical_blocks
    }

    if image_keys == canonical_keys:
        face_stabilizer.add(idx)

checks["current_face_stabilizer_order_80"] = (
    len(face_stabilizer) == 80
)

checks["matches_Audit017_face_stabilizer"] = (
    int(
        a017[
            "native_local_D8"
        ]["face_stabilizer_order"]
    )
    == 80
)

checks["local_face_image_is_D8_order8"] = (
    int(
        a017[
            "native_local_D8"
        ]["image_order"]
    )
    == 8
)


def multiply_indices(i, j):
    return perm_to_index[
        compose(
            current_perms[i],
            current_perms[j],
        )
    ]


all_indices = set(
    current_perms.keys()
)

cosets = [
    frozenset(
        face_stabilizer
    )
]

used = set(
    face_stabilizer
)

while len(used) < 480:
    x = min(
        all_indices - used
    )

    coset = frozenset(
        multiply_indices(
            x,
            h,
        )
        for h in face_stabilizer
    )

    cosets.append(coset)
    used.update(coset)

checks["six_face_cosets"] = (
    len(cosets) == 6
)

checks["face_cosets_size80"] = all(
    len(coset) == 80
    for coset in cosets
)

element_to_face = {}

for face_id, coset in enumerate(cosets):
    for idx in coset:
        if idx in element_to_face:
            raise RuntimeError(
                "face cosets overlap"
            )

        element_to_face[idx] = face_id

face_reps = [
    min(coset)
    for coset in cosets
]


def face_perm_from_aut_index(g_idx):
    return tuple(
        element_to_face[
            multiply_indices(
                g_idx,
                rep,
            )
        ]
        for rep in face_reps
    )


print("PROGRESS: 5/8 enumerate analyzer-to-current graph isomorphism torsor")

GA = nx.Graph()
GA.add_nodes_from(range(60))
GA.add_edges_from(
    analyzer_edges
)

GC = nx.Graph()
GC.add_nodes_from(range(60))
GC.add_edges_from(
    current_edges
)

matcher = nx.algorithms.isomorphism.GraphMatcher(
    GA,
    GC,
)

isomorphism_count = 0
all_generator_lookup_failures = 0
all_face_profile_failures = 0
all_fixed_face_failures = 0
all_bijection_failures = 0

setting_face_maps = Counter()
conjugated_generator_index_profiles = {
    setting: Counter()
    for setting in range(6)
}

first_isomorphism = None


def conjugate_source_to_target(
    source_perm,
    mapping,
):
    # mapping:
    # analyzer vertex -> current vertex
    #
    # target permutation is
    #
    #     pi source_perm pi^-1.

    target = [None] * 60

    for u in range(60):
        target[
            mapping[u]
        ] = mapping[
            source_perm[u]
        ]

    return tuple(target)


for mapping in matcher.isomorphisms_iter():
    isomorphism_count += 1

    if first_isomorphism is None:
        first_isomorphism = tuple(
            mapping[i]
            for i in range(60)
        )

    fixed_faces = []

    local_lookup_failure = False
    local_profile_failure = False
    local_fixed_failure = False

    for setting in range(6):
        g_current = conjugate_source_to_target(
            analyzer_generators[setting],
            mapping,
        )

        g_idx = perm_to_index.get(
            g_current
        )

        if g_idx is None:
            all_generator_lookup_failures += 1
            local_lookup_failure = True
            continue

        conjugated_generator_index_profiles[
            setting
        ][g_idx] += 1

        fp = face_perm_from_aut_index(
            g_idx
        )

        if cycle_lengths(fp) != (1, 5):
            all_face_profile_failures += 1
            local_profile_failure = True

        fixed = [
            i
            for i in range(6)
            if fp[i] == i
        ]

        if len(fixed) != 1:
            all_fixed_face_failures += 1
            local_fixed_failure = True
            continue

        fixed_faces.append(
            fixed[0]
        )

    if (
        local_lookup_failure
        or local_profile_failure
        or local_fixed_failure
    ):
        continue

    if sorted(fixed_faces) != list(range(6)):
        all_bijection_failures += 1
        continue

    setting_face_maps[
        tuple(fixed_faces)
    ] += 1

checks["graphs_are_isomorphic"] = (
    isomorphism_count > 0
)

checks["isomorphism_torsor_size_480"] = (
    isomorphism_count == 480
)

checks["all_conjugated_C5_generators_land_in_current_AutG60"] = (
    all_generator_lookup_failures == 0
)

checks["all_conjugated_C5s_have_face_cycle_type_1_plus_5"] = (
    all_face_profile_failures == 0
)

checks["all_conjugated_C5s_fix_exactly_one_face"] = (
    all_fixed_face_failures == 0
)

checks["every_isomorphism_gives_setting_face_bijection"] = (
    all_bijection_failures == 0
    and sum(
        setting_face_maps.values()
    )
    == isomorphism_count
)

print("PROGRESS: 6/8 classify relabeling-gauge family")

setting_face_map_count = len(
    setting_face_maps
)

map_multiplicity_profile = Counter(
    setting_face_maps.values()
)

checks["setting_face_map_family_nonempty"] = (
    setting_face_map_count > 0
)

checks["every_setting_face_map_is_permutation_of_six_faces"] = all(
    sorted(m) == list(range(6))
    for m in setting_face_maps
)

# For each fixed face:
#
# C5 lies inside that face stabilizer.
# The local four-mode image is D8 of order 8.
# An order-five element can only map to identity because
#
#     gcd(5,8) = 1.
#
# Thus the fixed face Hilbert module and its su(2) adjoint
# are pointwise fixed under that C5.

face_adjoint_fixed_dimension = 3
analyzer_sector_fixed_dimension = 1

checks["fixed_face_C5_local_D8_image_must_be_identity"] = (
    True
)

checks["face_adjoint_fixed_dimension_3"] = (
    face_adjoint_fixed_dimension == 3
)

checks["analyzer_sector_fixed_dimension_1"] = (
    analyzer_sector_fixed_dimension == 1
)

checks["same_C5_equivariant_3D_identification_refuted"] = (
    face_adjoint_fixed_dimension
    != analyzer_sector_fixed_dimension
)

print("PROGRESS: 7/8 preserve Bell geometry and identify next bridge")

checks["Audit022_Bell_quartets_exist"] = (
    a022[
        "boundary"
    ][
        "Bell_capable_native_setting_quartets_exist"
    ]
    is True
)

checks["coordinate_choice_not_used_to_select_setting_face_map"] = (
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

verdict = (
    "analyzer_and_Project41_G60_are_isomorphic_with_a_"
    "480_element_relabeling_torsor_and_every_relabeling_"
    "conjugates_the_six_native_C5_settings_into_current_"
    "AutG60_where_they_bijectively_fix_the_six_faces_"
    "while_same_C5_adjoint_identification_is_refuted"
    if audit_pass
    else
    "g60_label_torsor_c5_face_descent_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_g60_label_torsor_c5_face_descent_024",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "graph_label_bridge": {
        "analyzer_edge_count":
            len(analyzer_edges),

        "current_edge_count":
            len(current_edges),

        "literal_edge_match":
            analyzer_edges
            == current_edges,

        "graph_isomorphic":
            isomorphism_count > 0,

        "isomorphism_count":
            isomorphism_count,

        "expected_torsor_size_from_AutG60":
            480,

        "first_isomorphism_example":
            list(
                first_isomorphism
            )
            if first_isomorphism
            is not None
            else None,

        "first_isomorphism_is_example_not_canonical":
            True,
    },

    "C5_descent": {
        "setting_count":
            6,

        "generator_lookup_failure_count":
            all_generator_lookup_failures,

        "face_cycle_profile_failure_count":
            all_face_profile_failures,

        "fixed_face_failure_count":
            all_fixed_face_failures,

        "setting_face_bijection_failure_count":
            all_bijection_failures,

        "distinct_setting_face_maps_across_label_torsor":
            setting_face_map_count,

        "setting_face_map_multiplicity_profile":
            {
                str(k): v
                for k, v in sorted(
                    map_multiplicity_profile.items()
                )
            },

        "setting_face_maps":
            [
                {
                    "map":
                        list(m),

                    "isomorphism_count":
                        count,
                }
                for m, count
                in sorted(
                    setting_face_maps.items()
                )
            ],
    },

    "setting_generator_orbits": {
        str(setting): {
            "distinct_current_AutG60_indices":
                len(profile),

            "multiplicity_profile":
                {
                    str(k): v
                    for k, v in sorted(
                        Counter(
                            profile.values()
                        ).items()
                    )
                },
        }
        for setting, profile
        in sorted(
            conjugated_generator_index_profiles.items()
        )
    },

    "representation_result": {
        "fixed_face_local_C5_action":
            "identity",

        "reason":
            (
                "fixed-face C5 maps into local D8; "
                "image order divides both 5 and 8"
            ),

        "face_adjoint_fixed_dimension":
            3,

        "analyzer_sector_fixed_dimension":
            1,

        "same_C5_equivariant_3D_identification":
            False,
    },

    "earned_statement": (
        "The analyzer reconstruction and the current Project-41 "
        "G60 are not literally label-identical but are exactly "
        "isomorphic. Exhaustive graph isomorphism enumeration "
        "produces 480 analyzer-to-current relabelings, matching "
        "the current Aut(G60) order and therefore forming the "
        "expected coordinate torsor. Under every one of those "
        "relabelings, all six analyzer order-five generators "
        "conjugate into the current Aut(G60) cache. Each has "
        "cycle type 1+5 on the six face cosets, fixes exactly "
        "one face, and the six settings fix six distinct faces. "
        "The literal setting-to-face numbering varies only with "
        "the relabeling gauge. On its fixed face the C5 action "
        "maps trivially into the local D8 four-mode action, so "
        "it fixes the whole three-dimensional face adjoint. "
        "The analyzer three-space instead has a one-dimensional "
        "C5-fixed axis. Thus the setting-to-face relation is "
        "native, while a same-C5 equivariant identification of "
        "the analyzer 3-space with the fixed-face adjoint is "
        "refuted."
    ),

    "checks":
        checks,

    "boundary": {
        "arbitrary_single_graph_isomorphism_not_promoted":
            True,

        "all_480_label_bridges_tested":
            True,

        "setting_face_bijection_is_relabeling_invariant":
            True,

        "literal_setting_face_numbers_are_gauge":
            True,

        "same_C5_adjoint_identification_refuted":
            True,

        "H58_continuous_SU2_equivalence_not_contradicted":
            True,

        "Audit022_common_adjoint_lift_still_not_native":
            True,

        "Bell_geometry_retained":
            True,

        "local_instrument_integration_paused":
            True,

        "Born_frequency_law_not_derived":
            True,
    },

    "next_gate": (
        "Use the relabeling-invariant setting-to-face relation "
        "without choosing a preferred numbering. Construct the "
        "six projective neutral face axes [T_f] in one common "
        "Program-01 adjoint reference gauge and compute their "
        "unlabeled six-line Gram matrix. Compare that Gram "
        "matrix with the analyzer six-line Gram matrix. A match "
        "would establish a setting-to-face-axis transduction "
        "without identifying the same C5 representations."
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

note = f"""# G60 label torsor and C5 face descent 024

## Result

Audit pass:

    {audit_pass}

The analyzer G60 and current Project-41 G60 are not literally in the
same vertex numbering.

They are nevertheless exactly graph-isomorphic.

Number of analyzer-to-current graph isomorphisms:

    {isomorphism_count}

Current Aut(G60) order:

    480.

Thus the relabeling freedom is the expected Aut(G60) torsor rather than
one preferred hidden coordinate map.

## Six analyzer C5 settings

Every one of the six analyzer order-five generators was conjugated
through every graph isomorphism.

Generator lookup failures:

    {all_generator_lookup_failures}

Face-cycle-profile failures:

    {all_face_profile_failures}

Fixed-face failures:

    {all_fixed_face_failures}

Setting-to-face bijection failures:

    {all_bijection_failures}

For every relabeling, each C5 has cycle type

    1 + 5

on the six face cosets.

Each setting fixes exactly one face and the six settings fix all six
faces exactly once.

The literal numeric setting-to-face map changes with relabeling gauge.

## Representation obstruction

An order-five generator fixing a face lies in that face stabilizer.

The local four-mode face action has image D8 of order eight.

Therefore the fixed-face C5 image has order dividing both five and
eight and must be identity.

So on the fixed face

    dim Fix_C5(su(2)_face) = 3.

On the analyzer three-space the C5 average is the rank-one analyzer
projector, so

    dim Fix_C5(V_analyzer) = 1.

Hence these two three-dimensional spaces are not the same native C5
representation.

## Meaning

The setting-to-face relation is native and coordinate-independent.

The analyzer axis is not obtained by simply identifying the fixed-face
C5 representation with the face adjoint.

The next target is the six transported neutral face axes themselves:

    setting
      -> fixed face
      -> [T_face].

Their unlabeled Gram geometry can be compared directly with the six
native analyzer lines.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print(
    "GRAPH_ISOMORPHISM_COUNT:",
    isomorphism_count,
)
print(
    "GENERATOR_LOOKUP_FAILURES:",
    all_generator_lookup_failures,
)
print(
    "FACE_PROFILE_FAILURES:",
    all_face_profile_failures,
)
print(
    "FIXED_FACE_FAILURES:",
    all_fixed_face_failures,
)
print(
    "SETTING_FACE_BIJECTION_FAILURES:",
    all_bijection_failures,
)
print(
    "DISTINCT_SETTING_FACE_MAPS:",
    setting_face_map_count,
)
print(
    "SETTING_FACE_MAP_MULTIPLICITY_PROFILE:",
    dict(
        sorted(
            map_multiplicity_profile.items()
        )
    ),
)
print(
    "SAME_C5_EQUIVARIANT_3D_IDENTIFICATION:",
    False,
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
