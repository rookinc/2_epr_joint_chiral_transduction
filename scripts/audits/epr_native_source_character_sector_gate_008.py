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

A002 = (
    HERE
    / "artifacts/json"
    / "epr_native_deck_face_chiral_transduction_002.v1.json"
)

A006 = (
    HERE
    / "artifacts/json"
    / "epr_canonical_joint_projector_identity_006.v1.json"
)

A007 = (
    HERE
    / "artifacts/json"
    / "epr_canonical_joint_projector_gauge_descent_007.v1.json"
)

AUT = (
    P41
    / "artifacts/json"
    / "native_g60_automorphism_deck_cache.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_native_source_character_sector_gate_008.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_native_source_character_sector_gate_008.md"
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


print("PROGRESS: 1/7 load sealed interfaces")

a002 = load(A002)
a006 = load(A006)
a007 = load(A007)
aut = load(AUT)

checks = {}

checks["audit_002_passed"] = (
    a002["audit_pass"] is True
)

checks["audit_006_passed"] = (
    a006["audit_pass"] is True
)

checks["audit_007_passed"] = (
    a007["audit_pass"] is True
)

checks["canonical_boundary_projector_rank_one"] = (
    a006["projectors"]["complex_rank"] == 1
)

print("PROGRESS: 2/7 reconstruct native a-kernel G1800")

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

checks["native_a_is_involution"] = all(
    a_perm[a_perm[i]] == i
    for i in range(60)
)

checks["native_a_is_fixed_point_free"] = all(
    a_perm[i] != i
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
    return joint_index[
        diag_rep(u, v)
    ]


checks["native_g1800_class_count_1800"] = (
    len(joint_reps) == 1800
)

print("PROGRESS: 3/7 construct source symmetry V4")

X = {}
T = {}

for cid, (u, v) in enumerate(joint_reps):
    X[cid] = joint_class(v, u)

    T[cid] = joint_class(
        a_perm[u],
        v,
    )


def compose(p, q, i):
    return p[q[i]]


checks["X_is_involution"] = all(
    X[X[i]] == i
    for i in X
)

checks["tau_is_involution"] = all(
    T[T[i]] == i
    for i in T
)

checks["X_and_tau_commute"] = all(
    X[T[i]] == T[X[i]]
    for i in X
)

XT = {
    i: X[T[i]]
    for i in X
}

checks["Xtau_is_involution"] = all(
    XT[XT[i]] == i
    for i in XT
)

print("PROGRESS: 4/7 count fixed sets and V4 orbits")

fixed = {
    "identity":
        len(joint_reps),

    "X":
        sum(
            X[i] == i
            for i in X
        ),

    "tau":
        sum(
            T[i] == i
            for i in T
        ),

    "Xtau":
        sum(
            XT[i] == i
            for i in XT
        ),
}

checks["identity_fixed_count_1800"] = (
    fixed["identity"] == 1800
)

checks["X_fixed_count_60"] = (
    fixed["X"] == 60
)

checks["tau_fixed_count_0"] = (
    fixed["tau"] == 0
)

checks["Xtau_fixed_count_0"] = (
    fixed["Xtau"] == 0
)

seen = set()
orbit_sizes = []

for i in range(len(joint_reps)):
    if i in seen:
        continue

    orbit = {
        i,
        X[i],
        T[i],
        XT[i],
    }

    seen.update(orbit)
    orbit_sizes.append(len(orbit))

orbit_profile = Counter(
    orbit_sizes
)

checks["V4_orbit_profile_exact"] = (
    orbit_profile
    == Counter({
        2: 30,
        4: 435,
    })
)

checks["V4_orbit_count_465"] = (
    len(orbit_sizes) == 465
)

print("PROGRESS: 5/7 compute exact character multiplicities")

# Characters are labeled by their eigenvalues under
#
#   X    -> sx
#   tau  -> st
#
# where sx,st are +/-1.
#
# For a finite permutation representation,
#
# m_chi = (1/4) sum_g chi(g) Fix(g).

character_rows = []

for sx in (1, -1):
    for st in (1, -1):
        numerator = (
            fixed["identity"]
            + sx * fixed["X"]
            + st * fixed["tau"]
            + sx * st * fixed["Xtau"]
        )

        checks[
            f"character_integrality_X{sx}_tau{st}"
        ] = (
            numerator % 4 == 0
        )

        multiplicity = (
            numerator // 4
        )

        character_rows.append({
            "X_character": sx,
            "tau_character": st,
            "multiplicity": multiplicity,
        })

character_map = {
    (
        row["X_character"],
        row["tau_character"],
    ):
        row["multiplicity"]
    for row in character_rows
}

checks["X_even_tau_even_dim_465"] = (
    character_map[(1, 1)] == 465
)

checks["X_even_tau_odd_dim_465"] = (
    character_map[(1, -1)] == 465
)

checks["X_odd_tau_even_dim_435"] = (
    character_map[(-1, 1)] == 435
)

checks["X_odd_tau_odd_dim_435"] = (
    character_map[(-1, -1)] == 435
)

checks["character_dimensions_sum_1800"] = (
    sum(
        row["multiplicity"]
        for row in character_rows
    )
    == 1800
)

print("PROGRESS: 6/7 test raw native preparation sufficiency")

X_odd_dimension = (
    character_map[(-1, 1)]
    + character_map[(-1, -1)]
)

checks["native_X_odd_sector_dimension_870"] = (
    X_odd_dimension == 870
)

smallest_character_dimension = min(
    row["multiplicity"]
    for row in character_rows
)

checks["smallest_native_character_sector_dimension_435"] = (
    smallest_character_dimension == 435
)

checks["no_native_source_character_sector_is_rank_one"] = all(
    row["multiplicity"] != 1
    for row in character_rows
)

checks["raw_source_symmetry_does_not_prepare_boundary_line"] = (
    X_odd_dimension > 1
    and smallest_character_dimension > 1
    and checks[
        "canonical_boundary_projector_rank_one"
    ]
)

print("PROGRESS: 7/7 classify")

failed_checks = [
    name
    for name, passed in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "native_g1800_exchange_and_chiral_torsor_characters_"
    "define_large_source_sectors_not_the_rank_one_EPR_line"
    if audit_pass
    else
    "native_source_character_sector_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_native_source_character_sector_gate_008",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "native_source_carrier": {
        "carrier":
            "(G60 x G60)/diag(a)",

        "state_count":
            len(joint_reps),

        "symmetry_group":
            "<X,tau_a> ~= V4",

        "X":
            "factor exchange [u,v] -> [v,u]",

        "tau_a":
            (
                "relative lift exchange "
                "[u,v] -> [a u,v] = [u,a v]"
            ),
    },

    "fixed_counts": fixed,

    "V4_orbit_profile": {
        str(k): v
        for k, v in sorted(
            orbit_profile.items()
        )
    },

    "character_sectors":
        character_rows,

    "source_sector_dimensions": {
        "factor_exchange_odd_total":
            X_odd_dimension,

        "smallest_full_V4_character_sector":
            smallest_character_dimension,

        "boundary_canonical_line_complex_rank":
            1,
    },

    "earned_statement": (
        "The actual 1800-state common carrier has a canonical "
        "source-side V4 generated by factor exchange X and the "
        "native chiral relative-lift involution tau_a. Its "
        "permutation representation decomposes into character "
        "multiplicities 465,465,435,435. Factor-exchange "
        "oddness therefore selects an 870-dimensional source "
        "sector, not the rank-one boundary EPR line. Even "
        "specifying both X and tau_a characters leaves a "
        "435-dimensional sector. Thus discrete source parity "
        "alone cannot execute the canonical preparation found "
        "in Audits 006-007."
    ),

    "checks":
        checks,

    "boundary": {
        "native_source_symmetry_sector_exists":
            True,

        "native_exchange_odd_sector_is_setting_independent":
            True,

        "source_sector_not_identified_with_rank_one_boundary_line":
            True,

        "no_boundary_transduction_map_assumed":
            True,

        "group_algebra_sector_not_declared_physical_state":
            True,

        "no_Born_rule":
            True,

        "no_probability_table":
            True,

        "no_CHSH_claim":
            True,
    },

    "next_gate": (
        "Construct the explicit common-domain boundary "
        "transduction map from native G1800 preparation data "
        "into the two Program-01 local face doublets. The next "
        "question is whether that map sends an intrinsically "
        "specified native source sector onto the canonical "
        "rank-one joint line. Do not add another parity label; "
        "the raw source character census proves that discrete "
        "V4 character selection alone is insufficient."
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

note = f"""# Native source character-sector gate 008

## Result

Audit pass:

    {audit_pass}

The native a-kernel common carrier carries two commuting involutions:

    X
        factor exchange

    tau_a
        native relative-lift exchange.

Together they generate a source-side V4.

Their fixed counts are

    identity: {fixed["identity"]}
    X:        {fixed["X"]}
    tau_a:    {fixed["tau"]}
    X tau_a:  {fixed["Xtau"]}.

The resulting V4 orbit profile is

    {dict(sorted(orbit_profile.items()))}.

## Character census

The four source character-sector dimensions are

"""

for row in character_rows:
    note += (
        "    X={:+d}, tau={:+d}: {}\n".format(
            row["X_character"],
            row["tau_character"],
            row["multiplicity"],
        )
    )

note += f"""
Thus the full factor-exchange odd source sector has dimension

    {X_odd_dimension}.

Even after also fixing tau_a character, the smallest relevant source
sector has dimension

    {smallest_character_dimension}.

The canonical two-face EPR projector from Audit 006 has complex rank

    1.

## Consequence

Source-side exchange oddness is native and setting independent.

It is not enough to prepare the canonical boundary line.

Nor is adding the tau_a character enough.

The missing object is therefore not another discrete source parity.
It is the actual common-domain boundary transduction from native G1800
preparation data into the two functional Program-01 faces.

That map must explain how a large native source sector can present or
close onto the canonical rank-one joint line.

## Boundary

The source character sectors are representation-theoretic objects.

They are not declared physical states or probabilities.

No Born rule, trial weighting, no-signaling theorem, or CHSH value is
introduced here.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print("FIXED_COUNTS:", fixed)
print(
    "V4_ORBIT_PROFILE:",
    dict(sorted(orbit_profile.items())),
)
print(
    "CHARACTER_SECTORS:",
    character_rows,
)
print(
    "X_ODD_SOURCE_DIM:",
    X_odd_dimension,
)
print(
    "SMALLEST_V4_CHARACTER_DIM:",
    smallest_character_dimension,
)
print(
    "BOUNDARY_CANONICAL_LINE_RANK:",
    1,
)
print(
    "DISCRETE_SOURCE_CHARACTER_SELECTION_SUFFICIENT:",
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
