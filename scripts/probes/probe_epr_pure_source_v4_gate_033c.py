#!/usr/bin/env python3

from collections import Counter, defaultdict
from pathlib import Path
import json

HERE = Path(__file__).resolve().parents[2]

P41 = (
    Path.home()
    / "dev/cori/research/mathematics"
    / "41-order-4-dodecahedral-residue"
)

ART = HERE / "artifacts/json"

AUT = (
    P41
    / "artifacts/json"
    / "native_g60_automorphism_deck_cache.v1.json"
)

A008 = (
    ART
    / "epr_native_source_character_sector_gate_008.v1.json"
)

A009 = (
    ART
    / "epr_native_face_boundary_incidence_map_009.v1.json"
)

A016 = (
    ART
    / "epr_literal_phase_decorated_incidence_table_016.v1.json"
)

A018 = (
    ART
    / "epr_reference_gauge_projective_line_map_018.v1.json"
)

A019 = (
    ART
    / "epr_joint_preparation_wedge_gauge_stress_019.v1.json"
)


def load(path):
    return json.loads(
        path.read_text(
            encoding="utf-8",
        )
    )


aut = load(AUT)
a008 = load(A008)
a009 = load(A009)
a016 = load(A016)
a018 = load(A018)
a019 = load(A019)

print("== 033C PURE SOURCE V4 GATE ==")

print()
print("== LOAD CHECKS ==")

for label, data in (
    ("A008", a008),
    ("A009", a009),
    ("A016", a016),
    ("A018", a018),
    ("A019", a019),
):
    print(
        label,
        "AUDIT_PASS:",
        data.get("audit_pass"),
    )


# ------------------------------------------------------------
# Reconstruct the same native a-kernel G1800 carrier used
# in Audits 008 and 019.
# ------------------------------------------------------------

a_index = int(
    a009[
        "native_a_action"
    ][
        "automorphism_index"
    ]
)

aut_rows = {
    int(
        row[
            "automorphism_index"
        ]
    ):
        row
    for row in aut[
        "measurements"
    ][
        "automorphism_rows"
    ]
}

a = tuple(
    int(v)
    for v in aut_rows[
        a_index
    ][
        "permutation"
    ]
)

assert len(a) == 60

assert all(
    a[a[i]] == i
    for i in range(60)
)

assert all(
    a[i] != i
    for i in range(60)
)


def diag_rep(u, v):
    return min(
        (u, v),
        (a[u], a[v]),
    )


g1800_reps = sorted({
    diag_rep(u, v)
    for u in range(60)
    for v in range(60)
})

g1800_index = {
    rep: i
    for i, rep in enumerate(
        g1800_reps
    )
}


def state_class(u, v):
    return g1800_index[
        diag_rep(u, v)
    ]


assert len(g1800_reps) == 1800

print()
print(
    "G1800_STATE_COUNT:",
    len(g1800_reps),
)


# ------------------------------------------------------------
# Reconstruct source-side V4 exactly as Audit 008.
# ------------------------------------------------------------

X = {}
T = {}

for sid, (u, v) in enumerate(
    g1800_reps
):
    X[sid] = state_class(
        v,
        u,
    )

    T[sid] = state_class(
        a[u],
        v,
    )

XT = {
    sid:
        X[
            T[sid]
        ]
    for sid in range(
        len(g1800_reps)
    )
}

assert all(
    X[X[i]] == i
    for i in X
)

assert all(
    T[T[i]] == i
    for i in T
)

assert all(
    X[T[i]]
    ==
    T[X[i]]
    for i in X
)

assert all(
    XT[XT[i]] == i
    for i in XT
)

actions = {
    "identity":
        {
            i: i
            for i in range(
                len(g1800_reps)
            )
        },

    "X":
        X,

    "tau":
        T,

    "Xtau":
        XT,
}


# ------------------------------------------------------------
# Reconstruct the one-wing line label from Audits 016/018.
# ------------------------------------------------------------

decorated_rows = a016[
    "decorated_incidence_domain"
][
    "rows"
]

ell_rows = a018[
    "ell"
][
    "rows"
]

ell_lookup = {}

for row in ell_rows:
    key = (
        int(
            row[
                "phase_state_id"
            ]
        ),
        int(
            row[
                "signed_block_index"
            ]
        ),
    )

    line = int(
        row[
            "reference_line_id"
        ]
    )

    if key in ell_lookup:
        assert (
            ell_lookup[key]
            == line
        )

    ell_lookup[key] = line

assert len(ell_lookup) == 48


def row_line(row):
    return ell_lookup[
        (
            int(
                row[
                    "phase_state_id"
                ]
            ),
            int(
                row[
                    "signed_block_index"
                ]
            ),
        )
    ]


rows_by_g60 = defaultdict(list)

for row in decorated_rows:
    u = int(
        row[
            "g60_state"
        ]
    )

    rows_by_g60[u].append(
        row
    )

one_wing_profile = Counter(
    len(
        rows_by_g60[u]
    )
    for u in range(60)
)

print()
print(
    "ONE_WING_ROWS_PER_G60_STATE:",
    dict(
        sorted(
            one_wing_profile.items()
        )
    ),
)

assert one_wing_profile == Counter({
    4: 60,
})


# ------------------------------------------------------------
# Count the raw 32 decorated presentations over each G1800
# state. Since Audit 019 proved branch invariance under the
# diagonal-a quotient, quotient counts are exactly half.
# ------------------------------------------------------------

raw_branch = defaultdict(
    Counter
)

for u in range(60):
    for v in range(60):
        sid = state_class(
            u,
            v,
        )

        for row_a in rows_by_g60[u]:
            line_a = row_line(
                row_a
            )

            for row_b in rows_by_g60[v]:
                line_b = row_line(
                    row_b
                )

                branch = (
                    "zero_wedge"
                    if line_a
                    == line_b
                    else
                    "nonzero_wedge"
                )

                raw_branch[
                    sid
                ][
                    branch
                ] += 1


raw_total_profile = Counter(
    sum(
        profile.values()
    )
    for profile in raw_branch.values()
)

print()
print(
    "RAW_PRESENTATIONS_PER_G1800_STATE:",
    dict(
        sorted(
            raw_total_profile.items()
        )
    ),
)

assert raw_total_profile == Counter({
    32: 1800,
})


quotient_branch = {}

for sid in range(
    len(g1800_reps)
):
    z = raw_branch[
        sid
    ][
        "zero_wedge"
    ]

    n = raw_branch[
        sid
    ][
        "nonzero_wedge"
    ]

    assert z % 2 == 0
    assert n % 2 == 0

    quotient_branch[sid] = (
        z // 2,
        n // 2,
    )


profile_hist = Counter(
    quotient_branch.values()
)

print()
print(
    "RECONSTRUCTED_QUOTIENT_BRANCH_PROFILE:"
)

for key, count in sorted(
    profile_hist.items()
):
    print(
        " ",
        key,
        "->",
        count,
    )


sealed_profile = {
    tuple(
        int(v)
        for v in key.strip(
            "()"
        ).split(",")
    ):
        int(count)
    for key, count in a019[
        "wedge_classifier"
    ][
        "per_G1800_state_branch_profile"
    ].items()
}

profile_matches_019 = (
    profile_hist
    ==
    Counter(
        sealed_profile
    )
)

print()
print(
    "PROFILE_MATCHES_AUDIT019:",
    profile_matches_019,
)

assert profile_matches_019


# ------------------------------------------------------------
# Identify the two extreme 16-state sets.
# ------------------------------------------------------------

pure_nonzero = {
    sid
    for sid, profile
    in quotient_branch.items()
    if profile == (
        0,
        16,
    )
}

pure_zero = {
    sid
    for sid, profile
    in quotient_branch.items()
    if profile == (
        16,
        0,
    )
}

print()
print(
    "PURE_NONZERO_STATE_COUNT:",
    len(
        pure_nonzero
    ),
)

print(
    "PURE_ZERO_STATE_COUNT:",
    len(
        pure_zero
    ),
)

assert len(
    pure_nonzero
) == 16

assert len(
    pure_zero
) == 16

print()
print(
    "PURE_NONZERO_STATE_IDS:",
    sorted(
        pure_nonzero
    ),
)

print(
    "PURE_ZERO_STATE_IDS:",
    sorted(
        pure_zero
    ),
)

print()
print(
    "PURE_NONZERO_REPRESENTATIVES:"
)

for sid in sorted(
    pure_nonzero
):
    print(
        " ",
        sid,
        g1800_reps[
            sid
        ],
    )

print()
print(
    "PURE_ZERO_REPRESENTATIVES:"
)

for sid in sorted(
    pure_zero
):
    print(
        " ",
        sid,
        g1800_reps[
            sid
        ],
    )


# ------------------------------------------------------------
# Test each native source action on the two extreme sets.
# ------------------------------------------------------------

def image_set(
    subset,
    action,
):
    return {
        action[sid]
        for sid in subset
    }


print()
print(
    "== SOURCE ACTION ON EXTREME SETS =="
)

action_summary = {}

for name, action in actions.items():
    nz_image = image_set(
        pure_nonzero,
        action,
    )

    z_image = image_set(
        pure_zero,
        action,
    )

    summary = {
        "NZ_to_NZ":
            nz_image
            == pure_nonzero,

        "NZ_to_ZERO":
            nz_image
            == pure_zero,

        "ZERO_to_ZERO":
            z_image
            == pure_zero,

        "ZERO_to_NZ":
            z_image
            == pure_nonzero,

        "NZ_image_size":
            len(
                nz_image
            ),

        "ZERO_image_size":
            len(
                z_image
            ),
    }

    action_summary[
        name
    ] = summary

    print(
        name,
        summary,
    )


# ------------------------------------------------------------
# V4 orbit decomposition.
# ------------------------------------------------------------

def v4_orbit(sid):
    return frozenset({
        sid,
        X[sid],
        T[sid],
        XT[sid],
    })


def subset_orbits(
    subset,
):
    seen = set()
    orbits = []

    for sid in sorted(
        subset
    ):
        if sid in seen:
            continue

        orb = v4_orbit(
            sid
        )

        orbits.append(
            orb
        )

        seen.update(
            orb
        )

    return orbits


nz_orbits = subset_orbits(
    pure_nonzero
)

z_orbits = subset_orbits(
    pure_zero
)

nz_closed = all(
    orb <= pure_nonzero
    for orb in nz_orbits
)

z_closed = all(
    orb <= pure_zero
    for orb in z_orbits
)

nz_orbit_profile = Counter(
    len(orb)
    for orb in nz_orbits
)

z_orbit_profile = Counter(
    len(orb)
    for orb in z_orbits
)

print()
print(
    "PURE_NONZERO_V4_CLOSED:",
    nz_closed,
)

print(
    "PURE_ZERO_V4_CLOSED:",
    z_closed,
)

print(
    "PURE_NONZERO_V4_ORBIT_PROFILE:",
    dict(
        sorted(
            nz_orbit_profile.items()
        )
    ),
)

print(
    "PURE_ZERO_V4_ORBIT_PROFILE:",
    dict(
        sorted(
            z_orbit_profile.items()
        )
    ),
)

print()
print(
    "PURE_NONZERO_ORBITS:"
)

for orb in nz_orbits:
    print(
        " ",
        sorted(
            orb
        ),
    )

print()
print(
    "PURE_ZERO_ORBITS:"
)

for orb in z_orbits:
    print(
        " ",
        sorted(
            orb
        ),
    )


# ------------------------------------------------------------
# Compare with full-carrier source V4 orbit profile.
# ------------------------------------------------------------

seen = set()
full_orbits = []

for sid in range(
    len(g1800_reps)
):
    if sid in seen:
        continue

    orb = v4_orbit(
        sid
    )

    full_orbits.append(
        orb
    )

    seen.update(
        orb
    )

full_profile = Counter(
    len(orb)
    for orb in full_orbits
)

sealed_full_profile = {
    int(k):
        int(v)
    for k, v in a008[
        "V4_orbit_profile"
    ].items()
}

print()
print(
    "FULL_V4_ORBIT_PROFILE:",
    dict(
        sorted(
            full_profile.items()
        )
    ),
)

print(
    "FULL_PROFILE_MATCHES_AUDIT008:",
    dict(
        sorted(
            full_profile.items()
        )
    )
    ==
    sealed_full_profile,
)


# ------------------------------------------------------------
# Classification.
# ------------------------------------------------------------

exchange_actions = [
    name
    for name, row
    in action_summary.items()
    if (
        row[
            "NZ_to_ZERO"
        ]
        and row[
            "ZERO_to_NZ"
        ]
    )
]

preserve_actions = [
    name
    for name, row
    in action_summary.items()
    if (
        row[
            "NZ_to_NZ"
        ]
        and row[
            "ZERO_to_ZERO"
        ]
    )
]

print()
print(
    "== 033C CLASSIFICATION =="
)

print(
    "PURE_NONZERO_INTRINSIC_PREDICATE:",
    "all 16 decorated quotient presentations have nonzero wedge",
)

print(
    "PURE_NONZERO_SETTING_INDEPENDENT:",
    True,
)

print(
    "PURE_NONZERO_GAUGE_COVARIANT:",
    True,
)

print(
    "PURE_NONZERO_V4_INVARIANT:",
    nz_closed,
)

print(
    "PURE_ZERO_V4_INVARIANT:",
    z_closed,
)

print(
    "SOURCE_ACTIONS_PRESERVING_BOTH:",
    preserve_actions,
)

print(
    "SOURCE_ACTIONS_EXCHANGING_EXTREMES:",
    exchange_actions,
)

if nz_closed:
    print(
        "SOURCE_CANDIDATE_STATUS:",
        "native_V4_invariant_pure_nonzero_support_candidate",
    )
else:
    print(
        "SOURCE_CANDIDATE_STATUS:",
        "favorable_subset_not_closed_under_native_source_V4",
    )

print()
print(
    "BOUNDARY:"
)

print(
    "  This probe tests native support/invariance only."
)

print(
    "  It does not derive a preparation mechanism or probability law."
)

print(
    "  A V4-invariant pure-nonzero set is a source-support candidate,"
)

print(
    "  not yet a proof that the physical/native source selects it."
)
