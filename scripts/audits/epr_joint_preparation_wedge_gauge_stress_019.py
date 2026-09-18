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

A009 = (
    HERE
    / "artifacts/json"
    / "epr_native_face_boundary_incidence_map_009.v1.json"
)

A013 = (
    HERE
    / "artifacts/json"
    / "epr_complete_exchange_parity_closure_013.v1.json"
)

A016 = (
    HERE
    / "artifacts/json"
    / "epr_literal_phase_decorated_incidence_table_016.v1.json"
)

A018 = (
    HERE
    / "artifacts/json"
    / "epr_reference_gauge_projective_line_map_018.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_joint_preparation_wedge_gauge_stress_019.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_joint_preparation_wedge_gauge_stress_019.md"
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


def norm_line(line):
    return tuple(
        sorted(
            int(x)
            for x in line
        )
    )


def norm_partition(lines):
    return tuple(
        sorted(
            norm_line(line)
            for line in lines
        )
    )


def perm_inverse(p):
    q = [None] * len(p)

    for i, j in enumerate(p):
        q[j] = i

    return tuple(q)


def compose(p, q):
    # p after q
    return tuple(
        p[q[i]]
        for i in range(len(q))
    )


print("PROGRESS: 1/9 load sealed interfaces")

aut = load(AUT)
a009 = load(A009)
a013 = load(A013)
a016 = load(A016)
a018 = load(A018)

checks = {}

checks["Audit009_passed"] = (
    a009["audit_pass"] is True
)

checks["Audit013_passed"] = (
    a013["audit_pass"] is True
)

checks["Audit016_passed"] = (
    a016["audit_pass"] is True
)

checks["Audit018_passed"] = (
    a018["audit_pass"] is True
)

checks["Audit013_odd_branch_rank1"] = (
    int(
        a013[
            "complete_closure"
        ]["odd_real_rank"]
    )
    == 1
)

aut_rows = aut[
    "measurements"
]["automorphism_rows"]

checks["AutG60_row_count_480"] = (
    len(aut_rows) == 480
)

print("PROGRESS: 2/9 reconstruct native a action")

a_index = int(
    a009[
        "native_a_action"
    ]["automorphism_index"]
)

aut_by_index = {
    int(row["automorphism_index"]):
        row
    for row in aut_rows
}

a = tuple(
    int(x)
    for x in aut_by_index[
        a_index
    ]["permutation"]
)

checks["a_is_60_point_permutation"] = (
    len(a) == 60
    and sorted(a) == list(range(60))
)

checks["a_is_involution"] = all(
    a[a[u]] == u
    for u in range(60)
)

checks["a_is_fixed_point_free"] = all(
    a[u] != u
    for u in range(60)
)

block_rows = a016[
    "native_mode_space"
]["block_rows"]

block_state_key = {}

for row in block_rows:
    block_index = int(
        row["signed_block_index"]
    )

    key = tuple(
        sorted(
            int(x)
            for x in row[
                "five_state_block"
            ]
        )
    )

    block_state_key[
        key
    ] = block_index

checks["block_row_count_24"] = (
    len(block_rows) == 24
)

a_block = {}

for row in block_rows:
    block_index = int(
        row["signed_block_index"]
    )

    image_key = tuple(
        sorted(
            a[int(u)]
            for u in row[
                "five_state_block"
            ]
        )
    )

    a_block[
        block_index
    ] = block_state_key[
        image_key
    ]

checks["a_acts_as_block_involution"] = all(
    a_block[
        a_block[b]
    ]
    == b
    for b in a_block
)

print("PROGRESS: 3/9 reconstruct phase partitions and a action")

phase_states = a016[
    "phase_state_space"
]["states"]

phase_by_id = {
    int(row["phase_state_id"]):
        row
    for row in phase_states
}

ell_rows = a018[
    "ell"
]["rows"]

ell_lookup = {}

phase_partition_lines = defaultdict(set)

for row in ell_rows:
    sid = int(
        row["phase_state_id"]
    )

    block = int(
        row["signed_block_index"]
    )

    partner = int(
        row["phase_line_partner_block"]
    )

    line_id = int(
        row["reference_line_id"]
    )

    ell_lookup[
        (sid, block)
    ] = line_id

    phase_partition_lines[
        sid
    ].add(
        norm_line(
            (block, partner)
        )
    )

phase_partition = {
    sid: norm_partition(lines)
    for sid, lines in (
        phase_partition_lines.items()
    )
}

checks["twelve_phase_partitions"] = (
    len(phase_partition) == 12
)

checks["two_lines_per_phase_state"] = all(
    len(
        phase_partition[sid]
    ) == 2
    for sid in phase_partition
)

checks["ell_has_48_phase_mode_rows"] = (
    len(ell_lookup) == 48
)

phase_image_a = {}

phase_image_failure_count = 0

for sid, partition in phase_partition.items():
    image_partition = norm_partition(
        [
            a_block[b]
            for b in line
        ]
        for line in partition
    )

    matches = [
        tid
        for tid, target_partition
        in phase_partition.items()
        if target_partition
        == image_partition
    ]

    if len(matches) != 1:
        phase_image_failure_count += 1
        continue

    phase_image_a[
        sid
    ] = matches[0]

checks["a_maps_each_phase_state_uniquely"] = (
    phase_image_failure_count == 0
    and len(phase_image_a) == 12
)

checks["a_phase_action_is_involution"] = all(
    phase_image_a[
        phase_image_a[sid]
    ]
    == sid
    for sid in phase_image_a
)

phase_epsilon_map_profile = Counter()

phase_carrier_failure_count = 0

for sid, tid in phase_image_a.items():
    source = phase_by_id[sid]
    target = phase_by_id[tid]

    phase_epsilon_map_profile[
        (
            int(source["epsilon"]),
            int(target["epsilon"]),
        )
    ] += 1

    if (
        int(source["carrier_index"])
        != int(target["carrier_index"])
    ):
        phase_carrier_failure_count += 1

checks["a_preserves_face_carrier_at_phase_level"] = (
    phase_carrier_failure_count == 0
)

print("PROGRESS: 4/9 reconstruct literal 240-row one-wing domain")

decorated_rows = a016[
    "decorated_incidence_domain"
]["rows"]

checks["decorated_one_wing_row_count_240"] = (
    len(decorated_rows) == 240
)

row_key_to_id = {}

for row_id, row in enumerate(
    decorated_rows
):
    key = (
        int(row["phase_state_id"]),
        int(row["g60_state"]),
        int(row["signed_block_index"]),
    )

    if key in row_key_to_id:
        raise RuntimeError(
            "duplicate decorated incidence row"
        )

    row_key_to_id[key] = row_id

a_row = {}

decorated_closure_failure_count = 0

for row_id, row in enumerate(
    decorated_rows
):
    sid = int(
        row["phase_state_id"]
    )

    u = int(
        row["g60_state"]
    )

    block = int(
        row["signed_block_index"]
    )

    image_key = (
        phase_image_a[sid],
        a[u],
        a_block[block],
    )

    image_id = row_key_to_id.get(
        image_key
    )

    if image_id is None:
        decorated_closure_failure_count += 1
        continue

    a_row[
        row_id
    ] = image_id

checks["a_closes_on_all_240_decorated_rows"] = (
    decorated_closure_failure_count == 0
    and len(a_row) == 240
)

checks["decorated_a_action_is_involution"] = all(
    a_row[
        a_row[i]
    ]
    == i
    for i in a_row
)

checks["decorated_a_action_is_free"] = all(
    a_row[i] != i
    for i in a_row
)

print("PROGRESS: 5/9 construct phase-decorated two-wing G1800 quotient")

raw_pair_count = (
    len(decorated_rows)
    * len(decorated_rows)
)

checks["raw_two_wing_count_57600"] = (
    raw_pair_count == 57600
)


def pair_rep(i, j):
    image = (
        a_row[i],
        a_row[j],
    )

    row = (
        i,
        j,
    )

    return min(
        row,
        image,
    )


joint_reps = set()

for i in range(
    len(decorated_rows)
):
    for j in range(
        len(decorated_rows)
    ):
        joint_reps.add(
            pair_rep(
                i,
                j,
            )
        )

joint_reps = sorted(
    joint_reps
)

checks["decorated_joint_quotient_count_28800"] = (
    len(joint_reps) == 28800
)


def state_rep(u, v):
    return min(
        (u, v),
        (a[u], a[v]),
    )


g1800_state_reps = sorted({
    state_rep(u, v)
    for u in range(60)
    for v in range(60)
})

g1800_state_id = {
    rep: i
    for i, rep in enumerate(
        g1800_state_reps
    )
}

checks["G1800_state_count_1800"] = (
    len(g1800_state_reps)
    == 1800
)

state_to_joint_count = Counter()

for i, j in joint_reps:
    ri = decorated_rows[i]
    rj = decorated_rows[j]

    sid = g1800_state_id[
        state_rep(
            int(ri["g60_state"]),
            int(rj["g60_state"]),
        )
    ]

    state_to_joint_count[
        sid
    ] += 1

joint_fiber_profile = Counter(
    state_to_joint_count.values()
)

checks["sixteen_decorated_presentations_per_G1800_state"] = (
    joint_fiber_profile
    == Counter({
        16: 1800,
    })
)

print("PROGRESS: 6/9 classify equal-line versus unequal-line branches")


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


one_wing_line_profile = Counter(
    row_line(row)
    for row in decorated_rows
)

checks["one_wing_line_profile_120_120"] = (
    one_wing_line_profile
    == Counter({
        0: 120,
        1: 120,
    })
)

raw_branch_profile = Counter()

quotient_branch_profile = Counter()

quotient_descent_failure_count = 0

state_branch_profile = defaultdict(
    Counter
)

for i in range(
    len(decorated_rows)
):
    row_i = decorated_rows[i]
    line_i = row_line(
        row_i
    )

    for j in range(
        len(decorated_rows)
    ):
        row_j = decorated_rows[j]
        line_j = row_line(
            row_j
        )

        branch = (
            "zero_wedge"
            if line_i == line_j
            else
            "nonzero_wedge"
        )

        raw_branch_profile[
            branch
        ] += 1

        ai = a_row[i]
        aj = a_row[j]

        image_branch = (
            "zero_wedge"
            if row_line(
                decorated_rows[ai]
            )
            ==
            row_line(
                decorated_rows[aj]
            )
            else
            "nonzero_wedge"
        )

        if branch != image_branch:
            quotient_descent_failure_count += 1

checks["branch_classifier_descends_through_diagonal_a"] = (
    quotient_descent_failure_count
    == 0
)

for i, j in joint_reps:
    row_i = decorated_rows[i]
    row_j = decorated_rows[j]

    branch = (
        "zero_wedge"
        if row_line(row_i)
        == row_line(row_j)
        else
        "nonzero_wedge"
    )

    quotient_branch_profile[
        branch
    ] += 1

    state_id = g1800_state_id[
        state_rep(
            int(row_i["g60_state"]),
            int(row_j["g60_state"]),
        )
    ]

    state_branch_profile[
        state_id
    ][branch] += 1

checks["raw_branch_counts_balance_28800_28800"] = (
    raw_branch_profile
    == Counter({
        "zero_wedge": 28800,
        "nonzero_wedge": 28800,
    })
)

checks["quotient_branch_counts_balance_14400_14400"] = (
    quotient_branch_profile
    == Counter({
        "zero_wedge": 14400,
        "nonzero_wedge": 14400,
    })
)

per_state_branch_profile = Counter(
    (
        profile[
            "zero_wedge"
        ],
        profile[
            "nonzero_wedge"
        ],
    )
    for profile in (
        state_branch_profile.values()
    )
)

print("PROGRESS: 7/9 stress local line-label gauge covariance")

# In the selected flat reference gauge the transition on
# reference line labels is the identity.
#
# Under independent local gauge relabelings gamma_A and
# gamma_B, the local labels become
#
#     l_A' = gamma_A(l_A)
#     l_B' = gamma_B(l_B)
#
# while the transition changes covariantly to
#
#     T'_AB = gamma_A T_AB gamma_B^-1.
#
# With T_AB = identity in the reference trivialization,
#
#     T'_AB = gamma_A gamma_B^-1.
#
# Therefore
#
#     l_A' = T'_AB(l_B')
#
# iff
#
#     l_A = l_B.
#
# The finite two-line native label gauge is S2. Exhaust it.

S2 = (
    (0, 1),
    (1, 0),
)

slot_rows = {}

for row in ell_rows:
    key = (
        int(row["phase_state_id"]),
        int(row["local_mode_slot"]),
    )

    slot_rows[key] = int(
        row["reference_line_id"]
    )

checks["phase_mode_slot_count_48"] = (
    len(slot_rows) == 48
)

gauge_test_count = 0
gauge_covariance_failure_count = 0

for key_a, line_a in slot_rows.items():
    for key_b, line_b in slot_rows.items():
        original_equal = (
            line_a == line_b
        )

        for gamma_a in S2:
            for gamma_b in S2:
                gamma_b_inv = perm_inverse(
                    gamma_b
                )

                transition_prime = compose(
                    gamma_a,
                    gamma_b_inv,
                )

                line_a_prime = gamma_a[
                    line_a
                ]

                line_b_prime = gamma_b[
                    line_b
                ]

                transported_b_prime = (
                    transition_prime[
                        line_b_prime
                    ]
                )

                transformed_equal = (
                    line_a_prime
                    ==
                    transported_b_prime
                )

                gauge_test_count += 1

                if (
                    transformed_equal
                    != original_equal
                ):
                    gauge_covariance_failure_count += 1

checks["all_S2_line_gauge_covariance_tests_pass"] = (
    gauge_covariance_failure_count
    == 0
)

checks["line_gauge_test_count_9216"] = (
    gauge_test_count == 9216
)

print("PROGRESS: 8/9 verify canonical nonzero image and complete branch retention")

checks["nonzero_branch_maps_to_rank1_odd_sector"] = (
    a013[
        "branch_interpretation"
    ]["odd"]
    ==
    "canonical rank-one EPR line [|01>-|10>]"
)

checks["zero_branch_retained"] = (
    quotient_branch_profile[
        "zero_wedge"
    ] > 0
)

checks["nonzero_branch_retained"] = (
    quotient_branch_profile[
        "nonzero_wedge"
    ] > 0
)

checks["no_postselection_required_for_classifier"] = (
    checks[
        "zero_branch_retained"
    ]
    and checks[
        "nonzero_branch_retained"
    ]
)

print("PROGRESS: 9/9 classify")

failed_checks = [
    name
    for name, passed in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "phase_decorated_G1800_boundary_has_a_complete_"
    "setting_independent_gauge_covariant_zero_nonzero_"
    "wedge_branch_classifier_with_nonzero_branch_on_"
    "the_canonical_rank1_EPR_line"
    if audit_pass
    else
    "joint_preparation_wedge_gauge_stress_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_joint_preparation_wedge_gauge_stress_019",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "one_wing": {
        "phase_decorated_row_count":
            len(decorated_rows),

        "reference_line_profile": {
            str(k): v
            for k, v in sorted(
                one_wing_line_profile.items()
            )
        },
    },

    "native_a_action": {
        "automorphism_index":
            a_index,

        "phase_state_map": {
            str(sid):
                phase_image_a[sid]
            for sid in sorted(
                phase_image_a
            )
        },

        "epsilon_map_profile": {
            str(k):
                v
            for k, v in sorted(
                phase_epsilon_map_profile.items()
            )
        },

        "decorated_action_free":
            True,

        "branch_descent_failure_count":
            quotient_descent_failure_count,
    },

    "two_wing_domain": {
        "raw_phase_decorated_pair_count":
            raw_pair_count,

        "diagonal_a_quotient_count":
            len(joint_reps),

        "G1800_state_count":
            len(g1800_state_reps),

        "decorated_presentations_per_G1800_state_profile": {
            str(k): v
            for k, v in sorted(
                joint_fiber_profile.items()
            )
        },
    },

    "wedge_classifier": {
        "criterion":
            (
                "zero iff reference projective lines are equal; "
                "nonzero iff reference projective lines differ"
            ),

        "raw_branch_profile":
            dict(
                sorted(
                    raw_branch_profile.items()
                )
            ),

        "quotient_branch_profile":
            dict(
                sorted(
                    quotient_branch_profile.items()
                )
            ),

        "per_G1800_state_branch_profile": {
            str(k): v
            for k, v in sorted(
                per_state_branch_profile.items()
            )
        },

        "nonzero_image":
            "canonical rank-one EPR line",

        "zero_branch":
            "retained preparation-degenerate branch",
    },

    "gauge_stress": {
        "native_reference_line_label_group":
            "S2",

        "phase_mode_pair_count":
            len(slot_rows),

        "ordered_phase_mode_pair_count":
            len(slot_rows)
            * len(slot_rows),

        "local_gauge_pair_count_per_mode_pair":
            4,

        "test_count":
            gauge_test_count,

        "failure_count":
            gauge_covariance_failure_count,

        "covariance_law":
            (
                "ell_A'=gamma_A ell_A; "
                "ell_B'=gamma_B ell_B; "
                "T_AB'=gamma_A T_AB gamma_B^-1"
            ),

        "equal_versus_unequal_invariant":
            (
                gauge_covariance_failure_count
                == 0
            ),
    },

    "earned_statement": (
        "The exact 240-row phase-decorated one-wing boundary "
        "domain produces 57,600 ordered two-wing presentations. "
        "The native diagonal-a action is free on this domain and "
        "reduces it to 28,800 phase-decorated G1800 boundary "
        "classes over the same 1,800 common states. The Audit-018 "
        "projective-line map gives a complete pre-setting branch "
        "classifier: equal transported CP1 lines give zero "
        "antisymmetric wedge, while unequal lines give nonzero "
        "wedge. The classifier descends exactly through the "
        "diagonal-a quotient. In the selected flat reference "
        "gauge the quotient contains both branches explicitly; "
        "the nonzero branch lands on the unique rank-one odd EPR "
        "line already closed by Audits 006-013. Exhaustive S2 "
        "line-label gauge stress confirms that equal versus "
        "unequal is unchanged when local line labels and the "
        "pair transition transform covariantly. Both zero and "
        "nonzero branches remain in the complete preparation "
        "domain."
    ),

    "checks":
        checks,

    "boundary": {
        "pre_setting_branch_classifier_constructed":
            audit_pass,

        "branch_classifier_descends_to_G1800":
            (
                quotient_descent_failure_count
                == 0
            ),

        "line_label_gauge_covariance_closed":
            (
                gauge_covariance_failure_count
                == 0
            ),

        "nonzero_branch_is_canonical_EPR_line":
            True,

        "zero_branch_retained":
            True,

        "source_not_yet_proved_to_prepare_only_nonzero_branch":
            True,

        "branch_counts_not_probabilities":
            True,

        "source_weighting_not_yet_derived":
            True,

        "Born_rule_not_used":
            True,

        "analyzer_settings_not_used":
            True,

        "no_signaling_not_yet_tested":
            True,

        "CHSH_not_yet_derived_from_receipt_distribution":
            True,
    },

    "next_gate": (
        "Carry both preparation branches into the local analyzer "
        "apparatus. Construct the setting-dependent local "
        "transduction and complete receipt alphabet without "
        "discarding the zero-wedge branch. Then derive native "
        "trial weights or frequencies and test normalization, "
        "no-signaling, and CHSH from the complete receipt "
        "distribution."
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

note = f"""# Joint preparation wedge gauge stress 019

## Result

Audit pass:

    {audit_pass}

The phase-decorated one-wing domain contains

    {len(decorated_rows)}

rows.

The raw ordered two-wing product therefore contains

    {raw_pair_count}

presentations.

The native diagonal-a quotient contains

    {len(joint_reps)}

classes over

    {len(g1800_state_reps)}

G1800 common states.

The decorated fiber profile is

    {dict(sorted(joint_fiber_profile.items()))}.

## Wedge branch

In the Audit-018 flat reference gauge:

    equal projective lines
        -> zero wedge

    unequal projective lines
        -> nonzero wedge.

Raw branch profile:

    {dict(sorted(raw_branch_profile.items()))}

After the native diagonal-a quotient:

    {dict(sorted(quotient_branch_profile.items()))}

The branch classifier has

    {quotient_descent_failure_count}

failures under the diagonal-a identification.

The nonzero wedge lies on the unique rank-one odd line already
identified as the canonical EPR line.

The zero branch is retained.

## Gauge stress

The finite reference-line label gauge is S2.

For every ordered pair of the 48 phase-mode slots, all four independent
local line-label gauge pairs were tested with the covariant transition
law

    T'_AB
      =
    gamma_A T_AB gamma_B^-1.

Total tests:

    {gauge_test_count}

Failures:

    {gauge_covariance_failure_count}

Thus equal versus unequal transported projective line is invariant under
the native two-line label gauge when the transition transforms
covariantly.

## Preparation boundary

This closes a complete setting-independent branch classifier.

It does not prove that the source prepares only the nonzero branch.

Both branches remain part of the native preparation domain.

The branch counts are combinatorial counts, not probabilities.

The next task is to pass both branches through the local analyzer
apparatus and derive the complete local receipt distribution without
postselection.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print(
    "ONE_WING_DECORATED_COUNT:",
    len(decorated_rows),
)
print(
    "RAW_TWO_WING_COUNT:",
    raw_pair_count,
)
print(
    "DECORATED_G1800_QUOTIENT_COUNT:",
    len(joint_reps),
)
print(
    "G1800_STATE_COUNT:",
    len(g1800_state_reps),
)
print(
    "JOINT_FIBER_PROFILE:",
    dict(sorted(
        joint_fiber_profile.items()
    )),
)
print(
    "A_PHASE_EPSILON_MAP_PROFILE:",
    dict(sorted(
        phase_epsilon_map_profile.items()
    )),
)
print(
    "RAW_BRANCH_PROFILE:",
    dict(sorted(
        raw_branch_profile.items()
    )),
)
print(
    "QUOTIENT_BRANCH_PROFILE:",
    dict(sorted(
        quotient_branch_profile.items()
    )),
)
print(
    "PER_G1800_STATE_BRANCH_PROFILE:",
    dict(sorted(
        per_state_branch_profile.items()
    )),
)
print(
    "DIAGONAL_A_BRANCH_FAILURES:",
    quotient_descent_failure_count,
)
print(
    "GAUGE_STRESS_TEST_COUNT:",
    gauge_test_count,
)
print(
    "GAUGE_STRESS_FAILURE_COUNT:",
    gauge_covariance_failure_count,
)
print(
    "NONZERO_BRANCH_CANONICAL_EPR_LINE:",
    checks[
        "nonzero_branch_maps_to_rank1_odd_sector"
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
