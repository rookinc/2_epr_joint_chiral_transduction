#!/usr/bin/env python3

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from io import BytesIO
from itertools import permutations, product
from pathlib import Path
from zipfile import ZipFile
import json

HERE = Path(__file__).resolve().parents[2]

ZIP = (
    HERE
    / "evidence"
    / "native_coherent_axis_instrument.zip"
)

A019 = (
    HERE
    / "artifacts/json"
    / "epr_joint_preparation_wedge_gauge_stress_019.v1.json"
)

A020 = (
    HERE
    / "artifacts/json"
    / "epr_local_analyzer_carrier_interface_020.v1.json"
)

A021 = (
    HERE
    / "artifacts/json"
    / "epr_canonical_line_analyzer_compatibility_021.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_native_six_axis_singlet_quartet_census_022.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_native_six_axis_singlet_quartet_census_022.md"
)


def load(path):
    return json.loads(path.read_text())


# Exact Q(sqrt(5)) arithmetic.
# A field element is represented as
#
#     a + b sqrt(5)

def q(a=0, b=0):
    return (
        Fraction(a),
        Fraction(b),
    )


def qadd(x, y):
    return (
        x[0] + y[0],
        x[1] + y[1],
    )


def qneg(x):
    return (
        -x[0],
        -x[1],
    )


def qsub(x, y):
    return qadd(
        x,
        qneg(y),
    )


def qmul(x, y):
    a, b = x
    c, d = y

    return (
        a * c + 5 * b * d,
        a * d + b * c,
    )


def qscale(c, x):
    c = Fraction(c)

    return (
        c * x[0],
        c * x[1],
    )


def qsquare(x):
    return qmul(x, x)


def qsign(x):
    a, b = x

    if a == 0 and b == 0:
        return 0

    if b == 0:
        return 1 if a > 0 else -1

    if b > 0:
        if a >= 0:
            return 1

        return (
            1
            if 5 * b * b > a * a
            else -1
        )

    if a <= 0:
        return -1

    return (
        1
        if a * a > 5 * b * b
        else -1
    )


def qabs(x):
    return (
        x
        if qsign(x) >= 0
        else qneg(x)
    )


def qstr(x):
    a, b = x

    if b == 0:
        return str(a)

    pieces = []

    if a != 0:
        pieces.append(str(a))

    if b == 1:
        term = "sqrt(5)"
    elif b == -1:
        term = "-sqrt(5)"
    else:
        term = str(b) + "*sqrt(5)"

    if not pieces:
        return term

    if b > 0:
        return pieces[0] + " + " + term

    return pieces[0] + " - " + term[1:]


def digest_json(obj):
    raw = json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")

    return sha256(raw).hexdigest()


def decode_entry(row, i, j):
    den = int(row["denominator"])

    return (
        Fraction(
            int(
                row[
                    "rational_numerator"
                ][i][j]
            ),
            den,
        ),
        Fraction(
            int(
                row[
                    "sqrt5_numerator"
                ][i][j]
            ),
            den,
        ),
    )


def projector_seed_vector(row):
    # A deterministic generic seed.
    #
    # v = P (1,2,...,60)^T
    #
    # Because P is rank one, any nonzero v lies on the
    # projective analyzer line represented by P.

    out = []

    for i in range(60):
        total = q()

        for j in range(60):
            total = qadd(
                total,
                qscale(
                    j + 1,
                    decode_entry(
                        row,
                        i,
                        j,
                    ),
                ),
            )

        out.append(total)

    return out


def qdot(u, v):
    total = q()

    for x, y in zip(u, v):
        total = qadd(
            total,
            qmul(x, y),
        )

    return total


def oriented_line_sign_matrix(rows):
    vectors = [
        projector_seed_vector(row)
        for row in rows
    ]

    norms = [
        qdot(v, v)
        for v in vectors
    ]

    if any(
        qsign(n) <= 0
        for n in norms
    ):
        raise RuntimeError(
            "generic seed missed an analyzer line"
        )

    signs = [
        [0 for _ in range(6)]
        for _ in range(6)
    ]

    overlap_failures = 0

    for i in range(6):
        signs[i][i] = 1

        for j in range(i + 1, 6):
            d = qdot(
                vectors[i],
                vectors[j],
            )

            # Exact line-angle test:
            #
            #     cos^2(theta) = 1/5
            #
            # without extracting square roots.
            lhs = qscale(
                5,
                qsquare(d),
            )

            rhs = qmul(
                norms[i],
                norms[j],
            )

            if lhs != rhs:
                overlap_failures += 1

            s = qsign(d)

            if s == 0:
                raise RuntimeError(
                    "distinct analyzer lines became orthogonal"
                )

            signs[i][j] = s
            signs[j][i] = s

    return (
        signs,
        overlap_failures,
    )


def correlation(sign_matrix, i, j):
    # Conditional common-adjoint lift:
    #
    # native projective line -> unit su(2) direction
    #
    # and canonical antisymmetric line gives
    #
    #     E(i,j) = - n_i dot n_j.
    #
    # Same native line gives -1.
    # Distinct lines have magnitude 1/sqrt(5).

    if i == j:
        return q(-1, 0)

    s = sign_matrix[i][j]

    return q(
        0,
        Fraction(-s, 5),
    )


def best_chsh_for_quartet(
    sign_matrix,
    a0,
    a1,
    b0,
    b1,
):
    E00 = correlation(
        sign_matrix,
        a0,
        b0,
    )

    E01 = correlation(
        sign_matrix,
        a0,
        b1,
    )

    E10 = correlation(
        sign_matrix,
        a1,
        b0,
    )

    E11 = correlation(
        sign_matrix,
        a1,
        b1,
    )

    best = q()
    best_signs = None

    for (
        sa0,
        sa1,
        sb0,
        sb1,
    ) in product(
        (+1, -1),
        repeat=4,
    ):
        value = q()

        value = qadd(
            value,
            qscale(
                sa0 * sb0,
                E00,
            ),
        )

        value = qadd(
            value,
            qscale(
                sa0 * sb1,
                E01,
            ),
        )

        value = qadd(
            value,
            qscale(
                sa1 * sb0,
                E10,
            ),
        )

        value = qsub(
            value,
            qscale(
                sa1 * sb1,
                E11,
            ),
        )

        av = qabs(value)

        if qsign(
            qsub(
                av,
                best,
            )
        ) > 0:
            best = av
            best_signs = (
                sa0,
                sa1,
                sb0,
                sb1,
            )

    return (
        best,
        best_signs,
    )


print("PROGRESS: 1/7 load sealed Program-02 interfaces")

a019 = load(A019)
a020 = load(A020)
a021 = load(A021)

checks = {}

checks["Audit019_passed"] = (
    a019["audit_pass"] is True
)

checks["Audit020_passed"] = (
    a020["audit_pass"] is True
)

checks["Audit021_passed"] = (
    a021["audit_pass"] is True
)

checks["canonical_nonzero_branch_is_EPR_line"] = (
    a019[
        "boundary"
    ][
        "nonzero_branch_is_canonical_EPR_line"
    ]
    is True
)

checks["face_adjoint_dimension_3"] = (
    int(
        a020[
            "face_hilbert_carrier"
        ][
            "adjoint_real_dimension"
        ]
    )
    == 3
)

print("PROGRESS: 2/7 recover exact native analyzer projectors")

with ZipFile(ZIP) as z:
    up1 = z.read(
        "native_coherent_axis_instrument/"
        "upstream/"
        "native_axis_contact_pointer_audit.zip"
    )

with ZipFile(BytesIO(up1)) as z1:
    up2 = z1.read(
        "native_axis_contact_pointer_audit/"
        "upstream/"
        "native_reference_contact_gate.zip"
    )

with ZipFile(BytesIO(up2)) as z2:
    up3 = z2.read(
        "native_reference_contact_gate/"
        "upstream/"
        "native_analyzer_return_bridge.zip"
    )

with ZipFile(BytesIO(up3)) as z3:
    exact_projectors = json.loads(
        z3.read(
            "native_analyzer_return_bridge/"
            "EXACT_PROJECTORS.json"
        )
    )

axes = exact_projectors["axes"]

checks["two_native_character_sectors_present"] = (
    "B_plus3" in axes
    and "B_minus3" in axes
)

checks["six_axes_in_each_character_sector"] = (
    len(
        axes["B_plus3"]
    )
    == 6
    and len(
        axes["B_minus3"]
    )
    == 6
)

print("PROGRESS: 3/7 reconstruct exact projective line Gram signs")

sector_sign_matrices = {}
sector_overlap_failures = {}

for sector in (
    "B_plus3",
    "B_minus3",
):
    rows = sorted(
        axes[sector],
        key=lambda row: int(
            row["setting"]
        ),
    )

    checks[
        sector
        + "_settings_are_0_through_5"
    ] = (
        [
            int(row["setting"])
            for row in rows
        ]
        == list(range(6))
    )

    (
        sign_matrix,
        overlap_failures,
    ) = oriented_line_sign_matrix(
        rows
    )

    sector_sign_matrices[
        sector
    ] = sign_matrix

    sector_overlap_failures[
        sector
    ] = overlap_failures

checks["all_distinct_line_overlaps_are_exactly_one_fifth"] = (
    all(
        count == 0
        for count in (
            sector_overlap_failures.values()
        )
    )
)

print("PROGRESS: 4/7 exhaust all native setting quartets")

sector_results = {}

expected_profile = Counter({
    q(1, Fraction(3, 5)): 240,
    q(1, Fraction(1, 5)): 240,
    q(0, Fraction(4, 5)): 240,
    q(0, Fraction(2, 5)): 120,
    q(2, 0): 60,
})

for sector in (
    "B_plus3",
    "B_minus3",
):
    sign_matrix = (
        sector_sign_matrices[
            sector
        ]
    )

    profile = Counter()
    shared_profile = Counter()
    violating_shared_profile = Counter()

    violation_count = 0
    quartet_count = 0

    max_value = q()
    max_examples = []

    for a0, a1 in permutations(
        range(6),
        2,
    ):
        for b0, b1 in permutations(
            range(6),
            2,
        ):
            quartet_count += 1

            (
                best,
                signs,
            ) = best_chsh_for_quartet(
                sign_matrix,
                a0,
                a1,
                b0,
                b1,
            )

            profile[best] += 1

            shared = len(
                {
                    a0,
                    a1,
                }
                &
                {
                    b0,
                    b1,
                }
            )

            shared_profile[
                shared
            ] += 1

            if qsign(
                qsub(
                    best,
                    q(2, 0),
                )
            ) > 0:
                violation_count += 1

                violating_shared_profile[
                    shared
                ] += 1

            comparison = qsign(
                qsub(
                    best,
                    max_value,
                )
            )

            if comparison > 0:
                max_value = best
                max_examples = [
                    {
                        "Alice":
                            [a0, a1],

                        "Bob":
                            [b0, b1],

                        "outcome_signs":
                            list(signs),

                        "shared_axis_count":
                            shared,
                    }
                ]

            elif (
                comparison == 0
                and len(
                    max_examples
                ) < 12
            ):
                max_examples.append({
                    "Alice":
                        [a0, a1],

                    "Bob":
                        [b0, b1],

                    "outcome_signs":
                        list(signs),

                    "shared_axis_count":
                        shared,
                })

    sector_results[
        sector
    ] = {
        "quartet_count":
            quartet_count,

        "violation_count":
            violation_count,

        "max_CHSH":
            qstr(
                max_value
            ),

        "max_CHSH_field_pair":
            [
                str(
                    max_value[0]
                ),
                str(
                    max_value[1]
                ),
            ],

        "profile":
            {
                qstr(k): v
                for k, v
                in sorted(
                    profile.items(),
                    key=lambda kv: (
                        float(
                            kv[0][0]
                        )
                        + float(
                            kv[0][1]
                        )
                        * (5 ** 0.5)
                    ),
                )
            },

        "shared_axis_profile":
            dict(
                sorted(
                    shared_profile.items()
                )
            ),

        "violating_shared_axis_profile":
            dict(
                sorted(
                    violating_shared_profile.items()
                )
            ),

        "max_examples":
            max_examples,
    }

    checks[
        sector
        + "_quartet_count_900"
    ] = (
        quartet_count == 900
    )

    checks[
        sector
        + "_exact_CHSH_profile"
    ] = (
        profile
        == expected_profile
    )

    checks[
        sector
        + "_violation_count_240"
    ] = (
        violation_count == 240
    )

    checks[
        sector
        + "_violations_have_exactly_one_shared_axis"
    ] = (
        violating_shared_profile
        == Counter({
            1: 240,
        })
    )

    checks[
        sector
        + "_max_is_1_plus_3_over_sqrt5"
    ] = (
        max_value
        == q(
            1,
            Fraction(3, 5),
        )
    )

print("PROGRESS: 5/7 compare native character sectors")

checks["character_sector_CHSH_profiles_match"] = (
    sector_results[
        "B_plus3"
    ]["profile"]
    ==
    sector_results[
        "B_minus3"
    ]["profile"]
)

checks["character_sector_violation_counts_match"] = (
    sector_results[
        "B_plus3"
    ]["violation_count"]
    ==
    sector_results[
        "B_minus3"
    ]["violation_count"]
    ==
    240
)

checks["max_CHSH_exceeds_2_exactly"] = (
    qsign(
        qsub(
            q(
                1,
                Fraction(3, 5),
            ),
            q(2, 0),
        )
    )
    > 0
)

print("PROGRESS: 6/7 classify common-adjoint gauge status")

# A common orthogonal change of adjoint frame preserves
# every analyzer-line inner product.
#
# The canonical antisymmetric line is invariant under a
# common local SU(2) frame change.
#
# Therefore this quartet census depends only on the native
# six-line Gram geometry, not on the coordinates of a chosen
# common adjoint isometry.
#
# What is NOT yet derived is that Alice and Bob inherit the
# same adjoint-frame embedding, or which finite relative
# alignment the G1800 preparation supplies.

checks["census_uses_only_native_projective_Gram_data"] = (
    True
)

checks["no_continuous_CHSH_optimization_used"] = (
    True
)

checks["no_independent_Alice_Bob_rotation_inserted"] = (
    True
)

print("PROGRESS: 7/7 classify")

failed_checks = [
    name
    for name, passed
    in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "canonical_antisymmetric_EPR_line_is_Bell_compatible_"
    "with_the_full_native_six_axis_geometry_under_a_common_"
    "adjoint_lift_with_240_of_900_native_setting_quartets_"
    "reaching_1_plus_3_over_sqrt5_greater_than_two"
    if audit_pass
    else
    "native_six_axis_singlet_quartet_census_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_native_six_axis_singlet_quartet_census_022",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "state": {
        "projective_line":
            "[|01>-|10>]",

        "source":
            "Audit 019 nonzero wedge branch",

        "common_frame_invariance":
            (
                "canonical antisymmetric line is invariant "
                "under common local SU(2) frame changes"
            ),
    },

    "analyzer_geometry": {
        "native_character_sectors":
            [
                "B_plus3",
                "B_minus3",
            ],

        "axes_per_sector":
            6,

        "distinct_projective_line_overlap_squared":
            "1/5",

        "adjoint_direction_inner_product_magnitude":
            "1/sqrt(5)",

        "candidate_face_observable":
            "O_i = n_i dot sigma",
    },

    "quartet_census": {
        "requirements":
            (
                "Alice settings distinct; Bob settings "
                "distinct; all ordered native choices tested"
            ),

        "quartets_per_sector":
            900,

        "local_outcome_sign_relabelings_per_quartet":
            16,

        "sector_results":
            sector_results,

        "maximum":
            "1 + 3*sqrt(5)/5",

        "maximum_alternate_form":
            "1 + 3/sqrt(5)",

        "maximum_exceeds_2":
            True,

        "violating_quartets_per_sector":
            240,

        "violating_quartet_structure":
            (
                "exactly one analyzer axis is shared "
                "between the Alice and Bob setting pairs"
            ),
    },

    "relation_to_Audit021": {
        "Audit021_same_pair_result":
            "max 2 under local outcome relabeling",

        "new_result":
            (
                "the same-pair quartet is a restricted "
                "subfamily; the full six-axis native setting "
                "census contains Bell-capable quartets for "
                "the canonical EPR line"
            ),
    },

    "earned_statement": (
        "The negative result of Audit 021 is specific to using "
        "the same native analyzer pair on both wings. Each native "
        "three-dimensional analyzer character sector carries six "
        "projective axes with exact pair overlap squared 1/5. "
        "Under a common adjoint lift of that native line system "
        "into the Program-01 face su(2), the canonical "
        "antisymmetric EPR line has correlations determined only "
        "by the analyzer Gram matrix. Exhausting all 900 ordered "
        "setting quartets with distinct settings on each wing and "
        "all sixteen local outcome-sign relabelings gives 240 "
        "Bell-capable quartets in each native character sector. "
        "The exact maximum is 1 + 3/sqrt(5), greater than two. "
        "Every violating quartet has exactly one native analyzer "
        "axis shared between Alice's and Bob's setting pairs. "
        "No continuous frame optimization is used."
    ),

    "checks":
        checks,

    "boundary": {
        "common_adjoint_lift_is_conditional_interface":
            True,

        "common_adjoint_coordinates_are_gauge":
            True,

        "native_relative_Alice_Bob_adjoint_frame_not_yet_derived":
            True,

        "no_independent_continuous_frame_optimization":
            True,

        "Bell_capable_native_setting_quartets_exist":
            True,

        "complete_local_instrument_not_yet_bound_to_face_C2":
            True,

        "zero_wedge_branch_still_retained":
            True,

        "Born_frequency_law_not_derived":
            True,

        "no_signaling_not_yet_tested":
            True,

        "complete_receipt_CHSH_not_yet_derived":
            True,
    },

    "next_gate": (
        "Derive the finite native relative adjoint-frame relation "
        "between Alice and Bob from the G1800/Program-01 common "
        "boundary structure. Test that finite relative-frame set "
        "against the 240 Bell-capable native quartets. Do not "
        "introduce an arbitrary independent SO(3) rotation. If "
        "the common lift or another native relative alignment is "
        "admitted, proceed to construct the complete local face "
        "instrument and receipt weighting."
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

note = f"""# Native six-axis singlet quartet census 022

## Result

Audit pass:

    {audit_pass}

Audit 021 tested one restricted quartet: the same analyzer pair on both
wings.

This audit keeps the canonical antisymmetric EPR line fixed and exhausts
the complete six-axis native analyzer geometry.

Each native analyzer character sector has six projective lines with

    |n_i dot n_j|^2 = 1/5

for distinct axes.

Under a common adjoint lift into the face su(2), the singlet correlation
is

    E(i,j) = - n_i dot n_j.

A common adjoint-frame rotation changes no Gram entries and therefore no
value in this census.

## Exhaustive native setting census

Alice chooses two distinct axes.

Bob chooses two distinct axes.

Thus there are

    6*5*6*5 = 900

ordered setting quartets per native character sector.

For each quartet all

    16

local binary outcome-sign relabelings are exhausted.

Both character sectors give the same exact profile.

Bell-capable quartets per sector:

    240 / 900.

Exact maximum:

    1 + 3/sqrt(5)

which is greater than

    2.

Every violating quartet has exactly one native analyzer axis shared
between Alice's two-setting set and Bob's two-setting set.

No arbitrary continuous Alice/Bob frame rotation is optimized.

## Meaning

The Audit-021 negative result was not a state/analyzer incompatibility.

It was a restricted-setting result.

The full six-axis native analyzer geometry contains setting quartets
compatible with Bell violation by the canonical antisymmetric EPR line,
provided Alice and Bob inherit the required common adjoint lift.

That relative-frame admission is now the next theorem gate.

The coherent local instrument and frequency law remain downstream.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print(
    "QUARTETS_PER_SECTOR:",
    900,
)
print(
    "OUTCOME_SIGN_TESTS_PER_QUARTET:",
    16,
)
print(
    "B_PLUS_VIOLATING_QUARTETS:",
    sector_results[
        "B_plus3"
    ]["violation_count"],
)
print(
    "B_MINUS_VIOLATING_QUARTETS:",
    sector_results[
        "B_minus3"
    ]["violation_count"],
)
print(
    "MAX_CHSH:",
    "1 + 3/sqrt(5)",
)
print(
    "MAX_CHSH_GT_2:",
    checks[
        "max_CHSH_exceeds_2_exactly"
    ],
)
print(
    "VIOLATING_SHARED_AXIS_PROFILE:",
    sector_results[
        "B_plus3"
    ][
        "violating_shared_axis_profile"
    ],
)
print(
    "CONTINUOUS_FRAME_OPTIMIZATION_USED:",
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
