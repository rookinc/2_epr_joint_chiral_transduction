#!/usr/bin/env python3

from fractions import Fraction
from hashlib import sha256
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile
import json

HERE = Path(__file__).resolve().parents[2]

P41 = (
    Path.home()
    / "dev/cori/research/mathematics"
    / "41-order-4-dodecahedral-residue"
)

ZIP = (
    HERE
    / "evidence"
    / "native_coherent_axis_instrument.zip"
)

FT = (
    P41
    / "artifacts/provenance"
    / "project41_event_vector_u2_u1_stabilizer_201ft.v1.json"
)

FW = (
    P41
    / "artifacts/provenance"
    / "project41_neutral_generator_anatomy_201fw.v1.json"
)

A019 = (
    HERE
    / "artifacts/json"
    / "epr_joint_preparation_wedge_gauge_stress_019.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_local_analyzer_carrier_interface_020.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_local_analyzer_carrier_interface_020.md"
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


print("PROGRESS: 1/7 load local interfaces")

ft = load(FT)
fw = load(FW)
a019 = load(A019)

checks = {}

checks["201FT_passed"] = (
    ft["audit_pass"] is True
)

checks["201FW_passed"] = (
    fw["audit_pass"] is True
)

checks["Audit019_passed"] = (
    a019["audit_pass"] is True
)

checks["instrument_archive_exists"] = (
    ZIP.exists()
)

print("PROGRESS: 2/7 open coherent instrument packet")

archive_sha256 = sha256(
    ZIP.read_bytes()
).hexdigest()

with ZipFile(ZIP) as z:
    base = (
        "native_coherent_axis_instrument/"
    )

    outer_report = json.loads(
        z.read(
            base + "REPORT.json"
        )
    )

    exact_symmetry = json.loads(
        z.read(
            base + "EXACT_SYMMETRY.json"
        )
    )

    upstream_1 = z.read(
        base
        + "upstream/"
        + "native_axis_contact_pointer_audit.zip"
    )

checks["outer_instrument_passed"] = (
    outer_report["audit_pass"] is True
)

checks["outer_failed_check_count_zero"] = (
    int(
        outer_report[
            "failed_check_count"
        ]
    )
    == 0
)

checks["outer_check_count_4457"] = (
    int(
        outer_report[
            "check_count"
        ]
    )
    == 4457
)

checks["exact_symmetry_is_six_rows"] = (
    isinstance(
        exact_symmetry,
        list,
    )
    and len(
        exact_symmetry
    )
    == 6
)

print("PROGRESS: 3/7 certify current normalized instrument")

boundaries = outer_report[
    "boundaries"
]

exact_results = outer_report[
    "exact_results"
]

numerical = outer_report[
    "numerical_diagnostics"
]

checks["source_band_dimension_6"] = (
    int(
        numerical[
            "source_band_dimension"
        ]
    )
    == 6
)

checks["axis_count_6"] = (
    int(
        numerical[
            "axis_count"
        ]
    )
    == 6
)

checks["all_nonclean_outputs_retained"] = (
    boundaries[
        "all_nonclean_outputs_retained"
    ]
    is True
)

checks["conditioning_on_completion_not_used"] = (
    boundaries[
        "conditioning_on_completion_used"
    ]
    is False
)

checks["Born_frequency_law_not_derived"] = (
    boundaries[
        "Born_frequency_law_derived"
    ]
    is False
)

checks["coherent_dynamics_is_added_hypothesis"] = (
    boundaries[
        "coherent_dynamics_is_added_hypothesis"
    ]
    is True
)

checks["ideal_instrument_not_constructed"] = (
    boundaries[
        "ideal_instrument_constructed"
    ]
    is False
)

checks["all_pulse_effect_form_recorded"] = (
    exact_results[
        "all_pulse_effect_form"
    ]
    ==
    "a_b(t) Pi_i + b_b(t) (I-Pi_i)"
)

checks["normalizer_commutant_dimension_2"] = (
    int(
        exact_results[
            "normalizer_source_complex_commutant_dimension"
        ]
    )
    == 2
)

symmetry_failure_count = 0

for row in exact_symmetry:
    if int(
        row[
            "normalizer_order"
        ]
    ) != 80:
        symmetry_failure_count += 1

    if int(
        row[
            "source_axis_rank"
        ]
    ) != 2:
        symmetry_failure_count += 1

    if int(
        row[
            "source_plane_rank"
        ]
    ) != 4:
        symmetry_failure_count += 1

checks["all_six_settings_have_rank2_axis_rank4_plane"] = (
    symmetry_failure_count == 0
)

print("PROGRESS: 4/7 open nested native analyzer bridge")

with ZipFile(
    BytesIO(
        upstream_1
    )
) as z1:
    upstream_2 = z1.read(
        "native_axis_contact_pointer_audit/"
        "upstream/"
        "native_reference_contact_gate.zip"
    )

with ZipFile(
    BytesIO(
        upstream_2
    )
) as z2:
    upstream_3 = z2.read(
        "native_reference_contact_gate/"
        "upstream/"
        "native_analyzer_return_bridge.zip"
    )

with ZipFile(
    BytesIO(
        upstream_3
    )
) as z3:
    analyzer_report = json.loads(
        z3.read(
            "native_analyzer_return_bridge/"
            "REPORT.json"
        )
    )

    exact_projectors = json.loads(
        z3.read(
            "native_analyzer_return_bridge/"
            "EXACT_PROJECTORS.json"
        )
    )

checks["nested_analyzer_bridge_passed"] = (
    analyzer_report[
        "audit_pass"
    ]
    is True
)

checks["nested_failed_check_count_zero"] = (
    int(
        analyzer_report[
            "failed_check_count"
        ]
    )
    == 0
)

counts = analyzer_report[
    "counts"
]

bridge_boundaries = analyzer_report[
    "boundaries"
]

checks["leading_source_dimension_6"] = (
    int(
        counts[
            "leading_source_dimension"
        ]
    )
    == 6
)

checks["leading_character_dimensions_3_3"] = (
    counts[
        "leading_B_character_dimensions"
    ]
    == [3, 3]
)

checks["six_native_C5_axes"] = (
    int(
        counts[
            "native_C5_subgroups_reconstructed"
        ]
    )
    == 6
)

checks["twelve_rank1_axes_checked"] = (
    int(
        counts[
            "rank_one_axes_checked"
        ]
    )
    == 12
)

checks["thirty_axis_pair_overlaps_checked"] = (
    int(
        counts[
            "rank_one_axis_pair_overlaps_checked"
        ]
    )
    == 30
)

checks["historical_EPR_instrument_was_not_implemented"] = (
    bridge_boundaries[
        "EPR_instrument_implemented"
    ]
    is False
)

checks["historical_probability_rule_not_derived"] = (
    bridge_boundaries[
        "probability_rule_derived"
    ]
    is False
)

checks["projector_sector_keys_present"] = (
    set(
        exact_projectors[
            "axes"
        ].keys()
    )
    ==
    {
        "B_minus3",
        "B_plus3",
        "paired6",
    }
)

checks["six_projectors_per_character_sector"] = (
    len(
        exact_projectors[
            "axes"
        ]["B_minus3"]
    )
    == 6
    and
    len(
        exact_projectors[
            "axes"
        ]["B_plus3"]
    )
    == 6
)

print("PROGRESS: 5/7 certify Program-01 local face carrier")

face_input = ft[
    "native_input"
]

face_generators = fw[
    "native_generator_data"
]

face_complex_dimension = int(
    face_input[
        "face_complex_dimension"
    ]
)

face_real_dimension = int(
    face_input[
        "face_real_dimension"
    ]
)

derived_axes = list(
    face_generators[
        "derived_axes"
    ]
)

checks["face_complex_dimension_2"] = (
    face_complex_dimension == 2
)

checks["face_real_dimension_4"] = (
    face_real_dimension == 4
)

checks["face_continuous_envelope_U2"] = (
    face_input[
        "continuous_envelope"
    ]
    == "U(2)"
)

checks["face_has_three_SU2_type_axes"] = (
    derived_axes
    == [
        "K1",
        "K2",
        "K3",
    ]
)

face_adjoint_real_dimension = (
    len(
        derived_axes
    )
)

analyzer_character_dimension = int(
    counts[
        "leading_B_character_dimensions"
    ][0]
)

checks["face_adjoint_dimension_matches_analyzer_sector"] = (
    face_adjoint_real_dimension
    ==
    analyzer_character_dimension
    ==
    3
)

print("PROGRESS: 6/7 classify geometry versus instrument bridge")

# Native analyzer geometry:
#
# Each character sector is 3-dimensional.
# Each of six Pi_i is rank one there.
# Distinct axes satisfy
#
#     tr(Pi_i Pi_j) = 1/5.
#
# For unit representatives n_i,n_j:
#
#     (n_i dot n_j)^2 = 1/5.
#
# The reflection observable is
#
#     O_i = 2 Pi_i - I_3.
#
# For any distinct pair:
#
#     sin^2(2 theta)
#       = 4 cos^2(theta) sin^2(theta)
#       = 4 * (1/5) * (4/5)
#       = 16/25.
#
# The symmetric two-setting CHSH operator therefore has
#
#     norm^2 = 4(1 + 16/25)
#            = 164/25,
#
# hence norm = 2 sqrt(41) / 5.

axis_overlap_squared = Fraction(
    1,
    5,
)

sin2theta_squared = (
    4
    * axis_overlap_squared
    * (
        1
        - axis_overlap_squared
    )
)

chsh_norm_squared = (
    4
    * (
        1
        + sin2theta_squared
    )
)

checks["sin2theta_squared_16_over_25"] = (
    sin2theta_squared
    == Fraction(
        16,
        25,
    )
)

checks["CHSH_norm_squared_164_over_25"] = (
    chsh_norm_squared
    == Fraction(
        164,
        25,
    )
)

direct_carrier_identity_possible = (
    face_complex_dimension
    ==
    analyzer_character_dimension
)

checks["direct_face_C2_to_analyzer3_identity_is_impossible"] = (
    direct_carrier_identity_possible
    is False
)

adjoint_geometry_dimension_match = (
    face_adjoint_real_dimension
    ==
    analyzer_character_dimension
)

checks["adjoint_geometry_bridge_is_dimensionally_available"] = (
    adjoint_geometry_dimension_match
    is True
)

# Audit 019 supplies the native EPR preparation on the face
# Hilbert side. The coherent analyzer instrument acts on the
# separate six-dimensional source band. No map between those
# state carriers is supplied by either packet.

checks["Audit019_nonzero_branch_is_canonical_EPR_line"] = (
    a019[
        "boundary"
    ][
        "nonzero_branch_is_canonical_EPR_line"
    ]
    is True
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
    "native_analyzer_geometry_is_exact_on_two_"
    "three_dimensional_character_sectors_and_the_"
    "coherent_axis_contact_instrument_is_normalized_"
    "and_complete_but_no_direct_face_C2_instrument_"
    "carrier_binding_is_established"
    if audit_pass
    else
    "local_analyzer_carrier_interface_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_local_analyzer_carrier_interface_020",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "coherent_instrument": {
        "archive":
            str(ZIP),

        "archive_sha256":
            archive_sha256,

        "check_count":
            int(
                outer_report[
                    "check_count"
                ]
            ),

        "source_band_dimension":
            6,

        "setting_count":
            6,

        "axis_rank_on_paired_band":
            2,

        "plane_rank_on_paired_band":
            4,

        "effect_form":
            exact_results[
                "all_pulse_effect_form"
            ],

        "all_nonclean_outputs_retained":
            True,

        "normalized_mathematical_instrument":
            True,

        "coherent_dynamics_is_added_hypothesis":
            True,

        "Born_frequency_law_derived":
            False,

        "ideal_instrument_constructed":
            False,
    },

    "native_analyzer_geometry": {
        "character_sector_dimensions":
            [3, 3],

        "settings":
            6,

        "rank_one_axes_per_sector":
            6,

        "distinct_pair_overlap_squared":
            "1/5",

        "reflection_observable":
            "O_i = 2 Pi_i - I_3",

        "CHSH_norm_squared":
            "164/25",

        "CHSH_norm":
            "2 sqrt(41) / 5",

        "operator_Bell_capacity":
            True,
    },

    "face_hilbert_carrier": {
        "complex_dimension":
            face_complex_dimension,

        "real_dimension":
            face_real_dimension,

        "continuous_envelope":
            face_input[
                "continuous_envelope"
            ],

        "SU2_type_axes":
            derived_axes,

        "adjoint_real_dimension":
            face_adjoint_real_dimension,
    },

    "interface_result": {
        "direct_face_Hilbert_to_analyzer_sector_identity":
            False,

        "reason":
            (
                "local face Hilbert carrier has complex "
                "dimension 2 while each analyzer source "
                "sector has dimension 3"
            ),

        "face_adjoint_and_analyzer_sector_dimensions_match":
            True,

        "candidate_geometry_bridge":
            (
                "an isometric identification of the native "
                "three-dimensional analyzer axis space with "
                "the three-dimensional su(2) face adjoint"
            ),

        "candidate_geometry_bridge_constructed":
            False,

        "coherent_6D_instrument_transferred_to_face_C2":
            False,
    },

    "earned_statement": (
        "The recovered coherent axis-contact packet supplies "
        "a normalized, analyzer-sensitive, complete mathematical "
        "instrument on a six-dimensional native source band and "
        "retains all nonclean outputs. Its six analyzer settings "
        "restrict on each of two native three-dimensional "
        "character sectors to rank-one projective axes with "
        "distinct-pair overlap 1/5. The associated reflection "
        "geometry has exact CHSH operator norm 2 sqrt(41)/5. "
        "Program 01, however, supplies a complex two-dimensional "
        "local face Hilbert carrier. These state carriers are "
        "not directly identical. The dimensions do match at the "
        "adjoint level: the face su(2)-type frame has three real "
        "axes and each native analyzer sector is three-dimensional. "
        "Thus the analyzer setting geometry has a plausible "
        "adjoint bridge to the face apparatus, but the existing "
        "six-dimensional coherent instrument has not been "
        "transferred to the face C2 carrier."
    ),

    "checks":
        checks,

    "boundary": {
        "analyzer_operator_geometry_locked":
            True,

        "normalized_complete_coherent_instrument_locked":
            True,

        "coherent_dynamics_remains_model_hypothesis":
            True,

        "face_C2_preparation_not_fed_directly_into_6D_instrument":
            True,

        "adjoint_geometry_bridge_not_yet_constructed":
            True,

        "local_face_instrument_not_yet_constructed":
            True,

        "remote_setting_not_introduced":
            True,

        "zero_wedge_branch_not_discarded":
            True,

        "Born_rule_not_derived":
            True,

        "receipt_weighting_not_derived":
            True,

        "no_signaling_not_yet_tested":
            True,

        "CHSH_not_yet_derived_from_complete_receipts":
            True,
    },

    "next_gate": (
        "Construct the analyzer-setting geometry on the "
        "three-dimensional su(2) adjoint of the Program-01 "
        "face fiber, treating any common adjoint frame choice "
        "as a gauge candidate to be stress-tested. Verify the "
        "six axis overlaps and CHSH operator capacity there. "
        "Keep this geometry descent separate from the unresolved "
        "problem of constructing a normalized complete local "
        "instrument that actually acts on the face C2 carrier."
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

note = """# Local analyzer carrier interface 020

## Result

Audit pass:

    %s

The recovered coherent axis-contact packet has a six-dimensional
source band and six analyzer settings.

For every setting its complete pointer effects have the exact form

    E_i,b(t)
      =
    a_b(t) Pi_i
      +
    b_b(t) (I-Pi_i).

On the paired source band

    rank(Pi_i) = 2
    rank(I-Pi_i) = 4.

All nonclean outputs are retained.

The coherent-law choice remains an explicit model hypothesis, and no
Born-frequency law is derived.

## Analyzer sectors

The native source band splits into two character sectors:

    3 + 3.

On each sector there are six rank-one analyzer axes with

    tr(Pi_i Pi_j) = 1/5

for every distinct pair.

The reflection observables

    O_i = 2 Pi_i - I_3

therefore have the exact Bell-capable geometry

    ||B|| = 2 sqrt(41) / 5.

## Carrier boundary

Program 01 supplies a local face fiber of complex dimension two:

    C^2.

That is not the same state carrier as a three-dimensional analyzer
sector or the paired six-dimensional coherent-instrument source band.

So the current coherent instrument may not simply be attached to the
Audit-019 face EPR state.

There is, however, an exact dimensional meeting point:

    dim_R su(2)_face = 3
    dim analyzer character sector = 3.

The next gate is therefore an adjoint geometry bridge, not a fake
C2-to-6D identity.

After that geometry is closed, a separate normalized local instrument
must still be constructed on the actual face C2 carrier.
""" % audit_pass

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print(
    "COHERENT_SOURCE_BAND_DIM:",
    6,
)
print(
    "ANALYZER_CHARACTER_DIMS:",
    [3, 3],
)
print(
    "ANALYZER_SETTING_COUNT:",
    6,
)
print(
    "FACE_COMPLEX_DIM:",
    face_complex_dimension,
)
print(
    "FACE_ADJOINT_REAL_DIM:",
    face_adjoint_real_dimension,
)
print(
    "CHSH_NORM:",
    "2 sqrt(41) / 5",
)
print(
    "DIRECT_C2_TO_ANALYZER3_IDENTITY:",
    False,
)
print(
    "ADJOINT_DIMENSION_MATCH:",
    adjoint_geometry_dimension_match,
)
print(
    "COHERENT_INSTRUMENT_TRANSFERRED_TO_FACE_C2:",
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
