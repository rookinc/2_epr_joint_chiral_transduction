#!/usr/bin/env python3

from hashlib import sha256
from pathlib import Path
import json
import math

HERE = Path(__file__).resolve().parents[2]

A022 = (
    HERE
    / "artifacts/json"
    / "epr_native_six_axis_singlet_quartet_census_022.v1.json"
)

A027 = (
    HERE
    / "artifacts/json"
    / "epr_sampled_face_receipt_povm_descent_027.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_sampled_coherent_bell_capacity_obstruction_028.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_sampled_coherent_bell_capacity_obstruction_028.md"
)

TOL = 5.0e-12


def load(path):
    return json.loads(
        path.read_text()
    )


def digest_json(obj):
    raw = json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")

    return sha256(raw).hexdigest()


def calibration_id(row):
    return (
        str(row["model"])
        + "@"
        + str(
            row["dimensionless_pulse"]
        )
    )


print("PROGRESS: 1/7 load sealed interfaces")

a022 = load(A022)
a027 = load(A027)

checks = {}

checks["Audit022_passed"] = (
    a022["audit_pass"] is True
)

checks["Audit027_passed"] = (
    a027["audit_pass"] is True
)

checks["Audit027_face_POVM_constructed"] = (
    a027[
        "boundary"
    ][
        "sampled_face_qubit_POVM_constructed"
    ]
    is True
)

checks["Audit027_state_update_not_inserted"] = (
    a027[
        "boundary"
    ][
        "state_update_instrument_constructed"
    ]
    is False
)

print("PROGRESS: 2/7 recover binary observable coefficients")

calibrations = (
    a027[
        "calibrations"
    ]
)

checks["sixteen_calibrations"] = (
    len(calibrations) == 16
)

local_rows = []

for row in calibrations:
    b = float(
        row[
            "binary_observable"
        ]["bias"]
    )

    c = float(
        row[
            "binary_observable"
        ]["contrast"]
    )

    p0_axis = float(
        row[
            "pointer_0"
        ]["axis"]
    )

    p0_plane = float(
        row[
            "pointer_0"
        ]["plane"]
    )

    lambda_axis = (
        2.0 * p0_axis
        - 1.0
    )

    lambda_plane = (
        2.0 * p0_plane
        - 1.0
    )

    reconstructed_axis = (
        b + c
    )

    reconstructed_plane = (
        b - c
    )

    local_rows.append({
        "id":
            calibration_id(row),

        "model":
            row["model"],

        "dimensionless_pulse":
            float(
                row[
                    "dimensionless_pulse"
                ]
            ),

        "bias":
            b,

        "contrast":
            c,

        "observable_axis_eigenvalue":
            lambda_axis,

        "observable_plane_eigenvalue":
            lambda_plane,

        "reconstruction_error":
            max(
                abs(
                    lambda_axis
                    - reconstructed_axis
                ),
                abs(
                    lambda_plane
                    - reconstructed_plane
                ),
            ),

        "operator_norm":
            max(
                abs(
                    lambda_axis
                ),
                abs(
                    lambda_plane
                ),
            ),
    })

max_reconstruction_error = max(
    row[
        "reconstruction_error"
    ]
    for row in local_rows
)

max_local_operator_norm = max(
    row[
        "operator_norm"
    ]
    for row in local_rows
)

checks["bias_contrast_reconstruct_observable_eigenvalues"] = (
    max_reconstruction_error
    <= TOL
)

checks["all_local_binary_observables_are_contractions"] = (
    max_local_operator_norm
    <= 1.0 + TOL
)

print("PROGRESS: 3/7 derive state-independent CHSH bound")

# For
#
#   A_x = s_x (b_A I + c_A N_x)
#   B_y = t_y (b_B I + c_B M_y)
#
# with N_x^2=M_y^2=I and local outcome signs
# s_x,t_y in {+1,-1}, expand
#
#   C = A0 B0 + A0 B1 + A1 B0 - A1 B1.
#
# For every sign assignment:
#
# identity contribution:
#   norm <= 2 |b_A b_B|
#
# one Alice-local centered term survives:
#   norm <= 2 |c_A b_B|
#
# one Bob-local centered term survives:
#   norm <= 2 |b_A c_B|
#
# centered CHSH contribution:
#   norm <= 2 sqrt(2) |c_A c_B|
#
# by the universal Tsirelson operator inequality.
#
# Therefore
#
#   ||C|| <=
#     2|b_A b_B|
#     + 2|b_A c_B|
#     + 2|c_A b_B|
#     + 2sqrt(2)|c_A c_B|.
#
# This is state-independent and does not require a
# probability/frequency interpretation.

pair_rows = []

for alice in local_rows:
    for bob in local_rows:
        bA = alice["bias"]
        cA = alice["contrast"]
        bB = bob["bias"]
        cB = bob["contrast"]

        bound = (
            2.0
            * abs(
                bA * bB
            )
            +
            2.0
            * abs(
                bA * cB
            )
            +
            2.0
            * abs(
                cA * bB
            )
            +
            2.0
            * math.sqrt(2.0)
            * abs(
                cA * cB
            )
        )

        pair_rows.append({
            "Alice":
                alice["id"],

            "Bob":
                bob["id"],

            "Alice_bias":
                bA,

            "Alice_contrast":
                cA,

            "Bob_bias":
                bB,

            "Bob_contrast":
                cB,

            "state_independent_CHSH_upper_bound":
                bound,

            "exceeds_2":
                bound
                > 2.0 + TOL,
        })

checks["all_256_local_calibration_pairs_tested"] = (
    len(pair_rows)
    == 16 * 16
)

print("PROGRESS: 4/7 classify Bell-capacity surface")

max_pair = max(
    pair_rows,
    key=lambda row:
        row[
            "state_independent_CHSH_upper_bound"
        ],
)

max_bound = float(
    max_pair[
        "state_independent_CHSH_upper_bound"
    ]
)

violating_bound_rows = [
    row
    for row in pair_rows
    if row[
        "exceeds_2"
    ]
]

nonzero_pulse_pairs = [
    row
    for row in pair_rows
    if (
        float(
            row[
                "Alice"
            ].split("@")[1]
        )
        > 0.0
        and float(
            row[
                "Bob"
            ].split("@")[1]
        )
        > 0.0
    )
]

max_nonzero_pair = max(
    nonzero_pulse_pairs,
    key=lambda row:
        row[
            "state_independent_CHSH_upper_bound"
        ],
)

max_nonzero_bound = float(
    max_nonzero_pair[
        "state_independent_CHSH_upper_bound"
    ]
)

checks["no_sampled_pair_has_upper_bound_above_2"] = (
    len(
        violating_bound_rows
    )
    == 0
)

checks["global_maximum_bound_at_most_2"] = (
    max_bound
    <= 2.0 + TOL
)

checks["all_nonzero_pulse_pairs_strictly_below_2"] = (
    max_nonzero_bound
    < 2.0 - TOL
)

print("PROGRESS: 5/7 compare ideal geometry with sampled response")

ideal_singlet_max = (
    1.0
    + 3.0
    / math.sqrt(5.0)
)

ideal_unbiased_equal_visibility_threshold = math.sqrt(
    2.0
    / ideal_singlet_max
)

max_abs_contrast = float(
    a027[
        "binary_response_summary"
    ][
        "max_absolute_contrast"
    ]
)

checks["sampled_max_contrast_below_ideal_unbiased_threshold"] = (
    max_abs_contrast
    < ideal_unbiased_equal_visibility_threshold
)

print("PROGRESS: 6/7 preserve probability and state boundaries")

checks["no_joint_probability_law_used"] = (
    True
)

checks["no_Born_frequency_rule_used"] = (
    True
)

checks["no_state_update_rule_used"] = (
    True
)

checks["canonical_EPR_state_not_needed_for_obstruction"] = (
    True
)

checks["ideal_Bell_geometry_not_rejected"] = (
    a022[
        "boundary"
    ][
        "Bell_capable_native_setting_quartets_exist"
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
    "the_sampled_coherent_face_POVM_family_is_state_"
    "independently_CHSH_incapable_for_all_256_Alice_Bob_"
    "calibration_pairs_under_a_universal_operator_norm_"
    "bound_even_though_the_native_six_axis_geometry_itself_"
    "remains_Bell_capable"
    if audit_pass
    else
    "sampled_coherent_Bell_capacity_obstruction_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_sampled_coherent_bell_capacity_obstruction_028",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "local_binary_observable": {
        "formula":
            "A_i = bias I + contrast N_i",

        "N_i_squared":
            "I",

        "calibration_count":
            len(
                local_rows
            ),

        "max_reconstruction_error":
            max_reconstruction_error,

        "max_local_operator_norm":
            max_local_operator_norm,

        "max_absolute_contrast":
            max_abs_contrast,
    },

    "state_independent_bound": {
        "formula":
            (
                "2|b_A b_B| + 2|b_A c_B| "
                "+ 2|c_A b_B| "
                "+ 2 sqrt(2)|c_A c_B|"
            ),

        "derivation":
            (
                "identity term plus one surviving local centered "
                "term on each wing plus the universal Tsirelson "
                "bound on the centered CHSH operator"
            ),

        "Alice_Bob_calibration_pair_count":
            len(
                pair_rows
            ),

        "violating_upper_bound_count":
            len(
                violating_bound_rows
            ),

        "maximum_upper_bound":
            max_bound,

        "maximum_upper_bound_pair":
            max_pair,

        "maximum_nonzero_pulse_upper_bound":
            max_nonzero_bound,

        "maximum_nonzero_pulse_pair":
            max_nonzero_pair,
    },

    "ideal_geometry_reference": {
        "Audit022_ideal_singlet_max":
            ideal_singlet_max,

        "equal_unbiased_visibility_threshold_for_that_value":
            ideal_unbiased_equal_visibility_threshold,

        "sampled_max_absolute_contrast":
            max_abs_contrast,

        "comparison_is_reference_not_main_proof":
            True,
    },

    "earned_statement": (
        "The sampled coherent face POVMs from Audit 027 are "
        "locally normalized and positive, but their binary "
        "pointer observables are too weakly setting-dependent "
        "to support CHSH violation. Writing each local observable "
        "as bias times identity plus contrast times a Hermitian "
        "involution yields a state-independent CHSH operator "
        "upper bound. All 256 ordered Alice/Bob combinations of "
        "the sixteen sampled calibrations satisfy that bound at "
        "or below two, and every pair with nonzero pulse lies "
        "strictly below two. This obstruction does not use a "
        "joint probability law, Born-frequency interpretation, "
        "state-update rule, or the canonical EPR preparation. "
        "It therefore isolates the inherited sampled coherent "
        "response law itself as Bell-incapable while leaving the "
        "native six-axis Bell geometry from Audit 022 intact."
    ),

    "checks":
        checks,

    "boundary": {
        "sampled_coherent_response_Bell_capable":
            False,

        "obstruction_state_independent":
            True,

        "obstruction_probability_law_independent":
            True,

        "native_six_axis_geometry_Bell_capable":
            True,

        "canonical_EPR_preparation_rejected":
            False,

        "face_receipt_POVM_rejected":
            False,

        "sampled_scalar_response_strength_is_obstruction":
            True,

        "continuous_unsampled_controls_ruled_out":
            False,

        "other_dynamics_ruled_out":
            False,

        "time_dependent_or_weighted_controls_ruled_out":
            False,

        "Born_frequency_law_derived":
            False,

        "no_signaling_test_required_for_this_operator_obstruction":
            False,
    },

    "next_gate": (
        "Do not build joint probabilities from a sampled local "
        "response that is already Bell-incapable. Search the native "
        "local transduction/control space for a normalized response "
        "with substantially stronger analyzer contrast while retaining "
        "all clean/nonclean receipts and chi_deck covariance. Candidate "
        "directions include nonconstant controls, different native "
        "generators, or a genuinely derived weighting law. Any new "
        "response must first pass local normalization and then this "
        "state-independent Bell-capacity gate before joint EPR weighting."
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

note = f"""# Sampled coherent Bell-capacity obstruction 028

## Result

Audit pass:

    {audit_pass}

The local binary pointer observable is

    A_i = bias I + contrast N_i

with

    N_i^2 = I.

For arbitrary local outcome relabelings and arbitrary analyzer
directions, the CHSH operator obeys

    ||CHSH||
      <=
    2 |b_A b_B|
      +
    2 |b_A c_B|
      +
    2 |c_A b_B|
      +
    2 sqrt(2) |c_A c_B|.

This is state-independent.

It does not use a probability law.

## Sampled census

Local calibration count:

    {len(local_rows)}

Alice/Bob calibration pairs:

    {len(pair_rows)}

Upper bounds exceeding 2:

    {len(violating_bound_rows)}

Maximum upper bound:

    {max_bound}

Pair attaining it:

    {max_pair}

Maximum upper bound with nonzero pulse on both wings:

    {max_nonzero_bound}

Pair attaining it:

    {max_nonzero_pair}

## Interpretation

Audit 022 remains intact:

    the native six-axis geometry is Bell-capable.

Audit 027 also remains intact:

    the sampled coherent local receipt POVMs are positive and normalized.

The obstruction lies in the sampled response strength.

Maximum sampled analyzer contrast:

    {max_abs_contrast}

The inherited sampled coherent response therefore cannot transmit enough
setting dependence to support CHSH violation, regardless of joint state.

No Born-frequency law or state-update instrument is used in this
obstruction.

The next target is a stronger native local transduction law, not joint
probability bookkeeping for a response already below Bell capacity.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print(
    "CALIBRATION_COUNT:",
    len(local_rows),
)
print(
    "CALIBRATION_PAIR_COUNT:",
    len(pair_rows),
)
print(
    "MAX_LOCAL_OPERATOR_NORM:",
    max_local_operator_norm,
)
print(
    "MAX_ABS_CONTRAST:",
    max_abs_contrast,
)
print(
    "MAX_CHSH_UPPER_BOUND:",
    max_bound,
)
print(
    "MAX_CHSH_UPPER_BOUND_PAIR:",
    max_pair,
)
print(
    "MAX_NONZERO_PULSE_UPPER_BOUND:",
    max_nonzero_bound,
)
print(
    "MAX_NONZERO_PULSE_PAIR:",
    max_nonzero_pair,
)
print(
    "UPPER_BOUND_GT2_COUNT:",
    len(
        violating_bound_rows
    ),
)
print(
    "IDEAL_UNBIASED_EQUAL_VISIBILITY_THRESHOLD:",
    ideal_unbiased_equal_visibility_threshold,
)
print(
    "JOINT_PROBABILITY_LAW_USED:",
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
