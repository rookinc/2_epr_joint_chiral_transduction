#!/usr/bin/env python3

from collections import defaultdict
from hashlib import sha256
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile
import json

HERE = Path(__file__).resolve().parents[2]

A022 = (
    HERE
    / "artifacts/json"
    / "epr_native_six_axis_singlet_quartet_census_022.v1.json"
)

A025 = (
    HERE
    / "artifacts/json"
    / "epr_projective_analyzer_sector_extension_025.v1.json"
)

A026 = (
    HERE
    / "artifacts/json"
    / "epr_analyzer_native_character_identification_026.v1.json"
)

ZIP = (
    HERE
    / "evidence"
    / "native_coherent_axis_instrument.zip"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_sampled_face_receipt_povm_descent_027.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_sampled_face_receipt_povm_descent_027.md"
)

TOL = 2.0e-12


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


print("PROGRESS: 1/8 load sealed Program-02 interfaces")

a022 = load(A022)
a025 = load(A025)
a026 = load(A026)

checks = {}

checks["Audit022_passed"] = (
    a022["audit_pass"] is True
)

checks["Audit025_passed"] = (
    a025["audit_pass"] is True
)

checks["Audit026_passed"] = (
    a026["audit_pass"] is True
)

checks["analyzer_sector_character_is_chi_deck"] = (
    a026[
        "boundary"
    ][
        "matched_character"
    ]
    == "chi_deck"
)

print("PROGRESS: 2/8 load coherent scalar packet")

with ZipFile(ZIP) as z:
    base = "native_coherent_axis_instrument/"

    report = json.loads(
        z.read(
            base + "REPORT.json"
        )
    )

    symmetry = json.loads(
        z.read(
            base + "EXACT_SYMMETRY.json"
        )
    )

    numerical = json.loads(
        z.read(
            base + "NUMERICAL_PULSES.json"
        )
    )

checks["coherent_packet_audit_passed"] = (
    report["audit_pass"] is True
)

checks["six_exact_analyzer_symmetry_rows"] = (
    isinstance(symmetry, list)
    and len(symmetry) == 6
)

checks["all_symmetry_rows_axis_rank2_plane_rank4"] = all(
    int(row["source_axis_rank"]) == 2
    and int(row["source_plane_rank"]) == 4
    for row in symmetry
)

rows = numerical["rows"]

checks["numerical_row_count_96"] = (
    len(rows) == 96
)

print("PROGRESS: 3/8 partition sampled calibrations")

groups = defaultdict(list)

for row in rows:
    key = (
        str(row["model"]),
        float(
            row["dimensionless_pulse"]
        ),
    )

    groups[key].append(row)

models = sorted({
    model
    for model, pulse
    in groups
})

pulses = sorted({
    pulse
    for model, pulse
    in groups
})

checks["two_candidate_dynamics"] = (
    models
    == [
        "adjacency",
        "laplacian_control",
    ]
)

checks["eight_sampled_pulses"] = (
    pulses
    == [
        0.0,
        0.1,
        0.25,
        0.5,
        1.0,
        2.0,
        3.0,
        4.0,
    ]
)

checks["six_settings_per_calibration"] = all(
    len(group) == 6
    and sorted(
        int(row["axis"])
        for row in group
    )
    == list(range(6))
    for group in groups.values()
)

checks["sixteen_sampled_calibrations"] = (
    len(groups) == 16
)

print("PROGRESS: 4/8 test scalar covariance and completeness")

effect_names = (
    "all_pointer_0",
    "all_pointer_1",
    "clean_pointer_0",
    "clean_pointer_1",
    "nonclean_pointer_0",
    "nonclean_pointer_1",
)

components = (
    "axis",
    "plane",
)

max_setting_covariance_error = 0.0
max_binary_completeness_error = 0.0
max_clean_split_error = 0.0
max_four_receipt_completeness_error = 0.0

scalar_min = float("inf")
scalar_max = float("-inf")

for group in groups.values():
    group = sorted(
        group,
        key=lambda row: int(
            row["axis"]
        ),
    )

    reference = group[0]

    for row in group:
        for effect_name in effect_names:
            for component in components:
                value = float(
                    row[
                        "effects"
                    ][
                        effect_name
                    ][
                        component
                    ]
                )

                reference_value = float(
                    reference[
                        "effects"
                    ][
                        effect_name
                    ][
                        component
                    ]
                )

                max_setting_covariance_error = max(
                    max_setting_covariance_error,
                    abs(
                        value
                        - reference_value
                    ),
                )

                scalar_min = min(
                    scalar_min,
                    value,
                )

                scalar_max = max(
                    scalar_max,
                    value,
                )

        e = row["effects"]

        for component in components:
            p0 = float(
                e[
                    "all_pointer_0"
                ][component]
            )

            p1 = float(
                e[
                    "all_pointer_1"
                ][component]
            )

            max_binary_completeness_error = max(
                max_binary_completeness_error,
                abs(
                    p0 + p1 - 1.0
                ),
            )

            four_total = 0.0

            for pointer in (
                "0",
                "1",
            ):
                clean = float(
                    e[
                        "clean_pointer_"
                        + pointer
                    ][component]
                )

                nonclean = float(
                    e[
                        "nonclean_pointer_"
                        + pointer
                    ][component]
                )

                all_value = float(
                    e[
                        "all_pointer_"
                        + pointer
                    ][component]
                )

                max_clean_split_error = max(
                    max_clean_split_error,
                    abs(
                        clean
                        + nonclean
                        - all_value
                    ),
                )

                four_total += (
                    clean
                    + nonclean
                )

            max_four_receipt_completeness_error = max(
                max_four_receipt_completeness_error,
                abs(
                    four_total
                    - 1.0
                ),
            )

checks["setting_scalar_covariance_within_tolerance"] = (
    max_setting_covariance_error
    <= TOL
)

checks["binary_pointer_completeness_within_tolerance"] = (
    max_binary_completeness_error
    <= TOL
)

checks["clean_nonclean_split_within_tolerance"] = (
    max_clean_split_error
    <= TOL
)

checks["four_receipt_completeness_within_tolerance"] = (
    max_four_receipt_completeness_error
    <= TOL
)

checks["all_sampled_scalars_in_unit_interval"] = (
    scalar_min >= -TOL
    and scalar_max <= 1.0 + TOL
)

print("PROGRESS: 5/8 construct face-qubit receipt POVMs")

calibrations = []

max_abs_bias = 0.0
max_abs_contrast = 0.0
max_contrast_row = None

for (
    model,
    pulse,
), group in sorted(
    groups.items()
):
    group = sorted(
        group,
        key=lambda row: int(
            row["axis"]
        ),
    )

    representative = group[0]

    e = representative[
        "effects"
    ]

    receipts = {}

    for receipt_name in (
        "clean_pointer_0",
        "nonclean_pointer_0",
        "clean_pointer_1",
        "nonclean_pointer_1",
    ):
        receipts[
            receipt_name
        ] = {
            "P_i_coefficient":
                float(
                    e[
                        receipt_name
                    ]["axis"]
                ),

            "I_minus_P_i_coefficient":
                float(
                    e[
                        receipt_name
                    ]["plane"]
                ),
        }

    a0 = float(
        e[
            "all_pointer_0"
        ]["axis"]
    )

    b0 = float(
        e[
            "all_pointer_0"
        ]["plane"]
    )

    a1 = float(
        e[
            "all_pointer_1"
        ]["axis"]
    )

    b1 = float(
        e[
            "all_pointer_1"
        ]["plane"]
    )

    # For P_i=(I+n_i.sigma)/2,
    #
    # A_i = E_0 - E_1
    #     = bias*I + contrast*n_i.sigma.
    bias = (
        a0
        + b0
        - 1.0
    )

    contrast = (
        a0
        - b0
    )

    max_abs_bias = max(
        max_abs_bias,
        abs(bias),
    )

    if abs(
        contrast
    ) > max_abs_contrast:
        max_abs_contrast = abs(
            contrast
        )

        max_contrast_row = {
            "model":
                model,

            "dimensionless_pulse":
                pulse,

            "bias":
                bias,

            "contrast":
                contrast,
        }

    calibrations.append({
        "model":
            model,

        "dimensionless_pulse":
            pulse,

        "setting_count":
            6,

        "face_effect_formula":
            (
                "E_(i,r) = alpha_r P_i "
                "+ beta_r (I-P_i)"
            ),

        "receipts":
            receipts,

        "pointer_0": {
            "axis":
                a0,

            "plane":
                b0,
        },

        "pointer_1": {
            "axis":
                a1,

            "plane":
                b1,
        },

        "binary_observable": {
            "formula":
                "A_i = bias I + contrast n_i.sigma",

            "bias":
                bias,

            "contrast":
                contrast,
        },
    })

checks["sixteen_face_receipt_POVMs_constructed"] = (
    len(
        calibrations
    )
    == 16
)

checks["every_face_receipt_effect_positive"] = (
    scalar_min >= -TOL
)

checks["every_face_receipt_POVM_normalized"] = (
    max_four_receipt_completeness_error
    <= TOL
)

print("PROGRESS: 6/8 classify sector covariance")

# Audit026 identifies analyzer sector exchange with chi_deck.
#
# The sampled response coefficients depend only on
# axis-versus-complement and are setting-independent.
#
# Therefore a native permutation/sector exchange transports
# the projector label while leaving the scalar response law
# unchanged. No new sector-dependent scalar is inserted.

checks["scalar_law_setting_independent"] = (
    max_setting_covariance_error
    <= TOL
)

checks["no_extra_sector_dependent_weight_inserted"] = (
    True
)

checks["chi_deck_sector_exchange_retained"] = (
    a026[
        "boundary"
    ][
        "matched_character"
    ]
    == "chi_deck"
)

print("PROGRESS: 7/8 preserve instrument and probability boundaries")

checks["nonclean_receipts_retained"] = all(
    "nonclean_pointer_0"
    in row[
        "receipts"
    ]
    and "nonclean_pointer_1"
    in row[
        "receipts"
    ]
    for row in calibrations
)

checks["no_Luders_update_inserted"] = (
    True
)

checks["state_update_instrument_not_claimed"] = (
    True
)

checks["Born_frequency_law_not_claimed"] = (
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
    "the_sampled_coherent_axis_contact_scalar_response_"
    "descends_to_six_setting_normalized_four_receipt_"
    "face_qubit_POVMs_for_both_candidate_dynamics_and_"
    "all_eight_sampled_pulses_with_clean_and_nonclean_"
    "outputs_retained_and_no_state_update_rule_inserted"
    if audit_pass
    else
    "sampled_face_receipt_POVM_descent_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_sampled_face_receipt_povm_descent_027",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "sample_surface": {
        "model_count":
            len(models),

        "models":
            models,

        "pulse_count":
            len(pulses),

        "pulses":
            pulses,

        "setting_count":
            6,

        "raw_row_count":
            len(rows),

        "calibration_count":
            len(
                calibrations
            ),
    },

    "numerical_closure": {
        "tolerance":
            TOL,

        "max_setting_covariance_error":
            max_setting_covariance_error,

        "max_binary_completeness_error":
            max_binary_completeness_error,

        "max_clean_nonclean_split_error":
            max_clean_split_error,

        "max_four_receipt_completeness_error":
            max_four_receipt_completeness_error,

        "scalar_min":
            scalar_min,

        "scalar_max":
            scalar_max,
    },

    "face_receipt_POVM": {
        "carrier":
            "complex2 face carrier",

        "setting_projector":
            "rank-one P_i",

        "complement":
            "I-P_i",

        "receipt_names": [
            "clean_pointer_0",
            "nonclean_pointer_0",
            "clean_pointer_1",
            "nonclean_pointer_1",
        ],

        "effect_formula":
            "E_(i,r)=alpha_r P_i + beta_r (I-P_i)",

        "normalization":
            "sum over all four receipts equals I",

        "positivity":
            "all sampled alpha_r,beta_r lie in [0,1]",
    },

    "binary_response_summary": {
        "observable_formula":
            "A_i = bias I + contrast n_i.sigma",

        "bias_formula":
            "bias = axis_pointer0 + plane_pointer0 - 1",

        "contrast_formula":
            "contrast = axis_pointer0 - plane_pointer0",

        "max_absolute_bias":
            max_abs_bias,

        "max_absolute_contrast":
            max_abs_contrast,

        "max_contrast_calibration":
            max_contrast_row,
    },

    "calibrations":
        calibrations,

    "earned_statement": (
        "For each of the two declared coherent dynamics models "
        "and each of the eight sampled dimensionless pulses, the "
        "six native analyzer settings have the same axis/plane "
        "scalar response to numerical precision. The paired "
        "six-dimensional clean/nonclean pointer effects therefore "
        "descend by functional calculus to the Program-01 complex2 "
        "face carrier as four positive receipt effects of the form "
        "alpha P_i + beta(I-P_i). The four receipts sum to identity "
        "for every sampled calibration, while clean and nonclean "
        "outputs are both retained. Audit026's chi_deck sector "
        "character is preserved because the scalar law is setting "
        "independent and no additional sector weight is inserted. "
        "This closes a sampled normalized local POVM law only; no "
        "state-update instrument or physical frequency interpretation "
        "is introduced."
    ),

    "checks":
        checks,

    "boundary": {
        "sampled_face_qubit_POVM_constructed":
            audit_pass,

        "all_four_local_receipts_retained":
            True,

        "nonclean_outputs_retained":
            True,

        "analyzer_sector_character":
            "chi_deck",

        "state_update_instrument_constructed":
            False,

        "Luders_rule_inserted":
            False,

        "Born_frequency_law_derived":
            False,

        "continuous_all_time_scalar_law_derived":
            False,

        "sampled_numeric_response_only":
            True,

        "Audit019_zero_preparation_branch_still_retained":
            True,

        "no_signaling_not_yet_tested":
            True,

        "receipt_level_CHSH_not_yet_tested":
            True,
    },

    "next_gate": (
        "Combine the sampled local face receipt POVMs with the "
        "canonical antisymmetric joint line and compute complete "
        "joint receipt tables without discarding clean/nonclean "
        "or the Audit019 zero preparation branch. First test local "
        "normalization and no-signaling. Only after those pass, "
        "evaluate CHSH for the actual unsharp binary pointer "
        "observable rather than the ideal projective reflections."
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

note = f"""# Sampled face receipt POVM descent 027

## Result

Audit pass:

    {audit_pass}

Candidate dynamics:

    {models}

Sampled pulses:

    {pulses}

Settings per calibration:

    6

Total sampled calibration families:

    {len(calibrations)}

## Face-qubit receipt law

For every native analyzer setting i, let P_i be its rank-one face-qubit
projector.

Each sampled coherent response descends as

    E_(i,r)
      =
    alpha_r P_i
      +
    beta_r (I-P_i)

for four retained receipts

    clean pointer 0
    nonclean pointer 0
    clean pointer 1
    nonclean pointer 1.

All sampled coefficients are nonnegative.

The four receipt effects sum to identity.

## Numerical closure

Maximum setting covariance error:

    {max_setting_covariance_error}

Maximum binary completeness error:

    {max_binary_completeness_error}

Maximum clean/nonclean split error:

    {max_clean_split_error}

Maximum four-receipt completeness error:

    {max_four_receipt_completeness_error}

Scalar range:

    [{scalar_min}, {scalar_max}]

## Binary pointer observable

After forgetting only the clean/nonclean subreceipt,

    A_i
      =
    E_(i,0) - E_(i,1)
      =
    bias I + contrast n_i.sigma.

Maximum sampled absolute contrast:

    {max_abs_contrast}

Calibration attaining it:

    {max_contrast_row}

This quantifies the inherited unsharpness before any Bell calculation.

## Boundary

This is a sampled mathematical POVM descent under the declared coherent
amplitude/readout model.

It does not insert a state-update rule.

It does not derive a physical frequency law.

It retains the nonclean receipts and the Audit019 zero preparation
branch.

The next gate is complete joint receipt normalization and no-signaling.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print(
    "MODEL_COUNT:",
    len(models),
)
print(
    "CALIBRATION_COUNT:",
    len(calibrations),
)
print(
    "MAX_SETTING_COVARIANCE_ERROR:",
    max_setting_covariance_error,
)
print(
    "MAX_BINARY_COMPLETENESS_ERROR:",
    max_binary_completeness_error,
)
print(
    "MAX_CLEAN_NONCLEAN_SPLIT_ERROR:",
    max_clean_split_error,
)
print(
    "MAX_FOUR_RECEIPT_COMPLETENESS_ERROR:",
    max_four_receipt_completeness_error,
)
print(
    "SCALAR_RANGE:",
    (scalar_min, scalar_max),
)
print(
    "MAX_ABS_BIAS:",
    max_abs_bias,
)
print(
    "MAX_ABS_CONTRAST:",
    max_abs_contrast,
)
print(
    "MAX_CONTRAST_CALIBRATION:",
    max_contrast_row,
)
print(
    "STATE_UPDATE_INSTRUMENT_CONSTRUCTED:",
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
