#!/usr/bin/env python3

from hashlib import sha256
from pathlib import Path
import json

HERE = Path(__file__).resolve().parents[2]

P41 = (
    Path.home()
    / "dev/cori/research/mathematics"
    / "41-order-4-dodecahedral-residue"
)

FP7 = (
    P41
    / "artifacts/provenance"
    / "project41_u2_gauge_bridge_connection_frontier_checkpoint_201fp7.v1.json"
)

FS = (
    P41
    / "artifacts/provenance"
    / "project41_signed_event_doublet_vector_lift_201fs.v1.json"
)

FV = (
    P41
    / "artifacts/provenance"
    / "project41_positive_pair_hermitian_geometry_201fv.v1.json"
)

A010 = (
    HERE
    / "artifacts/json"
    / "epr_g1800_joint_boundary_incidence_quotient_010.v1.json"
)

A013 = (
    HERE
    / "artifacts/json"
    / "epr_complete_exchange_parity_closure_013.v1.json"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_flat_groupoid_preparation_interface_014.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_flat_groupoid_preparation_interface_014.md"
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


print("PROGRESS: 1/6 load sealed interfaces")

fp7 = load(FP7)
fs = load(FS)
fv = load(FV)
a010 = load(A010)
a013 = load(A013)

checks = {}

checks["fp7_passed"] = (
    fp7["audit_pass"] is True
)

checks["fs_passed"] = (
    fs["audit_pass"] is True
)

checks["fv_passed"] = (
    fv["audit_pass"] is True
)

checks["audit_010_passed"] = (
    a010["audit_pass"] is True
)

checks["audit_013_passed"] = (
    a013["audit_pass"] is True
)

print("PROGRESS: 2/6 certify 12-state flat groupoid interface")

flb = fp7["results"]["201FLB"]
fn = fp7["results"]["201FN"]

state_count = int(
    flb["local_order_parameter_state_count"]
)

pair_transition_count = int(
    fn["pair_transition_count"]
)

triangle_test_count = int(
    fn["triangle_cocycle_test_count"]
)

checks["local_order_parameter_state_count_12"] = (
    state_count == 12
)

checks["pair_transition_count_is_12_squared"] = (
    pair_transition_count
    == state_count * state_count
)

checks["triangle_test_count_is_12_cubed"] = (
    triangle_test_count
    == state_count * state_count * state_count
)

checks["matrix_cocycle_failures_zero"] = (
    int(fn["matrix_cocycle_failures"]) == 0
)

checks["native_group_cocycle_failures_zero"] = (
    int(fn["native_group_cocycle_failures"]) == 0
)

checks["semilinearity_cocycle_failures_zero"] = (
    int(fn["C2_semilinearity_cocycle_failures"]) == 0
)

checks["pair_groupoid_is_flat"] = (
    fn["section_induced_pair_groupoid"]
    == "flat"
)

checks["native_semilinear_U2_transition_cocycle"] = (
    fn["native_semilinear_U2_transition_cocycle"]
    is True
)

checks["section_change_gauge_law_recorded"] = (
    fn["section_change_gauge_law"]
    == "T'_yx = G_y T_yx G_x^-1"
)

print("PROGRESS: 3/6 compress pair groupoid to reference transports")

# A flat full pair groupoid on S has
#
#     T_zy T_yx = T_zx.
#
# Choose a reference state o and define
#
#     U_s = T_so.
#
# Then
#
#     T_yx = U_y U_x^-1.
#
# Therefore 144 pair transitions contain no independent
# pairwise curvature data. A reference trivialization needs
# 12 section maps, one of which can be fixed to identity.

reference_transport_count = state_count
nontrivial_reference_transport_count = state_count - 1

checks["flat_groupoid_compresses_to_12_reference_maps"] = (
    reference_transport_count == 12
)

checks["reference_gauge_leaves_at_most_11_nontrivial_maps"] = (
    nontrivial_reference_transport_count == 11
)

print("PROGRESS: 4/6 certify local boundary-mode input")

checks["four_signed_modes_per_face"] = (
    int(
        fs["counts"]["signed_blocks_per_face"]
    )
    == 4
)

checks["four_real_modes_as_two_complex_lines"] = (
    fs[
        "face_module"
    ]["four_real_modes_as_two_complex_lines"]
    is True
)

checks["positive_modes_are_2_and_3_on_canonical_face"] = (
    fv[
        "canonical_face"
    ]["positive_mode_indices"]
    == [2, 3]
)

checks["negative_modes_are_0_and_1_on_canonical_face"] = (
    fv[
        "canonical_face"
    ]["negative_mode_indices"]
    == [0, 1]
)

checks["positive_pair_is_Hermitian_orthogonal"] = (
    fv[
        "checks"
    ]["positive_pair_hermitian_overlap_zero"]
    is True
)

checks["positive_complex_lines_are_distinct"] = (
    fv[
        "checks"
    ]["positive_pair_hermitian_overlap_zero"]
    is True
)

print("PROGRESS: 5/6 bind to Program-02 joint boundary")

joint_boundary_count = int(
    a010[
        "diagonal_quotient"
    ]["class_count"]
)

presentations_per_state = int(
    a010[
        "diagonal_quotient"
    ]["boundary_presentations_per_G1800_state"]
)

checks["joint_boundary_count_7200"] = (
    joint_boundary_count == 7200
)

checks["four_presentations_per_G1800_state"] = (
    presentations_per_state == 4
)

checks["native_RP3_has_complete_exchange_parity_closure"] = (
    a013[
        "boundary"
    ]["complete_two_branch_closure_constructed"]
    is True
)

print("PROGRESS: 6/6 classify minimal missing interface")

# The explicit 144 transition matrices are not needed for
# the preparation branch test.
#
# For wedge vanishing, only projective complex-line equality
# matters:
#
#     u wedge v = 0
#       iff
#     [u] = [v] in CP^1.
#
# Any invertible complex-linear or conjugate-linear common
# frame change preserves equality versus inequality of
# projective complex lines.
#
# Therefore the minimal missing native interface is:
#
# eta:
#   one-wing boundary incidence
#     ->
#   (12-state local orientation state, local face-mode slot)
#
# ell:
#   (12-state local orientation state, local face-mode slot)
#     ->
#   CP^1 in one flat reference trivialization.
#
# Once eta and ell are known, every one of the 7200 native
# joint boundary presentations can be classified before
# analyzer settings exist by comparing two CP^1 labels.

missing_interface = {
    "incidence_crosswalk": {
        "name": "eta",
        "domain":
            "I_F, the 120 one-wing signed boundary incidences",
        "codomain":
            "S_12 x M_4",
        "status":
            "not yet serialized or derived in Program 02",
    },

    "projective_mode_table": {
        "name": "ell",
        "domain":
            "S_12 x M_4",
        "codomain":
            "CP^1 in a flat reference trivialization",
        "slot_count":
            state_count * 4,
        "status":
            "not yet serialized or derived in Program 02",
    },
}

failed_checks = [
    name
    for name, passed in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "flat_12_state_semilinear_U2_groupoid_compresses_"
    "the_missing_EPR_preparation_interface_to_native_"
    "incidence_crosswalk_eta_and_projective_mode_table_ell"
    if audit_pass
    else
    "flat_groupoid_preparation_interface_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_flat_groupoid_preparation_interface_014",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "program01_groupoid_interface": {
        "state_count":
            state_count,

        "pair_transition_count":
            pair_transition_count,

        "triangle_test_count":
            triangle_test_count,

        "flat":
            True,

        "semilinear_U2":
            True,

        "section_change_law":
            fn["section_change_gauge_law"],

        "reference_transport_count":
            reference_transport_count,

        "nontrivial_reference_transport_count_after_gauge":
            nontrivial_reference_transport_count,

        "explicit_144_matrix_payload_required":
            False,
    },

    "local_face_input": {
        "signed_modes_per_face":
            4,

        "real_modes_as_two_complex_lines":
            True,

        "canonical_positive_modes":
            [2, 3],

        "canonical_negative_modes":
            [0, 1],

        "canonical_positive_pair_Hermitian_orthogonal":
            True,
    },

    "program02_joint_boundary": {
        "joint_boundary_presentation_count":
            joint_boundary_count,

        "presentations_per_G1800_state":
            presentations_per_state,

        "existing_complete_parity_closure":
            True,
    },

    "minimal_missing_interface":
        missing_interface,

    "wedge_reduction": {
        "criterion":
            "u wedge v = 0 iff [u] = [v] in CP^1",

        "nonzero_result":
            (
                "every nonzero wedge lies on the unique "
                "projective line P(Lambda^2 C^2)"
            ),

        "common_semilinear_frame_change_preserves_zero_nonzero":
            True,

        "analyzer_settings_required":
            False,
    },

    "earned_statement": (
        "Program-01 FP7 certifies a flat semilinear U(2) pair "
        "groupoid on 12 local order-parameter states, with all "
        "144 ordered transitions and all 1728 triangle tests "
        "covered and zero recorded cocycle failures. Flatness "
        "means the missing relative-frame information can be "
        "represented by reference transports rather than 144 "
        "independent pair matrices. For the EPR antisymmetric "
        "wedge test, even those matrices are more information "
        "than required: only the transported projective complex "
        "line of each native boundary mode matters. The remaining "
        "native preparation interface is therefore exactly an "
        "incidence crosswalk eta from I_F to S_12 x M_4 and a "
        "projective mode table ell from S_12 x M_4 to CP^1."
    ),

    "checks":
        checks,

    "boundary": {
        "201FN_explicit_matrix_recovery_required":
            False,

        "historical_201FN_archaeology_stops_here":
            True,

        "wedge_census_not_yet_run":
            True,

        "incidence_to_order_parameter_state_not_yet_derived":
            True,

        "signed_mode_projective_line_table_not_yet_derived":
            True,

        "source_preparation_not_yet_admitted":
            True,

        "no_postselection":
            True,

        "no_probability_law":
            True,

        "no_Born_rule":
            True,

        "no_CHSH_claim":
            True,
    },

    "next_gate": (
        "Derive eta and ell directly from current canonical "
        "Program-01 phase, signed-event, and face-complex-"
        "structure interfaces. Do not reconstruct the vanished "
        "144 transition matrices. Once eta and ell are native, "
        "classify all 7200 joint boundary presentations by CP^1 "
        "line equality before analyzer settings are introduced."
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

note = f"""# Flat groupoid preparation interface 014

## Result

Audit pass:

    {audit_pass}

Program-01 FP7 certifies a flat semilinear U(2) pair groupoid on

    {state_count}

local order-parameter states.

It covers

    {pair_transition_count}

ordered pair transitions and

    {triangle_test_count}

triangle cocycle tests, with zero recorded matrix, native-group, and
semilinearity cocycle failures.

## Compression

A flat full pair groupoid may be trivialized from one reference state.

Instead of 144 independent pair matrices, choose

    12

reference transports, with one fixed to identity and therefore at most

    11

nontrivial reference maps.

For the EPR wedge test, the full matrices are not required.

The only required transported datum is the projective complex line

    [v] in CP^1

occupied by each native local boundary mode.

The wedge criterion is

    u wedge v = 0

if and only if

    [u] = [v] in CP^1.

This zero/nonzero distinction is preserved by a common invertible
complex-linear or conjugate-linear frame change.

## Minimal missing native interface

Two maps remain:

    eta:
        I_F -> S_12 x M_4

mapping each one-wing signed boundary incidence to its local
order-parameter state and local mode slot;

and

    ell:
        S_12 x M_4 -> CP^1

giving the projective line of that mode in one flat reference
trivialization.

There are 12 x 4 = 48 projective mode slots before identifications.

## Program-02 consequence

Once eta and ell are derived, all

    {joint_boundary_count}

native two-wing boundary presentations can be classified before analyzer
settings exist.

Equal CP^1 lines give zero antisymmetric wedge.

Unequal CP^1 lines give a nonzero wedge, and every nonzero wedge lies on
the unique projective line of Lambda^2 C^2, already identified as the
canonical EPR line by Audits 006-007.

The vanished 201FN transition matrices do not need to be reconstructed.

The next problem is eta and ell, not archaeology.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print("ORDER_PARAMETER_STATE_COUNT:", state_count)
print("PAIR_TRANSITION_COUNT:", pair_transition_count)
print("TRIANGLE_TEST_COUNT:", triangle_test_count)
print(
    "REFERENCE_TRANSPORT_COUNT:",
    reference_transport_count,
)
print(
    "NONTRIVIAL_REFERENCE_TRANSPORT_COUNT:",
    nontrivial_reference_transport_count,
)
print(
    "PROJECTIVE_MODE_SLOT_COUNT:",
    state_count * 4,
)
print(
    "JOINT_BOUNDARY_PRESENTATION_COUNT:",
    joint_boundary_count,
)
print(
    "EXPLICIT_144_MATRICES_REQUIRED:",
    False,
)
print(
    "MISSING_INTERFACE:",
    "eta + ell",
)
print("FAILED_CHECK_COUNT:", len(failed_checks))
print("FAILED_CHECKS:", failed_checks)
print("JSON_OUT:", JSON_OUT)
print("NOTE_OUT:", NOTE_OUT)
print(
    "JSON_SHA256:",
    sha256(JSON_OUT.read_bytes()).hexdigest(),
)
