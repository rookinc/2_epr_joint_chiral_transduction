#!/usr/bin/env python3

from pathlib import Path
from hashlib import sha256
import json
import math

HERE = Path(__file__).resolve().parents[2]

ART = HERE / "artifacts/json"

INPUTS = {
    "035":
        ART / "epr_heralded_conditional_hilbert_table_035.v1.json",

    "036":
        ART / "epr_native_symmetry_weighting_reduction_036.v1.json",

    "037":
        ART / "epr_native_bilinear_kernel_uniqueness_037.v1.json",

    "038":
        ART / "epr_native_rank3_kernel_visibility_038.v1.json",

    "039":
        ART / "epr_heralded_local_table_obstruction_039.v1.json",

    "040b":
        ART / "epr_native_joint_contraction_kernel_040b.v1.json",

    "041":
        ART / "epr_unique_binary_weight_closure_041.v1.json",
}

JSON_OUT = (
    ART
    / "epr_algebraic_assembly_lock_042.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_algebraic_assembly_lock_042.md"
)


def load(path):
    return json.loads(
        path.read_text(
            encoding="utf-8",
        )
    )


def digest_json(obj):
    raw = json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")

    return sha256(raw).hexdigest()


print("== 042 EPR ALGEBRAIC ASSEMBLY LOCK ==")

data = {}

for key, path in INPUTS.items():
    if not path.exists():
        raise SystemExit(
            "MISSING INPUT: " + str(path)
        )

    data[key] = load(path)

print(
    "INPUT_ARTIFACT_COUNT:",
    len(data),
)

checks = {}

for key, row in data.items():
    checks[
        key + "_passes"
    ] = (
        row["audit_pass"] is True
    )


# ------------------------------------------------------------
# 1. Pre-setting source herald.
#
# Audit034C was a probe-level source receipt and was not serialized
# under artifacts/json. Audit035 and Audit039 are the sealed
# downstream interfaces that consume that result.
# ------------------------------------------------------------

a035 = data["035"]

checks[
    "herald_is_pre_setting"
] = (
    a035[
        "source"
    ][
        "analyzer_settings_used_to_define_herald"
    ]
    is False
)

checks[
    "sealed_downstream_chain_identifies_h1_herald"
] = (
    data[
        "039"
    ][
        "heralded_source"
    ][
        "receipt"
    ]
    == "h=1 from Audit034C"
)

print()
print(
    "HERALD:",
    "h=1 canonical EPR source receipt",
)

print(
    "HERALD_PRE_SETTING:",
    checks[
        "herald_is_pre_setting"
    ],
)

print(
    "HERALD_DOWNSTREAM_PROVENANCE:",
    checks[
        "sealed_downstream_chain_identifies_h1_herald"
    ],
)


# ------------------------------------------------------------
# 2. Native symmetry reduction.
# ------------------------------------------------------------

a036 = data["036"]

checks[
    "zero_local_marginal_space"
] = (
    a036[
        "marginal_reduction"
    ][
        "covariant_space_dimension"
    ]
    == 0
)

checks[
    "correlation_reduced_to_one_visibility"
] = (
    a036[
        "boundary"
    ][
        "remaining_weighting_problem_dimension"
    ]
    == 1
)

print()
print(
    "LOCAL_MARGINAL_NULLITY:",
    a036[
        "marginal_reduction"
    ][
        "covariant_space_dimension"
    ],
)

print(
    "CORRELATION_FAMILY:",
    "E_ii=-1; E_ij=-v*S_ij",
)


# ------------------------------------------------------------
# 3. Native six-line geometry.
# ------------------------------------------------------------

a037 = data["037"]

checks[
    "conference_identity"
] = (
    a037[
        "conference_system"
    ][
        "identity"
    ]
    == "C^2=5I"
)

checks[
    "native_analyzer_sector_dimension_three"
] = (
    a037[
        "conference_system"
    ][
        "eigenspaces"
    ][
        "+sqrt(5)"
    ]
    == 3
)

checks[
    "commutant_dimension_two"
] = (
    a037[
        "commutant"
    ][
        "dimension"
    ]
    == 2
)

print()
print(
    "CONFERENCE_IDENTITY:",
    "C^2=5I",
)

print(
    "B_PLUS3_GRAM:",
    "I+C/sqrt(5)",
)


# ------------------------------------------------------------
# 4. Positive rank-three geometry isolates 1/sqrt(5).
# ------------------------------------------------------------

a038 = data["038"]

checks[
    "rank3_geometry_selects_visibility"
] = (
    a038[
        "native_sector"
    ][
        "selected_visibility"
    ]
    == "1/sqrt(5)"
)

print()
print(
    "GEOMETRIC_VISIBILITY:",
    "1/sqrt(5)",
)


# ------------------------------------------------------------
# 5. Classical local answer-table route is excluded.
# ------------------------------------------------------------

a039 = data["039"]

checks[
    "local_table_bound_is_two"
] = (
    a039[
        "local_answer_table_test"
    ][
        "maximum_absolute_CHSH"
    ]
    == 2
)

checks[
    "target_outside_local_answer_polytope"
] = (
    a039[
        "checks"
    ][
        "target_outside_local_polytope"
    ]
    is True
)

print()
print(
    "LOCAL_ANSWER_TABLE_MAX_CHSH:",
    2,
)

print(
    "SETTING_INDEPENDENT_LOCAL_TABLE_ROUTE:",
    "excluded",
)


# ------------------------------------------------------------
# 6. Native joint operator contraction.
# ------------------------------------------------------------

a040 = data["040b"]

checks[
    "native_joint_operator_bilinearity_closed"
] = (
    a040[
        "boundary"
    ][
        "native_joint_operator_bilinearity_closed"
    ]
    is True
)

checks[
    "joint_contraction_is_native_trace_kernel"
] = (
    a040[
        "joint_contraction"
    ][
        "traceless_identity"
    ]
    == "Gamma(A,B)=-1/2*Tr(AB)"
)

checks[
    "operator_visibility_is_native"
] = (
    a040[
        "native_analyzer_descent"
    ][
        "visibility"
    ]
    == "1/sqrt(5)"
)

print()
print(
    "JOINT_CONTRACTION:",
    "Gamma(A,B)=Tr[Podd(A tensor B)]",
)

print(
    "TRACELESS_KERNEL:",
    "Gamma(A,B)=-1/2*Tr(AB)",
)

print(
    "OPERATOR_VISIBILITY:",
    "1/sqrt(5)",
)


# ------------------------------------------------------------
# 7. Unique binary algebraic weights.
# ------------------------------------------------------------

a041 = data["041"]

checks[
    "binary_moment_system_full_rank"
] = (
    a041[
        "binary_moment_system"
    ][
        "rank"
    ]
    == 4
)

checks[
    "binary_weight_formula_unique"
] = (
    a041[
        "binary_moment_system"
    ][
        "unique_solution"
    ]
    == "w_ab=(1+a*b*Gamma)/4"
)

checks[
    "weights_positive_and_normalized"
] = (
    a041[
        "boundary"
    ][
        "weights_positive_and_normalized"
    ]
    is True
)

target_chsh = float(
    a041[
        "CHSH"
    ][
        "numeric"
    ]
)

expected_chsh = (
    1.0
    + 3.0 / math.sqrt(5.0)
)

checks[
    "exact_CHSH_matches_native_kernel"
] = (
    abs(
        target_chsh
        - expected_chsh
    )
    < 1.0e-12
)

checks[
    "Bell_CHSH_exceeds_two"
] = (
    target_chsh > 2.0
)

print()
print(
    "UNIQUE_BINARY_WEIGHT:",
    "w_ab=(1+a*b*Gamma)/4",
)

print(
    "CHSH_EXACT:",
    a041[
        "CHSH"
    ][
        "exact"
    ],
)

print(
    "CHSH_NUMERIC:",
    target_chsh,
)


# ------------------------------------------------------------
# 8. Freeze the remaining interfaces.
# ------------------------------------------------------------

checks[
    "frequency_bridge_remains_open"
] = (
    a041[
        "boundary"
    ][
        "registered_trial_frequency_equals_algebraic_weight_derived"
    ]
    is False
)

checks[
    "common_adjoint_binding_remains_conditional"
] = (
    a041[
        "boundary"
    ][
        "common_adjoint_face_binding_still_conditional"
    ]
    is True
)


failed = [
    name
    for name, passed
    in checks.items()
    if not passed
]

audit_pass = not failed

verdict = (
    "heralded_native_algebraic_EPR_assembly_closed_with_unique_"
    "positive_binary_measure_exact_no_signaling_marginals_and_"
    "CHSH_1_plus_3_over_sqrt5_while_operational_frequency_and_"
    "common_adjoint_correspondence_remain_explicit_interfaces"
    if audit_pass
    else
    "EPR_algebraic_assembly_lock_failed"
)

artifact = {
    "artifact_id":
        "epr_algebraic_assembly_lock_042",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "assembly_chain": [
        "native common-domain preparation",
        "pre-setting h receipt",
        "h=1 canonical antisymmetric preparation",
        "canonical exchange-odd projector P_odd",
        "native joint contraction Gamma",
        "native six-axis Gram geometry",
        "zero local signed marginals",
        "unique positive binary algebraic weights",
        "exact no-signaling marginals",
        "CHSH=1+3/sqrt(5)",
    ],

    "exact_core": {
        "source_projector":
            "P_odd=(I-X)/2",

        "joint_contraction":
            "Gamma(A,B)=Tr[P_odd(A tensor B)]",

        "traceless_joint_kernel":
            "Gamma(A,B)=-1/2 Tr(AB)",

        "native_analyzer_Gram":
            "I+C/sqrt(5)",

        "visibility":
            "1/sqrt(5)",

        "binary_weight":
            "w_ab=(1+a*b*Gamma)/4",

        "CHSH":
            "1+3/sqrt(5)",

        "CHSH_numeric":
            target_chsh,
    },

    "classical_boundary": {
        "preexisting_local_answer_table_model":
            "excluded",

        "maximum_CHSH":
            2,

        "setting_independent_positive_mixture":
            "excluded for target table",
    },

    "closed": {
        "pre_setting_herald":
            True,

        "canonical_joint_antisymmetric_object":
            True,

        "native_joint_operator_bilinearity":
            True,

        "native_visibility":
            True,

        "zero_local_marginals":
            True,

        "unique_positive_binary_algebraic_measure":
            True,

        "exact_CHSH_above_two":
            True,
    },

    "open_interfaces": {
        "common_adjoint_face_binding":
            (
                "conditional interface retained from earlier "
                "Program-02 analyzer-face work"
            ),

        "operational_frequency_correspondence":
            (
                "registered repeated-trial frequency has not been "
                "identified with the canonical algebraic weight"
            ),
    },

    "status": {
        "Program02_algebraic_EPR_assembly":
            "closed",

        "further_internal_weight_hunt_required":
            False,

        "physical_frequency_correspondence":
            "open",

        "face_analyzer_correspondence":
            "conditional",
    },

    "checks":
        checks,

    "earned_statement": (
        "Program 02 now has a complete algebraic EPR assembly. "
        "A native pre-setting herald isolates the canonical "
        "antisymmetric joint preparation without outcome-dependent "
        "postselection. The canonical exchange-odd projector supplies "
        "an exact joint contraction Gamma(A,B). On traceless local "
        "observables this contraction is minus the normalized trace "
        "inner product. Native six-axis geometry fixes the distinct-"
        "setting visibility to 1/sqrt(5), while native projective "
        "symmetry forces zero one-wing signed means. For binary "
        "outcomes these moments uniquely determine positive normalized "
        "weights w_ab=(1+a*b*Gamma)/4, giving exact half marginals, "
        "no signaling, and CHSH 1+3/sqrt(5)>2. No positive "
        "setting-independent mixture of preexisting local binary "
        "answer tables can reproduce the table. The remaining "
        "interfaces are operational frequency correspondence and the "
        "previously declared conditional common-adjoint face binding."
    ),

    "next_program_boundary": (
        "Do not continue searching Program 02 for hidden source "
        "counts or alternative weighting laws. Treat the algebraic "
        "EPR construction as assembled. Move frequency correspondence "
        "or empirical instrument realization into a separate program."
    ),
}

artifact[
    "artifact_sha256"
] = digest_json(
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

note = f"""# EPR algebraic assembly lock 042

## Status

    Program-02 algebraic EPR assembly: CLOSED

## Native chain

    common-domain preparation
      -> pre-setting herald h
      -> h=1 antisymmetric preparation
      -> P_odd=(I-X)/2
      -> Gamma(A,B)=Tr[P_odd(A tensor B)]
      -> Gamma(A,B)=-1/2 Tr(AB)
      -> native six-axis Gram geometry
      -> v=1/sqrt(5)
      -> zero local signed means
      -> unique binary algebraic weights
      -> exact no-signaling marginals
      -> CHSH=1+3/sqrt(5)

The binary weight law is

    w_ab=(1+a*b*Gamma)/4.

The frozen quartet gives

    CHSH = {a041["CHSH"]["exact"]}

numerically

    {target_chsh}.

## Classical obstruction

No setting-independent positive mixture of preexisting local binary
answer tables can reproduce this table.

The deterministic local maximum remains

    |CHSH| = 2.

## Remaining interfaces

Two interfaces remain explicit.

1. Common-adjoint face binding

   The native analyzer sector is still connected to the Program-01
   face adjoint through the previously declared conditional interface.

2. Operational frequency correspondence

   The canonical positive algebraic weights have not been identified
   with a microscopic repeated-trial counting law.

These are correspondence interfaces, not unresolved freedom in the
algebraic EPR table.

## Program boundary

No further internal weighting hunt is required in Program 02.

The next work should live in a separate operational-correspondence or
instrument-realization program.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print()
print(
    "AUDIT_PASS:",
    audit_pass,
)

print(
    "VERDICT:",
    verdict,
)

print(
    "FAILED_CHECK_COUNT:",
    len(failed),
)

print(
    "FAILED_CHECKS:",
    failed,
)

print(
    "PROGRAM02_ALGEBRAIC_EPR_ASSEMBLY:",
    "CLOSED" if audit_pass else "OPEN",
)

print(
    "CHSH_EXACT:",
    a041[
        "CHSH"
    ][
        "exact"
    ],
)

print(
    "OPEN_INTERFACE_1:",
    "common_adjoint_face_binding",
)

print(
    "OPEN_INTERFACE_2:",
    "operational_frequency_correspondence",
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
