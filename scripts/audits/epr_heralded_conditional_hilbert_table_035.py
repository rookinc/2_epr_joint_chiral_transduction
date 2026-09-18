#!/usr/bin/env python3

from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from math import sqrt
import contextlib
import io
import json
import runpy

HERE = Path(__file__).resolve().parents[2]

A022 = (
    HERE
    / "artifacts/json"
    / "epr_native_six_axis_singlet_quartet_census_022.v1.json"
)

S022 = (
    HERE
    / "scripts/audits"
    / "epr_native_six_axis_singlet_quartet_census_022.py"
)

S034C = (
    HERE
    / "scripts/probes"
    / "probe_epr_event_ready_source_receipt_034c.py"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_heralded_conditional_hilbert_table_035.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_heralded_conditional_hilbert_table_035.md"
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


# Exact quadratic-field value
#
#     a + b sqrt(5)
#
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


def qscale(c, x):
    c = Fraction(c)

    return (
        c * x[0],
        c * x[1],
    )


def qsub(x, y):
    return qadd(
        x,
        qscale(
            -1,
            y,
        ),
    )


def qnum(x):
    return (
        float(x[0])
        + float(x[1]) * sqrt(5.0)
    )


def fstr(x):
    if x.denominator == 1:
        return str(
            x.numerator
        )

    return (
        str(x.numerator)
        + "/"
        + str(x.denominator)
    )


def qstr(x):
    a, b = x

    if b == 0:
        return fstr(a)

    if a == 0:
        if b == 1:
            return "sqrt(5)"
        if b == -1:
            return "-sqrt(5)"

        return (
            fstr(b)
            + "*sqrt(5)"
        )

    sign = "+" if b > 0 else "-"

    bb = abs(b)

    if bb == 1:
        tail = "sqrt(5)"
    else:
        tail = (
            fstr(bb)
            + "*sqrt(5)"
        )

    return (
        fstr(a)
        + " "
        + sign
        + " "
        + tail
    )


def qeq(x, y):
    return x == y


print(
    "== 035 HERALDED CONDITIONAL HILBERT TABLE =="
)

a022 = load(A022)

checks = {}

checks["Audit022_passes"] = (
    a022["audit_pass"] is True
)

checks[
    "Audit022_common_adjoint_lift_is_conditional"
] = (
    a022[
        "boundary"
    ][
        "common_adjoint_lift_is_conditional_interface"
    ]
    is True
)

checks[
    "Audit022_Born_frequency_law_not_derived"
] = (
    a022[
        "boundary"
    ][
        "Born_frequency_law_not_derived"
    ]
    is True
)

print(
    "PROGRESS: 1/7 replay exact signed analyzer geometry"
)

sink = io.StringIO()

with contextlib.redirect_stdout(sink):
    ns022 = runpy.run_path(
        str(S022)
    )

sign_matrix = ns022[
    "sector_sign_matrices"
][
    "B_plus3"
]

checks[
    "signed_matrix_is_6_by_6"
] = (
    len(sign_matrix) == 6
    and all(
        len(row) == 6
        for row in sign_matrix
    )
)

checks[
    "signed_matrix_is_symmetric"
] = all(
    sign_matrix[i][j]
    == sign_matrix[j][i]
    for i in range(6)
    for j in range(6)
)

checks[
    "signed_matrix_diagonal_is_one"
] = all(
    sign_matrix[i][i] == 1
    for i in range(6)
)

print(
    "PROGRESS: 2/7 verify pre-setting herald"
)

with contextlib.redirect_stdout(sink):
    ns034 = runpy.run_path(
        str(S034C)
    )

checks[
    "Audit034C_passes"
] = (
    ns034[
        "audit_pass"
    ]
    is True
)

checks[
    "source_herald_is_pre_setting"
] = (
    ns034[
        "checks"
    ][
        "Audit019_pre_setting_classifier"
    ]
    is True
    and ns034[
        "checks"
    ][
        "Audit019_analyzer_settings_not_used"
    ]
    is True
)

checks[
    "source_herald_lands_on_canonical_EPR_line"
] = (
    ns034[
        "checks"
    ][
        "Audit019_nonzero_is_canonical_EPR_line"
    ]
    is True
)

print(
    "PROGRESS: 3/7 freeze one native Bell-capable quartet"
)

sector = "B_plus3"

alice_settings = (
    0,
    1,
)

bob_settings = (
    0,
    2,
)

outcome_signs = (
    1,
    1,
    1,
    1,
)

max_examples = a022[
    "quartet_census"
][
    "sector_results"
][
    sector
][
    "max_examples"
]

target_example = {
    "Alice": [
        0,
        1,
    ],
    "Bob": [
        0,
        2,
    ],
    "outcome_signs": [
        1,
        1,
        1,
        1,
    ],
    "shared_axis_count":
        1,
}

checks[
    "frozen_quartet_is_Audit022_max_example"
] = (
    target_example
    in max_examples
)

print(
    "SECTOR:",
    sector,
)

print(
    "ALICE_SETTINGS:",
    alice_settings,
)

print(
    "BOB_SETTINGS:",
    bob_settings,
)

print(
    "OUTCOME_SIGNS:",
    outcome_signs,
)

print()
print(
    "SIGNED_GRAM_ROWS:"
)

for row in sign_matrix:
    print(
        " ",
        row,
    )


def correlation(i, j):
    # Conditional common-adjoint singlet correlation.
    #
    # Same axis:
    #
    #     E = -1
    #
    # Distinct native analyzer lines:
    #
    #     E = -s / sqrt(5)
    #       = -(s/5) sqrt(5),
    #
    # where s is the exact oriented Gram sign.
    if i == j:
        return q(
            -1,
            0,
        )

    s = int(
        sign_matrix[i][j]
    )

    return q(
        0,
        Fraction(
            -s,
            5,
        ),
    )


print(
    "PROGRESS: 4/7 construct exact conditional outcome tables"
)

setting_pairs = (
    (
        "A0B0",
        alice_settings[0],
        bob_settings[0],
    ),
    (
        "A0B1",
        alice_settings[0],
        bob_settings[1],
    ),
    (
        "A1B0",
        alice_settings[1],
        bob_settings[0],
    ),
    (
        "A1B1",
        alice_settings[1],
        bob_settings[1],
    ),
)

tables = {}

all_probabilities_nonnegative = True
normalization_failures = 0
marginal_failures = 0

half = q(
    Fraction(1, 2),
    0,
)

one = q(
    1,
    0,
)

for label, i, j in setting_pairs:
    E = correlation(
        i,
        j,
    )

    outcomes = {}

    for a in (
        1,
        -1,
    ):
        for b in (
            1,
            -1,
        ):
            # Standard conditional singlet projective rule:
            #
            #     P(a,b|i,j,h=1)
            #       =
            #     (1 + a b E_ij) / 4.
            p = qscale(
                Fraction(
                    1,
                    4,
                ),
                qadd(
                    one,
                    qscale(
                        a * b,
                        E,
                    ),
                ),
            )

            outcomes[
                (
                    a,
                    b,
                )
            ] = p

            if qnum(p) < -1.0e-12:
                all_probabilities_nonnegative = False

    total = q()

    for p in outcomes.values():
        total = qadd(
            total,
            p,
        )

    if not qeq(
        total,
        one,
    ):
        normalization_failures += 1

    alice_marginals = {}

    for a in (
        1,
        -1,
    ):
        m = q()

        for b in (
            1,
            -1,
        ):
            m = qadd(
                m,
                outcomes[
                    (
                        a,
                        b,
                    )
                ],
            )

        alice_marginals[
            a
        ] = m

        if not qeq(
            m,
            half,
        ):
            marginal_failures += 1

    bob_marginals = {}

    for b in (
        1,
        -1,
    ):
        m = q()

        for a in (
            1,
            -1,
        ):
            m = qadd(
                m,
                outcomes[
                    (
                        a,
                        b,
                    )
                ],
            )

        bob_marginals[
            b
        ] = m

        if not qeq(
            m,
            half,
        ):
            marginal_failures += 1

    tables[
        label
    ] = {
        "Alice_setting":
            i,

        "Bob_setting":
            j,

        "Gram_sign":
            (
                1
                if i == j
                else int(
                    sign_matrix[i][j]
                )
            ),

        "correlation_exact":
            qstr(E),

        "correlation_numeric":
            qnum(E),

        "probabilities_exact": {
            "++":
                qstr(
                    outcomes[
                        (
                            1,
                            1,
                        )
                    ]
                ),

            "+-":
                qstr(
                    outcomes[
                        (
                            1,
                            -1,
                        )
                    ]
                ),

            "-+":
                qstr(
                    outcomes[
                        (
                            -1,
                            1,
                        )
                    ]
                ),

            "--":
                qstr(
                    outcomes[
                        (
                            -1,
                            -1,
                        )
                    ]
                ),
        },

        "probabilities_numeric": {
            "++":
                qnum(
                    outcomes[
                        (
                            1,
                            1,
                        )
                    ]
                ),

            "+-":
                qnum(
                    outcomes[
                        (
                            1,
                            -1,
                        )
                    ]
                ),

            "-+":
                qnum(
                    outcomes[
                        (
                            -1,
                            1,
                        )
                    ]
                ),

            "--":
                qnum(
                    outcomes[
                        (
                            -1,
                            -1,
                        )
                    ]
                ),
        },

        "Alice_marginals_exact": {
            "+1":
                qstr(
                    alice_marginals[
                        1
                    ]
                ),

            "-1":
                qstr(
                    alice_marginals[
                        -1
                    ]
                ),
        },

        "Bob_marginals_exact": {
            "+1":
                qstr(
                    bob_marginals[
                        1
                    ]
                ),

            "-1":
                qstr(
                    bob_marginals[
                        -1
                    ]
                ),
        },
    }

checks[
    "all_conditional_probabilities_nonnegative"
] = (
    all_probabilities_nonnegative
)

checks[
    "all_four_setting_tables_normalized"
] = (
    normalization_failures == 0
)

checks[
    "all_local_marginals_exactly_one_half"
] = (
    marginal_failures == 0
)

print()
print(
    "== CONDITIONAL HILBERT TABLES =="
)

for label, row in tables.items():
    print()
    print(
        label,
        "settings=",
        (
            row[
                "Alice_setting"
            ],
            row[
                "Bob_setting"
            ],
        ),
    )

    print(
        "  Gram_sign:",
        row[
            "Gram_sign"
        ],
    )

    print(
        "  E:",
        row[
            "correlation_exact"
        ],
    )

    print(
        "  P:",
        row[
            "probabilities_exact"
        ],
    )

    print(
        "  A marginal:",
        row[
            "Alice_marginals_exact"
        ],
    )

    print(
        "  B marginal:",
        row[
            "Bob_marginals_exact"
        ],
    )

print(
    "PROGRESS: 5/7 test exact no-signaling"
)

# Every local marginal is 1/2 for every remote setting.
checks[
    "conditional_table_no_signaling_exact"
] = (
    marginal_failures == 0
)

print(
    "PROGRESS: 6/7 evaluate exact CHSH"
)

E00 = correlation(
    alice_settings[0],
    bob_settings[0],
)

E01 = correlation(
    alice_settings[0],
    bob_settings[1],
)

E10 = correlation(
    alice_settings[1],
    bob_settings[0],
)

E11 = correlation(
    alice_settings[1],
    bob_settings[1],
)

raw_chsh = qsub(
    qadd(
        qadd(
            E00,
            E01,
        ),
        E10,
    ),
    E11,
)

if qnum(
    raw_chsh
) < 0:
    chsh = qscale(
        -1,
        raw_chsh,
    )
else:
    chsh = raw_chsh

expected_chsh = q(
    1,
    Fraction(
        3,
        5,
    ),
)

checks[
    "CHSH_exactly_1_plus_3_over_sqrt5"
] = qeq(
    chsh,
    expected_chsh,
)

checks[
    "CHSH_exceeds_2"
] = (
    qnum(chsh) > 2.0
)

# Any setting-independent Bell-local factorization with binary local
# responses obeys CHSH <= 2. Therefore this conditional probability
# table does not admit such a factorization.
checks[
    "standard_Bell_local_factorization_excluded_for_control_table"
] = (
    checks[
        "CHSH_exceeds_2"
    ]
)

print()
print(
    "CHSH_EXACT:",
    qstr(
        chsh
    ),
)

print(
    "CHSH_NUMERIC:",
    qnum(
        chsh
    ),
)

print(
    "CHSH_GT_2:",
    checks[
        "CHSH_exceeds_2"
    ],
)

print(
    "PROGRESS: 7/7 classify"
)

failed_checks = [
    name
    for name, passed
    in checks.items()
    if not passed
]

audit_pass = not failed_checks

verdict = (
    "the_034C_pre_setting_EPR_herald_combined_with_the_"
    "conditional_common_adjoint_singlet_rule_and_one_frozen_"
    "native_Audit022_quartet_gives_an_exact_normalized_"
    "no_signaling_probability_table_with_CHSH_1_plus_3_over_"
    "sqrt5"
    if audit_pass
    else
    "heralded_conditional_Hilbert_table_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_heralded_conditional_hilbert_table_035",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "source": {
        "pre_setting_herald":
            "Audit034C h=1",

        "herald_formula":
            "h=(1-relative_line_relation)/2",

        "heralded_state":
            "canonical antisymmetric line [|01>-|10>]",

        "analyzer_settings_used_to_define_herald":
            False,
    },

    "conditional_interface": {
        "sector":
            sector,

        "Alice_settings":
            list(
                alice_settings
            ),

        "Bob_settings":
            list(
                bob_settings
            ),

        "outcome_signs":
            list(
                outcome_signs
            ),

        "common_adjoint_lift":
            "conditional interface inherited from Audit022",

        "local_observable":
            "O_i=n_i dot sigma",

        "probability_rule":
            "P(a,b|x,y,h=1)=(1+a*b*E_xy)/4",

        "singlet_correlation_rule":
            "E_xy=-n_x dot n_y",
    },

    "signed_Gram_matrix_B_plus3":
        sign_matrix,

    "tables":
        tables,

    "CHSH": {
        "expression":
            "E00 + E01 + E10 - E11",

        "exact":
            qstr(
                chsh
            ),

        "alternate_exact":
            "1 + 3/sqrt(5)",

        "numeric":
            qnum(
                chsh
            ),

        "exceeds_2":
            checks[
                "CHSH_exceeds_2"
            ],
    },

    "checks":
        checks,

    "boundary": {
        "pre_setting_EPR_herald_closed":
            True,

        "conditional_Hilbert_probability_control":
            True,

        "common_adjoint_lift_is_conditional":
            True,

        "native_relative_Alice_Bob_adjoint_frame_derived":
            False,

        "complete_local_face_instrument_derived":
            False,

        "Born_frequency_law_derived":
            False,

        "native_receipt_frequency_law_derived":
            False,

        "attempt_frequency_law_derived":
            False,

        "physical_source_rate_derived":
            False,

        "control_table_no_signaling":
            checks[
                "conditional_table_no_signaling_exact"
            ],

        "control_table_CHSH_exceeds_2":
            checks[
                "CHSH_exceeds_2"
            ],

        "control_table_standard_Bell_local_factorization_excluded":
            checks[
                "standard_Bell_local_factorization_excluded_for_control_table"
            ],
    },

    "earned_statement": (
        "Audit 034C supplies a setting-independent herald h=1 whose "
        "nonzero preparation branch is the canonical antisymmetric EPR "
        "line. Freezing the native Audit-022 B_plus3 quartet "
        "Alice=(0,1), Bob=(0,2), with no outcome-sign relabeling, and "
        "supplying Audit 022's conditional common-adjoint singlet "
        "interface gives four exact binary joint probability tables. "
        "Every table is positive and normalized, all Alice and Bob "
        "marginals equal one half independently of the remote setting, "
        "and the exact CHSH value is 1+3/sqrt(5)>2. This is a "
        "conditional Hilbert probability control. It does not derive "
        "the finite native receipt-frequency law."
    ),

    "next_gate": (
        "Derive the four joint outcome weights from finite native "
        "receipt mechanics on the h=1 heralded preparation domain. "
        "The target table is now fixed exactly. Do not use the target "
        "probabilities as selector input. A native construction must "
        "reproduce normalization and remote-setting-independent "
        "marginals before its CHSH value is evaluated."
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

note = f"""# Heralded conditional Hilbert table 035

## Result

Audit pass:

    {audit_pass}

The pre-setting source herald is inherited from Audit 034C.

Frozen native analyzer quartet:

    sector: B_plus3
    Alice: (0,1)
    Bob:   (0,2)

The conditional common-adjoint singlet interface gives

    E(x,y) = - n_x dot n_y

and the conditional Hilbert probability control

    P(a,b|x,y,h=1)
      =
    (1 + a b E(x,y)) / 4.

Exact CHSH:

    {qstr(chsh)}

Equivalent form:

    1 + 3/sqrt(5)

Numeric value:

    {qnum(chsh)}

All four joint tables are normalized.

Every Alice marginal is exactly

    1/2

for both Bob settings.

Every Bob marginal is exactly

    1/2

for both Alice settings.

Thus the conditional control table is exactly no-signaling.

Because its CHSH value exceeds two, the completed control table does not
admit a setting-independent Bell-local factorization under the standard
CHSH assumptions.

## Boundary

This certificate does not derive the native outcome-frequency law.

The common face-adjoint lift remains a conditional interface.

No native Alice/Bob relative adjoint frame has been derived.

No complete local face instrument has been derived.

No Born frequency law is claimed as native.

The exact table is a fixed target for the next native receipt construction.

## Next gate

Derive the four joint outcome weights from finite native receipt mechanics
on the h=1 heralded preparation domain without using the target table as
selector input.
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
    len(
        failed_checks
    ),
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
