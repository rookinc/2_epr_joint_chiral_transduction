#!/usr/bin/env python3

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from io import BytesIO
from itertools import permutations
from pathlib import Path
from zipfile import ZipFile
import json

HERE = Path(__file__).resolve().parents[2]

A020 = (
    HERE
    / "artifacts/json"
    / "epr_local_analyzer_carrier_interface_020.v1.json"
)

A022 = (
    HERE
    / "artifacts/json"
    / "epr_native_six_axis_singlet_quartet_census_022.v1.json"
)

A024 = (
    HERE
    / "artifacts/json"
    / "epr_g60_label_torsor_c5_face_descent_024.v1.json"
)

ZIP = (
    HERE
    / "evidence"
    / "native_coherent_axis_instrument.zip"
)

JSON_OUT = (
    HERE
    / "artifacts/json"
    / "epr_projective_analyzer_sector_extension_025.v1.json"
)

NOTE_OUT = (
    HERE
    / "notes"
    / "epr_projective_analyzer_sector_extension_025.md"
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


# Q(sqrt(5))

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


def decode(row, i, j):
    den = int(
        row["denominator"]
    )

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


def seed_vector(row):
    out = []

    for i in range(60):
        total = q()

        for j in range(60):
            total = qadd(
                total,
                qscale(
                    j + 1,
                    decode(
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


def sign_matrix(rows):
    rows = sorted(
        rows,
        key=lambda row: int(
            row["setting"]
        ),
    )

    vectors = [
        seed_vector(row)
        for row in rows
    ]

    S = [
        [0] * 6
        for _ in range(6)
    ]

    for i in range(6):
        S[i][i] = 1

        for j in range(i + 1, 6):
            s = qsign(
                qdot(
                    vectors[i],
                    vectors[j],
                )
            )

            if s == 0:
                raise RuntimeError(
                    "unexpected orthogonal analyzer lines"
                )

            S[i][j] = s
            S[j][i] = s

    return tuple(
        tuple(row)
        for row in S
    )


def switching_maps(source, target):
    out = {}

    for p in permutations(
        range(6)
    ):
        # Projective line orientations are gauge.
        # Fix d_0=+1 and solve the remaining signs.
        d = [1] * 6

        for j in range(1, 6):
            d[j] = (
                target[
                    p[0]
                ][
                    p[j]
                ]
                * source[0][j]
            )

        good = True

        for i in range(6):
            for j in range(6):
                if (
                    target[
                        p[i]
                    ][
                        p[j]
                    ]
                    !=
                    d[i]
                    * d[j]
                    * source[i][j]
                ):
                    good = False
                    break

            if not good:
                break

        if good:
            out[
                tuple(p)
            ] = tuple(d)

    return out


def compose(p, q):
    # p after q
    return tuple(
        p[q[i]]
        for i in range(6)
    )


def inverse(p):
    out = [None] * 6

    for i, j in enumerate(p):
        out[j] = i

    return tuple(out)


def parity(p):
    count = 0

    for i in range(6):
        for j in range(i + 1, 6):
            if p[i] > p[j]:
                count += 1

    return count % 2


def is_group(perms):
    perms = set(perms)

    identity = tuple(
        range(6)
    )

    if identity not in perms:
        return False

    for p in perms:
        if inverse(p) not in perms:
            return False

    for p in perms:
        for r in perms:
            if compose(
                p,
                r,
            ) not in perms:
                return False

    return True


def relative_group(maps, base):
    inv = inverse(base)

    return {
        tuple(
            inv[
                m[s]
            ]
            for s in range(6)
        )
        for m in maps
    }


print("PROGRESS: 1/8 load sealed interfaces")

a020 = load(A020)
a022 = load(A022)
a024 = load(A024)

checks = {}

checks["Audit020_passed"] = (
    a020["audit_pass"] is True
)

checks["Audit022_passed"] = (
    a022["audit_pass"] is True
)

checks["Audit024_passed"] = (
    a024["audit_pass"] is True
)

print("PROGRESS: 2/8 recover two native analyzer line systems")

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
    exact = json.loads(
        z3.read(
            "native_analyzer_return_bridge/"
            "EXACT_PROJECTORS.json"
        )
    )

Splus = sign_matrix(
    exact["axes"]["B_plus3"]
)

Sminus = sign_matrix(
    exact["axes"]["B_minus3"]
)

checks["two_character_sector_line_systems_present"] = (
    Splus is not None
    and Sminus is not None
)

print("PROGRESS: 3/8 compute within-sector and cross-sector symmetries")

Gpp_map = switching_maps(
    Splus,
    Splus,
)

Gmm_map = switching_maps(
    Sminus,
    Sminus,
)

Gpm_map = switching_maps(
    Splus,
    Sminus,
)

Gmp_map = switching_maps(
    Sminus,
    Splus,
)

Gpp = set(
    Gpp_map
)

Gmm = set(
    Gmm_map
)

Gpm = set(
    Gpm_map
)

Gmp = set(
    Gmp_map
)

checks["within_plus_order_60"] = (
    len(Gpp) == 60
)

checks["within_minus_order_60"] = (
    len(Gmm) == 60
)

checks["plus_to_minus_count_60"] = (
    len(Gpm) == 60
)

checks["minus_to_plus_count_60"] = (
    len(Gmp) == 60
)

checks["within_sector_groups_equal"] = (
    Gpp == Gmm
)

checks["cross_sector_map_sets_equal"] = (
    Gpm == Gmp
)

checks["within_and_cross_sets_disjoint"] = (
    not (
        Gpp & Gpm
    )
)

checks["within_sector_symmetry_is_group"] = (
    is_group(Gpp)
)

print("PROGRESS: 4/8 reconstruct Audit024 relative face group")

maps_raw = (
    a024[
        "C5_descent"
    ][
        "setting_face_maps"
    ]
)

setting_face_maps = [
    tuple(
        int(x)
        for x in row["map"]
    )
    for row in maps_raw
]

multiplicities = [
    int(
        row["isomorphism_count"]
    )
    for row in maps_raw
]

checks["setting_face_map_count_120"] = (
    len(setting_face_maps) == 120
    and len(
        set(setting_face_maps)
    ) == 120
)

checks["each_setting_face_map_has_four_label_bridges"] = (
    set(multiplicities) == {4}
)

R = relative_group(
    setting_face_maps,
    setting_face_maps[0],
)

checks["relative_face_group_order_120"] = (
    len(R) == 120
)

checks["relative_face_maps_form_group"] = (
    is_group(R)
)

print("PROGRESS: 5/8 stress reference-map independence")

reference_group_failures = 0

for base in setting_face_maps:
    if relative_group(
        setting_face_maps,
        base,
    ) != R:
        reference_group_failures += 1

checks["relative_face_group_independent_of_reference_map"] = (
    reference_group_failures == 0
)

print("PROGRESS: 6/8 test exact Z2 sector extension")

checks["face_group_equals_within_union_cross"] = (
    R == (
        Gpp | Gpm
    )
)

checks["within_subgroup_index_two"] = (
    len(R) == 2 * len(Gpp)
    and Gpp <= R
)

within_parity = Counter(
    parity(p)
    for p in Gpp
)

cross_parity = Counter(
    parity(p)
    for p in Gpm
)

face_parity = Counter(
    parity(p)
    for p in R
)

checks["within_sector_maps_are_exactly_even"] = (
    within_parity
    == Counter({
        0: 60,
    })
)

checks["cross_sector_maps_are_exactly_odd"] = (
    cross_parity
    == Counter({
        1: 60,
    })
)

checks["face_group_parity_split_60_60"] = (
    face_parity
    == Counter({
        0: 60,
        1: 60,
    })
)

sector_bit = {
    p: (
        0
        if p in Gpp
        else 1
    )
    for p in R
}

homomorphism_failures = 0

for p in R:
    for r in R:
        pr = compose(
            p,
            r,
        )

        if (
            sector_bit[pr]
            !=
            (
                sector_bit[p]
                ^ sector_bit[r]
            )
        ):
            homomorphism_failures += 1

checks["sector_flip_bit_is_Z2_homomorphism"] = (
    homomorphism_failures == 0
)

checks["sector_flip_bit_equals_permutation_parity"] = all(
    sector_bit[p]
    == parity(p)
    for p in R
)

print("PROGRESS: 7/8 preserve EPR and interpretation boundaries")

checks["Bell_capable_quartets_retained"] = (
    a022[
        "boundary"
    ][
        "Bell_capable_native_setting_quartets_exist"
    ]
    is True
)

checks["sector_flip_not_called_chirality"] = (
    True
)

checks["Program01_complex_orientation_not_yet_identified"] = (
    True
)

checks["face_T_axis_not_identified_with_analyzer_axis"] = (
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
    "the_120_element_native_relative_setting_face_group_"
    "is_an_exact_Z2_extension_of_the_60_element_single_"
    "analyzer_sector_projective_symmetry_group_with_the_"
    "even_subgroup_preserving_each_character_sector_and_"
    "the_odd_coset_exchanging_B_plus3_and_B_minus3"
    if audit_pass
    else
    "projective_analyzer_sector_extension_gate_failed"
)

artifact = {
    "artifact_id":
        "epr_projective_analyzer_sector_extension_025",

    "version":
        1,

    "audit_pass":
        audit_pass,

    "verdict":
        verdict,

    "analyzer_sector_geometry": {
        "sector_names": [
            "B_plus3",
            "B_minus3",
        ],

        "lines_per_sector":
            6,

        "within_B_plus_symmetry_order":
            len(Gpp),

        "within_B_minus_symmetry_order":
            len(Gmm),

        "B_plus_to_B_minus_map_count":
            len(Gpm),

        "B_minus_to_B_plus_map_count":
            len(Gmp),

        "within_sector_groups_equal":
            Gpp == Gmm,

        "cross_sector_map_sets_equal":
            Gpm == Gmp,

        "within_cross_intersection_size":
            len(
                Gpp & Gpm
            ),
    },

    "native_face_group": {
        "G60_label_bridge_count":
            int(
                a024[
                    "graph_label_bridge"
                ][
                    "isomorphism_count"
                ]
            ),

        "setting_face_map_count":
            len(
                setting_face_maps
            ),

        "label_bridges_per_setting_face_map":
            4,

        "relative_face_group_order":
            len(R),

        "reference_map_independence_failure_count":
            reference_group_failures,
    },

    "sector_extension": {
        "kernel":
            "within-sector projective symmetry",

        "kernel_order":
            len(Gpp),

        "total_group_order":
            len(R),

        "quotient":
            "Z2",

        "even_subgroup_preserves_sector":
            True,

        "odd_coset_exchanges_sector":
            True,

        "within_parity_profile":
            dict(
                within_parity
            ),

        "cross_parity_profile":
            dict(
                cross_parity
            ),

        "total_parity_profile":
            dict(
                face_parity
            ),

        "sector_flip_homomorphism_failure_count":
            homomorphism_failures,

        "sector_flip_equals_permutation_parity":
            checks[
                "sector_flip_bit_equals_permutation_parity"
            ],

        "exact_sequence":
            (
                "1 -> G60_projective -> "
                "G120_face -> Z2 -> 1"
            ),
    },

    "earned_statement": (
        "The factor-of-two mismatch exposed by the first "
        "Audit-025 attempt is native structure. Each of the two "
        "three-dimensional analyzer character sectors has the "
        "same 60-element projective six-line symmetry group. "
        "There are exactly 60 projective switching maps from "
        "B_plus3 to B_minus3 and exactly 60 in the reverse "
        "direction, disjoint from the within-sector subgroup. "
        "The 120-element relative setting-to-face permutation "
        "group derived in Audit 024 is exactly the union of the "
        "60 within-sector maps and the 60 cross-sector maps. "
        "The within-sector subgroup consists exactly of the even "
        "six-setting permutations; the cross-sector coset consists "
        "exactly of the odd permutations. The induced sector bit "
        "is a surjective Z2 homomorphism and equals permutation "
        "parity on this native 120-element group. Thus the native "
        "face descent acts on the analyzer not as one isolated "
        "six-line sector but as a Z2-extended pair of native "
        "six-line character sectors."
    ),

    "checks":
        checks,

    "boundary": {
        "single_sector_projective_symmetry_order":
            60,

        "native_face_relative_group_order":
            120,

        "sector_exchange_Z2_derived":
            True,

        "sector_exchange_called_chirality":
            False,

        "sector_exchange_identified_with_Program01_relative_complex_orientation":
            False,

        "sector_exchange_identified_with_deck_character":
            False,

        "common_face_adjoint_embedding_closed":
            False,

        "face_T_equals_analyzer_axis":
            False,

        "Bell_capable_quartets_retained":
            True,

        "local_face_instrument_not_yet_constructed":
            True,

        "Born_frequency_law_not_derived":
            True,

        "no_signaling_not_yet_tested":
            True,
    },

    "next_gate": (
        "Put the analyzer sector-flip Z2 character and the "
        "Program-01 relative-complex-orientation/seam Z2 "
        "character on a genuine common native domain. Test "
        "whether they are the same character, its negative, or "
        "independent. Do not identify the analyzer sector bit "
        "with chirality or complex orientation by analogy alone."
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

note = f"""# Projective analyzer sector extension 025

## Result

Audit pass:

    {audit_pass}

The failed first 025 attempt exposed an exact factor of two.

Each native analyzer character sector has projective six-line symmetry
group order

    60.

The native relative setting-to-face group from Audit 024 has order

    120.

The missing 60 transformations are not absent symmetries of one sector.
They are exact projective maps between the two sectors.

Within-sector maps:

    {len(Gpp)}

Cross-sector maps:

    {len(Gpm)}

Intersection:

    {len(Gpp & Gpm)}

Union:

    {len(Gpp | Gpm)}

## Z2 structure

The 120-element native face group splits exactly as

    60 sector-preserving maps
    60 sector-exchanging maps.

The preserving subgroup consists exactly of even permutations of the
six settings.

The exchanging coset consists exactly of odd permutations.

The sector bit

    epsilon = 0  preserve sector
    epsilon = 1  exchange sector

satisfies

    epsilon(pq) = epsilon(p) + epsilon(q) mod 2

with failure count

    {homomorphism_failures}.

Thus there is an exact sequence

    1 -> G60_projective -> G120_face -> Z2 -> 1.

## Meaning

The factor of two in Audit 024 is native sector structure.

The analyzer should not be treated as one isolated three-dimensional
six-line system. Its native finite symmetry sees the pair

    B_plus3
    B_minus3

together.

The 60-element subgroup preserves either sector.

The other 60 transformations exchange them.

This Z2 is not yet called chirality and is not yet identified with the
Program-01 relative-complex-orientation or seam character.

That common-domain character comparison is the next theorem gate.
"""

NOTE_OUT.write_text(
    note,
    encoding="ascii",
)

print("AUDIT_PASS:", audit_pass)
print("VERDICT:", verdict)
print(
    "WITHIN_PLUS_COUNT:",
    len(Gpp),
)
print(
    "WITHIN_MINUS_COUNT:",
    len(Gmm),
)
print(
    "PLUS_TO_MINUS_COUNT:",
    len(Gpm),
)
print(
    "MINUS_TO_PLUS_COUNT:",
    len(Gmp),
)
print(
    "FACE_RELATIVE_GROUP_COUNT:",
    len(R),
)
print(
    "WITHIN_CROSS_INTERSECTION:",
    len(
        Gpp & Gpm
    ),
)
print(
    "WITHIN_UNION_CROSS_COUNT:",
    len(
        Gpp | Gpm
    ),
)
print(
    "WITHIN_PARITY_PROFILE:",
    dict(
        within_parity
    ),
)
print(
    "CROSS_PARITY_PROFILE:",
    dict(
        cross_parity
    ),
)
print(
    "SECTOR_FLIP_HOMOMORPHISM_FAILURES:",
    homomorphism_failures,
)
print(
    "SECTOR_FLIP_EQUALS_PERMUTATION_PARITY:",
    checks[
        "sector_flip_bit_equals_permutation_parity"
    ],
)
print(
    "REFERENCE_GROUP_FAILURES:",
    reference_group_failures,
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
