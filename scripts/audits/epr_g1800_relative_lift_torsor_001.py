#!/usr/bin/env python3

from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[2]
JSON_OUT = (
    ROOT
    / "artifacts"
    / "json"
    / "epr_g1800_relative_lift_torsor_001.v1.json"
)
NOTE_OUT = (
    ROOT
    / "notes"
    / "epr_g1800_relative_lift_torsor_001.md"
)

N = 60

# Canonical representative of the conjugacy class of free involutions
# on 60 points. Every free involution on 60 points consists of exactly
# 30 transpositions and is permutation-conjugate to this model.
def B(x):
    return x ^ 1


def diag_canonical(u, v):
    p0 = (u, v)
    p1 = (B(u), B(v))
    return min(p0, p1)


def q_local(x):
    return x // 2


def digest_json(obj):
    raw = json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("ascii")
    return sha256(raw).hexdigest()


def main():
    print("PROGRESS: 1/6 verify free involution")

    points = tuple(range(N))

    involution_ok = all(B(B(x)) == x for x in points)
    fixed_point_free = all(B(x) != x for x in points)

    b_orbits = {
        tuple(sorted((x, B(x))))
        for x in points
    }

    b_orbit_count = len(b_orbits)
    b_orbit_sizes = sorted({len(o) for o in b_orbits})

    print("PROGRESS: 2/6 construct diagonal quotient")

    classes = {}

    for u in points:
        for v in points:
            c = diag_canonical(u, v)
            classes.setdefault(c, len(classes))

    class_reps = tuple(sorted(classes))
    class_index = {
        rep: i
        for i, rep in enumerate(class_reps)
    }

    quotient_count = len(class_reps)

    def class_of(u, v):
        return class_index[diag_canonical(u, v)]

    print("PROGRESS: 3/6 construct local-address map")

    q_rows = {}
    fibers = defaultdict(list)

    for cid, (u, v) in enumerate(class_reps):
        qa = q_local(u)
        qb = q_local(v)
        q_rows[cid] = (qa, qb)
        fibers[(qa, qb)].append(cid)

    fiber_size_profile = Counter(
        len(rows)
        for rows in fibers.values()
    )

    print("PROGRESS: 4/6 construct residual joint involution")

    tau = {}

    tau_well_defined_failures = 0
    tau_same_from_other_factor_failures = 0

    for cid, (u, v) in enumerate(class_reps):
        left = class_of(B(u), v)
        right = class_of(u, B(v))

        if left != right:
            tau_same_from_other_factor_failures += 1

        # Check independence of the chosen representative.
        u2 = B(u)
        v2 = B(v)
        left_from_other_rep = class_of(B(u2), v2)

        if left != left_from_other_rep:
            tau_well_defined_failures += 1

        tau[cid] = left

    tau_involution_ok = all(
        tau[tau[cid]] == cid
        for cid in tau
    )

    tau_fixed_point_free = all(
        tau[cid] != cid
        for cid in tau
    )

    print("PROGRESS: 5/6 verify joint-only fiber law")

    tau_preserves_local_addresses = all(
        q_rows[tau[cid]] == q_rows[cid]
        for cid in tau
    )

    fiber_equals_tau_orbit_failures = 0

    for addr, rows in fibers.items():
        if len(rows) != 2:
            fiber_equals_tau_orbit_failures += 1
            continue

        a, b = sorted(rows)

        if tau[a] != b or tau[b] != a:
            fiber_equals_tau_orbit_failures += 1

    local_address_pair_count = len(fibers)

    local_pair_determines_joint_state = all(
        len(rows) == 1
        for rows in fibers.values()
    )

    print("PROGRESS: 6/6 classify")

    checks = {
        "B_is_involution": involution_ok,
        "B_is_fixed_point_free": fixed_point_free,
        "B_orbit_count_is_30": b_orbit_count == 30,
        "B_orbits_all_size_2": b_orbit_sizes == [2],
        "diagonal_quotient_count_is_1800": quotient_count == 1800,
        "local_address_pair_count_is_900":
            local_address_pair_count == 900,
        "every_local_pair_has_two_joint_lifts":
            fiber_size_profile == Counter({2: 900}),
        "tau_well_defined":
            tau_well_defined_failures == 0,
        "tau_left_equals_tau_right":
            tau_same_from_other_factor_failures == 0,
        "tau_is_involution":
            tau_involution_ok,
        "tau_is_fixed_point_free":
            tau_fixed_point_free,
        "tau_preserves_both_local_addresses":
            tau_preserves_local_addresses,
        "each_q_fiber_is_one_tau_orbit":
            fiber_equals_tau_orbit_failures == 0,
        "local_address_pair_does_not_determine_joint_state":
            not local_pair_determines_joint_state,
    }

    failed_checks = [
        name
        for name, passed in checks.items()
        if not passed
    ]

    audit_pass = not failed_checks

    artifact = {
        "artifact_id":
            "epr_g1800_relative_lift_torsor_001",
        "version": 1,
        "audit_pass": audit_pass,
        "verdict":
            (
                "g1800_diagonal_deck_quotient_has_exact_"
                "two_state_relative_lift_torsor"
                if audit_pass
                else
                "relative_lift_torsor_gate_failed"
            ),
        "proof_scope": {
            "finite_set_size": 60,
            "deck_action":
                "arbitrary free involution B on 60 points",
            "generality":
                (
                    "Every free involution on 60 points is "
                    "permutation-conjugate to the verified "
                    "30-transposition model."
                ),
            "native_binding":
                (
                    "Uses the previously earned native premise "
                    "that the selected G60 deck action B is free. "
                    "This audit does not replay that provenance."
                ),
        },
        "construction": {
            "local_cover":
                "G60 -> G30 = G60/<B>",
            "common_carrier":
                "G1800 = (G60 x G60)/diag(B)",
            "local_address_map":
                "q([u,v]) = ([u],[v])",
            "residual_joint_involution":
                "tau_AB([u,v]) = [Bu,v] = [u,Bv]",
        },
        "counts": {
            "G60_point_count": N,
            "G30_orbit_count": b_orbit_count,
            "product_point_count": N * N,
            "G1800_class_count": quotient_count,
            "local_address_pair_count":
                local_address_pair_count,
            "fiber_size_profile": {
                str(k): v
                for k, v in sorted(
                    fiber_size_profile.items()
                )
            },
        },
        "failures": {
            "tau_well_defined_failures":
                tau_well_defined_failures,
            "tau_same_from_other_factor_failures":
                tau_same_from_other_factor_failures,
            "fiber_equals_tau_orbit_failures":
                fiber_equals_tau_orbit_failures,
        },
        "checks": checks,
        "classification": {
            "locally_instantiated":
                [
                    "Alice G30 address",
                    "Bob G30 address",
                ],
            "joint_only":
                [
                    "two-state relative lift torsor over each "
                    "ordered local-address pair",
                    "residual swap tau_AB",
                ],
            "not_yet_established":
                [
                    "identification with Program-01 S1 sheet",
                    "identification with r1",
                    "Hilbert coherence",
                    "joint coherent preparation line",
                    "local measurement instrument",
                    "trial weighting law",
                    "Bell nonfactorizability",
                    "Bell violation",
                ],
        },
        "boundary": {
            "relative_torsor_is_not_a_Bell_resource": True,
            "relative_torsor_is_not_an_outcome": True,
            "tau_AB_is_not_identified_with_S1": True,
            "tau_AB_is_not_identified_with_r1": True,
            "no_probability_law_claim": True,
            "no_CHSH_claim": True,
            "no_physical_EPR_claim": True,
        },
        "next_gate": (
            "Determine whether the native local action induced "
            "by the selected deck involution B admits a "
            "Program-01 face-sheet transduction on each wing. "
            "Recover only that specific map if it is missing."
        ),
    }

    artifact["artifact_sha256"] = digest_json(artifact)

    JSON_OUT.parent.mkdir(parents=True, exist_ok=True)
    NOTE_OUT.parent.mkdir(parents=True, exist_ok=True)

    JSON_OUT.write_text(
        json.dumps(
            artifact,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="ascii",
    )

    note = """# EPR G1800 relative lift torsor 001

## Result

For any free involution B on a 60-point native carrier, form

    G30 = G60/<B>

and the diagonal quotient

    G1800 = (G60 x G60)/diag(B).

The natural local-address map is

    q([u,v]) = ([u],[v]).

It has 900 ordered local-address pairs and every fiber has exactly
two joint states.

The residual operation

    tau_AB([u,v]) = [Bu,v] = [u,Bv]

is well defined, fixed-point free, and swaps exactly the two states
inside each q-fiber.

Therefore the pair of local G30 addresses does not determine the
joint G1800 state.

The extra structure is a two-state relative lift torsor.

It is deliberately not named as an absolute hidden bit because naming
its two elements requires a basepoint choice.

## Native premise

The audit uses only the previously earned native premise that the
selected G60 deck operation B is a free involution.

It does not replay the historical derivation of that native deck action.

Every free involution on 60 points has cycle profile 2^30 and is
permutation-conjugate to the representative checked by the verifier.
The quotient and fiber theorem therefore depends only on the free
involution structure.

## EPR meaning

This is the first exact upstairs/downstairs separation.

Downstairs, the preparation exposes two local quotient addresses.

Upstairs, two distinct common states remain over the same ordered pair
of local addresses.

The residual swap tau_AB is invisible to both local quotient maps.

This proves that the common state contains relational information not
contained in the ordered pair of local addresses alone.

It does not prove Bell nonlocality.

A classical correlated source can carry exactly this kind of extra
joint datum.

## Boundary

This audit does not identify

    tau_AB = S1

and does not identify

    tau_AB = r1.

It does not yet prove that either local lift torsor is the Program-01
connection-sheet torsor.

It does not construct a coherent Hilbert preparation.

It does not define local outcomes, trial weights, probabilities,
no-signaling statistics, or CHSH.

The relative torsor is not yet a Bell resource and is not an outcome.

## Next gate

The next construction question is specific:

    How does the selected native G60 deck involution B act on the
    Program-01 local face bundle?

We need to determine whether the local lift ambiguity presented by
G1800 can be transduced into the already sealed Program-01
connection/event apparatus.

If that map is already available, use it.

If it is missing, recover or construct exactly that map and no broader
historical census.
"""

    NOTE_OUT.write_text(note, encoding="ascii")

    print("AUDIT_PASS:", audit_pass)
    print("VERDICT:", artifact["verdict"])
    print("G30_ORBIT_COUNT:", b_orbit_count)
    print("G1800_CLASS_COUNT:", quotient_count)
    print(
        "LOCAL_ADDRESS_PAIR_COUNT:",
        local_address_pair_count,
    )
    print(
        "FIBER_SIZE_PROFILE:",
        dict(sorted(fiber_size_profile.items())),
    )
    print(
        "TAU_FIXED_POINT_FREE:",
        tau_fixed_point_free,
    )
    print(
        "TAU_PRESERVES_LOCAL_ADDRESSES:",
        tau_preserves_local_addresses,
    )
    print(
        "LOCAL_PAIR_DETERMINES_JOINT_STATE:",
        local_pair_determines_joint_state,
    )
    print("FAILED_CHECK_COUNT:", len(failed_checks))
    print("FAILED_CHECKS:", failed_checks)
    print("JSON_OUT:", JSON_OUT)
    print("NOTE_OUT:", NOTE_OUT)
    print("JSON_SHA256:", sha256(JSON_OUT.read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
