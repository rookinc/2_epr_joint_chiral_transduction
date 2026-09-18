#!/usr/bin/env python3

from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parents[2]

ZIP = (
    HERE
    / "evidence"
    / "native_coherent_axis_instrument.zip"
)

TMAX = 20.0
DT = 0.5

with TemporaryDirectory(prefix="epr030_") as td0:
    td = Path(td0)

    with ZipFile(ZIP) as z:
        z.extractall(td)

    root = (
        td
        / "native_coherent_axis_instrument"
    )

    verify = (
        root
        / "verify.py"
    ).read_text(
        encoding="utf-8",
    )

    lines = verify.splitlines()

    cut = None

    for i, line in enumerate(lines):
        if (
            "piv=V.T@field_num(Pi)@V;"
            in line
            and "ip=np.eye(6)-piv"
            in line
        ):
            cut = i
            break

    if cut is None:
        raise SystemExit(
            "FAIL: projector construction not found"
        )

    prefix = "\n".join(
        lines[:cut + 1]
    )

    injected = r'''
        # 030 minimal two-stage local coherent transducer.
        #
        # We use only axis 0 because exact axis covariance
        # was already established upstream.
        #
        # Two native constant-control generators:
        #
        #   G = adjacency
        #   L = degree - adjacency
        #
        # We test both noncommuting orderings:
        #
        #   LG = exp(-i t2 L) exp(-i t1 G)
        #   GL = exp(-i t2 G) exp(-i t1 L)
        #
        # The goal is not global optimization.
        #
        # The goal is to find whether this minimal structured
        # control family can escape the Audit-028
        # state-independent Bell-capacity obstruction.

        from scipy.spatial import ConvexHull

        TMAX = __TMAX__
        DT = __DT__

        times = np.arange(
            0.0,
            TMAX + 0.5 * DT,
            DT,
        )

        print(
            "GRID:",
            len(times),
            "x",
            len(times),
            "TMAX=",
            TMAX,
            "DT=",
            DT,
            flush=True,
        )

        Gplus = Gp.astype(
            np.complex128
        )

        Gminus = Gm.astype(
            np.complex128
        )

        Lplus = (
            Deg - Gp
        ).astype(
            np.complex128
        )

        Lminus = (
            Deg - Gm
        ).astype(
            np.complex128
        )

        # Remove one common scalar from each plus/minus
        # generator pair. This changes both parity branches
        # only by a common global phase and therefore leaves
        # pointer effects unchanged.

        def centered_pair(
            plus,
            minus,
        ):
            n = plus.shape[0]

            mp = float(
                np.real(
                    plus.diagonal().sum()
                )
                / n
            )

            mm = float(
                np.real(
                    minus.diagonal().sum()
                )
                / n
            )

            mu = (
                mp + mm
            ) / 2.0

            ident = sparse.identity(
                n,
                format="csr",
                dtype=np.complex128,
            )

            return (
                plus - mu * ident,
                minus - mu * ident,
                mu,
                abs(mp - mm),
            )

        Gplus, Gminus, muG, gapG = (
            centered_pair(
                Gplus,
                Gminus,
            )
        )

        Lplus, Lminus, muL, gapL = (
            centered_pair(
                Lplus,
                Lminus,
            )
        )

        print(
            "COMMON_SHIFT_G:",
            muG,
            "TRACE_MEAN_GAP:",
            gapG,
            flush=True,
        )

        print(
            "COMMON_SHIFT_L:",
            muL,
            "TRACE_MEAN_GAP:",
            gapL,
            flush=True,
        )

        def evolve_grid(
            H,
            initial,
        ):
            if len(times) == 1:
                return np.asarray([
                    initial
                ])

            return expm_multiply(
                (-1j) * H,
                initial,
                start=0.0,
                stop=float(times[-1]),
                num=len(times),
                endpoint=True,
            )

        # Stage-one grids.
        #
        # These are reused for every second-stage duration.

        print(
            "PRECOMPUTE_STAGE1_G",
            flush=True,
        )

        G1_plus = evolve_grid(
            Gplus,
            W.astype(complex),
        )

        G1_minus = evolve_grid(
            Gminus,
            W.astype(complex),
        )

        print(
            "PRECOMPUTE_STAGE1_L",
            flush=True,
        )

        L1_plus = evolve_grid(
            Lplus,
            W.astype(complex),
        )

        L1_minus = evolve_grid(
            Lminus,
            W.astype(complex),
        )

        rows = []

        max_partition_error = 0.0
        max_norm_excess = 0.0

        def record(
            order,
            t1,
            t2,
            vp,
            vm,
        ):
            nonlocal max_partition_error
            nonlocal max_norm_excess

            y0 = (
                vp + vm
            ) / 2.0

            eff = (
                y0.conj().T
                @ y0
            )

            aa = float(
                np.trace(
                    piv @ eff
                ).real
                / 2.0
            )

            bb = float(
                np.trace(
                    ip @ eff
                ).real
                / 4.0
            )

            if not (
                np.isfinite(aa)
                and np.isfinite(bb)
            ):
                raise RuntimeError(
                    "nonfinite effect "
                    + order
                    + " "
                    + repr(t1)
                    + " "
                    + repr(t2)
                )

            err = maxabs(
                eff
                - aa * piv
                - bb * ip
            )

            max_partition_error = max(
                max_partition_error,
                err,
            )

            bias = (
                aa + bb - 1.0
            )

            contrast = (
                aa - bb
            )

            lam_axis = (
                bias + contrast
            )

            lam_plane = (
                bias - contrast
            )

            norm = max(
                abs(lam_axis),
                abs(lam_plane),
            )

            max_norm_excess = max(
                max_norm_excess,
                norm - 1.0,
            )

            p = abs(bias)
            q = abs(contrast)

            self_bound = 2.0 * (
                p * p
                + 2.0 * p * q
                + np.sqrt(2.0)
                * q * q
            )

            rows.append({
                "order":
                    order,
                "t1":
                    float(t1),
                "t2":
                    float(t2),
                "pointer0_axis":
                    aa,
                "pointer0_plane":
                    bb,
                "bias":
                    bias,
                "contrast":
                    contrast,
                "abs_bias":
                    p,
                "abs_contrast":
                    q,
                "operator_norm":
                    norm,
                "self_CHSH_upper_bound":
                    float(
                        self_bound
                    ),
            })

        # --------------------------------------------------
        # ORDER LG:
        #
        # first G for t1,
        # then L for t2.
        # --------------------------------------------------

        print(
            "SCAN_ORDER: LG",
            flush=True,
        )

        for i, t1 in enumerate(times):
            vp0 = G1_plus[i]
            vm0 = G1_minus[i]

            vp_grid = evolve_grid(
                Lplus,
                vp0,
            )

            vm_grid = evolve_grid(
                Lminus,
                vm0,
            )

            for j, t2 in enumerate(times):
                record(
                    "LG",
                    t1,
                    t2,
                    vp_grid[j],
                    vm_grid[j],
                )

            print(
                "PROGRESS_LG:",
                i + 1,
                "/",
                len(times),
                flush=True,
            )

        # --------------------------------------------------
        # ORDER GL:
        #
        # first L for t1,
        # then G for t2.
        # --------------------------------------------------

        print(
            "SCAN_ORDER: GL",
            flush=True,
        )

        for i, t1 in enumerate(times):
            vp0 = L1_plus[i]
            vm0 = L1_minus[i]

            vp_grid = evolve_grid(
                Gplus,
                vp0,
            )

            vm_grid = evolve_grid(
                Gminus,
                vm0,
            )

            for j, t2 in enumerate(times):
                record(
                    "GL",
                    t1,
                    t2,
                    vp_grid[j],
                    vm_grid[j],
                )

            print(
                "PROGRESS_GL:",
                i + 1,
                "/",
                len(times),
                flush=True,
            )

        # Bell capacity depends only on
        #
        #   p = |bias|
        #   q = |contrast|.
        #
        # Its maximum over the finite cloud is attained on
        # the convex hull of the cloud in the (p,q) plane.

        points = np.asarray([
            [
                row["abs_bias"],
                row["abs_contrast"],
            ]
            for row in rows
        ])

        rounded = np.round(
            points,
            13,
        )

        unique_points, unique_ix = np.unique(
            rounded,
            axis=0,
            return_index=True,
        )

        if len(unique_points) <= 2:
            hull_local_ix = np.arange(
                len(unique_points)
            )
        else:
            try:
                hull = ConvexHull(
                    unique_points
                )

                hull_local_ix = (
                    hull.vertices
                )

            except Exception:
                hull_local_ix = np.arange(
                    len(unique_points)
                )

        hull_row_ix = [
            int(
                unique_ix[i]
            )
            for i in hull_local_ix
        ]

        best_bound = -1.0
        best_a = None
        best_b = None

        for ia in hull_row_ix:
            A = rows[ia]

            pA = A["abs_bias"]
            qA = A["abs_contrast"]

            for ib in hull_row_ix:
                B = rows[ib]

                pB = B["abs_bias"]
                qB = B["abs_contrast"]

                bound = 2.0 * (
                    pA * pB
                    + pA * qB
                    + qA * pB
                    + np.sqrt(2.0)
                    * qA * qB
                )

                if bound > best_bound:
                    best_bound = float(
                        bound
                    )

                    best_a = A
                    best_b = B

        max_contrast = max(
            rows,
            key=lambda r:
                r["abs_contrast"],
        )

        max_norm = max(
            rows,
            key=lambda r:
                r["operator_norm"],
        )

        max_self = max(
            rows,
            key=lambda r:
                r[
                    "self_CHSH_upper_bound"
                ],
        )

        nontrivial_rows = [
            r
            for r in rows
            if (
                abs(
                    r["t1"]
                )
                > 1e-12
                or abs(
                    r["t2"]
                )
                > 1e-12
            )
        ]

        max_nontrivial_self = max(
            nontrivial_rows,
            key=lambda r:
                r[
                    "self_CHSH_upper_bound"
                ],
        )

        gate_pass = (
            best_bound
            > 2.0
            + 5.0e-10
        )

        result = {
            "probe":
                "epr_two_stage_bell_capacity_030",

            "grid": {
                "t_min":
                    0.0,
                "t_max":
                    TMAX,
                "dt":
                    DT,
                "count_per_axis":
                    len(times),
                "candidate_count":
                    len(rows),
                "orders": [
                    "LG",
                    "GL",
                ],
            },

            "control": {
                "G":
                    "adjacency",
                "L":
                    "degree_minus_adjacency",
                "LG":
                    "exp(-i t2 L) exp(-i t1 G)",
                "GL":
                    "exp(-i t2 G) exp(-i t1 L)",
            },

            "numerical_checks": {
                "max_axis_plane_partition_error":
                    max_partition_error,
                "max_operator_norm_excess":
                    max_norm_excess,
            },

            "max_abs_contrast":
                max_contrast,

            "max_operator_norm":
                max_norm,

            "max_self_CHSH_upper_bound":
                max_self,

            "max_nontrivial_self_CHSH_upper_bound":
                max_nontrivial_self,

            "pairwise_capacity": {
                "maximum_upper_bound":
                    best_bound,
                "Alice":
                    best_a,
                "Bob":
                    best_b,
                "exceeds_2":
                    gate_pass,
                "hull_vertex_count":
                    len(
                        hull_row_ix
                    ),
            },

            "boundary": {
                "Bell_violation_proven":
                    False,
                "Bell_capacity_gate_passed":
                    gate_pass,
                "joint_state_used":
                    False,
                "joint_probability_law_used":
                    False,
                "Born_rule_derived":
                    False,
                "coarse_grid_only":
                    True,
                "global_two_stage_optimum_proven":
                    False,
            },
        }

        out = (
            ROOT
            / "SCAN_030.json"
        )

        out.write_text(
            json.dumps(
                result,
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )

        print()
        print("== 030 SUMMARY ==")

        print(
            "CANDIDATE_COUNT:",
            len(rows),
        )

        print(
            "MAX_PARTITION_ERROR:",
            max_partition_error,
        )

        print(
            "MAX_OPERATOR_NORM_EXCESS:",
            max_norm_excess,
        )

        print(
            "MAX_ABS_CONTRAST:",
            max_contrast,
        )

        print(
            "MAX_OPERATOR_NORM:",
            max_norm,
        )

        print(
            "MAX_SELF_CHSH_UPPER_BOUND:",
            max_self,
        )

        print(
            "MAX_NONTRIVIAL_SELF_CHSH_UPPER_BOUND:",
            max_nontrivial_self,
        )

        print(
            "PAIRWISE_CHSH_CAPACITY_MAX:",
            best_bound,
        )

        print(
            "PAIRWISE_CHSH_CAPACITY_GT2:",
            gate_pass,
        )

        print(
            "BEST_ALICE:",
            best_a,
        )

        print(
            "BEST_BOB:",
            best_b,
        )

        print(
            "HULL_VERTEX_COUNT:",
            len(
                hull_row_ix
            ),
        )

        return
'''

    injected = (
        injected
        .replace(
            "__TMAX__",
            repr(TMAX),
        )
        .replace(
            "__DT__",
            repr(DT),
        )
    )

    runner_source = (
        prefix
        + "\n"
        + injected
        + "\n\n"
        + "if __name__ == '__main__':\n"
        + "    main()\n"
    )

    runner = (
        root
        / "scan_030_runner.py"
    )

    runner.write_text(
        runner_source,
        encoding="utf-8",
    )

    print(
        "RUNNING 030 TWO-STAGE BELL-CAPACITY SEARCH",
        flush=True,
    )

    run = subprocess.run(
        [
            sys.executable,
            str(runner),
            "--skip-upstream-replay",
        ],
        cwd=root,
    )

    if run.returncode != 0:
        raise SystemExit(
            "FAIL: 030 runner returned "
            + str(run.returncode)
        )

    result = json.loads(
        (
            root
            / "SCAN_030.json"
        ).read_text(
            encoding="utf-8",
        )
    )

    print()
    print("== 030 RESULT JSON ==")

    print(
        json.dumps(
            result,
            indent=2,
            sort_keys=True,
        )
    )
