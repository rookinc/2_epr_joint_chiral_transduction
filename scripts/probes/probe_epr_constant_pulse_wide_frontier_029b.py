#!/usr/bin/env python3

from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile
import json
import subprocess
import sys
import textwrap

HERE = Path(__file__).resolve().parents[2]

ZIP = (
    HERE
    / "evidence"
    / "native_coherent_axis_instrument.zip"
)

TMAX = 200.0
DT = 0.02
CHUNK = 100

with TemporaryDirectory(prefix="epr029b_") as td0:
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

    marker = (
        "        piv=V.T@field_num(Pi);"
    )

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
            "FAIL: could not locate axis/plane "
            "projector construction"
        )

    prefix = "\n".join(
        lines[:cut + 1]
    )

    injected = r'''
        # 029B wide constant-pulse Bell-capacity frontier.
        #
        # Only the first analyzer axis is needed because the
        # packet has already proved exact axis covariance.
        #
        # Pointer-0 all-receipt effect:
        #
        #   E0 = aa Pi + bb (I-Pi)
        #
        # Binary observable:
        #
        #   A = 2 E0 - I
        #     = bias I + contrast N
        #
        # where N = 2 Pi - I, hence
        #
        #   bias = aa + bb - 1
        #   contrast = aa - bb.

        from scipy.spatial import ConvexHull

        TMAX = __TMAX__
        DT = __DT__
        CHUNK = __CHUNK__

        count = int(round(TMAX / DT)) + 1
        times = np.linspace(
            0.0,
            TMAX,
            count,
        )

        print(
            "SCAN_GRID:",
            "0..",
            TMAX,
            "DT=",
            DT,
            "COUNT=",
            len(times),
            flush=True,
        )

        scan_rows = []
        scan_max_partition_error = 0.0
        scan_max_norm_excess = 0.0

        models = [
            (
                "adjacency",
                Gp,
                Gm,
            ),
            (
                "laplacian_control",
                Deg - Gp,
                Deg - Gm,
            ),
        ]

        for model, plus, minus in models:
            # Numerical stabilization.
            #
            # Both parity generators are Hermitian. Removing the
            # SAME scalar multiple of identity from both changes
            # vp and vm only by one common global phase and therefore
            # leaves pointer effects exactly unchanged.
            n_state = plus.shape[0]

            mean_plus = float(
                np.real(
                    plus.diagonal().sum()
                )
                / n_state
            )

            mean_minus = float(
                np.real(
                    minus.diagonal().sum()
                )
                / n_state
            )

            common_mu = (
                mean_plus
                + mean_minus
            ) / 2.0

            ident = sparse.identity(
                n_state,
                format="csr",
                dtype=np.complex128,
            )

            plus_eff = (
                plus.astype(np.complex128)
                - common_mu * ident
            )

            minus_eff = (
                minus.astype(np.complex128)
                - common_mu * ident
            )

            vp_current = W.astype(complex)
            vm_current = W.astype(complex)

            print(
                "SCAN_COMMON_SHIFT:",
                model,
                common_mu,
                "TRACE_MEAN_GAP=",
                abs(
                    mean_plus
                    - mean_minus
                ),
                flush=True,
            )

            print(
                "SCAN_MODEL_START:",
                model,
                flush=True,
            )

            for start in range(
                0,
                len(times),
                CHUNK,
            ):
                stop = min(
                    len(times),
                    start + CHUNK,
                )

                ts = times[start:stop]

                # Advance exactly one grid step from the previous
                # chunk endpoint to this chunk's first time.
                if start > 0:
                    vp_current = expm_multiply(
                        (-1j * DT) * plus_eff,
                        vp_current,
                    )

                    vm_current = expm_multiply(
                        (-1j * DT) * minus_eff,
                        vm_current,
                    )

                if len(ts) == 1:
                    vp_block = np.asarray([
                        vp_current
                    ])

                    vm_block = np.asarray([
                        vm_current
                    ])

                else:
                    span = float(
                        (len(ts) - 1)
                        * DT
                    )

                    vp_block = expm_multiply(
                        (-1j) * plus_eff,
                        vp_current,
                        start=0.0,
                        stop=span,
                        num=len(ts),
                        endpoint=True,
                    )

                    vm_block = expm_multiply(
                        (-1j) * minus_eff,
                        vm_current,
                        start=0.0,
                        stop=span,
                        num=len(ts),
                        endpoint=True,
                    )

                vp_current = np.array(
                    vp_block[-1],
                    copy=True,
                )

                vm_current = np.array(
                    vm_block[-1],
                    copy=True,
                )

                for j, t in enumerate(ts):
                    vp = vp_block[j]
                    vm = vm_block[j]

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
                            "non-finite pointer effect "
                            + model
                            + " t="
                            + repr(float(t))
                        )

                    err = maxabs(
                        eff
                        - aa * piv
                        - bb * ip
                    )

                    scan_max_partition_error = max(
                        scan_max_partition_error,
                        err,
                    )

                    bias = (
                        aa
                        + bb
                        - 1.0
                    )

                    contrast = (
                        aa
                        - bb
                    )

                    lam_axis = (
                        bias
                        + contrast
                    )

                    lam_plane = (
                        bias
                        - contrast
                    )

                    operator_norm = max(
                        abs(lam_axis),
                        abs(lam_plane),
                    )

                    scan_max_norm_excess = max(
                        scan_max_norm_excess,
                        operator_norm - 1.0,
                    )

                    p = abs(bias)
                    q = abs(contrast)

                    self_bound = 2.0 * (
                        p * p
                        + 2.0 * p * q
                        + np.sqrt(2.0)
                        * q * q
                    )

                    scan_rows.append({
                        "model":
                            model,
                        "t":
                            float(t),
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
                            operator_norm,
                        "self_CHSH_upper_bound":
                            self_bound,
                    })

                print(
                    "SCAN_PROGRESS:",
                    model,
                    stop,
                    "/",
                    len(times),
                    flush=True,
                )

            print(
                "SCAN_MODEL_DONE:",
                model,
                flush=True,
            )

        # Validate the 028 coefficient bridge at the known
        # laplacian t=0.1 calibration.
        lap01 = min(
            (
                r
                for r in scan_rows
                if r["model"]
                == "laplacian_control"
            ),
            key=lambda r:
                abs(
                    r["t"]
                    - 0.1
                ),
        )

        expected_b = 0.98049150206636
        expected_c = (
            1.4876389319962158e-05
        )

        calibration_error = max(
            abs(
                lap01["bias"]
                - expected_b
            ),
            abs(
                lap01["contrast"]
                - expected_c
            ),
        )

        # The Bell-capacity function uses only
        #
        #   p = |bias|
        #   q = |contrast|.
        #
        # Since it is bilinear in the Alice/Bob (p,q)
        # pairs, its maximum over the sampled cloud occurs
        # on the convex hull of that cloud.
        points = np.asarray([
            [
                r["abs_bias"],
                r["abs_contrast"],
            ]
            for r in scan_rows
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
            ra = scan_rows[ia]

            pA = ra["abs_bias"]
            qA = ra["abs_contrast"]

            for ib in hull_row_ix:
                rb = scan_rows[ib]

                pB = rb["abs_bias"]
                qB = rb["abs_contrast"]

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

                    best_a = ra
                    best_b = rb

        max_contrast_row = max(
            scan_rows,
            key=lambda r:
                r["abs_contrast"],
        )

        max_self_row = max(
            scan_rows,
            key=lambda r:
                r[
                    "self_CHSH_upper_bound"
                ],
        )

        model_summary = {}

        for model in (
            "adjacency",
            "laplacian_control",
        ):
            rows = [
                r
                for r in scan_rows
                if r["model"] == model
            ]

            model_summary[model] = {
                "max_abs_contrast":
                    max(
                        r["abs_contrast"]
                        for r in rows
                    ),

                "max_abs_contrast_row":
                    max(
                        rows,
                        key=lambda r:
                            r[
                                "abs_contrast"
                            ],
                    ),

                "max_self_CHSH_upper_bound":
                    max(
                        r[
                            "self_CHSH_upper_bound"
                        ]
                        for r in rows
                    ),

                "max_self_CHSH_row":
                    max(
                        rows,
                        key=lambda r:
                            r[
                                "self_CHSH_upper_bound"
                            ],
                    ),
            }

        top_contrast = sorted(
            scan_rows,
            key=lambda r:
                r["abs_contrast"],
            reverse=True,
        )[:12]

        top_self = sorted(
            scan_rows,
            key=lambda r:
                r[
                    "self_CHSH_upper_bound"
                ],
            reverse=True,
        )[:12]

        result = {
            "probe":
                "epr_constant_pulse_wide_frontier_029b",

            "grid": {
                "t_min":
                    0.0,
                "t_max":
                    TMAX,
                "dt":
                    DT,
                "count_per_model":
                    len(times),
                "model_count":
                    2,
                "axis_count":
                    1,
            },

            "interface_check": {
                "known_laplacian_t_0_1":
                    lap01,
                "expected_bias":
                    expected_b,
                "expected_contrast":
                    expected_c,
                "max_error":
                    calibration_error,
            },

            "numerical_checks": {
                "max_axis_plane_partition_error":
                    scan_max_partition_error,
                "max_operator_norm_excess":
                    scan_max_norm_excess,
            },

            "model_summary":
                model_summary,

            "global_max_abs_contrast":
                max_contrast_row,

            "global_max_self_CHSH_upper_bound":
                max_self_row,

            "convex_hull_vertex_count":
                len(
                    hull_row_ix
                ),

            "sampled_pairwise_CHSH_capacity": {
                "maximum_upper_bound":
                    best_bound,
                "Alice":
                    best_a,
                "Bob":
                    best_b,
                "exceeds_2":
                    best_bound
                    > 2.0
                    + 5.0e-10,
            },

            "top_contrast_rows":
                top_contrast,

            "top_self_bound_rows":
                top_self,

            "boundary": {
                "continuous_all_t_no_go_proven":
                    False,
                "wide_grid_only":
                    True,
                "joint_state_used":
                    False,
                "Born_rule_derived":
                    False,
                "axis_covariance_inherited":
                    True,
            },
        }

        out = (
            ROOT
            / "SCAN_029B.json"
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
        print(
            "INTERFACE_CALIBRATION_ERROR:",
            calibration_error,
        )
        print(
            "MAX_PARTITION_ERROR:",
            scan_max_partition_error,
        )
        print(
            "MAX_OPERATOR_NORM_EXCESS:",
            scan_max_norm_excess,
        )
        print(
            "MAX_ABS_CONTRAST:",
            max_contrast_row,
        )
        print(
            "MAX_SELF_CHSH_UPPER_BOUND:",
            max_self_row,
        )
        print(
            "PAIRWISE_CHSH_CAPACITY_MAX:",
            best_bound,
        )
        print(
            "PAIRWISE_CHSH_CAPACITY_GT2:",
            best_bound
            > 2.0
            + 5.0e-10,
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
        .replace(
            "__CHUNK__",
            repr(CHUNK),
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
        / "scan_029b_runner.py"
    )

    runner.write_text(
        runner_source,
        encoding="utf-8",
    )

    print(
        "RUNNING WIDE CONSTANT-PULSE FRONTIER",
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
            "FAIL: scan runner returned "
            + str(run.returncode)
        )

    result = json.loads(
        (
            root
            / "SCAN_029B.json"
        ).read_text(
            encoding="utf-8",
        )
    )

    print()
    print("== 029B SUMMARY ==")
    print(
        json.dumps(
            result,
            indent=2,
            sort_keys=True,
        )
    )
