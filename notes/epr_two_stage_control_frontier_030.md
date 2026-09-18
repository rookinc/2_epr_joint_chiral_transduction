# Two-stage coherent control frontier 030

The minimal structured local control family was tested using the two
native inherited generators

    G = adjacency
    L = degree - adjacency

in both orders

    exp(-i t2 L) exp(-i t1 G)

and

    exp(-i t2 G) exp(-i t1 L).

A coarse 41 x 41 grid was tested for each order over

    0 <= t1,t2 <= 20
    dt = 0.5.

Candidate count:

    3362.

Numerical response decomposition remained valid to approximately

    3.3e-15,

with no local operator norm violation.

Maximum absolute analyzer contrast:

    0.42750422948149297

at

    GL
    t1 = 17.5
    t2 = 1.0.

That response had

    abs(bias)     = 0.02305723296257267
    abs(contrast) = 0.42750422948149297
    operator norm = 0.45056146244406564.

Maximum nontrivial self CHSH-capacity upper bound:

    1.8912276527405207

at

    LG
    t1 = 0.5
    t2 = 0.5.

No sampled two-stage Alice/Bob calibration pair exceeded the Audit-028
state-independent Bell-capacity bound.

The global maximum remained the trivial zero-control detector at 2.

Conclusion:

The minimal noncommuting adjacency/Laplacian sequence does not escape
the inherited bias/contrast/sharpness obstruction on the tested grid.

No finer G/L timing search is planned for EPR assembly unless a later
result specifically motivates it.

Next target:

    native orientation-sensitive weighting/readout

rather than additional pulse scheduling.
