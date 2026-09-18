# Constant-pulse Bell-capacity frontier 029

A wide numerical frontier scan tested adjacency and Laplacian constant
controls over

    0 <= t <= 200

with grid spacing

    dt = 0.02.

The inherited exact response form is

    A_i(t) = bias(t) I + contrast(t) N_i

with

    N_i^2 = I.

The scan reproduced the Audit-028 calibration at Laplacian t=0.1 to
approximately 7.1e-15 and maintained the exact axis/plane response form
to approximately 1.1e-14. No local observable exceeded operator norm 1.

The wide scan found substantially more analyzer contrast than the
original eight sampled pulses.

Maximum absolute contrast:

    0.431782125555502

at

    laplacian_control
    t = 158.08

with

    bias     = 0.0813468201188412
    contrast = 0.431782125555502
    norm     = 0.5131289456743432.

The strongest adjacency contrast found was

    abs(contrast) = 0.4264948787147872

at

    t = 192.0

with

    abs(bias) = 0.3409009412302978
    norm      = 0.767395819945085.

Thus the sampled-pulse obstruction was not merely caused by failure to
sample a high-contrast constant pulse.

However, the high-contrast constant-pulse responses lose local
sharpness.

Using the Audit-028 state-independent Bell-capacity bound

    2|b_A b_B|
      + 2|b_A c_B|
      + 2|c_A b_B|
      + 2 sqrt(2)|c_A c_B|,

no scanned Alice/Bob constant-pulse pair exceeded 2.

The global maximum was the trivial zero-pulse detector:

    bias = 1
    contrast = 0
    Bell-capacity upper bound = 2.

Therefore the scanned constant-pulse family exhibits a
bias/contrast/sharpness obstruction.

This is a wide numerical frontier result, not a proof over all real
pulse durations.

For EPR assembly, no further constant-pulse optimization is planned
unless a later construction specifically requires it.

Next target:

    minimal structured local control that preserves local sharpness
    while transducing native analyzer orientation into the binary
    receipt.
