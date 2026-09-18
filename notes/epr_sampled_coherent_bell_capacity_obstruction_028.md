# Sampled coherent Bell-capacity obstruction 028

## Result

Audit pass:

    True

The local binary pointer observable is

    A_i = bias I + contrast N_i

with

    N_i^2 = I.

For arbitrary local outcome relabelings and arbitrary analyzer
directions, the CHSH operator obeys

    ||CHSH||
      <=
    2 |b_A b_B|
      +
    2 |b_A c_B|
      +
    2 |c_A b_B|
      +
    2 sqrt(2) |c_A c_B|.

This is state-independent.

It does not use a probability law.

## Sampled census

Local calibration count:

    16

Alice/Bob calibration pairs:

    256

Upper bounds exceeding 2:

    0

Maximum upper bound:

    2.0

Pair attaining it:

    {'Alice': 'adjacency@0.0', 'Bob': 'adjacency@0.0', 'Alice_bias': 1.0, 'Alice_contrast': 0.0, 'Bob_bias': 1.0, 'Bob_contrast': 0.0, 'state_independent_CHSH_upper_bound': 2.0, 'exceeds_2': False}

Maximum upper bound with nonzero pulse on both wings:

    1.9227855165678829

Pair attaining it:

    {'Alice': 'laplacian_control@0.1', 'Bob': 'laplacian_control@0.1', 'Alice_bias': 0.98049150206636, 'Alice_contrast': 1.4876389319962158e-05, 'Bob_bias': 0.98049150206636, 'Bob_contrast': 1.4876389319962158e-05, 'state_independent_CHSH_upper_bound': 1.9227855165678829, 'exceeds_2': False}

## Interpretation

Audit 022 remains intact:

    the native six-axis geometry is Bell-capable.

Audit 027 also remains intact:

    the sampled coherent local receipt POVMs are positive and normalized.

The obstruction lies in the sampled response strength.

Maximum sampled analyzer contrast:

    0.07296147126218

The inherited sampled coherent response therefore cannot transmit enough
setting dependence to support CHSH violation, regardless of joint state.

No Born-frequency law or state-update instrument is used in this
obstruction.

The next target is a stronger native local transduction law, not joint
probability bookkeeping for a response already below Bell capacity.
