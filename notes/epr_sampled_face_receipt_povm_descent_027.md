# Sampled face receipt POVM descent 027

## Result

Audit pass:

    True

Candidate dynamics:

    ['adjacency', 'laplacian_control']

Sampled pulses:

    [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 3.0, 4.0]

Settings per calibration:

    6

Total sampled calibration families:

    16

## Face-qubit receipt law

For every native analyzer setting i, let P_i be its rank-one face-qubit
projector.

Each sampled coherent response descends as

    E_(i,r)
      =
    alpha_r P_i
      +
    beta_r (I-P_i)

for four retained receipts

    clean pointer 0
    nonclean pointer 0
    clean pointer 1
    nonclean pointer 1.

All sampled coefficients are nonnegative.

The four receipt effects sum to identity.

## Numerical closure

Maximum setting covariance error:

    9.992007221626409e-15

Maximum binary completeness error:

    9.992007221626409e-15

Maximum clean/nonclean split error:

    1.0103029524088925e-14

Maximum four-receipt completeness error:

    9.992007221626409e-15

Scalar range:

    [0.0, 1.0]

## Binary pointer observable

After forgetting only the clean/nonclean subreceipt,

    A_i
      =
    E_(i,0) - E_(i,1)
      =
    bias I + contrast n_i.sigma.

Maximum sampled absolute contrast:

    0.07296147126218

Calibration attaining it:

    {'model': 'adjacency', 'dimensionless_pulse': 2.0, 'bias': -0.13365275831507994, 'contrast': 0.07296147126218}

This quantifies the inherited unsharpness before any Bell calculation.

## Boundary

This is a sampled mathematical POVM descent under the declared coherent
amplitude/readout model.

It does not insert a state-update rule.

It does not derive a physical frequency law.

It retains the nonclean receipts and the Audit019 zero preparation
branch.

The next gate is complete joint receipt normalization and no-signaling.
