# TDOA-STFT-Spatial-Model-Mismatch
# Analysis of TDOA-Induced Spatial Model Mismatch in STFT-Domain Microphone Arrays

MATLAB code accompanying the manuscript:

"Analysis of TDOA-Induced Spatial Model Mismatch in STFT-Domain Microphone Arrays"

## Requirements

- MATLAB
- Signal Processing Toolbox

## Scripts

- `stft_parameter_sweep.m`
  Generates the controlled TDOA and STFT window-length experiments
  corresponding to Figs. 1 and 2.

- `tdoa_error_test.m`
  Evaluates robustness to TDOA estimation error and generates Fig. 3.

- `baseline_compare.m`
  Compares the four STFT spatial representations and generates Fig. 4.

- `real_waveform_validation.m`
  Performs the real-waveform validation reported in Table 2.

- `frac_delay.m`
  Implements fractional-sample delay filtering used in the real-waveform experiment.

## Data

The ultrasonic recordings used in the real-waveform validation were
provided as teaching materials for a laboratory exercise at the
University of New South Wales.

The recordings are not redistributed because public redistribution
rights have not been established.