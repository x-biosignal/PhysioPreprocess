# Apply Artifact Subspace Reconstruction (ASR)

Cleans a signal with a calibrated ASR model
([`asrCalibrate`](https://x-biosignal.github.io/PhysioPreprocess/reference/asrCalibrate.md)).
For each overlapping sliding window it eigendecomposes the window
covariance, flags principal directions whose RMS exceeds the calibrated
threshold, and reconstructs those directions from the retained ones via
the clean-covariance mixing matrix; windows are recombined by
Hann-weighted overlap-add.

## Usage

``` r
asrProcess(
  x,
  calibration,
  window_len = NULL,
  step = NULL,
  assay_name = NULL,
  output_assay = "asr"
)
```

## Arguments

- x:

  A PhysioExperiment object with a 2D (time x channels) assay.

- calibration:

  An `"asr_calibration"` from
  [`asrCalibrate`](https://x-biosignal.github.io/PhysioPreprocess/reference/asrCalibrate.md).

- window_len:

  Sliding-window length in seconds (default: the calibration value).

- step:

  Window step in seconds (default: half the window length).

- assay_name:

  Input assay (default: `defaultAssay(x)`).

- output_assay:

  Name for the cleaned assay (default: "asr").

## Value

The PhysioExperiment with the cleaned signal in `output_assay`; the
removed-variance fraction and rejected-component fraction are stored in
`metadata(x)$asr` and logged in provenance.

## See also

[`asrCalibrate`](https://x-biosignal.github.io/PhysioPreprocess/reference/asrCalibrate.md),
[`cleanRawdata`](https://x-biosignal.github.io/PhysioPreprocess/reference/cleanRawdata.md)
