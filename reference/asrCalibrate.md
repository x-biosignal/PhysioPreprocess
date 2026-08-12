# Calibrate Artifact Subspace Reconstruction (ASR)

Learns the clean-baseline channel covariance and the per-component RMS
thresholds that
[`asrProcess`](https://x-biosignal.github.io/PhysioPreprocess/reference/asrProcess.md)
uses to detect and reconstruct artifact subspaces (Mullen et al., 2015).
Calibration data should be a relatively clean segment; supply
`calib_window` to restrict it.

## Usage

``` r
asrCalibrate(
  x,
  cutoff = 20,
  calib_window = NULL,
  window_len = 0.5,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with a 2D (time x channels) assay.

- cutoff:

  Rejection cutoff in robust standard deviations; a principal direction
  whose window RMS exceeds `median + cutoff * MAD` of the calibration
  distribution is treated as artifact (default: 20). `Inf` makes
  [`asrProcess`](https://x-biosignal.github.io/PhysioPreprocess/reference/asrProcess.md)
  an identity.

- calib_window:

  Optional numeric `c(start_sec, end_sec)` selecting the calibration
  segment. If `NULL` (default) the whole signal is used.

- window_len:

  Sliding-window length in seconds (default: 0.5).

- assay_name:

  Input assay (default: `defaultAssay(x)`).

## Value

An object of class `"asr_calibration"`: a list with the clean covariance
`M`, its square-root mixing matrix, the threshold matrix, per-component
thresholds, and settings.

## References

Mullen, T. R., et al. (2015). "Real-time neuroimaging and cognitive
monitoring using wearable dry EEG." *IEEE Transactions on Biomedical
Engineering*, 62(11), 2553-2567.
[doi:10.1109/TBME.2015.2481482](https://doi.org/10.1109/TBME.2015.2481482)

## See also

[`asrProcess`](https://x-biosignal.github.io/PhysioPreprocess/reference/asrProcess.md),
[`cleanRawdata`](https://x-biosignal.github.io/PhysioPreprocess/reference/cleanRawdata.md)
