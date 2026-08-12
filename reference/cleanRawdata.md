# Automated raw-data cleaning (EEGLAB clean_rawdata-style)

Orchestrates a standard cleaning pipeline mirroring EEGLAB's
`clean_rawdata`: flat-channel removal, high-pass drift removal,
correlation-based bad-channel interpolation, and Artifact Subspace
Reconstruction. The removed-variance fraction and any repaired channels
are recorded in provenance.

## Usage

``` r
cleanRawdata(
  x,
  flatline = TRUE,
  highpass = c(0.25, 0.75),
  channel_crit = 0.8,
  asr_cutoff = 20,
  window_len = 0.5,
  calib_window = NULL,
  assay_name = NULL,
  output_assay = "clean"
)
```

## Arguments

- x:

  A PhysioExperiment object with a 2D (time x channels) assay.

- flatline:

  Logical; detect and interpolate (near-)flat channels (default:
  `TRUE`).

- highpass:

  Numeric `c(low, high)` transition band in Hz for the drift-removal
  high-pass; the cutoff is the band midpoint (default `c(0.25, 0.75)`).
  `NULL` skips high-pass filtering.

- channel_crit:

  Minimum mean absolute correlation with the other channels for a
  channel to be kept; channels below this are interpolated (default:
  0.8). `NULL` skips the correlation criterion.

- asr_cutoff:

  ASR rejection cutoff in robust SDs (default: 20). `Inf` or `NULL`
  skips ASR.

- window_len:

  ASR sliding-window length in seconds (default: 0.5).

- calib_window:

  Optional `c(start_sec, end_sec)` calibration segment for ASR.

- assay_name:

  Input assay (default: `defaultAssay(x)`).

- output_assay:

  Name for the cleaned assay (default: "clean").

## Value

The PhysioExperiment with the cleaned signal in `output_assay`;
`metadata(x)$clean_rawdata` records the interpolated channels and the
removed-variance fraction, also logged in provenance.

## References

Mullen, T. R., et al. (2015). "Real-time neuroimaging and cognitive
monitoring using wearable dry EEG." *IEEE Transactions on Biomedical
Engineering*, 62(11), 2553-2567.
[doi:10.1109/TBME.2015.2481482](https://doi.org/10.1109/TBME.2015.2481482)

## See also

[`asrCalibrate`](https://x-biosignal.github.io/PhysioPreprocess/reference/asrCalibrate.md),
[`asrProcess`](https://x-biosignal.github.io/PhysioPreprocess/reference/asrProcess.md)
