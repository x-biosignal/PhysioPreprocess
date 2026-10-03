# Resampling operations for PhysioExperiment

Functions for resampling signal data to different sampling rates.
Resample signal data

## Usage

``` r
resample(
  x,
  target_rate,
  method = c("linear", "spline", "fft"),
  assay_name = NULL,
  output_assay = "resampled"
)
```

## Arguments

- x:

  A PhysioExperiment object.

- target_rate:

  Target sampling rate in Hz.

- method:

  Resampling method: "linear" (default), "spline", or "fft".

- assay_name:

  Optional assay name. If NULL, uses the default assay.

- output_assay:

  Name for the output assay. Default is "resampled".

## Value

A new `PhysioExperiment` object with sampling rate set to `target_rate`.
The time dimension is adjusted to match the new rate while preserving
the signal duration. Column data and metadata are carried over from the
input.

## Details

Resamples the signal data to a target sampling rate using interpolation.
Supports linear interpolation, spline interpolation, and FFT-based
resampling (zero-padding or truncation in the frequency domain).

## References

Crochiere, R.E. & Rabiner, L.R. (1983). "Multirate Digital Signal
Processing." Prentice Hall.

## See also

[`decimate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/decimate.md)
for integer-factor downsampling with anti-aliasing,
[`interpolate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/interpolate.md)
for integer-factor upsampling,
[`setAssaySamplingRate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/setAssaySamplingRate.md)
for per-assay rate tracking.

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(500 * 2), nrow = 500, ncol = 2)),
  samplingRate = 100
)
pe_50 <- resample(pe, target_rate = 50)
samplingRate(pe_50)
#> [1] 50
```
