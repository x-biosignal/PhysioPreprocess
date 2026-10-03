# Interpolate signal

Upsamples by an integer factor with interpolation. This is a convenience
wrapper around
[`resample()`](https://x-biosignal.github.io/PhysioPreprocess/reference/resample.md)
with `target_rate = sr * factor`.

## Usage

``` r
interpolate(
  x,
  factor,
  method = c("linear", "spline"),
  output_assay = "interpolated"
)
```

## Arguments

- x:

  A PhysioExperiment object.

- factor:

  Integer interpolation factor.

- method:

  Interpolation method: "linear" or "spline".

- output_assay:

  Name for the output assay.

## Value

A new `PhysioExperiment` object with sampling rate equal to the original
rate multiplied by `factor`. The time dimension is increased
accordingly.

## References

Crochiere, R.E. & Rabiner, L.R. (1983). "Multirate Digital Signal
Processing." Prentice Hall.

## See also

[`resample()`](https://x-biosignal.github.io/PhysioPreprocess/reference/resample.md)
for arbitrary-rate resampling,
[`decimate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/decimate.md)
for integer-factor downsampling.

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(200 * 2), nrow = 200, ncol = 2)),
  samplingRate = 100
)
pe_up <- interpolate(pe, factor = 2)
samplingRate(pe_up)
#> [1] 200
```
