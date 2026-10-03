# Decimate signal

Downsamples by an integer factor with anti-aliasing filter. An
80%-Nyquist lowpass Butterworth filter is applied before decimation to
prevent aliasing.

## Usage

``` r
decimate(x, factor, filter_order = 8L, output_assay = "decimated")
```

## Arguments

- x:

  A PhysioExperiment object.

- factor:

  Integer decimation factor.

- filter_order:

  Order of the anti-aliasing lowpass filter.

- output_assay:

  Name for the output assay.

## Value

A new `PhysioExperiment` object with sampling rate equal to the original
rate divided by `factor`. The time dimension is reduced accordingly.
Column data and metadata are carried over.

## References

Crochiere, R.E. & Rabiner, L.R. (1983). "Multirate Digital Signal
Processing." Prentice Hall.

## See also

[`resample()`](https://x-biosignal.github.io/PhysioPreprocess/reference/resample.md)
for arbitrary-rate resampling,
[`interpolate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/interpolate.md)
for integer-factor upsampling,
[`butterworthFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/butterworthFilter.md)
for the anti-aliasing filter used internally.

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(400 * 2), nrow = 400, ncol = 2)),
  samplingRate = 100
)
pe_dec <- decimate(pe, factor = 2)
samplingRate(pe_dec)
#> [1] 50
```
