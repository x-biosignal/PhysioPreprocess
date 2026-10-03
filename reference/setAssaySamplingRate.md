# Set a per-assay sampling rate (moved to PhysioCore)

**Deprecated**: delegates to
[`PhysioExperiment::setAssaySamplingRate()`](https://x-biosignal.r-universe.dev/PhysioExperiment/reference/setAssaySamplingRate.html).
New code should use PhysioExperiment directly.

## Usage

``` r
setAssaySamplingRate(x, assay_name, rate)
```

## Arguments

- x:

  A PhysioExperiment object.

- assay_name:

  Name of the assay.

- rate:

  Sampling rate for the assay in Hz.

## Value

The updated `PhysioExperiment`.

## See also

[`setAssaySamplingRate`](https://x-biosignal.r-universe.dev/PhysioExperiment/reference/setAssaySamplingRate.html)

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(100 * 2), nrow = 100, ncol = 2)),
  samplingRate = 100
)
pe <- setAssaySamplingRate(pe, "raw", 100)
assaySamplingRates(pe)
#> raw 
#> 100 
```
