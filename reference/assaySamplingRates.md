# Per-assay sampling rates (moved to PhysioCore)

**Deprecated**: these accessors now live in PhysioExperiment. The
version here delegates to
[`PhysioExperiment::assaySamplingRates()`](https://x-biosignal.r-universe.dev/PhysioExperiment/reference/assaySamplingRates.html)
for backward compatibility; new code should call the PhysioExperiment
function directly, or use the canonical multi-rate container
[`PhysioExperiment::MultiPhysioExperiment()`](https://x-biosignal.r-universe.dev/PhysioExperiment/reference/MultiPhysioExperiment.html)
when streams differ in length.

## Usage

``` r
assaySamplingRates(x)
```

## Arguments

- x:

  A PhysioExperiment object.

## Value

A named numeric vector of per-assay sampling rates (Hz).

## See also

[`assaySamplingRates`](https://x-biosignal.r-universe.dev/PhysioExperiment/reference/assaySamplingRates.html),
[`MultiPhysioExperiment`](https://x-biosignal.r-universe.dev/PhysioExperiment/reference/MultiPhysioExperiment.html)

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(100 * 2), nrow = 100, ncol = 2)),
  samplingRate = 100
)
assaySamplingRates(pe)
#> raw 
#> 100 
```
