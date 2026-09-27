# Set a per-assay sampling rate (moved to PhysioCore)

**Deprecated**: delegates to
[`PhysioCore::setAssaySamplingRate()`](https://x-biosignal.r-universe.dev/PhysioExperiment/reference/setAssaySamplingRate.html).
New code should use PhysioCore directly.

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
