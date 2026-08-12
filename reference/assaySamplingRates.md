# Per-assay sampling rates (moved to PhysioCore)

**Deprecated**: these accessors now live in PhysioCore. The version here
delegates to
[`PhysioCore::assaySamplingRates()`](https://x-biosignal.r-universe.dev/PhysioCore/reference/assaySamplingRates.html)
for backward compatibility; new code should call the PhysioCore function
directly, or use the canonical multi-rate container
[`PhysioCore::MultiRatePhysioExperiment()`](https://x-biosignal.r-universe.dev/PhysioCore/reference/MultiRatePhysioExperiment.html)
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

[`assaySamplingRates`](https://x-biosignal.r-universe.dev/PhysioCore/reference/assaySamplingRates.html),
[`MultiRatePhysioExperiment`](https://x-biosignal.r-universe.dev/PhysioCore/reference/MultiRatePhysioExperiment.html)
