# Baseline correction

Subtracts baseline from epochs.

## Usage

``` r
baselineCorrect(
  x,
  baseline = c(-0.2, 0),
  method = c("mean", "median"),
  output_assay = "baseline_corrected"
)
```

## Arguments

- x:

  An epoched PhysioExperiment object.

- baseline:

  Numeric vector of length 2 (tmin, tmax) for baseline period.

- method:

  Correction method: "mean" or "median".

- output_assay:

  Name for the output assay.

## Value

Modified PhysioExperiment with baseline-corrected data.
