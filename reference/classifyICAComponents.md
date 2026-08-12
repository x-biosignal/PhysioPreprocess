# Classify ICA components as brain or artifact

Automatically identifies artifact ICA components using statistical
criteria.

## Usage

``` r
classifyICAComponents(
  x,
  method = c("autocorrelation", "kurtosis", "frequency"),
  eog_channels = NULL,
  threshold = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with ICA decomposition (from
  [`icaDecompose`](https://x-biosignal.github.io/PhysioPreprocess/reference/icaDecompose.md)).

- method:

  Classification method: "autocorrelation", "kurtosis", or "frequency".

- eog_channels:

  Integer vector of EOG channel indices for correlation-based detection.
  If provided, also computes EOG correlation scores.

- threshold:

  Classification threshold. If NULL, uses method-specific defaults:
  autocorrelation=0.9, kurtosis=3.0, frequency=0.5.

## Value

A list with:

- `artifact_indices`: Integer vector of artifact component indices

- `scores`: Named numeric vector of scores per component

- `threshold`: The threshold used

- `method`: The method used
