# Detect bad channels

Identifies channels with abnormal characteristics.

## Usage

``` r
detectBadChannels(
  x,
  method = c("zscore", "correlation", "flatline"),
  threshold = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object.

- method:

  Detection method: "zscore", "correlation", or "flatline".

- threshold:

  Threshold for detection (depends on method).

## Value

Integer vector of bad channel indices.
