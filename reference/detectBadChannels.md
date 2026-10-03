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

## Examples

``` r
m <- matrix(rnorm(300 * 3), nrow = 300, ncol = 3)
m[, 2] <- m[, 2] * 10  # channel 2 is an amplitude outlier
pe <- PhysioExperiment(assays = list(raw = m), samplingRate = 100)
detectBadChannels(pe, method = "zscore")
#> [1] 2
```
