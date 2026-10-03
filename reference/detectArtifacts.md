# Detect artifacts in continuous data

Scans continuous data for artifacts using amplitude or gradient
criteria.

## Usage

``` r
detectArtifacts(
  x,
  method = c("amplitude", "gradient", "joint"),
  threshold = NULL,
  window_sec = 1
)
```

## Arguments

- x:

  A PhysioExperiment object.

- method:

  Detection method: "amplitude", "gradient", or "joint".

- threshold:

  Detection threshold. If NULL, computed as median + 5 \* MAD.

- window_sec:

  Detection window in seconds (default: 1.0).

## Value

A data.frame with columns: onset (seconds), offset (seconds), channel
(index), type (detection method).

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(500 * 2), nrow = 500, ncol = 2)),
  samplingRate = 100
)
arts <- detectArtifacts(pe, method = "amplitude")
head(arts)
#>   onset offset channel      type
#> 1  3.36   3.37       1 amplitude
```
