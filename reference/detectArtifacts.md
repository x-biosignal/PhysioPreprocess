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
