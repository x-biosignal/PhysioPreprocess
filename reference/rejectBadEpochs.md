# Reject bad epochs

Identifies and optionally removes epochs with artifacts.

## Usage

``` r
rejectBadEpochs(
  x,
  threshold = 100,
  method = c("amplitude", "gradient", "variance"),
  remove = TRUE
)
```

## Arguments

- x:

  An epoched PhysioExperiment object.

- threshold:

  Amplitude threshold for rejection.

- method:

  Detection method: "amplitude", "gradient", or "variance".

- remove:

  If TRUE, removes bad epochs. If FALSE, returns indices only.

## Value

If remove=TRUE, modified object. If remove=FALSE, indices of bad epochs.
