# Clean data using an artifact removal pipeline

Runs a configurable sequence of artifact removal steps.

## Usage

``` r
cleanData(
  x,
  steps = c("bad_channels", "ica"),
  bad_channel_method = "zscore",
  bad_channel_threshold = NULL,
  interpolation_method = "average",
  ica_method = "fastica",
  ica_classify_method = "kurtosis",
  output_assay = "cleaned"
)
```

## Arguments

- x:

  A PhysioExperiment object.

- steps:

  Character vector of steps to run. Options: "bad_channels", "ica".

- bad_channel_method:

  Method for
  [`detectBadChannels`](https://x-biosignal.github.io/PhysioPreprocess/reference/detectBadChannels.md).

- bad_channel_threshold:

  Threshold for bad channel detection.

- interpolation_method:

  Method for
  [`interpolateBadChannels`](https://x-biosignal.github.io/PhysioPreprocess/reference/interpolateBadChannels.md).

- ica_method:

  Method for
  [`icaDecompose`](https://x-biosignal.github.io/PhysioPreprocess/reference/icaDecompose.md).

- ica_classify_method:

  Method for
  [`classifyICAComponents`](https://x-biosignal.github.io/PhysioPreprocess/reference/classifyICAComponents.md).

- output_assay:

  Name for the output assay.

## Value

Modified PhysioExperiment with cleaned data.
