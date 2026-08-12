# Interpolate bad channels

Replaces bad channels with interpolated values from neighboring
channels.

## Usage

``` r
interpolateBadChannels(
  x,
  bad_channels,
  method = c("average", "spline"),
  output_assay = "interpolated"
)
```

## Arguments

- x:

  A PhysioExperiment object.

- bad_channels:

  Integer vector of channel indices to interpolate.

- method:

  Interpolation method: "average" or "spline".

- output_assay:

  Name for the output assay.

## Value

Modified PhysioExperiment with interpolated channels.
