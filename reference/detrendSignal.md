# Detrend signal

Removes linear or polynomial trends from the signal.

## Usage

``` r
detrendSignal(x, type = c("linear", "constant"), output_assay = "detrended")
```

## Arguments

- x:

  A `PhysioExperiment` object.

- type:

  Type of detrending: "linear" or "constant" (mean removal).

- output_assay:

  Name for the output assay. Default is "detrended".

## Value

A `PhysioExperiment` object with a new assay named `output_assay`
containing detrended data. For "constant" type, the channel mean is
subtracted. For "linear" type, a least-squares linear fit is removed.

## References

Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems." 2nd ed.
Prentice Hall.

## See also

[`detrendSignals()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detrendSignals.md)
for the alternative detrending implementation with polynomial support,
[`butterworthFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/butterworthFilter.md)
for highpass filtering as an alternative to detrending.
