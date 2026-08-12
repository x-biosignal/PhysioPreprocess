# Moving average filter

Applies a moving average filter along the first dimension (time axis) of
the default assay.

## Usage

``` r
filterSignals(x, window = 5L, na.rm = FALSE, output_assay = "filtered")
```

## Arguments

- x:

  A `PhysioExperiment` object.

- window:

  Integer window length for the moving average.

- na.rm:

  Logical. If TRUE, NA values are ignored in the filter computation.

- output_assay:

  Name for the output assay. Default is "filtered".

## Value

A `PhysioExperiment` object with a new assay named `output_assay`
containing the moving-average-filtered data. Dimensions match the input
assay. Edge values where the full window cannot be applied are set to
`NA`.

## References

Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems." 2nd ed.
Prentice Hall.

## See also

[`butterworthFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/butterworthFilter.md)
for IIR filtering,
[`firFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/firFilter.md)
for FIR filtering,
[`notchFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/notchFilter.md)
for power line noise removal,
[`detrendSignal()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detrendSignal.md)
for trend removal.
