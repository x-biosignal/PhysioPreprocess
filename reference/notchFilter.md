# Notch filter (power line noise removal)

Applies a notch filter to remove power line noise (50 Hz or 60 Hz) and
optionally its harmonics. Implemented as a Butterworth bandstop filter
with zero-phase filtering.

## Usage

``` r
notchFilter(
  x,
  freq = 50,
  bandwidth = 2,
  harmonics = 1L,
  output_assay = "filtered"
)
```

## Arguments

- x:

  A `PhysioExperiment` object.

- freq:

  Center frequency to remove in Hz. Default is 50 (European power line).

- bandwidth:

  Bandwidth of the notch in Hz. Default is 2.

- harmonics:

  Number of harmonics to remove. Default is 1 (only fundamental).

- output_assay:

  Name for the output assay. Default is "filtered".

## Value

A `PhysioExperiment` object with a new assay named `output_assay`
containing the notch-filtered data. Harmonics above the Nyquist
frequency are skipped with a warning.

## References

Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems." 2nd ed.
Prentice Hall.

## See also

[`butterworthFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/butterworthFilter.md)
for general Butterworth filtering,
[`firFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/firFilter.md)
for FIR filtering,
[`filterSignals()`](https://x-biosignal.github.io/PhysioPreprocess/reference/filterSignals.md)
for moving average filtering.

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(1000 * 4), nrow = 1000)),
  samplingRate = 250
)

# Remove 50 Hz power line noise (Europe/Asia)
pe <- notchFilter(pe, freq = 50)

# Remove 60 Hz and harmonics (Americas)
pe <- notchFilter(pe, freq = 60, harmonics = 2)
```
