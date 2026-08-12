# Remove baseline

Subtracts baseline from a specified time window.

## Usage

``` r
removeBaseline(
  pe,
  baseline_start,
  baseline_end,
  assay_name = NULL,
  output_assay = "baseline_corrected"
)
```

## Arguments

- pe:

  A PhysioExperiment object.

- baseline_start:

  Start time of baseline window in seconds.

- baseline_end:

  End time of baseline window in seconds.

- assay_name:

  Name of the assay to baseline correct.

- output_assay:

  Name for the baseline-corrected output assay.

## Value

PhysioExperiment with baseline-corrected data.

## References

Oppenheim AV, Willsky AS, Nawab SH (1997). "Signals and Systems." 2nd
ed. Prentice Hall.

## See also

[`detrendSignals()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detrendSignals.md)
for trend removal,
[`baselineCorrect()`](https://x-biosignal.github.io/PhysioPreprocess/reference/baselineCorrect.md)
for epoch-based baseline correction,
[`filterSignals()`](https://x-biosignal.github.io/PhysioPreprocess/reference/filterSignals.md)
for moving average smoothing.

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(1000) + 5, nrow = 100, ncol = 10)),
  colData = S4Vectors::DataFrame(label = paste0("Ch", 1:10)),
  samplingRate = 100
)
pe_corrected <- removeBaseline(pe, baseline_start = 0, baseline_end = 0.2)
```
