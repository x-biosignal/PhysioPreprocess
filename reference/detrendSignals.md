# Signal Detrending Functions for PhysioExperiment

Functions for removing trends from physiological signals. Detrend
signals

## Usage

``` r
detrendSignals(
  pe,
  method = c("linear", "mean", "polynomial"),
  order = 2,
  assay_name = NULL,
  output_assay = "detrended"
)
```

## Arguments

- pe:

  A PhysioExperiment object.

- method:

  Detrending method: "linear", "mean", or "polynomial".

- order:

  Polynomial order for method="polynomial".

- assay_name:

  Name of the assay to detrend.

- output_assay:

  Name for the detrended output assay.

## Value

PhysioExperiment with detrended data in output_assay.

## Details

Removes linear or polynomial trends from signals.

## References

Oppenheim AV, Willsky AS, Nawab SH (1997). "Signals and Systems." 2nd
ed. Prentice Hall.

## See also

[`removeBaseline()`](https://x-biosignal.github.io/PhysioPreprocess/reference/removeBaseline.md)
for baseline subtraction over a time window,
[`detrendSignal()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detrendSignal.md)
for the alternative linear/constant detrending implementation,
[`butterworthFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/butterworthFilter.md)
for highpass filtering as an alternative to detrending.

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(1000) + 1:100/10, nrow = 100, ncol = 10)),
  colData = S4Vectors::DataFrame(label = paste0("Ch", 1:10)),
  samplingRate = 256
)
pe_detrended <- detrendSignals(pe, method = "linear")
```
