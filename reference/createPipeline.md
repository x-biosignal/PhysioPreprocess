# Preprocessing Pipeline for PhysioExperiment

Functions for creating and applying preprocessing pipelines. Create
preprocessing pipeline

## Usage

``` r
createPipeline(...)
```

## Arguments

- ...:

  Preprocessing steps as named function calls.

## Value

A preprocessing pipeline object (list).

## Details

Creates a preprocessing pipeline specification that can be applied to
data.

## References

Oppenheim AV, Willsky AS, Nawab SH (1997). "Signals and Systems." 2nd
ed. Prentice Hall.

## See also

[`applyPipeline()`](https://x-biosignal.github.io/PhysioPreprocess/reference/applyPipeline.md)
for executing a pipeline on data,
[`filterSignals()`](https://x-biosignal.github.io/PhysioPreprocess/reference/filterSignals.md)
and
[`detrendSignals()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detrendSignals.md)
for individual preprocessing steps that can be included in a pipeline.

## Examples

``` r
pipeline <- createPipeline(
  filter = list(fn = "filterSignals", lowcut = 1, highcut = 40),
  detrend = list(fn = "detrendSignals", method = "linear")
)
```
