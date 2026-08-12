# Apply preprocessing pipeline

Applies a preprocessing pipeline to a PhysioExperiment object.

## Usage

``` r
applyPipeline(pe, pipeline, verbose = FALSE)
```

## Arguments

- pe:

  A PhysioExperiment object.

- pipeline:

  A pipeline object created with createPipeline().

- verbose:

  Logical. If TRUE, prints progress messages.

## Value

PhysioExperiment with all pipeline steps applied.

## References

Oppenheim AV, Willsky AS, Nawab SH (1997). "Signals and Systems." 2nd
ed. Prentice Hall.

## See also

[`createPipeline()`](https://x-biosignal.github.io/PhysioPreprocess/reference/createPipeline.md)
for defining the pipeline steps,
[`butterworthFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/butterworthFilter.md)
and
[`resample()`](https://x-biosignal.github.io/PhysioPreprocess/reference/resample.md)
for common preprocessing operations used within pipelines.

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(1000), nrow = 100, ncol = 10)),
  colData = S4Vectors::DataFrame(label = paste0("Ch", 1:10)),
  samplingRate = 256
)
pipeline <- createPipeline(
  detrend = list(fn = "detrendSignals", method = "linear")
)
pe_processed <- applyPipeline(pe, pipeline)
```
