# Remove ICA components

Reconstructs signals after removing specified ICA components (e.g.,
artifacts).

## Usage

``` r
removeICAComponents(
  pe,
  ica_result,
  remove_components,
  assay_name = NULL,
  output_assay = "ica_cleaned"
)
```

## Arguments

- pe:

  A PhysioExperiment object.

- ica_result:

  Result from runICA().

- remove_components:

  Integer vector of component indices to remove.

- assay_name:

  Name of the assay to reconstruct.

- output_assay:

  Name for the cleaned output assay.

## Value

PhysioExperiment with cleaned data in output_assay.

## References

Hyvarinen A, Oja E (2000). "Independent component analysis: algorithms
and applications." *Neural Networks*, 13(4-5), 411-430.

## See also

[`runICA()`](https://x-biosignal.github.io/PhysioPreprocess/reference/runICA.md)
for performing the ICA decomposition,
[`icaRemove()`](https://x-biosignal.github.io/PhysioPreprocess/reference/icaRemove.md)
for the alternative component removal implementation,
[`detectBadChannels()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detectBadChannels.md)
for channel-level artifact detection.

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(1000), nrow = 100, ncol = 10)),
  colData = S4Vectors::DataFrame(label = paste0("Ch", 1:10)),
  samplingRate = 256
)
# \donttest{
if (requireNamespace("fastICA", quietly = TRUE)) {
  ica_result <- runICA(pe, n_components = 5)
  pe_clean <- removeICAComponents(pe, ica_result, remove_components = c(1, 3))
}
# }
```
