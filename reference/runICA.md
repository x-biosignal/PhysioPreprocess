# ICA-Based Artifact Removal for PhysioExperiment

Functions for Independent Component Analysis (ICA) decomposition and
artifact removal using the fastICA algorithm. Perform ICA decomposition

## Usage

``` r
runICA(
  pe,
  n_components = NULL,
  assay_name = NULL,
  method = "fastica",
  max_iter = 200,
  ...
)
```

## Arguments

- pe:

  A PhysioExperiment object.

- n_components:

  Number of components to extract. If NULL, uses all channels.

- assay_name:

  Name of the assay to decompose.

- method:

  ICA algorithm: "fastica" (default).

- max_iter:

  Maximum number of iterations.

- ...:

  Additional arguments passed to fastICA::fastICA.

## Value

A list with components:

- S:

  Source matrix (time x components)

- A:

  Mixing matrix (channels x components)

- W:

  Unmixing matrix (components x channels)

## Details

Decomposes signals into independent components using the fastICA
algorithm.

## References

Hyvarinen A, Oja E (2000). "Independent component analysis: algorithms
and applications." *Neural Networks*, 13(4-5), 411-430.

## See also

[`removeICAComponents()`](https://x-biosignal.github.io/PhysioPreprocess/reference/removeICAComponents.md)
for reconstructing signals after removing artifact components,
[`icaDecompose()`](https://x-biosignal.github.io/PhysioPreprocess/reference/icaDecompose.md)
for the built-in ICA implementation,
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
}
# }
```
