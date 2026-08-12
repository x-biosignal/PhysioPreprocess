# ICA and Artifact Removal for PhysioExperiment

Functions for Independent Component Analysis (ICA) and artifact removal
from physiological signals.

## Usage

``` r
icaDecompose(
  x,
  n_components = NULL,
  method = c("fastica", "jade"),
  max_iter = 200L,
  tol = 1e-04
)
```

## Arguments

- x:

  A PhysioExperiment object.

- n_components:

  Number of components to extract. If NULL, uses number of channels.

- method:

  ICA method: "fastica" (default) or "jade".

- max_iter:

  Maximum iterations for convergence.

- tol:

  Tolerance for convergence.

## Value

A list with four elements:

- components:

  The independent components as a matrix or 3D array (time by component
  by samples for 3D data).

- mixing:

  The mixing matrix (channels x components).

- unmixing:

  The unmixing matrix (components x channels).

- object:

  The input `PhysioExperiment` with the ICA decomposition (components,
  mixing/unmixing matrices, column means) stored in `metadata(x)$ica`.
  The components are held in metadata rather than as an assay because
  their column count is `n_components`, which need not equal the
  object's channel count.

## References

Hyvarinen, A. & Oja, E. (2000). "Independent component analysis:
algorithms and applications." Neural Networks, 13(4-5), 411-430.
doi:10.1016/S0893-6080(00)00026-5 Perform ICA decomposition

Decomposes the signal into independent components using FastICA
algorithm.

Hyvarinen, A. & Oja, E. (2000). "Independent component analysis:
algorithms and applications." Neural Networks, 13(4-5), 411-430.
doi:10.1016/S0893-6080(00)00026-5

## See also

[`icaRemove()`](https://x-biosignal.github.io/PhysioPreprocess/reference/icaRemove.md)
for removing specific components after decomposition,
[`runICA()`](https://x-biosignal.github.io/PhysioPreprocess/reference/runICA.md)
for an alternative ICA implementation using the fastICA package,
[`detectBadChannels()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detectBadChannels.md)
for channel-level artifact detection.
