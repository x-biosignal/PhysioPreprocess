# Remove ICA components

Removes specified ICA components and reconstructs the signal.

## Usage

``` r
icaRemove(x, components, output_assay = "ica_cleaned")
```

## Arguments

- x:

  A PhysioExperiment object with ICA decomposition.

- components:

  Integer vector of component indices to remove.

- output_assay:

  Name for the output assay.

## Value

Modified PhysioExperiment with cleaned signal.
