# Get current reference

Returns the current reference electrode information.

## Usage

``` r
getCurrentReference(x)
```

## Arguments

- x:

  A PhysioExperiment object.

## Value

A character string describing the current reference (e.g., "average",
"Cz", "M1+M2"), or `NULL` if no reference has been set.

## References

Nunez, P.L. & Srinivasan, R. (2006). "Electric Fields of the Brain." 2nd
ed. Oxford University Press.

## See also

[`rereference()`](https://x-biosignal.github.io/PhysioPreprocess/reference/rereference.md)
for changing the reference,
[`isAverageReferenced()`](https://x-biosignal.github.io/PhysioPreprocess/reference/isAverageReferenced.md)
for checking average reference status.

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(100), nrow = 10, ncol = 10)),
  samplingRate = 100
)
pe <- setReference(pe, "Cz")
getCurrentReference(pe)  # "Cz"
#> [1] "Cz"
```
