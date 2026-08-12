# Check if data is average referenced

Check if data is average referenced

## Usage

``` r
isAverageReferenced(x)
```

## Arguments

- x:

  A PhysioExperiment object.

## Value

A logical scalar: `TRUE` if the reference metadata is set to "average",
`FALSE` otherwise.

## References

Nunez, P.L. & Srinivasan, R. (2006). "Electric Fields of the Brain." 2nd
ed. Oxford University Press.

## See also

[`rereference()`](https://x-biosignal.github.io/PhysioPreprocess/reference/rereference.md)
for changing the reference,
[`getCurrentReference()`](https://x-biosignal.github.io/PhysioPreprocess/reference/getCurrentReference.md)
for querying the current reference.

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(100), nrow = 10, ncol = 10)),
  samplingRate = 100
)
pe <- rereference(pe, ref_type = "average")
isAverageReferenced(pe)  # TRUE
#> [1] TRUE
```
