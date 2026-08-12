# Zero-phase SOS filtering with state-initialised edge handling

Forward-backward (zero-phase) SOS filtering that mirrors
`scipy.signal.sosfiltfilt`: it odd-reflects the signal at both edges,
initialises the delay-line to steady state (scaled by the boundary
sample) to suppress start-up transients, filters forward and then
backward. The result has zero phase distortion but is non-causal (both
directions are used). For strictly causal / real-time-equivalent
filtering, use
[`sosfilt`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfilt.md)
or
[`butterworthFilter`](https://x-biosignal.github.io/PhysioPreprocess/reference/butterworthFilter.md)
with `causal = TRUE`.

## Usage

``` r
sosfiltfilt(sos, x, padlen = NULL)
```

## Arguments

- sos:

  SOS matrix (`n_sections` x 6).

- x:

  Numeric input vector.

- padlen:

  Number of samples of odd reflection padding at each edge. Defaults to
  `min(3 * (2 * n_sections + 1), length(x) - 1)`.

## Value

The zero-phase filtered numeric vector, same length as `x`.

## See also

[`sosfilt`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfilt.md),
[`sosDesign`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosDesign.md)

## Examples

``` r
sos <- sosDesign(high = 40, type = "low", sr = 250, order = 4)
y <- sosfiltfilt(sos, rnorm(1000))
```
