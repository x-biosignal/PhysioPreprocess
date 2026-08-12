# One-pass causal SOS filtering with optional carried state

Filters `x` through a cascade of second-order sections, forward in time,
mirroring `scipy.signal.sosfilt`. Each section uses the DF2T biquad
recursion. Optionally carries the two-sample-per-section state so a
signal can be processed in consecutive chunks identically to processing
it whole.

## Usage

``` r
sosfilt(sos, x, zi = NULL)
```

## Arguments

- sos:

  SOS matrix (`n_sections` x 6), e.g. from
  [`sosDesign`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosDesign.md).

- x:

  Numeric input vector.

- zi:

  Optional initial state, a matrix `n_sections` x 2. When supplied the
  return value is a list carrying the final state; when `NULL` (default)
  state starts at zero and a plain vector is returned. Use
  [`sosfiltInit`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfiltInit.md)
  for a transient-free warm start.

## Value

If `zi` is `NULL`, the filtered numeric vector. Otherwise a list with
`y` (filtered vector) and `zi` (final `n_sections` x 2 state matrix).

## See also

[`sosDesign`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosDesign.md),
[`sosfiltInit`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfiltInit.md),
[`sosfiltfilt`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfiltfilt.md)

## Examples

``` r
sos <- sosDesign(high = 40, type = "low", sr = 250, order = 4)
y <- sosfilt(sos, rnorm(1000))
# chunked with carried state == whole-signal filtering
x <- rnorm(1000)
zi0 <- matrix(0, nrow(sos), 2)
r1 <- sosfilt(sos, x[1:500], zi = zi0)
r2 <- sosfilt(sos, x[501:1000], zi = r1$zi)
max(abs(c(r1$y, r2$y) - sosfilt(sos, x)))
#> [1] 0
```
