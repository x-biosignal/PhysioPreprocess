# One-pass causal IIR filtering (Direct Form II Transposed)

Applies a rational transfer function defined by numerator `b` and
denominator `a` to `x`, forward in time, mirroring
`scipy.signal.lfilter`. Optionally carries the delay-line state so a
signal can be filtered in consecutive chunks with a result identical to
filtering it whole.

## Usage

``` r
lfilter(b, a, x, zi = NULL)
```

## Arguments

- b:

  Numerator (feed-forward) coefficients.

- a:

  Denominator (feedback) coefficients (`a[1]` normalised internally if
  not 1).

- x:

  Numeric input vector.

- zi:

  Optional initial delay-line state, length
  `max(length(a), length(b)) - 1`. When supplied, the return value is a
  list carrying the final state; when `NULL` (default) the state starts
  at zero and a plain filtered vector is returned. Use
  [`lfilterInit`](https://x-biosignal.github.io/PhysioPreprocess/reference/lfilterInit.md)
  for a transient-free warm start.

## Value

If `zi` is `NULL`, the filtered numeric vector. Otherwise a list with
`y` (filtered vector) and `zi` (final state).

## See also

[`lfilterInit`](https://x-biosignal.github.io/PhysioPreprocess/reference/lfilterInit.md),
[`sosfilt`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfilt.md)

## Examples

``` r
ba <- signal::butter(4, 0.2, "low")
# single call
y <- lfilter(ba$b, ba$a, rnorm(500))
# streamed in two chunks with carried state
x <- rnorm(500)
r1 <- lfilter(ba$b, ba$a, x[1:250], zi = numeric(4))
r2 <- lfilter(ba$b, ba$a, x[251:500], zi = r1$zi)
identical(length(c(r1$y, r2$y)), length(x))
#> [1] TRUE
```
