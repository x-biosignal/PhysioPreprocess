# Stateful streaming filter object

Constructs a mutable filter object that carries its delay-line state
across successive `$apply()` calls, so a signal fed in arbitrary chunks
produces exactly the same output as filtering it whole. The identical
math is used offline and online, which lets a real-time streaming
backend reuse this object without re-deriving coefficients.

## Usage

``` r
StreamFilter(sos = NULL, b = NULL, a = NULL, warmup = TRUE)
```

## Arguments

- sos:

  SOS matrix (`n_sections` x 6), e.g. from
  [`sosDesign`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosDesign.md).
  Provide either `sos` or both `b` and `a`.

- b, a:

  Transfer-function coefficients (converted to a single SOS section only
  when they describe a biquad; otherwise supply `sos`).

- warmup:

  Logical; if `TRUE` (default) the first `$apply()` after construction
  or `$reset()` initialises the state to steady state scaled by the
  first sample, suppressing the start-up transient. If `FALSE`, the
  state starts at zero.

## Value

An object of class `"StreamFilter"` (an environment).

## Details

The returned object is an environment (reference semantics) with
methods:

- `$apply(x)`:

  Filter numeric vector `x`, advancing and storing the internal state;
  returns the filtered vector.

- `$reset()`:

  Clear the carried state so the next `$apply()` starts fresh (cold or
  warm start per `warmup`).

- `$state()`:

  Return the current internal state matrix.

## See also

[`sosfilt`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfilt.md),
[`sosDesign`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosDesign.md)

## Examples

``` r
sos <- sosDesign(high = 30, type = "low", sr = 250, order = 4)
sf <- StreamFilter(sos)
x <- rnorm(1000)
y_chunked <- c(sf$apply(x[1:400]), sf$apply(x[401:1000]))
sf$reset()
y_whole <- sf$apply(x)
max(abs(y_chunked - y_whole))
#> [1] 0
```
