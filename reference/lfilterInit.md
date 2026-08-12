# Steady-state initial conditions for a one-pass IIR filter

Computes the filter delay-line state that, when scaled by the first
input sample and passed to
[`lfilter`](https://x-biosignal.github.io/PhysioPreprocess/reference/lfilter.md),
makes the filter's step response start already in steady state (no
start-up transient). This is the direct analogue of
`scipy.signal.lfilter_zi`.

## Usage

``` r
lfilterInit(b, a)
```

## Arguments

- b:

  Numerator (feed-forward) coefficients.

- a:

  Denominator (feedback) coefficients; `a[1]` need not be 1
  (coefficients are normalised internally).

## Value

A numeric vector of length `max(length(a), length(b)) - 1`, the initial
delay-line state (unscaled; multiply by the first sample level).

## See also

[`lfilter`](https://x-biosignal.github.io/PhysioPreprocess/reference/lfilter.md),
[`sosfiltInit`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfiltInit.md)

## Examples

``` r
ba <- signal::butter(2, 0.2, "low")
zi <- lfilterInit(ba$b, ba$a)
```
