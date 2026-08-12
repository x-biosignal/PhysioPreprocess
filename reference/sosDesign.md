# Design a Butterworth filter in second-order-sections (SOS) form

Builds the cascaded second-order-sections representation of a
Butterworth filter directly from the analog prototype poles (bilinear
transform), which is numerically stable for all orders. This is the
design entry point for the causal / stateful filtering functions
([`sosfilt`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfilt.md),
[`sosfiltfilt`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfiltfilt.md),
[`StreamFilter`](https://x-biosignal.github.io/PhysioPreprocess/reference/StreamFilter.md)).

## Usage

``` r
sosDesign(
  low = NULL,
  high = NULL,
  order = 4L,
  type = c("pass", "low", "high", "stop"),
  sr
)
```

## Arguments

- low:

  Lower cutoff frequency in Hz. Required for `"high"` and
  `"pass"`/`"stop"`.

- high:

  Upper cutoff frequency in Hz. Required for `"low"` and
  `"pass"`/`"stop"`.

- order:

  Filter order (per band). Default `4`.

- type:

  Filter type: `"low"`, `"high"`, `"pass"` (bandpass), or `"stop"`
  (bandstop).

- sr:

  Sampling rate in Hz.

## Value

A numeric matrix with one second-order section per row and six columns
`c(b0, b1, b2, a0, a1, a2)`.

## See also

[`sosfilt`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfilt.md),
[`sosfiltfilt`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfiltfilt.md),
[`StreamFilter`](https://x-biosignal.github.io/PhysioPreprocess/reference/StreamFilter.md)

## Examples

``` r
sos <- sosDesign(high = 40, type = "low", sr = 250, order = 4)
y <- sosfilt(sos, rnorm(1000))
```
