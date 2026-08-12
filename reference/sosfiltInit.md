# Steady-state initial conditions for a cascaded SOS filter

Per-section analogue of
[`lfilterInit`](https://x-biosignal.github.io/PhysioPreprocess/reference/lfilterInit.md)
for second-order-sections form; equivalent to `scipy.signal.sosfilt_zi`.
Each section's steady state is scaled by the DC gain of all preceding
sections so that the whole cascade starts in steady state when the
returned matrix (scaled by the first sample) is handed to
[`sosfilt`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfilt.md).

## Usage

``` r
sosfiltInit(sos)
```

## Arguments

- sos:

  SOS matrix (`n_sections` x 6).

## Value

A numeric matrix (`n_sections` x 2) of per-section initial states
(unscaled).

## See also

[`sosfilt`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfilt.md),
[`lfilterInit`](https://x-biosignal.github.io/PhysioPreprocess/reference/lfilterInit.md)

## Examples

``` r
sos <- sosDesign(high = 40, type = "low", sr = 250)
zi <- sosfiltInit(sos)
```
