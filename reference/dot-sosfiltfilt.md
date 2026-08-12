# Zero-phase SOS filtering (forward-backward)

Applies the SOS filter forward, then reverses the result and applies the
filter again to achieve zero-phase distortion. This doubles the
effective filter order.

## Usage

``` r
.sosfiltfilt(sos, x)
```

## Arguments

- sos:

  SOS matrix (n_sections x 6).

- x:

  Numeric vector to filter.

## Value

Filtered numeric vector with zero phase distortion.
