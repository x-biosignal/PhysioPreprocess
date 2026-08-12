# Apply SOS filter to a signal (single pass)

Cascades second-order sections to filter a signal. Each section is
applied sequentially using
[`signal::filter`](https://rdrr.io/pkg/signal/man/filter.html).

## Usage

``` r
.sosfilt(sos, x)
```

## Arguments

- sos:

  SOS matrix (n_sections x 6).

- x:

  Numeric vector to filter.

## Value

Filtered numeric vector.
