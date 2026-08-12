# Design Butterworth filter directly in SOS form

Avoids polynomial root-finding by computing SOS sections directly from
the Butterworth analog prototype poles via bilinear transform. This is
numerically stable for all filter orders.

## Usage

``` r
.buttersos(n, W, type = c("low", "high", "pass", "stop"))
```

## Arguments

- n:

  Filter order.

- W:

  Normalized frequency (0-1, Nyquist = 1). Scalar for low/high, length-2
  vector for pass/stop.

- type:

  Filter type: "low", "high", "pass", or "stop".

## Value

SOS matrix (n_sections x 6).
