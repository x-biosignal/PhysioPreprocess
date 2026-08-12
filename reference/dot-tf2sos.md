# Advanced signal filtering functions

This file provides advanced filtering operations including Butterworth,
FIR, and notch filters for physiological signal processing. Convert
transfer function to second-order sections

## Usage

``` r
.tf2sos(b, a)
```

## Arguments

- b:

  Numerator polynomial coefficients.

- a:

  Denominator polynomial coefficients.

## Value

A matrix with one row per section and 6 columns: b0, b1, b2, a0, a1, a2.

## Details

Converts filter coefficients from transfer function (b, a) form to
cascaded second-order sections (SOS) form. SOS form is numerically
stable for higher filter orders.
