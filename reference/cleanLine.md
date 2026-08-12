# Remove line noise with CleanLine

Removes power-line noise (and optional harmonics) by sliding-window
sinusoidal regression with a Thomson F-test (the CleanLine / Chronux
approach; Mitra & Bokil, 2008). In each Hann-overlapped window each
channel is regressed onto the line sine/cosine pair; when the fit is
significant (`p < p_thresh`) the fitted sinusoid is subtracted. This
attenuates the line deeply without a wide spectral notch and works on
any number of channels.

## Usage

``` r
cleanLine(
  x,
  line_freq = 50,
  bandwidth = 2,
  harmonics = 1,
  p_thresh = 0.01,
  window_sec = 4,
  assay_name = NULL,
  output_assay = "cleanline"
)
```

## Arguments

- x:

  A PhysioExperiment object with a 2D (time x channels) assay.

- line_freq:

  Line frequency in Hz (default: 50).

- bandwidth:

  Reserved spectral bandwidth in Hz around the line (default: 2).

- harmonics:

  Number of harmonics to remove, including the fundamental (default: 1,
  i.e. `line_freq` only; 3 removes 50/100/150 Hz).

- p_thresh:

  F-test significance threshold for subtracting the line (default:
  0.01).

- window_sec:

  Sliding-window length in seconds (default: 4).

- assay_name:

  Input assay (default: `defaultAssay(x)`).

- output_assay:

  Name for the cleaned assay (default: "cleanline").

## Value

The PhysioExperiment with the cleaned signal in `output_assay`; the
removed line-power fraction is stored in `metadata(x)$line_noise` and
logged in provenance.

## References

Mitra, P., & Bokil, H. (2008). *Observed Brain Dynamics*. Oxford
University Press.

## See also

[`zapLine`](https://x-biosignal.github.io/PhysioPreprocess/reference/zapLine.md),
[`notchFilter`](https://x-biosignal.github.io/PhysioPreprocess/reference/notchFilter.md)
(the fast fixed notch)
