# Remove line noise with ZapLine

Removes power-line noise while preserving the broadband spectrum using
ZapLine (de Cheveigne, 2020): the line is first suppressed by a
comb-notch smoother, then a DSS/PCA of the notched-out residual isolates
the line-noise spatial subspace, and only its top `nremove` components
are removed from the data. This avoids the wide spectral notch of a
fixed filter and needs at least two channels (line noise is treated as a
spatial component).

## Usage

``` r
zapLine(
  x,
  line_freq = 50,
  nremove = 1,
  nfft = NULL,
  assay_name = NULL,
  output_assay = "zapline"
)
```

## Arguments

- x:

  A PhysioExperiment object with a 2D (time x channels) assay.

- line_freq:

  Line frequency in Hz (default: 50).

- nremove:

  Number of line-noise components to remove (default: 1).

- nfft:

  Optional analysis FFT length controlling the line-band width for the
  DSS bias (default: `NULL` = 1 Hz half-band).

- assay_name:

  Input assay (default: `defaultAssay(x)`).

- output_assay:

  Name for the cleaned assay (default: "zapline").

## Value

The PhysioExperiment with the cleaned signal in `output_assay`; the
removed line-power fraction is stored in `metadata(x)$line_noise` and
logged in provenance.

## References

de Cheveigne, A. (2020). "ZapLine: A simple and effective method to
remove power line artifacts." *NeuroImage*, 207, 116356.
[doi:10.1016/j.neuroimage.2019.116356](https://doi.org/10.1016/j.neuroimage.2019.116356)

## See also

[`cleanLine`](https://x-biosignal.github.io/PhysioPreprocess/reference/cleanLine.md),
[`notchFilter`](https://x-biosignal.github.io/PhysioPreprocess/reference/notchFilter.md)
(the fast fixed notch)

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(400 * 3), nrow = 400, ncol = 3)),
  samplingRate = 200
)
pe <- zapLine(pe, line_freq = 50)
```
