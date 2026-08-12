# Butterworth filter

Applies a Butterworth filter (lowpass, highpass, bandpass, or bandstop)
along the time axis of the specified assay.

## Usage

``` r
butterworthFilter(
  x,
  low = NULL,
  high = NULL,
  order = 4L,
  type = c("pass", "low", "high", "stop"),
  use_sos = TRUE,
  causal = FALSE,
  output_assay = "filtered"
)
```

## Arguments

- x:

  A `PhysioExperiment` object.

- low:

  Lower cutoff frequency in Hz. Required for highpass and bandpass.

- high:

  Upper cutoff frequency in Hz. Required for lowpass and bandpass.

- order:

  Filter order. Default is 4.

- type:

  Filter type: "low", "high", "pass" (bandpass), or "stop" (bandstop).

- use_sos:

  Logical. If TRUE (default), use second-order sections (SOS) form for
  improved numerical stability, especially for high filter orders. If
  FALSE, use the traditional transfer function (ba) form.

- causal:

  Logical. If FALSE (default), apply zero-phase forward-backward
  filtering
  ([`signal::filtfilt()`](https://rdrr.io/pkg/signal/man/filtfilt.html)
  / SOS `filtfilt`), which is non-causal and symmetrically smears sharp
  transitions in time. If TRUE, apply a single forward (causal) pass via
  [`sosfilt()`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfilt.md)
  with a steady-state warm start
  ([`sosfiltInit()`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfiltInit.md)
  scaled by the first sample), producing a real-time- equivalent result
  with a causal group delay and no acausal pre-ringing. Causal filtering
  always uses SOS form for numerical stability.

- output_assay:

  Name for the output assay. Default is "filtered".

## Value

The input object with a new assay containing filtered data.

## Details

By default, second-order sections (SOS) form is used for filtering,
which provides numerical stability for higher filter orders (\> 8). The
traditional transfer function (ba) form can exhibit coefficient
quantization errors at high orders, leading to unstable filters. Set
`use_sos = FALSE` to revert to the legacy ba-form behavior.

## See also

[`sosfilt()`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfilt.md),
[`StreamFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/StreamFilter.md)
for the underlying causal/stateful filtering primitives.

## Examples

``` r
# Create example EEG data
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(1000 * 4), nrow = 1000)),
  samplingRate = 250
)

# Bandpass filter (1-40 Hz) - common for EEG
pe <- butterworthFilter(pe, low = 1, high = 40, type = "pass")

# Lowpass filter (30 Hz)
pe <- butterworthFilter(pe, high = 30, type = "low",
                        output_assay = "lowpass")

# Highpass filter (0.5 Hz) to remove DC drift
pe <- butterworthFilter(pe, low = 0.5, type = "high",
                        output_assay = "highpass")

# High-order filter with SOS (numerically stable)
pe <- butterworthFilter(pe, low = 1, high = 40, type = "pass",
                        order = 10)

# Legacy ba-form filtering
pe <- butterworthFilter(pe, low = 1, high = 40, type = "pass",
                        use_sos = FALSE, output_assay = "filtered_ba")
```
