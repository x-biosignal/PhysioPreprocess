# FIR filter

Applies a Finite Impulse Response (FIR) filter along the time axis. Uses
zero-phase forward-backward filtering via
[`signal::filtfilt()`](https://rdrr.io/pkg/signal/man/filtfilt.html) to
avoid phase distortion.

## Usage

``` r
firFilter(
  x,
  low = NULL,
  high = NULL,
  order = 100L,
  type = c("pass", "low", "high", "stop"),
  window = "hamming",
  causal = FALSE,
  output_assay = "filtered"
)
```

## Arguments

- x:

  A `PhysioExperiment` object.

- low:

  Lower cutoff frequency in Hz.

- high:

  Upper cutoff frequency in Hz.

- order:

  Filter order (number of taps - 1). Default is 100.

- type:

  Filter type: "low", "high", "pass" (bandpass), or "stop" (bandstop).

- window:

  Window function for FIR design. Default is "hamming".

- causal:

  Logical. If FALSE (default), apply zero-phase forward-backward
  filtering
  ([`signal::filtfilt()`](https://rdrr.io/pkg/signal/man/filtfilt.html)),
  which is non-causal. If TRUE, apply a single forward (causal) pass via
  [`lfilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/lfilter.md)
  with a steady-state warm start, producing a real-time-equivalent
  result (linear-phase group delay of `order / 2` samples, no acausal
  pre-ringing).

- output_assay:

  Name for the output assay. Default is "filtered".

## Value

A `PhysioExperiment` object with a new assay named `output_assay`
containing the FIR-filtered data. Dimensions match the input assay.

## References

Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems." 2nd ed.
Prentice Hall.

## See also

[`butterworthFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/butterworthFilter.md)
for IIR filtering,
[`notchFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/notchFilter.md)
for power line noise removal,
[`filterSignals()`](https://x-biosignal.github.io/PhysioPreprocess/reference/filterSignals.md)
for moving average filtering,
[`lfilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/lfilter.md)
for the underlying causal filtering primitive.
