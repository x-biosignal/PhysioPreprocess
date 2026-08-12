# Signal Preprocessing Guide

## Introduction

PhysioPreprocess provides a comprehensive set of tools for preprocessing
physiological signals stored in `PhysioExperiment` objects. This
vignette covers the core preprocessing operations: filtering,
resampling, re-referencing, and detrending.

## Setup

``` r

library(PhysioCore)
library(PhysioPreprocess)
```

## Creating Example Data

We begin by constructing a `PhysioExperiment` object with simulated
EEG-like data. The object holds a time-by-channels matrix as its primary
assay along with channel metadata and a sampling rate.

``` r

set.seed(42)
n_time <- 2560
n_channels <- 10
sr <- 256

# Simulate 10 seconds of data at 256 Hz
raw_data <- matrix(rnorm(n_time * n_channels), nrow = n_time, ncol = n_channels)

pe <- PhysioExperiment(
  assays = list(raw = raw_data),
  colData = S4Vectors::DataFrame(
    label = c("Fp1", "Fp2", "F3", "F4", "C3", "C4", "P3", "P4", "O1", "O2"),
    type = rep("EEG", n_channels)
  ),
  samplingRate = sr
)
pe
```

## Filtering

PhysioPreprocess offers several filtering functions for removing
unwanted frequency components from recorded signals.

### Moving Average Filter

The simplest filter is a moving average, which smooths the signal by
averaging over a sliding window. This is useful for quick noise
reduction but offers limited frequency selectivity.

``` r

pe_ma <- filterSignals(pe, window = 5)
```

### Butterworth Filter

The Butterworth filter is an infinite impulse response (IIR) filter
widely used in physiological signal processing. It provides a maximally
flat magnitude response in the passband and supports lowpass, highpass,
bandpass, and bandstop configurations. Zero-phase forward-backward
filtering is applied automatically via
[`signal::filtfilt()`](https://rdrr.io/pkg/signal/man/filtfilt.html).

``` r

# Bandpass filter between 1 and 40 Hz (common for EEG)
pe_bp <- butterworthFilter(pe, low = 1, high = 40, type = "pass")

# Highpass filter at 0.5 Hz to remove slow drift
pe_hp <- butterworthFilter(pe, low = 0.5, type = "high",
                           output_assay = "highpass")

# Lowpass filter at 30 Hz
pe_lp <- butterworthFilter(pe, high = 30, type = "low",
                           output_assay = "lowpass")
```

### FIR Filter

Finite impulse response (FIR) filters offer linear phase and guaranteed
stability. They require a higher order than Butterworth filters for the
same transition bandwidth but are more predictable in their frequency
response.

``` r

pe_fir <- firFilter(pe, low = 1, high = 40, order = 100, type = "pass")
```

### Notch Filter

Power line interference (50 Hz in Europe/Asia, 60 Hz in the Americas) is
a common artifact in electrophysiological recordings. The notch filter
removes the fundamental frequency and optionally its harmonics.

``` r

# Remove 50 Hz power line noise
pe_notch <- notchFilter(pe, freq = 50)

# Remove 60 Hz and its second harmonic
pe_notch60 <- notchFilter(pe, freq = 60, harmonics = 2)
```

## Resampling

Changing the sampling rate is sometimes necessary when combining
datasets recorded at different rates or when reducing data size for
storage.

### Arbitrary Resampling

The
[`resample()`](https://x-biosignal.github.io/PhysioPreprocess/reference/resample.md)
function converts data to any target sampling rate. Three interpolation
methods are available: linear interpolation (fast, adequate for most
cases), spline interpolation (smoother), and FFT-based resampling
(preserves frequency content).

``` r

# Downsample from 256 Hz to 128 Hz
pe_128 <- resample(pe, target_rate = 128, method = "linear")

# Upsample with spline interpolation
pe_512 <- resample(pe, target_rate = 512, method = "spline")
```

### Integer-Factor Decimation

When the downsampling factor is an integer,
[`decimate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/decimate.md)
applies an anti-aliasing lowpass filter before subsampling to prevent
spectral aliasing.

``` r

# Decimate by a factor of 2 (256 Hz -> 128 Hz)
pe_dec <- decimate(pe, factor = 2)
```

### Integer-Factor Interpolation

For integer-factor upsampling,
[`interpolate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/interpolate.md)
inserts new samples using the specified interpolation method.

``` r

# Upsample by a factor of 2 (256 Hz -> 512 Hz)
pe_up <- interpolate(pe, factor = 2, method = "spline")
```

## Re-Referencing

Re-referencing changes the reference electrode for EEG data. The choice
of reference affects spatial interpretation and is an important
preprocessing decision.

### Average Reference

The average reference subtracts the mean of all channels at each time
point. This is the most common choice for high-density EEG because it
approximates a reference-free measure.

``` r

pe_avg <- rereference(pe, ref_type = "average")
getCurrentReference(pe_avg)
```

### Single Channel Reference

A single electrode can serve as the reference. Common choices include
Cz, the nose, or a mastoid electrode.

``` r

pe_c3 <- rereference(pe, ref_type = "channel", ref_channels = "C3")
```

### Linked Reference

Averaging two electrodes (e.g., linked mastoids) provides a symmetric
reference.

``` r

# Assuming M1 and M2 channels exist
# pe_linked <- rereference(pe, ref_type = "channels",
#                          ref_channels = c("M1", "M2"))
```

### Checking the Reference

After re-referencing, the reference information is stored in the object
metadata and can be queried.

``` r

getCurrentReference(pe_avg)
isAverageReferenced(pe_avg)
```

## Detrending

Slow drifts and DC offsets are common in electrophysiological
recordings. Detrending removes these trends without altering faster
signal components.

### Linear and Mean Detrending

The
[`detrendSignals()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detrendSignals.md)
function supports mean removal, linear detrending, and polynomial
detrending.

``` r

# Remove the mean (DC offset)
pe_mean <- detrendSignals(pe, method = "mean")

# Remove a linear trend
pe_linear <- detrendSignals(pe, method = "linear")

# Remove a quadratic trend
pe_poly <- detrendSignals(pe, method = "polynomial", order = 2)
```

### Alternative Detrending

The
[`detrendSignal()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detrendSignal.md)
function in `filters-advanced.R` provides a simpler interface for linear
and constant detrending.

``` r

pe_dt <- detrendSignal(pe, type = "linear")
```

### Baseline Removal

When a known baseline period exists (e.g., a pre-stimulus interval), you
can subtract the baseline mean from the entire signal using
[`removeBaseline()`](https://x-biosignal.github.io/PhysioPreprocess/reference/removeBaseline.md).

``` r

pe_bl <- removeBaseline(pe, baseline_start = 0, baseline_end = 0.2)
```

## Preprocessing Pipelines

For reproducible analysis, PhysioPreprocess supports pipelines that
chain multiple preprocessing steps together. Pipelines are defined as a
sequence of named steps, each specifying a function and its arguments.

### Creating a Pipeline

``` r

pipeline <- createPipeline(
  detrend = list(fn = "detrendSignals", method = "linear"),
  filter  = list(fn = "butterworthFilter", low = 1, high = 40, type = "pass"),
  notch   = list(fn = "notchFilter", freq = 50)
)
print(pipeline)
```

### Applying a Pipeline

``` r

pe_processed <- applyPipeline(pe, pipeline, verbose = TRUE)
```

## Summary

This vignette demonstrated the core preprocessing operations available
in PhysioPreprocess:

- **Filtering**: Moving average, Butterworth, FIR, and notch filters
- **Resampling**: Arbitrary-rate resampling, decimation, and
  interpolation
- **Re-referencing**: Average, single-channel, and multi-channel
  references
- **Detrending**: Mean removal, linear/polynomial detrending, and
  baseline subtraction
- **Pipelines**: Reproducible, chainable preprocessing workflows

For artifact removal using Independent Component Analysis, see the
companion vignette
[`vignette("ica-artifact-removal", package = "PhysioPreprocess")`](https://x-biosignal.github.io/PhysioPreprocess/articles/ica-artifact-removal.md).

## References

- Oppenheim AV, Willsky AS, Nawab SH (1997). *Signals and Systems*. 2nd
  ed. Prentice Hall.
- Nunez PL, Srinivasan R (2006). *Electric Fields of the Brain*. 2nd ed.
  Oxford University Press.
- Crochiere RE, Rabiner LR (1983). *Multirate Digital Signal
  Processing*. Prentice Hall.
