# PhysioPreprocess ![PhysioPreprocess logo](reference/figures/logo.png)

**Preprocessing Functions for PhysioExperiment Objects**

PhysioPreprocess provides a comprehensive preprocessing toolkit for
multi-modal physiological signal data. Built on top of PhysioCore, it
delivers 30 exported functions covering digital filtering (Butterworth,
FIR, notch), resampling and interpolation, ICA-based artifact removal,
bad channel detection, baseline correction, re-referencing, and
reproducible preprocessing pipelines – all operating directly on
`PhysioExperiment` objects.

## Installation

You can install PhysioPreprocess from
[r-universe](https://x-biosignal.r-universe.dev):

``` r

# the containers build on Bioconductor, so its repositories are needed too
install.packages("BiocManager", repos = "https://cloud.r-project.org")
install.packages("PhysioPreprocess",
  repos = c("https://x-biosignal.r-universe.dev", BiocManager::repositories()))
```

Or install the development version from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("x-biosignal/PhysioPreprocess")
```

## Quick Start

``` r

library(PhysioPreprocess)

# Create sample EEG data (4 seconds, 4 channels, 250 Hz)
signal_matrix <- matrix(rnorm(1000 * 4), nrow = 1000, ncol = 4)
colnames(signal_matrix) <- c("Fz", "Cz", "Pz", "Oz")

pe <- PhysioExperiment(
  assays = list(raw = signal_matrix),
  samplingRate = 250
)

# Apply a 1-40 Hz bandpass Butterworth filter
pe <- filterSignals(pe, lowcut = 1, highcut = 40, order = 4)

# Remove 50 Hz line noise with a notch filter
pe <- notchFilter(pe, freq = 50)

# ICA-based artifact removal
ica_result <- runICA(pe, n_components = 4)
pe <- removeICAComponents(pe, ica_result, components = c(1))

# Re-reference to average
pe <- rereference(pe, ref = "average")

# Check reference status
isAverageReferenced(pe)  # TRUE
```

## Features

### Digital Filters

A full suite of frequency-domain filters for physiological signals:

- **Butterworth:**
  [`filterSignals()`](https://x-biosignal.github.io/PhysioPreprocess/reference/filterSignals.md),
  [`butterworthFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/butterworthFilter.md)
  – bandpass, lowpass, and highpass IIR filters with configurable order
- **FIR:**
  [`firFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/firFilter.md)
  – finite impulse response filters for linear-phase filtering
- **Notch:**
  [`notchFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/notchFilter.md)
  – remove power line interference (50/60 Hz) and harmonics
- **Detrending:**
  [`detrendSignal()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detrendSignal.md),
  [`detrendSignals()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detrendSignals.md)
  – remove linear or polynomial trends from signals

### Resampling

Flexible sample rate conversion and multi-rate signal support:

- **Rate conversion:**
  [`resample()`](https://x-biosignal.github.io/PhysioPreprocess/reference/resample.md),
  [`decimate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/decimate.md),
  [`interpolate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/interpolate.md)
  – change sampling rates with anti-aliasing
- **Multi-rate support:**
  [`assaySamplingRates()`](https://x-biosignal.github.io/PhysioPreprocess/reference/assaySamplingRates.md),
  [`setAssaySamplingRate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/setAssaySamplingRate.md)
  – manage per-assay sampling rates for mixed-rate recordings

### Artifact Handling

Automated and semi-automated artifact detection and removal:

- **ICA decomposition:**
  [`icaDecompose()`](https://x-biosignal.github.io/PhysioPreprocess/reference/icaDecompose.md),
  [`runICA()`](https://x-biosignal.github.io/PhysioPreprocess/reference/runICA.md)
  – decompose signals into independent components for artifact
  identification
- **ICA removal:**
  [`icaRemove()`](https://x-biosignal.github.io/PhysioPreprocess/reference/icaRemove.md),
  [`removeICAComponents()`](https://x-biosignal.github.io/PhysioPreprocess/reference/removeICAComponents.md)
  – remove artifact components (e.g., eye blinks, muscle activity) and
  reconstruct clean signals
- **Bad channels:**
  [`detectBadChannels()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detectBadChannels.md),
  [`interpolateBadChannels()`](https://x-biosignal.github.io/PhysioPreprocess/reference/interpolateBadChannels.md)
  – identify noisy or flat channels and interpolate from neighbors
- **Epoch rejection:**
  [`rejectBadEpochs()`](https://x-biosignal.github.io/PhysioPreprocess/reference/rejectBadEpochs.md)
  – reject epochs exceeding amplitude or variance thresholds

### Baseline Correction

Remove baseline drift and DC offsets:

- [`baselineCorrect()`](https://x-biosignal.github.io/PhysioPreprocess/reference/baselineCorrect.md)
  – subtract mean of a specified baseline window from each epoch
- [`removeBaseline()`](https://x-biosignal.github.io/PhysioPreprocess/reference/removeBaseline.md)
  – flexible baseline removal with configurable time windows and methods

### Re-referencing

Electrode re-referencing for EEG and related modalities:

- [`rereference()`](https://x-biosignal.github.io/PhysioPreprocess/reference/rereference.md)
  – re-reference to average, specific electrode(s), or linked mastoids
- [`getCurrentReference()`](https://x-biosignal.github.io/PhysioPreprocess/reference/getCurrentReference.md)
  – query the current reference scheme
- [`isAverageReferenced()`](https://x-biosignal.github.io/PhysioPreprocess/reference/isAverageReferenced.md)
  – check whether average reference has been applied

### Pipeline API

Build reproducible, shareable preprocessing chains:

- [`createPipeline()`](https://x-biosignal.github.io/PhysioPreprocess/reference/createPipeline.md)
  – define an ordered sequence of preprocessing steps with parameters
- [`applyPipeline()`](https://x-biosignal.github.io/PhysioPreprocess/reference/applyPipeline.md)
  – apply a saved pipeline to new data for consistent processing across
  datasets

``` r

# Build a reusable preprocessing pipeline
pipeline <- createPipeline(
  list(filterSignals, lowcut = 1, highcut = 40),
  list(notchFilter, freq = 50),
  list(rereference, ref = "average")
)

# Apply to any PhysioExperiment
pe_clean <- applyPipeline(pe, pipeline)
```

## Dependencies

- **R** (\>= 4.2)
- **PhysioCore** – core data structures
- **signal** – DSP primitives
- **methods**, **SummarizedExperiment**, **S4Vectors**, **stats**
- **Suggests:** fastICA, testthat, knitr, rmarkdown

## PhysioExperiment Ecosystem

PhysioPreprocess is part of the PhysioExperiment ecosystem, a suite of R
packages for multi-modal physiological signal analysis:

| Package | Description |
|----|----|
| [PhysioCore](https://github.com/x-biosignal/PhysioCore) | Core data structures and accessors |
| [PhysioIO](https://github.com/x-biosignal/PhysioIO) | File I/O (EDF, HDF5, BIDS, CSV, MAT) |
| **PhysioPreprocess** | Signal preprocessing and artifact removal |
| [PhysioAnalysis](https://github.com/x-biosignal/PhysioAnalysis) | Spectral analysis, epoching, statistics, visualization |
| [PhysioMoCap](https://github.com/x-biosignal/PhysioMoCap) | Motion capture data processing |
| [PhysioOpenSim](https://github.com/x-biosignal/PhysioOpenSim) | OpenSim biomechanical modeling integration |

Visit the [r-universe page](https://x-biosignal.r-universe.dev) to
browse all available packages.

## License

MIT License. See
[LICENSE](https://x-biosignal.github.io/PhysioPreprocess/LICENSE) for
details.

## Author

Yusuke Matsui

## Governance & support

Part of the [Physio ecosystem](https://x-biosignal.r-universe.dev).
Community and policy documents (Code of Conduct, Contributing,
Governance, Support, Security, and Deprecation & lifecycle) live in the
[umbrella repository](https://github.com/x-biosignal/PhysioExperiment).
