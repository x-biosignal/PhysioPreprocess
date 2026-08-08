# PhysioPreprocess 0.3.1

## Provenance

- Filtering operations now record W3C-PROV provenance like the rest of the
  ecosystem: `butterworthFilter()`, `notchFilter()`, `firFilter()`,
  `filterSignals()`, and `detrendSignal()` append an activity entry (with their
  parameters, input/output assays, and package version) via `.recordProv()`, so
  a filtering step is now visible in `provenance()` and in run-DAG capture.
  Previously these operations were silent, leaving a gap in the provenance
  record. No filtering behaviour changes.

# PhysioPreprocess 0.3.0

## Performance

- The IIR filtering core now runs in compiled code (Rcpp). `sosfilt()`,
  `lfilter()`, and the zero-phase SOS path used by `butterworthFilter()`
  (`.sosfiltfilt`) call C++ kernels implementing the identical DF2T
  recursions; the previous pure-R loops are retained internally as
  reference implementations and covered by numerical-parity tests
  (agreement ~1e-12 or better). `butterworthFilter()` and `notchFilter()`
  additionally process plain time x channels matrices in a single C++
  call, parallelized across channels with OpenMP when available (thread
  count respects `OMP_NUM_THREADS` / `OMP_THREAD_LIMIT`). On a 64-channel,
  5-minute, 256 Hz recording, the default 1--40 Hz bandpass drops from
  ~3.2 s to ~0.07 s.
- `notchFilter()` and the ba-form (`use_sos = FALSE`) path of
  `butterworthFilter()` now run through a compiled exact port of
  `signal::filtfilt()` (identical zero-padding and zero-state recursion;
  agreement with `signal::filtfilt()` well within the golden tolerances),
  also matrix-vectorized and OpenMP-parallel. The 50 Hz notch on the same
  recording drops from ~1.4 s to ~0.08 s.

# PhysioPreprocess 0.2.0

Initial release as a standalone package in the x-biosignal ecosystem, split
out of the monolithic PhysioExperiment codebase. PhysioPreprocess provides
preprocessing operations for physiological signal data held in
`PhysioExperiment` objects (a `SummarizedExperiment` extension from
PhysioCore). All operations write their results into a new, named assay,
leaving the input data intact so that processing history is preserved.

## New Features

- Digital filtering along the time axis for 2D and 3D signal data:
  - `butterworthFilter()` for lowpass, highpass, bandpass, and bandstop
    zero-phase (forward-backward) filtering. Uses second-order-sections
    (SOS) form by default for numerical stability at high filter orders,
    with a fallback to the classic transfer-function form via `use_sos`.
  - `firFilter()` for zero-phase FIR filtering with configurable window.
  - `notchFilter()` for power-line (50/60 Hz) noise removal, including
    harmonics, with harmonics above Nyquist skipped and warned.
  - `filterSignals()` for moving-average smoothing and `detrendSignal()`
    for linear or constant trend removal.
- Resampling and rate management:
  - `resample()` to an arbitrary target rate via linear, spline, or
    FFT-based interpolation; `decimate()` for integer-factor downsampling
    with an anti-aliasing lowpass filter; `interpolate()` for
    integer-factor upsampling.
  - `assaySamplingRates()` and `setAssaySamplingRate()` to track per-assay
    sampling rates when assays are resampled independently.
- EEG re-referencing via `rereference()`, supporting common average,
  single-channel, multi-channel (e.g. linked mastoids), and a REST
  fallback, with optional removal of reference channels.
  `getCurrentReference()` and `isAverageReferenced()` query the current
  reference recorded in object metadata.
- Artifact detection and correction:
  - `detectBadChannels()` (z-score, correlation, and flatline criteria)
    and `interpolateBadChannels()` (average or position-weighted spline).
  - `detectArtifacts()` for continuous-data amplitude/gradient scanning
    with automatic median + 5*MAD thresholding and segment merging.
  - `rejectBadEpochs()` and `baselineCorrect()` for epoched data.
- Independent component analysis for artifact removal, with two backends:
  - A built-in FastICA implementation (`icaDecompose()` /
    `icaRemove()`) that stores components and mixing matrices in the
    object, plus `classifyICAComponents()` to flag artifact components by
    kurtosis, autocorrelation, or high-frequency power.
  - A `fastICA`-package backend (`runICA()` / `removeICAComponents()`).
  - `cleanData()` chains bad-channel repair and ICA artifact removal into
    a single call.
- Detrending and baseline utilities: `detrendSignals()` (linear, mean, or
  polynomial) and `removeBaseline()` (window-based subtraction).
- Composable pipelines: `createPipeline()` builds a `PhysioPipeline` from
  named steps and `applyPipeline()` runs them in sequence, with a `print()`
  method for inspection.

## Documentation

- Full roxygen2 documentation and runnable examples for all exported
  functions, with cross-references between related preprocessing steps.
- Re-exports the PhysioCore API so `PhysioExperiment` construction and
  accessors are available directly.
