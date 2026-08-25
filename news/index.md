# Changelog

## PhysioPreprocess 0.4.0

### Breaking changes

- [`icaDecompose()`](https://x-biosignal.github.io/PhysioPreprocess/reference/icaDecompose.md)
  now stores the ICA decomposition (components, mixing/unmixing
  matrices, column means) in `metadata(object)$ica` instead of adding an
  `"ica_components"` assay to the returned object. The old assay-based
  storage could not represent a component count different from the
  channel count: `icaDecompose(x, n_components = k)` with
  `k < n_channels` (dimensionality reduction — a core ICA use) errored
  with “all assays must have the same nrow and ncol”, and even
  `k == n_channels` failed on objects whose channels carried names
  (“colnames … not identical”). With the decomposition in metadata, both
  dimensionality reduction and named channels work. Code that read
  `assay(result$object, "ica_components")` should now read
  `metadata(result$object)$ica$components`; the `$components` element of
  the returned list is unchanged.

### Bug fixes

- [`icaRemove()`](https://x-biosignal.github.io/PhysioPreprocess/reference/icaRemove.md)
  reconstructed the signal with the mixing matrix in the wrong
  orientation (`ic %*% mixing` rather than `ic %*% t(mixing)`), which
  only had valid dimensions when `n_components == n_channels` and was
  not a faithful inverse even then. It now uses the correct orientation
  — removing no components round-trips the input to numerical precision
  (~1e-15) — and works for any `n_components <= n_channels`. It also no
  longer reads the (removed) `"ica_components"` assay, and it matches
  the source assay’s dimnames when writing the cleaned assay so
  named-channel objects no longer error.

## PhysioPreprocess 0.3.3

- [`resample()`](https://x-biosignal.github.io/PhysioPreprocess/reference/resample.md)
  now records a W3C-PROV provenance activity (it previously recorded
  none), and
  [`rereference()`](https://x-biosignal.github.io/PhysioPreprocess/reference/rereference.md)
  provenance was standardized from the ad-hoc `withProvenance` call to
  the shared `.recordProv` helper (still one `rereference` activity;
  excluded channels remain in metadata). No numeric change.

## PhysioPreprocess 0.3.2

- Relaxed the `lfilter` C++-vs-pure-R parity test tolerances from 1e-12
  to 1e-8 (the tolerance already used elsewhere in the file). The
  compiled and reference recursions are identical and agree bit-for-bit
  without multiply-add contraction; on an FMA-enabled toolchain
  (r-universe binaries) this order-8 direct-form filter with a 1 Hz band
  edge (near-unit-circle poles) amplifies ~1 ULP/op rounding to ~5e-9.
  Pure FP noise, not an algorithmic change; the well-conditioned SOS
  parity checks stay \<1e-12. No kernel code changed.

## PhysioPreprocess 0.3.1

### Provenance

- Filtering operations now record W3C-PROV provenance like the rest of
  the ecosystem:
  [`butterworthFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/butterworthFilter.md),
  [`notchFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/notchFilter.md),
  [`firFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/firFilter.md),
  [`filterSignals()`](https://x-biosignal.github.io/PhysioPreprocess/reference/filterSignals.md),
  and
  [`detrendSignal()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detrendSignal.md)
  append an activity entry (with their parameters, input/output assays,
  and package version) via `.recordProv()`, so a filtering step is now
  visible in
  [`provenance()`](https://x-biosignal.github.io/PhysioCore//reference/provenance.html)
  and in run-DAG capture. Previously these operations were silent,
  leaving a gap in the provenance record. No filtering behaviour
  changes.

## PhysioPreprocess 0.3.0

### Performance

- The IIR filtering core now runs in compiled code (Rcpp).
  [`sosfilt()`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfilt.md),
  [`lfilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/lfilter.md),
  and the zero-phase SOS path used by
  [`butterworthFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/butterworthFilter.md)
  (`.sosfiltfilt`) call C++ kernels implementing the identical DF2T
  recursions; the previous pure-R loops are retained internally as
  reference implementations and covered by numerical-parity tests
  (agreement ~1e-12 or better).
  [`butterworthFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/butterworthFilter.md)
  and
  [`notchFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/notchFilter.md)
  additionally process plain time x channels matrices in a single C++
  call, parallelized across channels with OpenMP when available (thread
  count respects `OMP_NUM_THREADS` / `OMP_THREAD_LIMIT`). On a
  64-channel, 5-minute, 256 Hz recording, the default 1–40 Hz bandpass
  drops from ~3.2 s to ~0.07 s.
- [`notchFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/notchFilter.md)
  and the ba-form (`use_sos = FALSE`) path of
  [`butterworthFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/butterworthFilter.md)
  now run through a compiled exact port of
  [`signal::filtfilt()`](https://rdrr.io/pkg/signal/man/filtfilt.html)
  (identical zero-padding and zero-state recursion; agreement with
  [`signal::filtfilt()`](https://rdrr.io/pkg/signal/man/filtfilt.html)
  well within the golden tolerances), also matrix-vectorized and
  OpenMP-parallel. The 50 Hz notch on the same recording drops from ~1.4
  s to ~0.08 s.

## PhysioPreprocess 0.2.0

Initial release as a standalone package in the x-biosignal ecosystem,
split out of the monolithic PhysioExperiment codebase. PhysioPreprocess
provides preprocessing operations for physiological signal data held in
`PhysioExperiment` objects (a `SummarizedExperiment` extension from
PhysioCore). All operations write their results into a new, named assay,
leaving the input data intact so that processing history is preserved.

### New Features

- Digital filtering along the time axis for 2D and 3D signal data:
  - [`butterworthFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/butterworthFilter.md)
    for lowpass, highpass, bandpass, and bandstop zero-phase
    (forward-backward) filtering. Uses second-order-sections (SOS) form
    by default for numerical stability at high filter orders, with a
    fallback to the classic transfer-function form via `use_sos`.
  - [`firFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/firFilter.md)
    for zero-phase FIR filtering with configurable window.
  - [`notchFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/notchFilter.md)
    for power-line (50/60 Hz) noise removal, including harmonics, with
    harmonics above Nyquist skipped and warned.
  - [`filterSignals()`](https://x-biosignal.github.io/PhysioPreprocess/reference/filterSignals.md)
    for moving-average smoothing and
    [`detrendSignal()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detrendSignal.md)
    for linear or constant trend removal.
- Resampling and rate management:
  - [`resample()`](https://x-biosignal.github.io/PhysioPreprocess/reference/resample.md)
    to an arbitrary target rate via linear, spline, or FFT-based
    interpolation;
    [`decimate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/decimate.md)
    for integer-factor downsampling with an anti-aliasing lowpass
    filter;
    [`interpolate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/interpolate.md)
    for integer-factor upsampling.
  - [`assaySamplingRates()`](https://x-biosignal.github.io/PhysioPreprocess/reference/assaySamplingRates.md)
    and
    [`setAssaySamplingRate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/setAssaySamplingRate.md)
    to track per-assay sampling rates when assays are resampled
    independently.
- EEG re-referencing via
  [`rereference()`](https://x-biosignal.github.io/PhysioPreprocess/reference/rereference.md),
  supporting common average, single-channel, multi-channel (e.g. linked
  mastoids), and a REST fallback, with optional removal of reference
  channels.
  [`getCurrentReference()`](https://x-biosignal.github.io/PhysioPreprocess/reference/getCurrentReference.md)
  and
  [`isAverageReferenced()`](https://x-biosignal.github.io/PhysioPreprocess/reference/isAverageReferenced.md)
  query the current reference recorded in object metadata.
- Artifact detection and correction:
  - [`detectBadChannels()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detectBadChannels.md)
    (z-score, correlation, and flatline criteria) and
    [`interpolateBadChannels()`](https://x-biosignal.github.io/PhysioPreprocess/reference/interpolateBadChannels.md)
    (average or position-weighted spline).
  - [`detectArtifacts()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detectArtifacts.md)
    for continuous-data amplitude/gradient scanning with automatic
    median + 5\*MAD thresholding and segment merging.
  - [`rejectBadEpochs()`](https://x-biosignal.github.io/PhysioPreprocess/reference/rejectBadEpochs.md)
    and
    [`baselineCorrect()`](https://x-biosignal.github.io/PhysioPreprocess/reference/baselineCorrect.md)
    for epoched data.
- Independent component analysis for artifact removal, with two
  backends:
  - A built-in FastICA implementation
    ([`icaDecompose()`](https://x-biosignal.github.io/PhysioPreprocess/reference/icaDecompose.md)
    /
    [`icaRemove()`](https://x-biosignal.github.io/PhysioPreprocess/reference/icaRemove.md))
    that stores components and mixing matrices in the object, plus
    [`classifyICAComponents()`](https://x-biosignal.github.io/PhysioPreprocess/reference/classifyICAComponents.md)
    to flag artifact components by kurtosis, autocorrelation, or
    high-frequency power.
  - A `fastICA`-package backend
    ([`runICA()`](https://x-biosignal.github.io/PhysioPreprocess/reference/runICA.md)
    /
    [`removeICAComponents()`](https://x-biosignal.github.io/PhysioPreprocess/reference/removeICAComponents.md)).
  - [`cleanData()`](https://x-biosignal.github.io/PhysioPreprocess/reference/cleanData.md)
    chains bad-channel repair and ICA artifact removal into a single
    call.
- Detrending and baseline utilities:
  [`detrendSignals()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detrendSignals.md)
  (linear, mean, or polynomial) and
  [`removeBaseline()`](https://x-biosignal.github.io/PhysioPreprocess/reference/removeBaseline.md)
  (window-based subtraction).
- Composable pipelines:
  [`createPipeline()`](https://x-biosignal.github.io/PhysioPreprocess/reference/createPipeline.md)
  builds a `PhysioPipeline` from named steps and
  [`applyPipeline()`](https://x-biosignal.github.io/PhysioPreprocess/reference/applyPipeline.md)
  runs them in sequence, with a
  [`print()`](https://rdrr.io/r/base/print.html) method for inspection.

### Documentation

- Full roxygen2 documentation and runnable examples for all exported
  functions, with cross-references between related preprocessing steps.
- Re-exports the PhysioCore API so `PhysioExperiment` construction and
  accessors are available directly.
