# Package index

## All functions

- [`StreamFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/StreamFilter.md)
  : Stateful streaming filter object
- [`applyPipeline()`](https://x-biosignal.github.io/PhysioPreprocess/reference/applyPipeline.md)
  : Apply preprocessing pipeline
- [`asrCalibrate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/asrCalibrate.md)
  : Calibrate Artifact Subspace Reconstruction (ASR)
- [`asrProcess()`](https://x-biosignal.github.io/PhysioPreprocess/reference/asrProcess.md)
  : Apply Artifact Subspace Reconstruction (ASR)
- [`assaySamplingRates()`](https://x-biosignal.github.io/PhysioPreprocess/reference/assaySamplingRates.md)
  : Per-assay sampling rates (moved to PhysioCore)
- [`baselineCorrect()`](https://x-biosignal.github.io/PhysioPreprocess/reference/baselineCorrect.md)
  : Baseline correction
- [`butterworthFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/butterworthFilter.md)
  : Butterworth filter
- [`classifyICAComponents()`](https://x-biosignal.github.io/PhysioPreprocess/reference/classifyICAComponents.md)
  : Classify ICA components as brain or artifact
- [`cleanData()`](https://x-biosignal.github.io/PhysioPreprocess/reference/cleanData.md)
  : Clean data using an artifact removal pipeline
- [`cleanLine()`](https://x-biosignal.github.io/PhysioPreprocess/reference/cleanLine.md)
  : Remove line noise with CleanLine
- [`cleanRawdata()`](https://x-biosignal.github.io/PhysioPreprocess/reference/cleanRawdata.md)
  : Automated raw-data cleaning (EEGLAB clean_rawdata-style)
- [`createPipeline()`](https://x-biosignal.github.io/PhysioPreprocess/reference/createPipeline.md)
  : Preprocessing Pipeline for PhysioExperiment
- [`decimate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/decimate.md)
  : Decimate signal
- [`detectArtifacts()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detectArtifacts.md)
  : Detect artifacts in continuous data
- [`detectBadChannels()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detectBadChannels.md)
  : Detect bad channels
- [`detrendSignal()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detrendSignal.md)
  : Detrend signal
- [`detrendSignals()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detrendSignals.md)
  : Signal Detrending Functions for PhysioExperiment
- [`filterSignals()`](https://x-biosignal.github.io/PhysioPreprocess/reference/filterSignals.md)
  : Moving average filter
- [`filters-causal`](https://x-biosignal.github.io/PhysioPreprocess/reference/filters-causal.md)
  : Causal and stateful IIR filtering primitives
- [`firFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/firFilter.md)
  : FIR filter
- [`getCurrentReference()`](https://x-biosignal.github.io/PhysioPreprocess/reference/getCurrentReference.md)
  : Get current reference
- [`icaDecompose()`](https://x-biosignal.github.io/PhysioPreprocess/reference/icaDecompose.md)
  : ICA and Artifact Removal for PhysioExperiment
- [`icaRemove()`](https://x-biosignal.github.io/PhysioPreprocess/reference/icaRemove.md)
  : Remove ICA components
- [`interpolate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/interpolate.md)
  : Interpolate signal
- [`interpolateBadChannels()`](https://x-biosignal.github.io/PhysioPreprocess/reference/interpolateBadChannels.md)
  : Interpolate bad channels
- [`isAverageReferenced()`](https://x-biosignal.github.io/PhysioPreprocess/reference/isAverageReferenced.md)
  : Check if data is average referenced
- [`lfilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/lfilter.md)
  : One-pass causal IIR filtering (Direct Form II Transposed)
- [`lfilterInit()`](https://x-biosignal.github.io/PhysioPreprocess/reference/lfilterInit.md)
  : Steady-state initial conditions for a one-pass IIR filter
- [`notchFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/notchFilter.md)
  : Notch filter (power line noise removal)
- [`print(`*`<PhysioPipeline>`*`)`](https://x-biosignal.github.io/PhysioPreprocess/reference/print.PhysioPipeline.md)
  : Print pipeline summary
- [`rejectBadEpochs()`](https://x-biosignal.github.io/PhysioPreprocess/reference/rejectBadEpochs.md)
  : Reject bad epochs
- [`removeBaseline()`](https://x-biosignal.github.io/PhysioPreprocess/reference/removeBaseline.md)
  : Remove baseline
- [`removeICAComponents()`](https://x-biosignal.github.io/PhysioPreprocess/reference/removeICAComponents.md)
  : Remove ICA components
- [`rereference()`](https://x-biosignal.github.io/PhysioPreprocess/reference/rereference.md)
  : Re-referencing Operations for EEG Data
- [`resample()`](https://x-biosignal.github.io/PhysioPreprocess/reference/resample.md)
  : Resampling operations for PhysioExperiment
- [`runICA()`](https://x-biosignal.github.io/PhysioPreprocess/reference/runICA.md)
  : ICA-Based Artifact Removal for PhysioExperiment
- [`setAssaySamplingRate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/setAssaySamplingRate.md)
  : Set a per-assay sampling rate (moved to PhysioCore)
- [`sosDesign()`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosDesign.md)
  : Design a Butterworth filter in second-order-sections (SOS) form
- [`sosfilt()`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfilt.md)
  : One-pass causal SOS filtering with optional carried state
- [`sosfiltInit()`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfiltInit.md)
  : Steady-state initial conditions for a cascaded SOS filter
- [`sosfiltfilt()`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfiltfilt.md)
  : Zero-phase SOS filtering with state-initialised edge handling
- [`zapLine()`](https://x-biosignal.github.io/PhysioPreprocess/reference/zapLine.md)
  : Remove line noise with ZapLine
