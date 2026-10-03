#' PhysioPreprocess: preprocessing for PhysioExperiment objects
#'
#' PhysioPreprocess is the signal-conditioning layer of the ecosystem. It
#' operates directly on [PhysioExperiment::PhysioExperiment] containers,
#' reading an input assay and writing the processed signal into a new named
#' assay so that every step is traceable and the raw data is never overwritten.
#'
#' @section Filtering:
#' [butterworthFilter()] (IIR, SOS or ba form), [firFilter()] (linear-phase
#' FIR), [notchFilter()] (fixed power-line notch), [filterSignals()] (moving
#' average), and [detrendSignal()] / [detrendSignals()] / [removeBaseline()]
#' for trend and offset removal.
#'
#' @section Causal and stateful filtering:
#' [sosDesign()], [sosfilt()], [sosfiltInit()], [sosfiltfilt()], [lfilter()],
#' [lfilterInit()], and [StreamFilter()] for block-wise, state-carrying
#' filtering (the primitives used by the streaming layer).
#'
#' @section Resampling:
#' [resample()] (arbitrary rate), [decimate()] (integer downsampling with
#' anti-aliasing), [interpolate()] (integer upsampling), with per-assay rate
#' tracking via [assaySamplingRates()] and [setAssaySamplingRate()].
#'
#' @section Artifact handling:
#' [detectBadChannels()] and [interpolateBadChannels()], [detectArtifacts()]
#' on continuous data, [rejectBadEpochs()] and [baselineCorrect()] on epoched
#' data, and the configurable [cleanData()] wrapper.
#'
#' @section ICA and ASR:
#' [runICA()] / [removeICAComponents()] and the earlier
#' [icaDecompose()] / [icaRemove()] / [classifyICAComponents()] interface for
#' independent-component artifact removal; [asrCalibrate()], [asrProcess()] and
#' [cleanRawdata()] for Artifact Subspace Reconstruction.
#'
#' @section Line-noise removal:
#' [zapLine()] (ZapLine / DSS) and [cleanLine()] (sinusoidal regression) for
#' deep power-line attenuation without a wide spectral notch.
#'
#' @section Re-referencing and pipelines:
#' [rereference()] with [getCurrentReference()] / [isAverageReferenced()], and
#' [createPipeline()] / [applyPipeline()] to compose reproducible, replayable
#' preprocessing chains.
#'
#' @section Where to go next:
#' The data container itself lives in \pkg{PhysioExperiment}. See the
#' \code{preprocessing-guide} and \code{ica-artifact-removal} vignettes for
#' worked, end-to-end examples.
#'
#' @keywords internal
"_PACKAGE"
