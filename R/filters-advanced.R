#' Advanced signal filtering functions
#'
#' This file provides advanced filtering operations including Butterworth,
#' FIR, and notch filters for physiological signal processing.

#' Butterworth filter
#'
#' Applies a Butterworth filter (lowpass, highpass, bandpass, or bandstop)
#' along the time axis of the specified assay. Uses zero-phase forward-backward
#' filtering via [signal::filtfilt()] to avoid phase distortion.
#'
#' @param x A `PhysioExperiment` object.
#' @param low Lower cutoff frequency in Hz. Required for highpass and bandpass.
#' @param high Upper cutoff frequency in Hz. Required for lowpass and bandpass.
#' @param order Filter order. Default is 4.
#' @param type Filter type: "low", "high", "pass" (bandpass), or "stop" (bandstop).
#' @param output_assay Name for the output assay. Default is "filtered".
#' @return A `PhysioExperiment` object with a new assay named `output_assay`
#'   containing the Butterworth-filtered data. Dimensions match the input assay.
#' @references Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems."
#'   2nd ed. Prentice Hall.
#' @seealso [firFilter()] for FIR filtering, [notchFilter()] for power line
#'   noise removal, [filterSignals()] for moving average filtering,
#'   [detrendSignal()] for trend removal.
#' @export
#' @examples
#' # Create example EEG data
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(1000 * 4), nrow = 1000)),
#'   samplingRate = 250
#' )
#'
#' # Bandpass filter (1-40 Hz) - common for EEG
#' pe <- butterworthFilter(pe, low = 1, high = 40, type = "pass")
#'
#' # Lowpass filter (30 Hz)
#' pe <- butterworthFilter(pe, high = 30, type = "low",
#'                         output_assay = "lowpass")
#'
#' # Highpass filter (0.5 Hz) to remove DC drift
#' pe <- butterworthFilter(pe, low = 0.5, type = "high",
#'                         output_assay = "highpass")
butterworthFilter <- function(x, low = NULL, high = NULL, order = 4L,
                               type = c("pass", "low", "high", "stop"),
                               output_assay = "filtered") {
  stopifnot(inherits(x, "PhysioExperiment"))
  type <- match.arg(type)


  sr <- samplingRate(x)
  if (is.na(sr) || sr <= 0) {
    stop("Valid sampling rate is required for Butterworth filter", call. = FALSE)
  }

  nyquist <- sr / 2

  # Determine filter type and frequencies
  if (type == "low") {
    if (is.null(high)) stop("'high' frequency required for lowpass filter", call. = FALSE)
    W <- high / nyquist
    ftype <- "low"
  } else if (type == "high") {
    if (is.null(low)) stop("'low' frequency required for highpass filter", call. = FALSE)
    W <- low / nyquist
    ftype <- "high"
  } else if (type == "pass") {
    if (is.null(low) || is.null(high)) {
      stop("Both 'low' and 'high' frequencies required for bandpass filter", call. = FALSE)
    }
    W <- c(low / nyquist, high / nyquist)
    ftype <- "pass"
  } else if (type == "stop") {
    if (is.null(low) || is.null(high)) {
      stop("Both 'low' and 'high' frequencies required for bandstop filter", call. = FALSE)
    }
    W <- c(low / nyquist, high / nyquist)
    ftype <- "stop"
  }

  # Validate frequency range

if (any(W <= 0) || any(W >= 1)) {
    stop("Filter frequencies must be between 0 and Nyquist frequency", call. = FALSE)
  }

  # Design filter
  bf <- signal::butter(n = order, W = W, type = ftype)

  assay_name <- defaultAssay(x)
  if (is.na(assay_name)) {
    stop("No assays available to filter", call. = FALSE)
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  dims <- dim(data)

  apply_butter <- function(vec) {
    signal::filtfilt(bf, vec)
  }

  filtered <- data
  if (length(dims) == 1) {
    filtered[] <- apply_butter(data)
  } else {
    filtered[] <- apply(data, seq_along(dims)[-1], apply_butter)
  }

  assays <- SummarizedExperiment::assays(x)
  assays[[output_assay]] <- filtered
  SummarizedExperiment::assays(x) <- assays
  x
}

#' FIR filter
#'
#' Applies a Finite Impulse Response (FIR) filter along the time axis.
#' Uses zero-phase forward-backward filtering via [signal::filtfilt()] to
#' avoid phase distortion.
#'
#' @param x A `PhysioExperiment` object.
#' @param low Lower cutoff frequency in Hz.
#' @param high Upper cutoff frequency in Hz.
#' @param order Filter order (number of taps - 1). Default is 100.
#' @param type Filter type: "low", "high", "pass" (bandpass), or "stop" (bandstop).
#' @param window Window function for FIR design. Default is "hamming".
#' @param output_assay Name for the output assay. Default is "filtered".
#' @return A `PhysioExperiment` object with a new assay named `output_assay`
#'   containing the FIR-filtered data. Dimensions match the input assay.
#' @references Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems."
#'   2nd ed. Prentice Hall.
#' @seealso [butterworthFilter()] for IIR filtering, [notchFilter()] for power
#'   line noise removal, [filterSignals()] for moving average filtering.
#' @export
firFilter <- function(x, low = NULL, high = NULL, order = 100L,
                      type = c("pass", "low", "high", "stop"),
                      window = "hamming", output_assay = "filtered") {
  stopifnot(inherits(x, "PhysioExperiment"))
  type <- match.arg(type)

  sr <- samplingRate(x)
  if (is.na(sr) || sr <= 0) {
    stop("Valid sampling rate is required for FIR filter", call. = FALSE)
  }

  nyquist <- sr / 2

  # Determine filter type and frequencies
  if (type == "low") {
    if (is.null(high)) stop("'high' frequency required for lowpass filter", call. = FALSE)
    W <- high / nyquist
    ftype <- "low"
  } else if (type == "high") {
    if (is.null(low)) stop("'low' frequency required for highpass filter", call. = FALSE)
    W <- low / nyquist
    ftype <- "high"
  } else if (type == "pass") {
    if (is.null(low) || is.null(high)) {
      stop("Both 'low' and 'high' frequencies required for bandpass filter", call. = FALSE)
    }
    W <- c(low / nyquist, high / nyquist)
    ftype <- "pass"
  } else if (type == "stop") {
    if (is.null(low) || is.null(high)) {
      stop("Both 'low' and 'high' frequencies required for bandstop filter", call. = FALSE)
    }
    W <- c(low / nyquist, high / nyquist)
    ftype <- "stop"
  }

  # Validate frequency range
  if (any(W <= 0) || any(W >= 1)) {
    stop("Filter frequencies must be between 0 and Nyquist frequency", call. = FALSE)
  }

  # Design FIR filter
  fir_coef <- signal::fir1(n = order, w = W, type = ftype)

  assay_name <- defaultAssay(x)
  if (is.na(assay_name)) {
    stop("No assays available to filter", call. = FALSE)
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  dims <- dim(data)

  apply_fir <- function(vec) {
    signal::filtfilt(fir_coef, 1, vec)
  }

  filtered <- data
  if (length(dims) == 1) {
    filtered[] <- apply_fir(data)
  } else {
    filtered[] <- apply(data, seq_along(dims)[-1], apply_fir)
  }

  assays <- SummarizedExperiment::assays(x)
  assays[[output_assay]] <- filtered
  SummarizedExperiment::assays(x) <- assays
  x
}

#' Notch filter (power line noise removal)
#'
#' Applies a notch filter to remove power line noise (50 Hz or 60 Hz) and
#' optionally its harmonics. Implemented as a Butterworth bandstop filter
#' with zero-phase filtering.
#'
#' @param x A `PhysioExperiment` object.
#' @param freq Center frequency to remove in Hz. Default is 50 (European power line).
#' @param bandwidth Bandwidth of the notch in Hz. Default is 2.
#' @param harmonics Number of harmonics to remove. Default is 1 (only fundamental).
#' @param output_assay Name for the output assay. Default is "filtered".
#' @return A `PhysioExperiment` object with a new assay named `output_assay`
#'   containing the notch-filtered data. Harmonics above the Nyquist frequency
#'   are skipped with a warning.
#' @references Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems."
#'   2nd ed. Prentice Hall.
#' @seealso [butterworthFilter()] for general Butterworth filtering,
#'   [firFilter()] for FIR filtering, [filterSignals()] for moving average
#'   filtering.
#' @export
#' @examples
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(1000 * 4), nrow = 1000)),
#'   samplingRate = 250
#' )
#'
#' # Remove 50 Hz power line noise (Europe/Asia)
#' pe <- notchFilter(pe, freq = 50)
#'
#' # Remove 60 Hz and harmonics (Americas)
#' pe <- notchFilter(pe, freq = 60, harmonics = 2)
notchFilter <- function(x, freq = 50, bandwidth = 2, harmonics = 1L,
                        output_assay = "filtered") {
  stopifnot(inherits(x, "PhysioExperiment"))

  sr <- samplingRate(x)
  if (is.na(sr) || sr <= 0) {
    stop("Valid sampling rate is required for notch filter", call. = FALSE)
  }

  nyquist <- sr / 2

  assay_name <- defaultAssay(x)
  if (is.na(assay_name)) {
    stop("No assays available to filter", call. = FALSE)
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  dims <- dim(data)
  filtered <- data

  # Apply notch filter for each harmonic
  for (h in seq_len(harmonics)) {
    center_freq <- freq * h

    if (center_freq >= nyquist) {
      warning(sprintf("Harmonic %d (%.1f Hz) exceeds Nyquist frequency; skipping",
                      h, center_freq), call. = FALSE)
      next
    }

    low <- (center_freq - bandwidth / 2) / nyquist
    high <- (center_freq + bandwidth / 2) / nyquist

    if (low <= 0) low <- 0.001
    if (high >= 1) high <- 0.999

    bf <- signal::butter(n = 4, W = c(low, high), type = "stop")

    apply_notch <- function(vec) {
      signal::filtfilt(bf, vec)
    }

    if (length(dims) == 1) {
      filtered[] <- apply_notch(filtered)
    } else {
      filtered[] <- apply(filtered, seq_along(dims)[-1], apply_notch)
    }
  }

  assays <- SummarizedExperiment::assays(x)
  assays[[output_assay]] <- filtered
  SummarizedExperiment::assays(x) <- assays
  x
}

#' Detrend signal
#'
#' Removes linear or polynomial trends from the signal.
#'
#' @param x A `PhysioExperiment` object.
#' @param type Type of detrending: "linear" or "constant" (mean removal).
#' @param output_assay Name for the output assay. Default is "detrended".
#' @return A `PhysioExperiment` object with a new assay named `output_assay`
#'   containing detrended data. For "constant" type, the channel mean is
#'   subtracted. For "linear" type, a least-squares linear fit is removed.
#' @references Oppenheim, A.V. & Willsky, A.S. (1997). "Signals and Systems."
#'   2nd ed. Prentice Hall.
#' @seealso [detrendSignals()] for the alternative detrending implementation
#'   with polynomial support, [butterworthFilter()] for highpass filtering
#'   as an alternative to detrending.
#' @export
detrendSignal <- function(x, type = c("linear", "constant"),
                          output_assay = "detrended") {
  stopifnot(inherits(x, "PhysioExperiment"))
  type <- match.arg(type)

  assay_name <- defaultAssay(x)
  if (is.na(assay_name)) {
    stop("No assays available to detrend", call. = FALSE)
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  dims <- dim(data)
  n <- dims[1]
  t <- seq_len(n)

  detrend_vec <- function(vec) {
    if (type == "constant") {
      vec - mean(vec, na.rm = TRUE)
    } else {
      fit <- stats::lm.fit(cbind(1, t), vec)
      vec - fit$fitted.values
    }
  }

  detrended <- data
  if (length(dims) == 1) {
    detrended[] <- detrend_vec(data)
  } else {
    detrended[] <- apply(data, seq_along(dims)[-1], detrend_vec)
  }

  assays <- SummarizedExperiment::assays(x)
  assays[[output_assay]] <- detrended
  SummarizedExperiment::assays(x) <- assays
  x
}
