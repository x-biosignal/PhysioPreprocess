#' Signal Detrending Functions for PhysioExperiment
#'
#' Functions for removing trends from physiological signals.

#' Detrend signals
#'
#' Removes linear or polynomial trends from signals.
#'
#' @param pe A PhysioExperiment object.
#' @param method Detrending method: "linear", "mean", or "polynomial".
#' @param order Polynomial order for method="polynomial".
#' @param assay_name Name of the assay to detrend.
#' @param output_assay Name for the detrended output assay.
#' @return PhysioExperiment with detrended data in output_assay.
#' @references Oppenheim AV, Willsky AS, Nawab SH (1997). "Signals and Systems."
#'   2nd ed. Prentice Hall.
#' @seealso [removeBaseline()] for baseline subtraction over a time window,
#'   [detrendSignal()] for the alternative linear/constant detrending
#'   implementation, [butterworthFilter()] for highpass filtering as an
#'   alternative to detrending.
#' @export
#' @examples
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(1000) + 1:100/10, nrow = 100, ncol = 10)),
#'   colData = S4Vectors::DataFrame(label = paste0("Ch", 1:10)),
#'   samplingRate = 256
#' )
#' pe_detrended <- detrendSignals(pe, method = "linear")
detrendSignals <- function(pe, method = c("linear", "mean", "polynomial"),
                           order = 2, assay_name = NULL,
                           output_assay = "detrended") {
  stopifnot(inherits(pe, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(assay_name)) {
    assay_name <- defaultAssay(pe)
  }

  data <- SummarizedExperiment::assay(pe, assay_name)

  # Ensure 2D
  if (length(dim(data)) > 2) {
    data <- data[, , 1]
  }

  n_samples <- nrow(data)
  n_channels <- ncol(data)
  time_vec <- seq_len(n_samples)

  detrended <- matrix(0, nrow = n_samples, ncol = n_channels)

  for (ch in seq_len(n_channels)) {
    signal <- data[, ch]

    if (method == "mean") {
      detrended[, ch] <- signal - mean(signal, na.rm = TRUE)
    } else if (method == "linear") {
      fit <- lm(signal ~ time_vec)
      detrended[, ch] <- residuals(fit)
    } else if (method == "polynomial") {
      fit <- lm(signal ~ poly(time_vec, order))
      detrended[, ch] <- residuals(fit)
    }
  }

  SummarizedExperiment::assay(pe, output_assay) <- detrended

  pe
}

#' Remove baseline
#'
#' Subtracts baseline from a specified time window.
#'
#' @param pe A PhysioExperiment object.
#' @param baseline_start Start time of baseline window in seconds.
#' @param baseline_end End time of baseline window in seconds.
#' @param assay_name Name of the assay to baseline correct.
#' @param output_assay Name for the baseline-corrected output assay.
#' @return PhysioExperiment with baseline-corrected data.
#' @references Oppenheim AV, Willsky AS, Nawab SH (1997). "Signals and Systems."
#'   2nd ed. Prentice Hall.
#' @seealso [detrendSignals()] for trend removal, [baselineCorrect()] for
#'   epoch-based baseline correction, [filterSignals()] for moving average
#'   smoothing.
#' @export
#' @examples
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(1000) + 5, nrow = 100, ncol = 10)),
#'   colData = S4Vectors::DataFrame(label = paste0("Ch", 1:10)),
#'   samplingRate = 100
#' )
#' pe_corrected <- removeBaseline(pe, baseline_start = 0, baseline_end = 0.2)
removeBaseline <- function(pe, baseline_start, baseline_end,
                           assay_name = NULL,
                           output_assay = "baseline_corrected") {
  stopifnot(inherits(pe, "PhysioExperiment"))

  if (is.null(assay_name)) {
    assay_name <- defaultAssay(pe)
  }

  sr <- samplingRate(pe)
  if (is.na(sr)) {
    stop("Valid sampling rate required for baseline removal", call. = FALSE)
  }

  data <- SummarizedExperiment::assay(pe, assay_name)

  # Ensure 2D
  if (length(dim(data)) > 2) {
    data <- data[, , 1]
  }

  # Convert time to samples
  start_sample <- max(1, floor(baseline_start * sr) + 1)
  end_sample <- min(nrow(data), ceiling(baseline_end * sr) + 1)

  # Calculate baseline mean for each channel
  baseline_means <- colMeans(data[start_sample:end_sample, , drop = FALSE],
                             na.rm = TRUE)

  # Subtract baseline
  corrected <- sweep(data, 2, baseline_means, "-")

  SummarizedExperiment::assay(pe, output_assay) <- corrected

  pe
}
