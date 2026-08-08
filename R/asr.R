# Artifact Subspace Reconstruction (ASR), after Mullen et al. (2015) and the
# EEGLAB clean_rawdata / clean_asr pipeline. ASR learns a clean-baseline
# covariance and per-component amplitude thresholds; in each sliding window it
# finds principal directions whose amplitude exceeds the clean threshold and
# reconstructs that "artifact subspace" from the retained directions using the
# clean-covariance mixing matrix.

# Moore-Penrose pseudoinverse via SVD.
.asr_ginv <- function(A, tol = 1e-9) {
  s <- svd(A)
  keep <- s$d > tol * max(s$d)
  if (!any(keep)) return(matrix(0, ncol(A), nrow(A)))
  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
}

# Symmetric positive-definite matrix square root.
.asr_msqrt <- function(S) {
  e <- eigen((S + t(S)) / 2, symmetric = TRUE)
  v <- pmax(e$values, 1e-12)
  e$vectors %*% diag(sqrt(v), length(v)) %*% t(e$vectors)
}

# Calibrate on a clean matrix (time x channels).
.asr_calibrate_matrix <- function(data, sr, cutoff, window_len) {
  data <- as.matrix(data)
  n <- nrow(data); m <- ncol(data)
  M <- stats::cov(data)
  M <- (M + t(M)) / 2
  V <- eigen(M, symmetric = TRUE)$vectors
  mixing <- .asr_msqrt(M)

  identity <- !is.finite(cutoff)
  if (identity) {
    return(list(M = M, mixing = mixing, T = NULL, V = V, tau = rep(Inf, m),
                sr = sr, window_len = window_len, cutoff = cutoff,
                identity = TRUE))
  }

  # Distribution of per-component RMS over sliding calibration windows.
  Y <- data %*% V
  L <- max(2L, as.integer(round(window_len * sr)))
  starts <- seq(1L, n - L + 1L, by = max(1L, L %/% 2L))
  rms_mat <- vapply(seq_len(m), function(j)
    vapply(starts, function(s) sqrt(mean(Y[s:(s + L - 1L), j]^2)), numeric(1)),
    numeric(length(starts)))
  mu <- apply(rms_mat, 2, stats::median)
  sg <- apply(rms_mat, 2, stats::mad)
  bad <- sg <= 0
  if (any(bad)) sg[bad] <- apply(rms_mat[, bad, drop = FALSE], 2, stats::sd)
  sg[sg <= 0] <- .Machine$double.eps
  tau <- mu + cutoff * sg
  Tmat <- V %*% diag(tau, m) %*% t(V)
  list(M = M, mixing = mixing, T = Tmat, V = V, tau = tau, sr = sr,
       window_len = window_len, cutoff = cutoff, identity = FALSE)
}

# Process a matrix (time x channels) with a calibration.
.asr_process_matrix <- function(data, cal, window_len, step) {
  data <- as.matrix(data)
  n <- nrow(data); m <- ncol(data)
  if (isTRUE(cal$identity)) {
    return(list(cleaned = data, removed_var = 0, frac_rejected = 0))
  }
  if (is.null(window_len)) window_len <- cal$window_len
  L <- max(2L, as.integer(round(window_len * cal$sr)))
  step <- if (is.null(step)) max(1L, L %/% 2L) else max(1L, as.integer(round(step * cal$sr)))
  if (L >= n) L <- n
  hann <- 0.5 - 0.5 * cos(2 * pi * (seq_len(L) - 1) / max(1L, L - 1))

  out <- matrix(0, n, m); wsum <- numeric(n)
  starts <- seq(1L, n - L + 1L, by = step)
  if (starts[length(starts)] != n - L + 1L) starts <- c(starts, n - L + 1L)
  n_rej <- 0; n_comp <- 0
  for (s in starts) {
    idx <- s:(s + L - 1L)
    W <- data[idx, , drop = FALSE]
    C <- stats::cov(W); C <- (C + t(C)) / 2
    eg <- eigen(C, symmetric = TRUE)
    U <- eg$vectors
    rms <- sqrt(pmax(eg$values, 0))
    thr <- sqrt(colSums((cal$T %*% U)^2))            # allowed RMS = ||T u_k||
    keep <- rms <= thr
    keep[is.na(keep)] <- TRUE
    n_rej <- n_rej + sum(!keep); n_comp <- n_comp + m
    if (all(keep)) {
      Wc <- W
    } else {
      Uk <- U[, keep, drop = FALSE]
      R <- cal$mixing %*% .asr_ginv(t(Uk) %*% cal$mixing) %*% t(Uk)
      Wc <- W %*% t(R)
    }
    out[idx, ] <- out[idx, ] + Wc * hann
    wsum[idx] <- wsum[idx] + hann
  }
  wsum[wsum < 1e-12] <- 1
  cleaned <- out / wsum
  list(cleaned = cleaned,
       removed_var = 1 - sum(cleaned^2) / sum(data^2),
       frac_rejected = n_rej / max(1L, n_comp))
}

#' Calibrate Artifact Subspace Reconstruction (ASR)
#'
#' Learns the clean-baseline channel covariance and the per-component RMS
#' thresholds that \code{\link{asrProcess}} uses to detect and reconstruct
#' artifact subspaces (Mullen et al., 2015). Calibration data should be a
#' relatively clean segment; supply \code{calib_window} to restrict it.
#'
#' @param x A PhysioExperiment object with a 2D (time x channels) assay.
#' @param cutoff Rejection cutoff in robust standard deviations; a principal
#'   direction whose window RMS exceeds \code{median + cutoff * MAD} of the
#'   calibration distribution is treated as artifact (default: 20). \code{Inf}
#'   makes \code{\link{asrProcess}} an identity.
#' @param calib_window Optional numeric \code{c(start_sec, end_sec)} selecting
#'   the calibration segment. If \code{NULL} (default) the whole signal is used.
#' @param window_len Sliding-window length in seconds (default: 0.5).
#' @param assay_name Input assay (default: \code{defaultAssay(x)}).
#' @return An object of class \code{"asr_calibration"}: a list with the clean
#'   covariance \code{M}, its square-root mixing matrix, the threshold matrix,
#'   per-component thresholds, and settings.
#' @references
#' Mullen, T. R., et al. (2015). "Real-time neuroimaging and cognitive
#' monitoring using wearable dry EEG." \emph{IEEE Transactions on Biomedical
#' Engineering}, 62(11), 2553-2567. \doi{10.1109/TBME.2015.2481482}
#' @seealso \code{\link{asrProcess}}, \code{\link{cleanRawdata}}
#' @export
asrCalibrate <- function(x, cutoff = 20, calib_window = NULL, window_len = 0.5,
                         assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  if (length(dim(data)) != 2) {
    stop("asrCalibrate requires a 2D (time x channels) assay.", call. = FALSE)
  }
  sr <- samplingRate(x)
  if (!is.null(calib_window)) {
    stopifnot(is.numeric(calib_window), length(calib_window) == 2)
    i0 <- max(1L, as.integer(floor(calib_window[1] * sr)) + 1L)
    i1 <- min(nrow(data), as.integer(floor(calib_window[2] * sr)) + 1L)
    data <- data[i0:i1, , drop = FALSE]
  }
  cal <- .asr_calibrate_matrix(data, sr, cutoff, window_len)
  structure(cal, class = "asr_calibration")
}

#' Apply Artifact Subspace Reconstruction (ASR)
#'
#' Cleans a signal with a calibrated ASR model (\code{\link{asrCalibrate}}). For
#' each overlapping sliding window it eigendecomposes the window covariance,
#' flags principal directions whose RMS exceeds the calibrated threshold, and
#' reconstructs those directions from the retained ones via the clean-covariance
#' mixing matrix; windows are recombined by Hann-weighted overlap-add.
#'
#' @param x A PhysioExperiment object with a 2D (time x channels) assay.
#' @param calibration An \code{"asr_calibration"} from \code{\link{asrCalibrate}}.
#' @param window_len Sliding-window length in seconds (default: the calibration
#'   value).
#' @param step Window step in seconds (default: half the window length).
#' @param assay_name Input assay (default: \code{defaultAssay(x)}).
#' @param output_assay Name for the cleaned assay (default: "asr").
#' @return The PhysioExperiment with the cleaned signal in \code{output_assay};
#'   the removed-variance fraction and rejected-component fraction are stored in
#'   \code{metadata(x)$asr} and logged in provenance.
#' @seealso \code{\link{asrCalibrate}}, \code{\link{cleanRawdata}}
#' @export
asrProcess <- function(x, calibration, window_len = NULL, step = NULL,
                       assay_name = NULL, output_assay = "asr") {
  stopifnot(inherits(x, "PhysioExperiment"),
            inherits(calibration, "asr_calibration"))
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  if (length(dim(data)) != 2) {
    stop("asrProcess requires a 2D (time x channels) assay.", call. = FALSE)
  }
  if (ncol(data) != ncol(calibration$M)) {
    stop("Channel count differs between the data and the ASR calibration.",
         call. = FALSE)
  }
  res <- .asr_process_matrix(data, calibration, window_len, step)

  SummarizedExperiment::assay(x, output_assay) <- res$cleaned
  meta <- S4Vectors::metadata(x)
  meta$asr <- list(removed_var = res$removed_var,
                   frac_rejected = res$frac_rejected,
                   cutoff = calibration$cutoff)
  S4Vectors::metadata(x) <- meta
  x <- PhysioCore::withProvenance(x, x, step = "asrProcess",
                                  params = list(cutoff = calibration$cutoff,
                                                removed_var = res$removed_var))
  x
}

#' Automated raw-data cleaning (EEGLAB clean_rawdata-style)
#'
#' Orchestrates a standard cleaning pipeline mirroring EEGLAB's
#' \code{clean_rawdata}: flat-channel removal, high-pass drift removal,
#' correlation-based bad-channel interpolation, and Artifact Subspace
#' Reconstruction. The removed-variance fraction and any repaired channels are
#' recorded in provenance.
#'
#' @param x A PhysioExperiment object with a 2D (time x channels) assay.
#' @param flatline Logical; detect and interpolate (near-)flat channels
#'   (default: \code{TRUE}).
#' @param highpass Numeric \code{c(low, high)} transition band in Hz for the
#'   drift-removal high-pass; the cutoff is the band midpoint (default
#'   \code{c(0.25, 0.75)}). \code{NULL} skips high-pass filtering.
#' @param channel_crit Minimum mean absolute correlation with the other channels
#'   for a channel to be kept; channels below this are interpolated
#'   (default: 0.8). \code{NULL} skips the correlation criterion.
#' @param asr_cutoff ASR rejection cutoff in robust SDs (default: 20). \code{Inf}
#'   or \code{NULL} skips ASR.
#' @param window_len ASR sliding-window length in seconds (default: 0.5).
#' @param calib_window Optional \code{c(start_sec, end_sec)} calibration segment
#'   for ASR.
#' @param assay_name Input assay (default: \code{defaultAssay(x)}).
#' @param output_assay Name for the cleaned assay (default: "clean").
#' @return The PhysioExperiment with the cleaned signal in \code{output_assay};
#'   \code{metadata(x)$clean_rawdata} records the interpolated channels and the
#'   removed-variance fraction, also logged in provenance.
#' @references
#' Mullen, T. R., et al. (2015). "Real-time neuroimaging and cognitive
#' monitoring using wearable dry EEG." \emph{IEEE Transactions on Biomedical
#' Engineering}, 62(11), 2553-2567. \doi{10.1109/TBME.2015.2481482}
#' @seealso \code{\link{asrCalibrate}}, \code{\link{asrProcess}}
#' @export
cleanRawdata <- function(x, flatline = TRUE, highpass = c(0.25, 0.75),
                         channel_crit = 0.8, asr_cutoff = 20, window_len = 0.5,
                         calib_window = NULL, assay_name = NULL,
                         output_assay = "clean") {
  stopifnot(inherits(x, "PhysioExperiment"))
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  if (length(dim(data)) != 2) {
    stop("cleanRawdata requires a 2D (time x channels) assay.", call. = FALSE)
  }
  sr <- samplingRate(x)
  m <- ncol(data)
  repaired <- integer(0)

  # Interpolate a bad channel by the mean of the good channels.
  interp <- function(mat, bad) {
    if (length(bad) == 0 || length(bad) >= ncol(mat)) return(mat)
    good <- setdiff(seq_len(ncol(mat)), bad)
    gmean <- rowMeans(mat[, good, drop = FALSE])
    for (b in bad) mat[, b] <- gmean
    mat
  }

  # (1) flat channels
  if (isTRUE(flatline)) {
    v <- apply(data, 2, stats::var)
    flat <- which(v < 1e-10 | v < 1e-8 * stats::median(v[v > 0]))
    if (length(flat)) { data <- interp(data, flat); repaired <- union(repaired, flat) }
  }

  # (2) high-pass drift removal (Butterworth), via a temporary object
  if (!is.null(highpass)) {
    stopifnot(is.numeric(highpass), length(highpass) == 2)
    cutoff <- mean(highpass)
    tmp <- PhysioCore::PhysioExperiment(
      assays = list(raw = data),
      colData = SummarizedExperiment::colData(x),
      samplingRate = sr)
    tmp <- butterworthFilter(tmp, low = cutoff, type = "high",
                             output_assay = "._hp")
    data <- SummarizedExperiment::assay(tmp, "._hp")
  }

  # (3) correlation-based bad channels
  if (!is.null(channel_crit)) {
    R <- suppressWarnings(stats::cor(data))
    diag(R) <- NA
    mean_cor <- colMeans(abs(R), na.rm = TRUE)
    bad <- which(mean_cor < channel_crit)
    if (length(bad)) { data <- interp(data, bad); repaired <- union(repaired, bad) }
  }

  # (4) ASR
  removed_var <- 0
  if (!is.null(asr_cutoff)) {
    cal_data <- data
    if (!is.null(calib_window)) {
      i0 <- max(1L, as.integer(floor(calib_window[1] * sr)) + 1L)
      i1 <- min(nrow(data), as.integer(floor(calib_window[2] * sr)) + 1L)
      cal_data <- data[i0:i1, , drop = FALSE]
    }
    cal <- .asr_calibrate_matrix(cal_data, sr, asr_cutoff, window_len)
    res <- .asr_process_matrix(data, cal, window_len, NULL)
    data <- res$cleaned
    removed_var <- res$removed_var
  }

  SummarizedExperiment::assay(x, output_assay) <- data
  meta <- S4Vectors::metadata(x)
  meta$clean_rawdata <- list(
    repaired_channels = repaired,
    removed_var = removed_var,
    asr_cutoff = asr_cutoff)
  S4Vectors::metadata(x) <- meta
  x <- PhysioCore::withProvenance(x, x, step = "cleanRawdata",
                                  params = list(repaired_channels = repaired,
                                                removed_var = removed_var,
                                                asr_cutoff = asr_cutoff))
  x
}
