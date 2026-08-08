# Spectral line-noise removal that, unlike a fixed notch, removes the line while
# preserving the broadband spectrum: ZapLine (de Cheveigne, 2020) via a
# DSS/PCA line subspace, and CleanLine (Mullen; Mitra & Bokil) via sliding-window
# sinusoidal regression with a Thomson F-test. The existing notchFilter() remains
# the fast option (a fixed notch is cheaper but attenuates a wider band).

# Moore-Penrose pseudoinverse via SVD.
.line_ginv <- function(A, tol = 1e-9) {
  s <- svd(A)
  keep <- s$d > tol * max(s$d)
  if (!any(keep)) return(matrix(0, ncol(A), nrow(A)))
  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
}

# ZapLine on a matrix (time x channels).
.zapline_matrix <- function(data, sr, line_freq, nremove, nfft) {
  data <- as.matrix(data)
  n <- nrow(data); m <- ncol(data)
  Lm <- max(2L, as.integer(round(sr / line_freq)))       # moving-average -> comb notch

  smooth1 <- function(v) {
    padded <- c(rev(v[seq_len(Lm)]), v, rev(v[(n - Lm + 1):n]))
    stats::filter(padded, rep(1 / Lm, Lm), sides = 2)[(Lm + 1):(Lm + n)]
  }
  xs <- apply(data, 2, function(v) as.numeric(smooth1(as.numeric(smooth1(v)))))
  resid <- data - xs                                     # line-rich residual

  if (m == 1L) {
    # No spatial subspace to separate; the comb-notched signal is the result.
    return(list(cleaned = xs, removed_power = sum(resid^2) / sum(data^2)))
  }

  # Bias the covariance toward the line band to order DSS components by
  # line-noise concentration.
  halfbw <- if (is.null(nfft)) 1 else max(0.5, sr / nfft)
  f <- (0:(n - 1)) * sr / n
  keep_bins <- abs(pmin(f, sr - f) - line_freq) <= halfbw
  resid_bp <- apply(resid, 2, function(v)
    Re(stats::fft(stats::fft(v) * keep_bins, inverse = TRUE)) / n)

  c0 <- stats::cov(resid); c0 <- (c0 + t(c0)) / 2
  c1 <- stats::cov(resid_bp); c1 <- (c1 + t(c1)) / 2
  e0 <- eigen(c0, symmetric = TRUE)
  W <- e0$vectors %*% diag(1 / sqrt(pmax(e0$values, 1e-12)), m)   # whitening
  cc <- t(W) %*% c1 %*% W
  todss <- W %*% eigen((cc + t(cc)) / 2, symmetric = TRUE)$vectors

  nremove <- min(as.integer(nremove), m - 1L)
  keep <- seq_len(max(1L, nremove))
  todss_inv <- .line_ginv(todss)
  comp <- resid %*% todss
  line_est <- comp[, keep, drop = FALSE] %*% todss_inv[keep, , drop = FALSE]
  cleaned <- data - line_est
  list(cleaned = cleaned, removed_power = sum(line_est^2) / sum(data^2))
}

# CleanLine on a matrix (time x channels).
.cleanline_matrix <- function(data, sr, line_freq, harmonics, p_thresh,
                              window_sec) {
  data <- as.matrix(data)
  n <- nrow(data); m <- ncol(data)
  L <- min(n, max(2L, as.integer(round(window_sec * sr))))
  step <- max(1L, L %/% 2L)
  hann <- 0.5 - 0.5 * cos(2 * pi * (seq_len(L) - 1) / max(1L, L - 1))
  freqs <- line_freq * seq_len(harmonics)
  k <- 2L * harmonics

  out <- matrix(0, n, m); wsum <- numeric(n)
  removed_power <- 0
  starts <- unique(c(seq(1L, n - L + 1L, by = step), n - L + 1L))
  for (s in starts) {
    idx <- s:(s + L - 1L)
    tt <- (seq_len(L) - 1) / sr
    design <- matrix(1, L, 1)
    for (fr in freqs) design <- cbind(design, cos(2 * pi * fr * tt),
                                      sin(2 * pi * fr * tt))
    dfres <- L - (k + 1L)
    for (ch in seq_len(m)) {
      seg <- data[idx, ch]
      fit <- stats::lm.fit(design, seg)
      b <- fit$coefficients
      line_part <- design[, -1, drop = FALSE] %*% b[-1]
      rss_full <- sum(fit$residuals^2)
      rss_red <- sum((seg - mean(seg))^2)
      Fstat <- ((rss_red - rss_full) / k) / (rss_full / dfres)
      pval <- stats::pf(Fstat, k, dfres, lower.tail = FALSE)
      if (is.finite(pval) && pval < p_thresh) {
        out[idx, ch] <- out[idx, ch] + (seg - line_part) * hann
        removed_power <- removed_power + sum((line_part * hann)^2)
      } else {
        out[idx, ch] <- out[idx, ch] + seg * hann
      }
    }
    wsum[idx] <- wsum[idx] + hann
  }
  wsum[wsum < 1e-12] <- 1
  list(cleaned = out / wsum, removed_power = removed_power / sum(data^2))
}

#' Remove line noise with ZapLine
#'
#' Removes power-line noise while preserving the broadband spectrum using
#' ZapLine (de Cheveigne, 2020): the line is first suppressed by a comb-notch
#' smoother, then a DSS/PCA of the notched-out residual isolates the line-noise
#' spatial subspace, and only its top \code{nremove} components are removed from
#' the data. This avoids the wide spectral notch of a fixed filter and needs at
#' least two channels (line noise is treated as a spatial component).
#'
#' @param x A PhysioExperiment object with a 2D (time x channels) assay.
#' @param line_freq Line frequency in Hz (default: 50).
#' @param nremove Number of line-noise components to remove (default: 1).
#' @param nfft Optional analysis FFT length controlling the line-band width for
#'   the DSS bias (default: \code{NULL} = 1 Hz half-band).
#' @param assay_name Input assay (default: \code{defaultAssay(x)}).
#' @param output_assay Name for the cleaned assay (default: "zapline").
#' @return The PhysioExperiment with the cleaned signal in \code{output_assay};
#'   the removed line-power fraction is stored in \code{metadata(x)$line_noise}
#'   and logged in provenance.
#' @references
#' de Cheveigne, A. (2020). "ZapLine: A simple and effective method to remove
#' power line artifacts." \emph{NeuroImage}, 207, 116356.
#' \doi{10.1016/j.neuroimage.2019.116356}
#' @seealso \code{\link{cleanLine}}, \code{\link{notchFilter}} (the fast fixed
#'   notch)
#' @export
zapLine <- function(x, line_freq = 50, nremove = 1, nfft = NULL,
                    assay_name = NULL, output_assay = "zapline") {
  stopifnot(inherits(x, "PhysioExperiment"))
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  if (length(dim(data)) != 2) {
    stop("zapLine requires a 2D (time x channels) assay.", call. = FALSE)
  }
  res <- .zapline_matrix(data, samplingRate(x), line_freq, nremove, nfft)
  x <- .store_line_result(x, res, output_assay, "zapLine",
                          list(line_freq = line_freq, nremove = nremove))
  x
}

#' Remove line noise with CleanLine
#'
#' Removes power-line noise (and optional harmonics) by sliding-window
#' sinusoidal regression with a Thomson F-test (the CleanLine / Chronux
#' approach; Mitra & Bokil, 2008). In each Hann-overlapped window each channel is
#' regressed onto the line sine/cosine pair; when the fit is significant
#' (\code{p < p_thresh}) the fitted sinusoid is subtracted. This attenuates the
#' line deeply without a wide spectral notch and works on any number of channels.
#'
#' @param x A PhysioExperiment object with a 2D (time x channels) assay.
#' @param line_freq Line frequency in Hz (default: 50).
#' @param bandwidth Reserved spectral bandwidth in Hz around the line
#'   (default: 2).
#' @param harmonics Number of harmonics to remove, including the fundamental
#'   (default: 1, i.e. \code{line_freq} only; 3 removes 50/100/150 Hz).
#' @param p_thresh F-test significance threshold for subtracting the line
#'   (default: 0.01).
#' @param window_sec Sliding-window length in seconds (default: 4).
#' @param assay_name Input assay (default: \code{defaultAssay(x)}).
#' @param output_assay Name for the cleaned assay (default: "cleanline").
#' @return The PhysioExperiment with the cleaned signal in \code{output_assay};
#'   the removed line-power fraction is stored in \code{metadata(x)$line_noise}
#'   and logged in provenance.
#' @references
#' Mitra, P., & Bokil, H. (2008). \emph{Observed Brain Dynamics}. Oxford
#' University Press.
#' @seealso \code{\link{zapLine}}, \code{\link{notchFilter}} (the fast fixed
#'   notch)
#' @export
cleanLine <- function(x, line_freq = 50, bandwidth = 2, harmonics = 1,
                      p_thresh = 0.01, window_sec = 4, assay_name = NULL,
                      output_assay = "cleanline") {
  stopifnot(inherits(x, "PhysioExperiment"))
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  if (length(dim(data)) != 2) {
    stop("cleanLine requires a 2D (time x channels) assay.", call. = FALSE)
  }
  res <- .cleanline_matrix(data, samplingRate(x), line_freq, harmonics,
                           p_thresh, window_sec)
  x <- .store_line_result(x, res, output_assay, "cleanLine",
                          list(line_freq = line_freq, harmonics = harmonics))
  x
}

# Store a line-removal result: assay, metadata, provenance.
.store_line_result <- function(x, res, output_assay, method, params) {
  SummarizedExperiment::assay(x, output_assay) <- res$cleaned
  meta <- S4Vectors::metadata(x)
  meta$line_noise <- c(list(method = method, removed_power = res$removed_power),
                       params)
  S4Vectors::metadata(x) <- meta
  PhysioCore::withProvenance(x, x, step = method,
                             params = c(params, list(removed_power = res$removed_power)))
}
