#' Causal and stateful IIR filtering primitives
#'
#' This file provides one-pass, causal, state-carrying IIR filtering that
#' mirrors \code{scipy.signal} (\code{lfilter}, \code{lfilter_zi},
#' \code{sosfilt}, \code{sosfilt_zi}, \code{sosfiltfilt}). Unlike the zero-phase
#' \code{\link{butterworthFilter}} default (which uses forward-backward
#' \code{filtfilt} and is inherently non-causal), these functions process a
#' signal strictly forward in time and can carry filter state across
#' consecutive chunks. This is the numerical core that a real-time streaming
#' backend (e.g. a future \code{PhysioStream}) reuses so that offline and
#' online filtering produce identical results.
#'
#' All recursions use the Direct-Form-II-Transposed (DF2T) structure, which is
#' the same structure \code{scipy.signal} and \code{signal::filter} use, so the
#' outputs agree to machine precision.
#'
#' @name filters-causal
#' @references
#' Oppenheim, A.V. & Schafer, R.W. (2010). \emph{Discrete-Time Signal
#' Processing}, 3rd ed. Prentice Hall.
#'
#' Virtanen, P. et al. (2020). SciPy 1.0. \emph{Nature Methods} 17, 261-272.
#' (reference implementation: \code{scipy.signal.sosfilt} /
#' \code{lfilter} / \code{lfilter_zi} / \code{sosfilt_zi}).
NULL

# ---------------------------------------------------------------------------
# Internal validation helpers
# ---------------------------------------------------------------------------

.check_sos <- function(sos) {
  if (is.null(dim(sos)) || length(dim(sos)) != 2L || ncol(sos) != 6L) {
    stop("`sos` must be a matrix with 6 columns (b0, b1, b2, a0, a1, a2).",
         call. = FALSE)
  }
  if (!is.numeric(sos)) {
    stop("`sos` must be a numeric matrix.", call. = FALSE)
  }
  if (any(sos[, 4] == 0)) {
    stop("Every SOS section must have a non-zero leading denominator (a0).",
         call. = FALSE)
  }
  invisible(sos)
}

# ---------------------------------------------------------------------------
# Filter design
# ---------------------------------------------------------------------------

#' Design a Butterworth filter in second-order-sections (SOS) form
#'
#' Builds the cascaded second-order-sections representation of a Butterworth
#' filter directly from the analog prototype poles (bilinear transform), which
#' is numerically stable for all orders. This is the design entry point for the
#' causal / stateful filtering functions (\code{\link{sosfilt}},
#' \code{\link{sosfiltfilt}}, \code{\link{StreamFilter}}).
#'
#' @param low Lower cutoff frequency in Hz. Required for \code{"high"} and
#'   \code{"pass"}/\code{"stop"}.
#' @param high Upper cutoff frequency in Hz. Required for \code{"low"} and
#'   \code{"pass"}/\code{"stop"}.
#' @param order Filter order (per band). Default \code{4}.
#' @param type Filter type: \code{"low"}, \code{"high"}, \code{"pass"}
#'   (bandpass), or \code{"stop"} (bandstop).
#' @param sr Sampling rate in Hz.
#' @return A numeric matrix with one second-order section per row and six
#'   columns \code{c(b0, b1, b2, a0, a1, a2)}.
#' @seealso \code{\link{sosfilt}}, \code{\link{sosfiltfilt}},
#'   \code{\link{StreamFilter}}
#' @export
#' @examples
#' sos <- sosDesign(high = 40, type = "low", sr = 250, order = 4)
#' y <- sosfilt(sos, rnorm(1000))
sosDesign <- function(low = NULL, high = NULL, order = 4L,
                      type = c("pass", "low", "high", "stop"), sr) {
  type <- match.arg(type)
  if (missing(sr) || !is.numeric(sr) || length(sr) != 1L || is.na(sr) || sr <= 0) {
    stop("`sr` (sampling rate in Hz) must be a positive number.", call. = FALSE)
  }
  order <- as.integer(order)
  if (is.na(order) || order < 1L) {
    stop("`order` must be a positive integer.", call. = FALSE)
  }
  nyquist <- sr / 2

  if (type == "low") {
    if (is.null(high)) stop("'high' frequency required for lowpass filter", call. = FALSE)
    W <- high / nyquist
  } else if (type == "high") {
    if (is.null(low)) stop("'low' frequency required for highpass filter", call. = FALSE)
    W <- low / nyquist
  } else { # pass / stop
    if (is.null(low) || is.null(high)) {
      stop(sprintf("Both 'low' and 'high' frequencies required for %s filter", type),
           call. = FALSE)
    }
    W <- c(low / nyquist, high / nyquist)
  }

  if (any(W <= 0) || any(W >= 1)) {
    stop("Filter frequencies must be strictly between 0 and the Nyquist frequency.",
         call. = FALSE)
  }
  if (length(W) == 2L && W[1] >= W[2]) {
    stop("'low' must be below 'high' for band filters.", call. = FALSE)
  }

  sos <- .buttersos(n = order, W = W, type = type)
  if (is.null(dim(sos))) sos <- matrix(sos, nrow = 1L)
  .check_sos(sos)
  sos
}

# ---------------------------------------------------------------------------
# Steady-state initial conditions (lfilter_zi / sosfilt_zi)
# ---------------------------------------------------------------------------

#' Steady-state initial conditions for a one-pass IIR filter
#'
#' Computes the filter delay-line state that, when scaled by the first input
#' sample and passed to \code{\link{lfilter}}, makes the filter's step response
#' start already in steady state (no start-up transient). This is the direct
#' analogue of \code{scipy.signal.lfilter_zi}.
#'
#' @param b Numerator (feed-forward) coefficients.
#' @param a Denominator (feedback) coefficients; \code{a[1]} need not be 1
#'   (coefficients are normalised internally).
#' @return A numeric vector of length \code{max(length(a), length(b)) - 1}, the
#'   initial delay-line state (unscaled; multiply by the first sample level).
#' @seealso \code{\link{lfilter}}, \code{\link{sosfiltInit}}
#' @export
#' @examples
#' ba <- signal::butter(2, 0.2, "low")
#' zi <- lfilterInit(ba$b, ba$a)
lfilterInit <- function(b, a) {
  b <- as.numeric(b)
  a <- as.numeric(a)
  # Drop leading zeros in a (keep at least one coefficient).
  while (length(a) > 1L && a[1] == 0) a <- a[-1]
  if (length(a) < 1L || a[1] == 0) {
    stop("`a` must have a non-zero leading coefficient.", call. = FALSE)
  }
  if (a[1] != 1) {
    b <- b / a[1]
    a <- a / a[1]
  }
  n <- max(length(a), length(b))
  if (length(a) < n) a <- c(a, rep(0, n - length(a)))
  if (length(b) < n) b <- c(b, rep(0, n - length(b)))
  K <- n - 1L
  if (K < 1L) return(numeric(0))

  # IminusA = I - companion(a)^T
  # companion(a)^T has first column -a[2:n] and a super-diagonal of ones.
  CT <- matrix(0, K, K)
  CT[, 1] <- -a[2:n]
  if (K >= 2L) {
    for (i in seq_len(K - 1L)) CT[i, i + 1L] <- 1
  }
  IminusA <- diag(K) - CT
  B <- b[2:n] - a[2:n] * b[1]
  as.numeric(solve(IminusA, B))
}

#' Steady-state initial conditions for a cascaded SOS filter
#'
#' Per-section analogue of \code{\link{lfilterInit}} for second-order-sections
#' form; equivalent to \code{scipy.signal.sosfilt_zi}. Each section's steady
#' state is scaled by the DC gain of all preceding sections so that the whole
#' cascade starts in steady state when the returned matrix (scaled by the first
#' sample) is handed to \code{\link{sosfilt}}.
#'
#' @param sos SOS matrix (\code{n_sections} x 6).
#' @return A numeric matrix (\code{n_sections} x 2) of per-section initial
#'   states (unscaled).
#' @seealso \code{\link{sosfilt}}, \code{\link{lfilterInit}}
#' @export
#' @examples
#' sos <- sosDesign(high = 40, type = "low", sr = 250)
#' zi <- sosfiltInit(sos)
sosfiltInit <- function(sos) {
  .check_sos(sos)
  n_sections <- nrow(sos)
  zi <- matrix(0, n_sections, 2L)
  scale <- 1
  for (s in seq_len(n_sections)) {
    b <- sos[s, 1:3]
    a <- sos[s, 4:6]
    zi[s, ] <- scale * lfilterInit(b, a)
    dc_den <- sum(a)
    if (dc_den == 0) {
      # DC gain undefined (pole at z = 1); leave running scale unchanged.
      next
    }
    scale <- scale * sum(b) / dc_den
  }
  zi
}

# ---------------------------------------------------------------------------
# Core one-pass filters
# ---------------------------------------------------------------------------

#' One-pass causal IIR filtering (Direct Form II Transposed)
#'
#' Applies a rational transfer function defined by numerator \code{b} and
#' denominator \code{a} to \code{x}, forward in time, mirroring
#' \code{scipy.signal.lfilter}. Optionally carries the delay-line state so a
#' signal can be filtered in consecutive chunks with a result identical to
#' filtering it whole.
#'
#' @param b Numerator (feed-forward) coefficients.
#' @param a Denominator (feedback) coefficients (\code{a[1]} normalised
#'   internally if not 1).
#' @param x Numeric input vector.
#' @param zi Optional initial delay-line state, length
#'   \code{max(length(a), length(b)) - 1}. When supplied, the return value is a
#'   list carrying the final state; when \code{NULL} (default) the state starts
#'   at zero and a plain filtered vector is returned. Use
#'   \code{\link{lfilterInit}} for a transient-free warm start.
#' @return If \code{zi} is \code{NULL}, the filtered numeric vector. Otherwise a
#'   list with \code{y} (filtered vector) and \code{zi} (final state).
#' @seealso \code{\link{lfilterInit}}, \code{\link{sosfilt}}
#' @export
#' @examples
#' ba <- signal::butter(4, 0.2, "low")
#' # single call
#' y <- lfilter(ba$b, ba$a, rnorm(500))
#' # streamed in two chunks with carried state
#' x <- rnorm(500)
#' r1 <- lfilter(ba$b, ba$a, x[1:250], zi = numeric(4))
#' r2 <- lfilter(ba$b, ba$a, x[251:500], zi = r1$zi)
#' identical(length(c(r1$y, r2$y)), length(x))
lfilter <- function(b, a, x, zi = NULL) {
  b <- as.numeric(b)
  a <- as.numeric(a)
  x <- as.numeric(x)
  if (length(a) < 1L || a[1] == 0) {
    stop("`a` must have a non-zero leading coefficient.", call. = FALSE)
  }
  if (a[1] != 1) {
    b <- b / a[1]
    a <- a / a[1]
  }
  n <- max(length(a), length(b))
  if (length(a) < n) a <- c(a, rep(0, n - length(a)))
  if (length(b) < n) b <- c(b, rep(0, n - length(b)))
  K <- n - 1L

  return_state <- !is.null(zi)
  if (return_state) {
    z <- as.numeric(zi)
    if (length(z) != K) {
      stop(sprintf("`zi` must have length %d (max(length(a), length(b)) - 1).", K),
           call. = FALSE)
    }
  } else {
    z <- numeric(K)
  }

  res <- cpp_lfilter(b, a, x, z)

  if (return_state) list(y = res$y, zi = res$zi) else res$y
}

# Pure-R reference implementation of the lfilter recursion (pre-0.3.0 loop),
# retained for numerical-parity tests against the C++ kernel.
.lfilter_r <- function(b, a, x, zi = NULL) {
  b <- as.numeric(b)
  a <- as.numeric(a)
  x <- as.numeric(x)
  if (a[1] != 1) {
    b <- b / a[1]
    a <- a / a[1]
  }
  n <- max(length(a), length(b))
  if (length(a) < n) a <- c(a, rep(0, n - length(a)))
  if (length(b) < n) b <- c(b, rep(0, n - length(b)))
  K <- n - 1L
  return_state <- !is.null(zi)
  z <- if (return_state) as.numeric(zi) else numeric(K)

  nx <- length(x)
  y <- numeric(nx)

  if (K == 0L) {
    # Pure gain (order 0): y = b0 * x
    y <- b[1] * x
  } else {
    for (m in seq_len(nx)) {
      xm <- x[m]
      ym <- b[1] * xm + z[1]
      y[m] <- ym
      if (K >= 2L) {
        for (i in seq_len(K - 1L)) {
          z[i] <- b[i + 1L] * xm + z[i + 1L] - a[i + 1L] * ym
        }
      }
      z[K] <- b[K + 1L] * xm - a[K + 1L] * ym
    }
  }

  if (return_state) list(y = y, zi = z) else y
}

#' One-pass causal SOS filtering with optional carried state
#'
#' Filters \code{x} through a cascade of second-order sections, forward in time,
#' mirroring \code{scipy.signal.sosfilt}. Each section uses the DF2T biquad
#' recursion. Optionally carries the two-sample-per-section state so a signal
#' can be processed in consecutive chunks identically to processing it whole.
#'
#' @param sos SOS matrix (\code{n_sections} x 6), e.g. from
#'   \code{\link{sosDesign}}.
#' @param x Numeric input vector.
#' @param zi Optional initial state, a matrix \code{n_sections} x 2. When
#'   supplied the return value is a list carrying the final state; when
#'   \code{NULL} (default) state starts at zero and a plain vector is returned.
#'   Use \code{\link{sosfiltInit}} for a transient-free warm start.
#' @return If \code{zi} is \code{NULL}, the filtered numeric vector. Otherwise a
#'   list with \code{y} (filtered vector) and \code{zi} (final
#'   \code{n_sections} x 2 state matrix).
#' @seealso \code{\link{sosDesign}}, \code{\link{sosfiltInit}},
#'   \code{\link{sosfiltfilt}}
#' @export
#' @examples
#' sos <- sosDesign(high = 40, type = "low", sr = 250, order = 4)
#' y <- sosfilt(sos, rnorm(1000))
#' # chunked with carried state == whole-signal filtering
#' x <- rnorm(1000)
#' zi0 <- matrix(0, nrow(sos), 2)
#' r1 <- sosfilt(sos, x[1:500], zi = zi0)
#' r2 <- sosfilt(sos, x[501:1000], zi = r1$zi)
#' max(abs(c(r1$y, r2$y) - sosfilt(sos, x)))
sosfilt <- function(sos, x, zi = NULL) {
  .check_sos(sos)
  x <- as.numeric(x)
  n_sections <- nrow(sos)

  return_state <- !is.null(zi)
  if (return_state) {
    if (is.null(dim(zi)) || nrow(zi) != n_sections || ncol(zi) != 2L) {
      stop(sprintf("`zi` must be a %d x 2 matrix.", n_sections), call. = FALSE)
    }
    state <- matrix(as.numeric(zi), n_sections, 2L)
  } else {
    state <- matrix(0, n_sections, 2L)
  }

  res <- cpp_sosfilt(sos, x, state)

  if (return_state) list(y = res$y, zi = res$zi) else res$y
}

# Pure-R reference implementation of the DF2T SOS cascade (pre-0.3.0 loop),
# retained for numerical-parity tests against the C++ kernel.
.sosfilt_r <- function(sos, x, zi = NULL) {
  n_sections <- nrow(sos)
  return_state <- !is.null(zi)
  state <- if (return_state) {
    matrix(as.numeric(zi), n_sections, 2L)
  } else {
    matrix(0, n_sections, 2L)
  }

  nx <- length(x)
  y <- as.numeric(x)
  for (s in seq_len(n_sections)) {
    a0 <- sos[s, 4]
    b0 <- sos[s, 1] / a0; b1 <- sos[s, 2] / a0; b2 <- sos[s, 3] / a0
    a1 <- sos[s, 5] / a0; a2 <- sos[s, 6] / a0
    z1 <- state[s, 1]; z2 <- state[s, 2]
    out <- numeric(nx)
    for (m in seq_len(nx)) {
      xm <- y[m]
      ym <- b0 * xm + z1
      z1 <- b1 * xm - a1 * ym + z2
      z2 <- b2 * xm - a2 * ym
      out[m] <- ym
    }
    state[s, 1] <- z1; state[s, 2] <- z2
    y <- out
  }

  if (return_state) list(y = y, zi = state) else y
}

#' Zero-phase SOS filtering with state-initialised edge handling
#'
#' Forward-backward (zero-phase) SOS filtering that mirrors
#' \code{scipy.signal.sosfiltfilt}: it odd-reflects the signal at both edges,
#' initialises the delay-line to steady state (scaled by the boundary sample) to
#' suppress start-up transients, filters forward and then backward. The result
#' has zero phase distortion but is non-causal (both directions are used). For
#' strictly causal / real-time-equivalent filtering, use \code{\link{sosfilt}}
#' or \code{\link{butterworthFilter}} with \code{causal = TRUE}.
#'
#' @param sos SOS matrix (\code{n_sections} x 6).
#' @param x Numeric input vector.
#' @param padlen Number of samples of odd reflection padding at each edge.
#'   Defaults to \code{min(3 * (2 * n_sections + 1), length(x) - 1)}.
#' @return The zero-phase filtered numeric vector, same length as \code{x}.
#' @seealso \code{\link{sosfilt}}, \code{\link{sosDesign}}
#' @export
#' @examples
#' sos <- sosDesign(high = 40, type = "low", sr = 250, order = 4)
#' y <- sosfiltfilt(sos, rnorm(1000))
sosfiltfilt <- function(sos, x, padlen = NULL) {
  .check_sos(sos)
  x <- as.numeric(x)
  n <- length(x)
  n_sections <- nrow(sos)

  if (is.null(padlen)) {
    padlen <- min(3L * (2L * n_sections + 1L), n - 1L)
  }
  padlen <- max(0L, as.integer(padlen))
  if (padlen >= n) padlen <- n - 1L

  # Odd (point-symmetric) reflection at both ends.
  if (padlen > 0L) {
    x_start <- 2 * x[1] - x[(padlen + 1L):2]
    x_end <- 2 * x[n] - x[(n - 1L):(n - padlen)]
    ext <- c(x_start, x, x_end)
  } else {
    ext <- x
  }

  zi <- sosfiltInit(sos)

  # Forward pass with steady-state warm start.
  fwd <- sosfilt(sos, ext, zi = zi * ext[1])
  # Backward pass with steady-state warm start.
  yrev <- rev(fwd$y)
  bwd <- sosfilt(sos, yrev, zi = zi * yrev[1])
  y <- rev(bwd$y)

  if (padlen > 0L) {
    y <- y[(padlen + 1L):(padlen + n)]
  }
  y
}

# ---------------------------------------------------------------------------
# StreamFilter reference object (mutable, carries state across apply() calls)
# ---------------------------------------------------------------------------

#' Stateful streaming filter object
#'
#' Constructs a mutable filter object that carries its delay-line state across
#' successive \code{$apply()} calls, so a signal fed in arbitrary chunks
#' produces exactly the same output as filtering it whole. The identical math is
#' used offline and online, which lets a real-time streaming backend reuse this
#' object without re-deriving coefficients.
#'
#' The returned object is an environment (reference semantics) with methods:
#' \describe{
#'   \item{\code{$apply(x)}}{Filter numeric vector \code{x}, advancing and
#'     storing the internal state; returns the filtered vector.}
#'   \item{\code{$reset()}}{Clear the carried state so the next \code{$apply()}
#'     starts fresh (cold or warm start per \code{warmup}).}
#'   \item{\code{$state()}}{Return the current internal state matrix.}
#' }
#'
#' @param sos SOS matrix (\code{n_sections} x 6), e.g. from
#'   \code{\link{sosDesign}}. Provide either \code{sos} or both \code{b} and
#'   \code{a}.
#' @param b,a Transfer-function coefficients (converted to a single SOS section
#'   only when they describe a biquad; otherwise supply \code{sos}).
#' @param warmup Logical; if \code{TRUE} (default) the first \code{$apply()}
#'   after construction or \code{$reset()} initialises the state to steady state
#'   scaled by the first sample, suppressing the start-up transient. If
#'   \code{FALSE}, the state starts at zero.
#' @return An object of class \code{"StreamFilter"} (an environment).
#' @seealso \code{\link{sosfilt}}, \code{\link{sosDesign}}
#' @export
#' @examples
#' sos <- sosDesign(high = 30, type = "low", sr = 250, order = 4)
#' sf <- StreamFilter(sos)
#' x <- rnorm(1000)
#' y_chunked <- c(sf$apply(x[1:400]), sf$apply(x[401:1000]))
#' sf$reset()
#' y_whole <- sf$apply(x)
#' max(abs(y_chunked - y_whole))
StreamFilter <- function(sos = NULL, b = NULL, a = NULL, warmup = TRUE) {
  if (is.null(sos)) {
    if (is.null(b) || is.null(a)) {
      stop("Provide either `sos` or both `b` and `a`.", call. = FALSE)
    }
    b <- as.numeric(b); a <- as.numeric(a)
    if (length(b) > 3L || length(a) > 3L) {
      stop("b/a construction supports biquads only (length <= 3); pass `sos` for higher orders.",
           call. = FALSE)
    }
    b <- c(b, rep(0, 3L - length(b)))
    a <- c(a, rep(0, 3L - length(a)))
    sos <- matrix(c(b, a), nrow = 1L)
  }
  .check_sos(sos)

  self <- new.env(parent = emptyenv())
  self$sos <- sos
  self$warmup <- isTRUE(warmup)
  self$steady <- sosfiltInit(sos)
  self$zi <- matrix(0, nrow(sos), 2L)
  self$primed <- FALSE  # TRUE once a warm/cold start has been applied

  self$reset <- function() {
    self$zi <- matrix(0, nrow(self$sos), 2L)
    self$primed <- FALSE
    invisible(self)
  }

  self$state <- function() self$zi

  self$apply <- function(x) {
    x <- as.numeric(x)
    if (length(x) == 0L) return(numeric(0))
    if (!self$primed) {
      if (self$warmup) {
        self$zi <- self$steady * x[1]
      }
      self$primed <- TRUE
    }
    res <- sosfilt(self$sos, x, zi = self$zi)
    self$zi <- res$zi
    res$y
  }

  class(self) <- "StreamFilter"
  self
}

#' @export
print.StreamFilter <- function(x, ...) {
  cat(sprintf("<StreamFilter: %d second-order section(s), warmup=%s, primed=%s>\n",
              nrow(x$sos), x$warmup, x$primed))
  invisible(x)
}
