library(testthat)
library(PhysioPreprocess)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_chirp <- function(N = 2000, sr = 250, f0 = 1, f1 = 60) {
  t <- seq(0, (N - 1) / sr, length.out = N)
  dur <- (N - 1) / sr
  sin(2 * pi * (f0 * t + (f1 - f0) / (2 * dur) * t^2))
}

make_pe_1ch <- function(x, sr = 250) {
  n <- length(x)
  PhysioExperiment(
    assays = S4Vectors::SimpleList(raw = matrix(x, ncol = 1)),
    rowData = S4Vectors::DataFrame(time_idx = seq_len(n)),
    colData = S4Vectors::DataFrame(label = "Ch1"),
    samplingRate = sr
  )
}

# ---------------------------------------------------------------------------
# AC1 -- causal sosfilt equals signal::filter() within 1e-9 on a chirp
# ---------------------------------------------------------------------------

test_that("sosfilt matches signal::filter on a chirp (< 1e-9)", {
  sr <- 250
  chirp <- make_chirp(sr = sr)

  # Cover every band type against the signal::butter reference so a design
  # error (e.g. mis-placed band-stop zeros) cannot ship green. Orders are kept
  # low (2, 4) where the SOS cascade and the reference ba-form realisation
  # agree to machine precision; higher orders diverge only because the ba
  # reference itself becomes ill-conditioned.
  cases <- list(
    list(type = "low",  args = list(high = 30)),
    list(type = "high", args = list(low = 5)),
    list(type = "pass", args = list(low = 8, high = 30)),
    list(type = "stop", args = list(low = 45, high = 55))
  )
  for (cs in cases) {
    for (ord in c(2L, 4L)) {
      W <- switch(cs$type,
                  low  = cs$args$high / (sr / 2),
                  high = cs$args$low / (sr / 2),
                  c(cs$args$low, cs$args$high) / (sr / 2))
      bf <- signal::butter(ord, W, cs$type)
      ref <- as.numeric(signal::filter(bf, chirp))
      sos <- do.call(sosDesign, c(cs$args, list(type = cs$type, sr = sr, order = ord)))
      got <- sosfilt(sos, chirp)
      expect_equal(length(got), length(chirp))
      expect_lt(max(abs(got - ref)), 1e-9)
    }
  }
})

test_that("lfilter matches signal::filter for the same (b, a)", {
  sr <- 250
  chirp <- make_chirp(sr = sr)
  bf <- signal::butter(4, 30 / (sr / 2), "low")
  ref <- as.numeric(signal::filter(bf, chirp))
  got <- lfilter(bf$b, bf$a, chirp)
  expect_lt(max(abs(got - ref)), 1e-9)
})

test_that("SOS-form band designs match the signal::butter reference (regression)", {
  # Regression guard for the .buttersos band-pass gain and band-stop zero
  # placement: both centres must be the geometric/pre-warped mean, not the
  # arithmetic mean of the digital band edges. Compare the whole magnitude
  # response of the SOS design to signal::butter over a frequency grid.
  sr <- 250
  freqs <- seq(1, 124, by = 1)
  z <- exp(1i * 2 * pi * freqs / sr)
  sos_mag <- function(sos, z) {
    H <- rep(1 + 0i, length(z))
    for (r in seq_len(nrow(sos))) {
      H <- H * (sos[r, 1] * z^-0 + sos[r, 2] * z^-1 + sos[r, 3] * z^-2) /
               (sos[r, 4] * z^-0 + sos[r, 5] * z^-1 + sos[r, 6] * z^-2)
    }
    abs(H)
  }
  ba_mag <- function(bf, z) {
    nb <- length(bf$b); na <- length(bf$a)
    num <- Reduce(`+`, lapply(seq_len(nb), function(k) bf$b[k] * z^-(k - 1)))
    den <- Reduce(`+`, lapply(seq_len(na), function(k) bf$a[k] * z^-(k - 1)))
    abs(num / den)
  }
  for (spec in list(
    list(type = "pass", low = 8,  high = 30),
    list(type = "stop", low = 45, high = 55)
  )) {
    W <- c(spec$low, spec$high) / (sr / 2)
    sos <- sosDesign(low = spec$low, high = spec$high, type = spec$type,
                     sr = sr, order = 4)
    bf <- signal::butter(4, W, spec$type)
    expect_lt(max(abs(sos_mag(sos, z) - ba_mag(bf, z))), 1e-3)
  }

  # Band-stop pass-band must be flat (~unit gain) well away from the notch.
  sos_stop <- sosDesign(low = 45, high = 55, type = "stop", sr = sr, order = 4)
  passband_hz <- c(5, 20, 35, 70, 90, 110)
  zp <- exp(1i * 2 * pi * passband_hz / sr)
  expect_lt(max(abs(sos_mag(sos_stop, zp) - 1)), 1e-2)
})

test_that("sosfilt matches lfilter of the convolved (b, a)", {
  sr <- 250
  chirp <- make_chirp(sr = sr, N = 800)
  sos <- sosDesign(low = 8, high = 30, type = "pass", sr = sr, order = 4)
  # Convolve section polynomials into overall (b, a).
  b <- 1; a <- 1
  for (s in seq_len(nrow(sos))) {
    b <- convolve(b, rev(sos[s, 1:3] / sos[s, 4]), type = "open")
    a <- convolve(a, rev(sos[s, 4:6] / sos[s, 4]), type = "open")
  }
  got_sos <- sosfilt(sos, chirp)
  got_lf <- lfilter(b, a, chirp)
  # Cascade vs direct form differ only by rounding for this order.
  expect_lt(max(abs(got_sos - got_lf)), 1e-6)
})

# ---------------------------------------------------------------------------
# AC2 -- state continuity: 2-chunk carried zi == whole-signal filtering
# ---------------------------------------------------------------------------

test_that("sosfilt state continuity: 2 chunks equal the whole (< 1e-10)", {
  sr <- 250
  chirp <- make_chirp(sr = sr)
  sos <- sosDesign(high = 30, type = "low", sr = sr, order = 4)
  N <- length(chirp)
  split <- 900L

  zi0 <- matrix(0, nrow(sos), 2)
  whole <- sosfilt(sos, chirp, zi = zi0)$y
  r1 <- sosfilt(sos, chirp[1:split], zi = zi0)
  r2 <- sosfilt(sos, chirp[(split + 1):N], zi = r1$zi)
  chunked <- c(r1$y, r2$y)

  expect_equal(max(abs(whole - chunked)), 0, tolerance = 1e-10)
  # zero-state list form equals the plain-vector form
  expect_equal(whole, sosfilt(sos, chirp), tolerance = 1e-12)
})

test_that("lfilter state continuity: 2 chunks equal the whole (< 1e-10)", {
  sr <- 250
  chirp <- make_chirp(sr = sr)
  bf <- signal::butter(4, 30 / (sr / 2), "low")
  K <- max(length(bf$b), length(bf$a)) - 1L
  N <- length(chirp)
  split <- 900L

  whole <- lfilter(bf$b, bf$a, chirp, zi = numeric(K))$y
  r1 <- lfilter(bf$b, bf$a, chirp[1:split], zi = numeric(K))
  r2 <- lfilter(bf$b, bf$a, chirp[(split + 1):N], zi = r1$zi)
  expect_lt(max(abs(whole - c(r1$y, r2$y))), 1e-10)
})

test_that("multi-chunk continuity holds for many small chunks", {
  sr <- 250
  chirp <- make_chirp(sr = sr, N = 1500)
  sos <- sosDesign(low = 5, high = 40, type = "pass", sr = sr, order = 6)
  whole <- sosfilt(sos, chirp)

  bounds <- c(0, 137, 400, 401, 900, 1200, 1500)
  zi <- matrix(0, nrow(sos), 2)
  pieces <- list()
  for (k in seq_len(length(bounds) - 1L)) {
    idx <- (bounds[k] + 1L):bounds[k + 1L]
    r <- sosfilt(sos, chirp[idx], zi = zi)
    pieces[[k]] <- r$y
    zi <- r$zi
  }
  expect_lt(max(abs(whole - unlist(pieces))), 1e-10)
})

# ---------------------------------------------------------------------------
# Steady-state warm start (lfilterInit / sosfiltInit)
# ---------------------------------------------------------------------------

test_that("sosfiltInit warm start removes the step start-up transient", {
  sr <- 250
  sos <- sosDesign(high = 30, type = "low", sr = sr, order = 4)
  step <- rep(1, 500)

  y_warm <- sosfilt(sos, step, zi = sosfiltInit(sos) * step[1])$y
  y_cold <- sosfilt(sos, step)

  # Unit DC gain lowpass: warm start sits at the steady value immediately.
  expect_lt(max(abs(y_warm - 1)), 1e-8)
  # Cold start must ramp up (transient present).
  expect_gt(abs(y_cold[1] - 1), 0.5)
})

test_that("lfilterInit has the correct length and warm-starts a step", {
  bf <- signal::butter(4, 0.2, "low")
  zi <- lfilterInit(bf$b, bf$a)
  expect_length(zi, max(length(bf$b), length(bf$a)) - 1L)

  step <- rep(2, 300)
  y <- lfilter(bf$b, bf$a, step, zi = zi * step[1])$y
  expect_lt(max(abs(y - 2)), 1e-8)
})

# ---------------------------------------------------------------------------
# sosfiltfilt -- zero phase
# ---------------------------------------------------------------------------

test_that("sosfiltfilt is zero-phase and length-preserving", {
  sr <- 250
  N <- 1000
  t <- seq(0, (N - 1) / sr, length.out = N)
  x <- sin(2 * pi * 8 * t) + 0.4 * sin(2 * pi * 55 * t)
  sos <- sosDesign(high = 20, type = "low", sr = sr, order = 4)

  y <- sosfiltfilt(sos, x)
  expect_length(y, N)
  # Zero phase: the passband 8 Hz sinusoid stays aligned (no phase lag).
  # Cross-correlation peak should be at lag 0.
  seg <- 300:700
  cc <- sapply(-5:5, function(lag) sum(x[seg] * y[seg + lag]))
  expect_equal(which.max(cc) - 6L, 0L)
})

# ---------------------------------------------------------------------------
# AC3 -- butterworthFilter(causal = TRUE)
# ---------------------------------------------------------------------------

test_that("butterworthFilter(causal=TRUE) preserves dims and has no pre-ringing", {
  sr <- 250
  n <- 400
  xstep <- c(rep(0, 200), rep(1, 200))
  pe <- make_pe_1ch(xstep, sr = sr)

  pe_c <- butterworthFilter(pe, high = 10, type = "low", causal = TRUE)
  pe_z <- butterworthFilter(pe, high = 10, type = "low", causal = FALSE)

  expect_s4_class(pe_c, "PhysioExperiment")
  expect_true("filtered" %in% SummarizedExperiment::assayNames(pe_c))
  expect_identical(dim(SummarizedExperiment::assay(pe_c, "filtered")),
                   dim(SummarizedExperiment::assay(pe, "raw")))

  yc <- as.numeric(SummarizedExperiment::assay(pe_c, "filtered"))
  yz <- as.numeric(SummarizedExperiment::assay(pe_z, "filtered"))

  # Causal filter: strictly no response before the transition (no pre-ringing).
  expect_lt(max(abs(yc[1:195])), 1e-9)
  # Zero-phase filtfilt smears the edge backward in time (acausal pre-ringing).
  expect_gt(max(abs(yz[180:200])), 1e-2)
})

test_that("butterworthFilter(causal=TRUE) is dimension-preserving on 3D assays", {
  sr <- 250
  n <- 500
  x <- make_chirp(N = n, sr = sr)
  arr <- array(c(x, x * 0.8), dim = c(n, 2, 1))
  pe <- PhysioExperiment(
    assays = S4Vectors::SimpleList(raw = arr),
    rowData = S4Vectors::DataFrame(time_idx = seq_len(n)),
    colData = S4Vectors::DataFrame(label = c("Ch1", "Ch2")),
    samplingRate = sr
  )
  out <- butterworthFilter(pe, low = 5, high = 40, type = "pass", causal = TRUE)
  fd <- SummarizedExperiment::assay(out, "filtered")
  expect_identical(dim(fd), dim(arr))
  expect_true(all(is.finite(fd)))

  # Per-channel result equals the direct sosfilt warm-start on that channel.
  sos <- sosDesign(low = 5, high = 40, type = "pass", sr = sr, order = 4)
  ref_ch1 <- sosfilt(sos, arr[, 1, 1], zi = sosfiltInit(sos) * arr[1, 1, 1])$y
  expect_lt(max(abs(fd[, 1, 1] - ref_ch1)), 1e-10)
})

test_that("firFilter(causal=TRUE) preserves dims and is finite", {
  sr <- 250
  x <- make_chirp(N = 600, sr = sr)
  pe <- make_pe_1ch(x, sr = sr)
  out <- firFilter(pe, high = 20, type = "low", order = 50, causal = TRUE)
  fd <- SummarizedExperiment::assay(out, "filtered")
  expect_identical(dim(fd), dim(SummarizedExperiment::assay(pe, "raw")))
  expect_true(all(is.finite(fd)))
})

# ---------------------------------------------------------------------------
# StreamFilter
# ---------------------------------------------------------------------------

test_that("StreamFilter chunked apply equals whole-signal filtering", {
  sr <- 250
  chirp <- make_chirp(sr = sr)
  sos <- sosDesign(high = 30, type = "low", sr = sr, order = 4)

  sf <- StreamFilter(sos)
  yc <- c(sf$apply(chirp[1:400]),
          sf$apply(chirp[401:1200]),
          sf$apply(chirp[1201:length(chirp)]))
  sf$reset()
  yw <- sf$apply(chirp)
  expect_equal(max(abs(yc - yw)), 0, tolerance = 1e-12)
})

test_that("StreamFilter warmup removes the start-up transient on a step", {
  sr <- 250
  sos <- sosDesign(high = 30, type = "low", sr = sr, order = 4)
  step <- rep(3, 400)

  sf_warm <- StreamFilter(sos, warmup = TRUE)
  y_warm <- sf_warm$apply(step)
  expect_lt(max(abs(y_warm - 3)), 1e-8)

  sf_cold <- StreamFilter(sos, warmup = FALSE)
  y_cold <- sf_cold$apply(step)
  expect_gt(abs(y_cold[1] - 3), 1)
})

test_that("StreamFilter can be built from biquad b/a and prints", {
  bf <- signal::butter(2, 0.2, "low")
  sf <- StreamFilter(b = bf$b, a = bf$a)
  expect_s3_class(sf, "StreamFilter")
  expect_output(print(sf), "StreamFilter")

  x <- make_chirp(N = 500)
  y <- sf$apply(x)
  ref <- sosfilt(matrix(c(bf$b, bf$a), nrow = 1), x,
                 zi = sosfiltInit(matrix(c(bf$b, bf$a), nrow = 1)) * x[1])$y
  expect_lt(max(abs(y - ref)), 1e-10)
})

test_that("StreamFilter reset clears carried state", {
  sos <- sosDesign(high = 30, type = "low", sr = 250, order = 4)
  sf <- StreamFilter(sos)
  x <- make_chirp(N = 300)
  sf$apply(x)
  expect_false(all(sf$state() == 0))
  sf$reset()
  expect_true(all(sf$state() == 0))
})

# ---------------------------------------------------------------------------
# Input validation / edge cases
# ---------------------------------------------------------------------------

test_that("sosDesign validates arguments", {
  expect_error(sosDesign(high = 30, type = "low", sr = -1), "positive")
  expect_error(sosDesign(type = "low", sr = 250), "high")
  expect_error(sosDesign(low = 5, type = "high", sr = 250, order = 0), "positive integer")
  expect_error(sosDesign(high = 200, type = "low", sr = 250), "Nyquist|between")
  expect_error(sosDesign(low = 40, high = 10, type = "pass", sr = 250), "below")
})

test_that("sosfilt validates sos shape and zi shape", {
  expect_error(sosfilt(matrix(0, 2, 5), rnorm(10)), "6 columns")
  sos <- sosDesign(high = 30, type = "low", sr = 250, order = 4)
  expect_error(sosfilt(sos, rnorm(10), zi = matrix(0, 1, 2)), "matrix")
})

test_that("lfilter validates a leading coefficient", {
  expect_error(lfilter(c(1, 1), c(0, 1), rnorm(5)), "non-zero leading")
})

test_that("empty input returns empty output", {
  sos <- sosDesign(high = 30, type = "low", sr = 250, order = 4)
  expect_length(sosfilt(sos, numeric(0)), 0)
  sf <- StreamFilter(sos)
  expect_length(sf$apply(numeric(0)), 0)
})

test_that("non-unit a0 in a section is normalised correctly", {
  sr <- 250
  x <- make_chirp(N = 500, sr = sr)
  sos <- sosDesign(high = 30, type = "low", sr = sr, order = 4)
  ref <- sosfilt(sos, x)
  # Scale one section's b and a by the same factor -> identical filter.
  sos2 <- sos
  sos2[1, ] <- sos2[1, ] * 3
  got <- sosfilt(sos2, x)
  expect_lt(max(abs(got - ref)), 1e-10)
})
