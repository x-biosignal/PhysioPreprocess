library(testthat)
library(PhysioPreprocess)

# Golden regression tests for PhysioPreprocess (WSF-08).
#
# Each test rebuilds the SAME deterministic input the generator used
# (data-raw/golden.R), computes the value via the PACKAGE function, and
# compares against the stored golden with expect_equal_golden(). The goldens
# were captured from INDEPENDENT references (signal::butter+filtfilt for the
# filters, base-R rowMeans subtraction for average reference); the ICA golden
# is a fixed-seed characterization (regression-guard only).
#
# If a golden .rds is absent, expect_equal_golden() skips (helper-golden.R).

# --- Shared deterministic filter input (mirrors data-raw/golden.R) ----------
make_filter_input <- function() {
  set.seed(42)
  n <- 1024L
  sr <- 250
  tt <- seq(0, (n - 1) / sr, length.out = n)
  sig <- sin(2 * pi * 10 * tt) + 0.5 * sin(2 * pi * 50 * tt) + 0.2 * rnorm(n)
  mat <- cbind(sig, 0.8 * sig + 0.1 * rnorm(n))
  list(mat = mat, n = n, sr = sr)
}

make_filter_pe <- function() {
  inp <- make_filter_input()
  PhysioExperiment(
    assays = S4Vectors::SimpleList(raw = array(inp$mat, dim = c(inp$n, 2L))),
    rowData = S4Vectors::DataFrame(t = seq_len(inp$n)),
    colData = S4Vectors::DataFrame(label = c("Ch1", "Ch2")),
    samplingRate = inp$sr
  )
}

# Package assays are unnamed matrices; strip any dimnames for a content-only
# comparison against the goldens (which were saved without dimnames).
undim <- function(m) {
  m <- as.matrix(m)
  dimnames(m) <- NULL
  m
}

test_that("butterworthFilter bandpass (ba-form) matches signal::butter+filtfilt", {
  pe <- make_filter_pe()
  out <- butterworthFilter(pe, low = 1, high = 40, order = 4,
                           type = "pass", use_sos = FALSE)
  actual <- undim(SummarizedExperiment::assay(out, "filtered"))
  expect_equal_golden(actual, "butter_bandpass_ba", tol = 1e-8)
})

test_that("butterworthFilter lowpass (ba-form) matches signal::butter+filtfilt", {
  pe <- make_filter_pe()
  out <- butterworthFilter(pe, high = 30, order = 4,
                           type = "low", use_sos = FALSE)
  actual <- undim(SummarizedExperiment::assay(out, "filtered"))
  expect_equal_golden(actual, "butter_lowpass_ba", tol = 1e-8)
})

test_that("butterworthFilter highpass (ba-form) matches signal::butter+filtfilt", {
  pe <- make_filter_pe()
  out <- butterworthFilter(pe, low = 0.5, order = 4,
                           type = "high", use_sos = FALSE)
  actual <- undim(SummarizedExperiment::assay(out, "filtered"))
  expect_equal_golden(actual, "butter_highpass_ba", tol = 1e-8)
})

test_that("notchFilter 50 Hz matches signal::butter(stop)+filtfilt", {
  pe <- make_filter_pe()
  out <- notchFilter(pe, freq = 50, bandwidth = 2, harmonics = 1L)
  actual <- undim(SummarizedExperiment::assay(out, "filtered"))
  expect_equal_golden(actual, "notch_50", tol = 1e-8)
})

test_that("rereference average matches base-R rowMeans subtraction", {
  set.seed(7)
  n <- 200L
  nch <- 8L
  mat <- matrix(rnorm(n * nch), n, nch)
  pe <- PhysioExperiment(
    assays = S4Vectors::SimpleList(raw = array(mat, dim = c(n, nch))),
    rowData = S4Vectors::DataFrame(t = seq_len(n)),
    colData = S4Vectors::DataFrame(label = paste0("Ch", seq_len(nch))),
    samplingRate = 250
  )
  out <- rereference(pe, ref_type = "average")
  actual <- undim(SummarizedExperiment::assay(out, "rereferenced"))
  expect_equal_golden(actual, "rereference_average", tol = 1e-8)
})

test_that("icaDecompose reproduces its fixed-seed characterization golden", {
  set.seed(123)
  n <- 200L
  s1 <- sin(2 * pi * 3 * seq_len(n) / 250)
  s2 <- sign(sin(2 * pi * 7 * seq_len(n) / 250))
  s3 <- rnorm(n)
  S <- cbind(s1, s2, s3)
  A_true <- matrix(c(1, 0.5, 0.3,
                     0.4, 1, 0.6,
                     0.2, 0.7, 1), nrow = 3, byrow = TRUE)
  X <- S %*% t(A_true)
  pe <- PhysioExperiment(
    assays = S4Vectors::SimpleList(raw = array(X, dim = c(n, 3L))),
    rowData = S4Vectors::DataFrame(t = seq_len(n)),
    colData = S4Vectors::DataFrame(label = paste0("Ch", 1:3)),
    samplingRate = 250
  )
  # Fix the seed of the FastICA random init immediately before decomposition,
  # exactly as the generator did (regression-guard, not an accuracy check).
  set.seed(99)
  actual <- undim(icaDecompose(pe, n_components = 3L)$components)
  expect_equal_golden(actual, "ica_characterization", tol = 1e-6)
})
