library(testthat)
library(PhysioPreprocess)

# Simple Welch PSD for attenuation checks.
.psd <- function(x, sr, nperseg = 2048) {
  n <- length(x); step <- nperseg %/% 2
  win <- 0.5 - 0.5 * cos(2 * pi * (seq_len(nperseg) - 1) / (nperseg - 1))
  starts <- seq(1, n - nperseg + 1, by = step)
  P <- 0
  for (s in starts) P <- P + Mod(stats::fft(x[s:(s + nperseg - 1)] * win))^2
  P <- P / length(starts)
  half <- seq_len(nperseg %/% 2 + 1)
  list(freq = (0:(nperseg - 1))[half] * sr / nperseg, psd = P[half])
}
.atten_db <- function(before, after, sr, f0) {
  pb <- .psd(before, sr); pa <- .psd(after, sr)
  bi <- which.min(abs(pb$freq - f0))
  10 * log10(pb$psd[bi] / pa$psd[bi])
}
.broadband_db <- function(before, after, sr, lines) {
  pb <- .psd(before, sr); pa <- .psd(after, sr)
  m <- vapply(pb$freq, function(f) all(abs(f - lines) > 2) & f > 1 & f < sr / 2 - 1,
              logical(1))
  max(abs(10 * log10(pb$psd[m]) - 10 * log10(pa$psd[m])))
}

make_line_pe <- function(seed = 1, sr = 500, m = 6, dur = 60,
                         amps = c(5, 2, 1)) {
  set.seed(seed)
  n <- sr * dur; t <- (0:(n - 1)) / sr
  S <- matrix(rnorm(n * 3), n, 3); Mix <- matrix(rnorm(m * 3), m, 3)
  X <- S %*% t(Mix) + matrix(rnorm(n * m, 0, 0.3), n, m)
  line <- amps[1] * sin(2 * pi * 50 * t) + amps[2] * sin(2 * pi * 100 * t) +
    amps[3] * sin(2 * pi * 150 * t)
  Xl <- X + outer(line, runif(m, 0.5, 1.5))
  pe <- function(mat) PhysioExperiment(
    assays = list(raw = mat),
    colData = S4Vectors::DataFrame(label = paste0("E", seq_len(m)),
                                   type = rep("EEG", m)),
    samplingRate = sr)
  list(X = X, Xl = Xl, pe_line = pe(Xl), pe_clean = pe(X), sr = sr, m = m)
}

test_that("cleanLine attenuates the 50 Hz line > 30 dB and preserves broadband < 1 dB", {
  d <- make_line_pe()
  res <- cleanLine(d$pe_line, line_freq = 50, harmonics = 3)
  cl <- SummarizedExperiment::assay(res, "cleanline")
  expect_equal(dim(cl), dim(d$Xl))
  expect_gt(.atten_db(d$Xl[, 1], cl[, 1], d$sr, 50), 30)
  expect_lt(.broadband_db(d$Xl[, 1], cl[, 1], d$sr, c(50, 100, 150)), 1)
})

test_that("cleanLine harmonics option removes 100 and 150 Hz when requested", {
  d <- make_line_pe()
  h3 <- SummarizedExperiment::assay(cleanLine(d$pe_line, 50, harmonics = 3), "cleanline")
  expect_gt(.atten_db(d$Xl[, 1], h3[, 1], d$sr, 100), 15)
  expect_gt(.atten_db(d$Xl[, 1], h3[, 1], d$sr, 150), 15)

  # with harmonics = 1 the 100 Hz component is left essentially untouched
  h1 <- SummarizedExperiment::assay(cleanLine(d$pe_line, 50, harmonics = 1), "cleanline")
  expect_gt(.atten_db(d$Xl[, 1], h1[, 1], d$sr, 50), 30)
  expect_lt(.atten_db(d$Xl[, 1], h1[, 1], d$sr, 100), 3)
})

test_that("cleanLine records removed line-power in metadata and provenance", {
  d <- make_line_pe()
  res <- cleanLine(d$pe_line, line_freq = 50, harmonics = 3)
  md <- S4Vectors::metadata(res)$line_noise
  expect_equal(md$method, "cleanLine")
  expect_gt(md$removed_power, 0)
  expect_true(any(PhysioCore::provenance(res)$step == "cleanLine"))
})

test_that("zapLine attenuates the 50 Hz line > 30 dB and preserves broadband < 1 dB", {
  d <- make_line_pe()
  res <- zapLine(d$pe_line, line_freq = 50, nremove = 1)
  z <- SummarizedExperiment::assay(res, "zapline")
  expect_equal(dim(z), dim(d$Xl))
  expect_gt(.atten_db(d$Xl[, 1], z[, 1], d$sr, 50), 30)
  expect_lt(.broadband_db(d$Xl[, 1], z[, 1], d$sr, c(50, 100, 150)), 1)
  expect_true(any(PhysioCore::provenance(res)$step == "zapLine"))
})

test_that("line-noise removers require a 2D assay and validate input", {
  d <- make_line_pe()
  expect_error(cleanLine(list()), "PhysioExperiment")
  expect_error(zapLine(list()), "PhysioExperiment")
  arr <- array(rnorm(100 * 4 * 2), dim = c(100, 4, 2))
  pe3 <- PhysioExperiment(
    assays = list(raw = arr),
    colData = S4Vectors::DataFrame(label = paste0("E", 1:4)),
    samplingRate = 100)
  expect_error(zapLine(pe3), "2D")
  expect_error(cleanLine(pe3), "2D")
})
