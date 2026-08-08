library(testthat)
library(PhysioPreprocess)

# Correlated multichannel EEG-like data (shared sources + small channel noise)
# with optional injected 10x-amplitude transient bursts.
make_asr_data <- function(seed = 1, sr = 100, m = 8, dur = 60, n_sources = 4,
                          burst = TRUE) {
  set.seed(seed)
  n <- sr * dur
  S <- matrix(rnorm(n * n_sources), n, n_sources)
  Mix <- matrix(rnorm(m * n_sources), m, n_sources)
  X <- S %*% t(Mix) + matrix(rnorm(n * m, 0, 0.1), n, m)
  crms <- sqrt(mean(X^2))
  Xb <- X; burst_idx <- integer(0)
  if (burst) {
    for (b in seq(sr * 35, sr * 55, by = sr * 5)) {
      ii <- b:(b + sr - 1); burst_idx <- c(burst_idx, ii)
      ch <- sample(m, 1)
      Xb[ii, ch] <- Xb[ii, ch] + 10 * crms * sin(2 * pi * 15 * seq_along(ii) / sr)
    }
  }
  pe <- function(mat) PhysioExperiment(
    assays = list(raw = mat),
    colData = S4Vectors::DataFrame(label = paste0("E", seq_len(m)),
                                   type = rep("EEG", m)),
    samplingRate = sr)
  list(X = X, Xb = Xb, pe_clean = pe(X), pe_burst = pe(Xb),
       burst_idx = burst_idx, sr = sr, m = m, pe = pe)
}

rms <- function(a) sqrt(mean(a^2))

test_that("asrCalibrate returns an SPD covariance and an asr_calibration object", {
  d <- make_asr_data()
  cal <- asrCalibrate(d$pe(d$X[seq_len(30 * d$sr), ]), cutoff = 20)
  expect_s3_class(cal, "asr_calibration")
  ev <- eigen(cal$M, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > 0))                              # SPD
  expect_equal(dim(cal$M), c(d$m, d$m))
})

test_that("ASR reduces injected artifacts > 80% and changes clean data < 5%", {
  d <- make_asr_data()
  cal <- asrCalibrate(d$pe(d$X[seq_len(30 * d$sr), ]), cutoff = 20)
  res <- asrProcess(d$pe_burst, cal)
  cl <- SummarizedExperiment::assay(res, "asr")

  # channel count / dimensions preserved
  expect_equal(dim(cl), dim(d$Xb))

  # injected artifact (Xb - X) reduced > 80% in the burst regions
  before <- rms(d$Xb[d$burst_idx, ] - d$X[d$burst_idx, ])
  after <- rms(cl[d$burst_idx, ] - d$X[d$burst_idx, ])
  expect_gt(1 - after / before, 0.8)

  # clean first half changes < 5%
  ci <- seq_len(30 * d$sr)
  change <- rms(cl[ci, ] - d$Xb[ci, ]) / rms(d$Xb[ci, ])
  expect_lt(change, 0.05)
})

test_that("ASR with cutoff = Inf is (near) identity on clean data (< 1%)", {
  d <- make_asr_data(burst = FALSE)
  cal <- asrCalibrate(d$pe(d$X[seq_len(30 * d$sr), ]), cutoff = Inf)
  res <- asrProcess(d$pe_clean, cal)
  cl <- SummarizedExperiment::assay(res, "asr")
  expect_lt(rms(cl - d$X) / rms(d$X), 0.01)
})

test_that("asrProcess records removed variance in metadata and provenance", {
  d <- make_asr_data()
  cal <- asrCalibrate(d$pe(d$X[seq_len(30 * d$sr), ]), cutoff = 20)
  res <- asrProcess(d$pe_burst, cal)
  md <- S4Vectors::metadata(res)$asr
  expect_true(is.numeric(md$removed_var))
  expect_true(md$removed_var >= 0 && md$removed_var <= 1)
  prov <- PhysioCore::provenance(res)
  expect_true(any(prov$step == "asrProcess"))
})

test_that("asrProcess errors on a channel-count mismatch", {
  d <- make_asr_data()
  cal <- asrCalibrate(d$pe(d$X[seq_len(30 * d$sr), ]), cutoff = 20)
  wrong <- PhysioExperiment(                             # 4 channels vs 8
    assays = list(raw = d$X[, 1:4]),
    colData = S4Vectors::DataFrame(label = paste0("E", 1:4),
                                   type = rep("EEG", 4)),
    samplingRate = d$sr)
  expect_error(asrProcess(wrong, cal), "Channel count")
})

test_that("cleanRawdata orchestrates the pipeline and preserves dimensions", {
  d <- make_asr_data()
  cr <- cleanRawdata(d$pe_burst, asr_cutoff = 20)
  clean <- SummarizedExperiment::assay(cr, "clean")
  expect_equal(dim(clean), dim(d$Xb))

  before <- rms(d$Xb[d$burst_idx, ] - d$X[d$burst_idx, ])
  after <- rms(clean[d$burst_idx, ] - d$X[d$burst_idx, ])
  expect_gt(1 - after / before, 0.8)

  md <- S4Vectors::metadata(cr)$clean_rawdata
  expect_true(is.numeric(md$removed_var))
  expect_true(any(PhysioCore::provenance(cr)$step == "cleanRawdata"))
})

test_that("cleanRawdata with asr_cutoff = NULL skips ASR", {
  d <- make_asr_data(burst = FALSE)
  cr <- cleanRawdata(d$pe_clean, asr_cutoff = NULL, channel_crit = NULL)
  expect_true("clean" %in% SummarizedExperiment::assayNames(cr))
  expect_equal(S4Vectors::metadata(cr)$clean_rawdata$removed_var, 0)
})
