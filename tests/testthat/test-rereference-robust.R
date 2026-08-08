library(testthat)
library(PhysioPreprocess)

# Multichannel matrix with a shared signal plus one injected high-amplitude channel.
make_bad_channel_pe <- function(seed = 1, nch = 19, nt = 3000, bad = 7,
                                artifact_sd = 200) {
  set.seed(seed)
  base <- matrix(rnorm(nt * nch, sd = 10), nt, nch) +
    20 * sin(2 * pi * 10 * seq_len(nt) / 500)
  base[, bad] <- base[, bad] + rnorm(nt, sd = artifact_sd)
  pe <- PhysioExperiment(
    assays = list(raw = base),
    colData = S4Vectors::DataFrame(label = paste0("E", seq_len(nch)),
                                   type = rep("EEG", nch)),
    samplingRate = 500)
  list(pe = pe, bad = bad, base = base)
}

rms <- function(m) sqrt(mean(m^2))

test_that("robust reference matches the clean average and beats plain average", {
  fx <- make_bad_channel_pe()
  good <- setdiff(seq_len(ncol(fx$base)), fx$bad)
  clean_avg <- rowMeans(fx$base[, good])
  clean_ref <- fx$base[, good] - clean_avg

  rob <- SummarizedExperiment::assay(
    rereference(fx$pe, ref_type = "robust"), "rereferenced")
  expect_lt(rms(rob[, good] - clean_ref) / rms(clean_ref), 0.01)

  plain <- SummarizedExperiment::assay(
    rereference(fx$pe, ref_type = "average"), "rereferenced")
  expect_gt(rms(plain[, good] - clean_ref), rms(rob[, good] - clean_ref))
  expect_gt(rms(plain[, good] - clean_ref), 200 / ncol(fx$base) * 0.5)
})

test_that("robust reference excludes the injected channel and records provenance", {
  fx <- make_bad_channel_pe()
  pr <- rereference(fx$pe, ref_type = "robust")
  expect_equal(S4Vectors::metadata(pr)$reference, "robust")
  expect_true(paste0("E", fx$bad) %in%
                S4Vectors::metadata(pr)$reference_excluded_channels)
  prov <- PhysioCore::provenance(pr)
  expect_true(any(prov$step == "rereference"))
})

test_that("robust reference is idempotent (< 1e-6)", {
  fx <- make_bad_channel_pe()
  rob <- SummarizedExperiment::assay(
    rereference(fx$pe, ref_type = "robust"), "rereferenced")
  pe2 <- fx$pe
  SummarizedExperiment::assay(pe2, "raw") <- rob
  rob2 <- SummarizedExperiment::assay(
    rereference(pe2, ref_type = "robust"), "rereferenced")
  expect_lt(max(abs(rob2 - rob)), 1e-6)
})

test_that("median reference is robust and idempotent", {
  fx <- make_bad_channel_pe()
  med <- SummarizedExperiment::assay(
    rereference(fx$pe, ref_type = "median"), "rereferenced")
  pe2 <- fx$pe
  SummarizedExperiment::assay(pe2, "raw") <- med
  med2 <- SummarizedExperiment::assay(
    rereference(pe2, ref_type = "median"), "rereferenced")
  expect_lt(max(abs(med2 - med)), 1e-6)
})

test_that("robust reference works on epoched (3D) data", {
  fx <- make_bad_channel_pe()
  d3 <- array(fx$base, dim = c(nrow(fx$base), ncol(fx$base), 2))
  pe3 <- PhysioExperiment(
    assays = list(raw = d3),
    colData = S4Vectors::DataFrame(label = paste0("E", seq_len(ncol(fx$base))),
                                   type = rep("EEG", ncol(fx$base))),
    samplingRate = 500)
  r3 <- rereference(pe3, ref_type = "robust")
  expect_true(paste0("E", fx$bad) %in%
                S4Vectors::metadata(r3)$reference_excluded_channels)
  expect_equal(dim(SummarizedExperiment::assay(r3, "rereferenced")), dim(d3))
})

test_that("existing reference types still work", {
  fx <- make_bad_channel_pe()
  av <- rereference(fx$pe, ref_type = "average")
  expect_true("rereferenced" %in% SummarizedExperiment::assayNames(av))
})
