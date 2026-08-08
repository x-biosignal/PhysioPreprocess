library(testthat)
library(PhysioPreprocess)

make_example <- function(n = 1000, sr = 500) {
  # Create test signal with known frequency components
  t <- seq(0, (n - 1) / sr, length.out = n)
  # 10 Hz signal + 50 Hz noise
  signal <- sin(2 * pi * 10 * t) + 0.5 * sin(2 * pi * 50 * t)
  # Create 2-channel data
  signal_data <- cbind(signal, signal * 0.9)
  assays <- S4Vectors::SimpleList(raw = array(signal_data, dim = c(n, 2, 1)))
  # rowData must have n rows (matching dim[1] = time points)
  rowData <- S4Vectors::DataFrame(time_idx = seq_len(n))
  # colData must have 2 rows (matching dim[2] = channels)
  colData <- S4Vectors::DataFrame(
    label = c("Ch1", "Ch2"),
    sensor_type = rep("EEG", 2)
  )
  PhysioExperiment(assays, rowData, colData, samplingRate = sr)
}

test_that("filterSignals works with different window sizes", {
  x <- make_example()

  y <- filterSignals(x, window = 5)
  expect_true("filtered" %in% SummarizedExperiment::assayNames(y))

  y <- filterSignals(x, window = 10, output_assay = "smoothed")
  expect_true("smoothed" %in% SummarizedExperiment::assayNames(y))
})

test_that("filterSignals handles edge cases", {
  x <- make_example(n = 3)

  # Window larger than data - should warn
  expect_warning(filterSignals(x, window = 10))

  # Invalid window
  expect_error(filterSignals(x, window = 0))
})

test_that("butterworthFilter applies bandpass correctly", {
  x <- make_example()

  # Bandpass 5-15 Hz should keep 10 Hz, remove 50 Hz
  y <- butterworthFilter(x, low = 5, high = 15, type = "pass")
  expect_true("filtered" %in% SummarizedExperiment::assayNames(y))

  filtered_data <- SummarizedExperiment::assay(y, "filtered")
  expect_true(all(is.finite(filtered_data)))
})

test_that("butterworthFilter validates parameters", {
  x <- make_example()

  # Missing frequency for lowpass
  expect_error(butterworthFilter(x, type = "low"))

  # Missing frequency for highpass
  expect_error(butterworthFilter(x, type = "high"))

  # Invalid frequency (above Nyquist)
  expect_error(butterworthFilter(x, high = 300, type = "low"))
})

test_that("firFilter works", {
  x <- make_example()

  y <- firFilter(x, low = 5, high = 30, order = 50, type = "pass")
  expect_true("filtered" %in% SummarizedExperiment::assayNames(y))
})

test_that("notchFilter removes line noise", {
  x <- make_example()

  # Remove 50 Hz
  y <- notchFilter(x, freq = 50, bandwidth = 2)
  expect_true("filtered" %in% SummarizedExperiment::assayNames(y))

  # Multiple harmonics
  y <- notchFilter(x, freq = 50, harmonics = 2)
  expect_true("filtered" %in% SummarizedExperiment::assayNames(y))
})

test_that("detrendSignal removes linear trends", {
  # Create signal with linear trend
  n <- 100
  t <- seq_len(n)
  signal <- 2 * t + rnorm(n)
  assays <- S4Vectors::SimpleList(raw = array(signal, dim = c(n, 1, 1)))
  rowData <- S4Vectors::DataFrame(time_idx = seq_len(n))
  colData <- S4Vectors::DataFrame(label = "Ch1")
  x <- PhysioExperiment(assays, rowData, colData, samplingRate = 100)

  y <- detrendSignal(x, type = "linear")
  expect_true("detrended" %in% SummarizedExperiment::assayNames(y))

  detrended <- SummarizedExperiment::assay(y, "detrended")
  # Mean should be near zero after detrending
  expect_true(abs(mean(detrended)) < 1)
})

test_that("detrendSignal removes constant (mean)", {
  n <- 100
  signal <- rnorm(n) + 100
  assays <- S4Vectors::SimpleList(raw = array(signal, dim = c(n, 1, 1)))
  rowData <- S4Vectors::DataFrame(time_idx = seq_len(n))
  colData <- S4Vectors::DataFrame(label = "Ch1")
  x <- PhysioExperiment(assays, rowData, colData, samplingRate = 100)

  y <- detrendSignal(x, type = "constant")
  detrended <- SummarizedExperiment::assay(y, "detrended")
  expect_true(abs(mean(detrended)) < 1e-10)
})
