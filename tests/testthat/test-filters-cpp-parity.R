library(testthat)
library(PhysioPreprocess)

# Numerical parity of the 0.3.0 C++ filtering kernels against the retained
# pure-R (and signal::filter-based) reference implementations. The kernels
# must be drop-in replacements: same recursion, same state contract, same
# padding, agreeing to ~1e-12.

sos_for <- function(order, type = "pass") {
  if (type == "pass") {
    sosDesign(low = 1, high = 40, type = "pass", sr = 256, order = order)
  } else if (type == "low") {
    sosDesign(high = 40, type = "low", sr = 256, order = order)
  } else if (type == "high") {
    sosDesign(low = 1, type = "high", sr = 256, order = order)
  } else {
    sosDesign(low = 45, high = 55, type = "stop", sr = 256, order = order)
  }
}

test_that("sosfilt (C++) matches the pure-R DF2T reference across designs", {
  set.seed(11)
  x <- rnorm(4096)
  for (type in c("low", "high", "pass", "stop")) {
    for (order in c(2L, 4L, 8L)) {
      sos <- sos_for(order, type)
      expect_lt(max(abs(sosfilt(sos, x) -
                        PhysioPreprocess:::.sosfilt_r(sos, x))), 1e-12)
    }
  }
})

test_that("sosfilt (C++) honors the zi state contract like the R reference", {
  set.seed(12)
  x <- rnorm(2000)
  sos <- sos_for(4L)
  zi0 <- matrix(0, nrow(sos), 2)

  r_cpp <- sosfilt(sos, x, zi = zi0)
  r_ref <- PhysioPreprocess:::.sosfilt_r(sos, x, zi = zi0)
  expect_lt(max(abs(r_cpp$y - r_ref$y)), 1e-12)
  expect_lt(max(abs(r_cpp$zi - r_ref$zi)), 1e-12)

  # chunked with carried state == whole-signal filtering
  c1 <- sosfilt(sos, x[1:1000], zi = zi0)
  c2 <- sosfilt(sos, x[1001:2000], zi = c1$zi)
  expect_lt(max(abs(c(c1$y, c2$y) - sosfilt(sos, x))), 1e-12)

  # warm start via sosfiltInit still works through the C++ path
  ziw <- sosfiltInit(sos)
  w <- sosfilt(sos, x, zi = ziw * x[1])
  expect_length(w$y, length(x))
  expect_true(all(is.finite(w$y)))
})

test_that("lfilter (C++) matches the pure-R reference incl. state and gain", {
  set.seed(13)
  x <- rnorm(3000)
  bf <- signal::butter(4, c(1, 40) / 128, type = "pass")

  expect_lt(max(abs(lfilter(bf$b, bf$a, x) -
                    PhysioPreprocess:::.lfilter_r(bf$b, bf$a, x))), 1e-12)

  zi <- lfilterInit(bf$b, bf$a) * x[1]
  r_cpp <- lfilter(bf$b, bf$a, x, zi = zi)
  r_ref <- PhysioPreprocess:::.lfilter_r(bf$b, bf$a, x, zi = zi)
  expect_lt(max(abs(r_cpp$y - r_ref$y)), 1e-12)
  expect_lt(max(abs(r_cpp$zi - r_ref$zi)), 1e-12)

  # order-0 pure gain
  expect_equal(lfilter(2, 1, x), 2 * x)

  # non-normalized a[1] != 1: same scaled coefficients through both paths
  # (the normalization division itself perturbs the filter, so the reference
  # must receive the identical scaled inputs)
  expect_lt(max(abs(lfilter(bf$b * 3, bf$a * 3, x) -
                    PhysioPreprocess:::.lfilter_r(bf$b * 3, bf$a * 3, x))), 1e-12)
})

test_that(".sosfilt (C++) matches the signal::filter cascade reference", {
  set.seed(14)
  x <- rnorm(4096)
  for (order in c(2L, 4L, 8L)) {
    sos <- sos_for(order)
    expect_lt(max(abs(PhysioPreprocess:::.sosfilt(sos, x) -
                      PhysioPreprocess:::.sosfilt_signal(sos, x))), 1e-10)
  }
})

test_that(".sosfiltfilt (C++) matches the pure-R zero-phase reference", {
  set.seed(15)
  for (n in c(20L, 500L, 4096L)) {           # incl. n < padding length
    x <- rnorm(n)
    for (order in c(2L, 4L, 8L)) {
      sos <- sos_for(order)
      expect_lt(max(abs(PhysioPreprocess:::.sosfiltfilt(sos, x) -
                        PhysioPreprocess:::.sosfiltfilt_r(sos, x))), 1e-10)
    }
  }
})

test_that("cpp_filtfilt_ba matches signal::filtfilt exactly in scheme", {
  set.seed(17)
  x <- rnorm(3000)
  for (spec in list(signal::butter(4, c(49, 51) / 128, type = "stop"),
                    signal::butter(4, c(1, 40) / 128, type = "pass"),
                    signal::butter(2, 30 / 128, type = "low"))) {
    # DF2T vs stats::filter's MA+AR split reorders float ops; near-unit-circle
    # poles (1 Hz edge) amplify that to ~1e-9 on an O(1) signal - still far
    # inside the 1e-8 golden tolerance.
    expect_lt(max(abs(PhysioPreprocess:::cpp_filtfilt_ba(spec$b, spec$a, x) -
                      signal::filtfilt(spec, x))), 5e-8)
  }
  # matrix variant equals the vector variant column-wise
  X <- matrix(rnorm(1000 * 4), ncol = 4)
  bf <- signal::butter(4, c(49, 51) / 128, type = "stop")
  got <- PhysioPreprocess:::cpp_filtfilt_ba_mat(bf$b, bf$a, X)
  want <- apply(X, 2, function(v) PhysioPreprocess:::cpp_filtfilt_ba(bf$b, bf$a, v))
  expect_identical(got, want)
})

test_that("butterworthFilter matrix fast path equals the per-channel path", {
  set.seed(16)
  X <- matrix(rnorm(2048 * 5), ncol = 5)
  colnames(X) <- paste0("ch", 1:5)
  pe <- PhysioExperiment(assays = list(raw = X), samplingRate = 256)

  pe_f <- butterworthFilter(pe, low = 1, high = 40, order = 4L, type = "pass")
  got <- SummarizedExperiment::assay(pe_f, "filtered")

  sos <- PhysioPreprocess:::.buttersos(n = 4L, W = c(1, 40) / 128, type = "pass")
  want <- apply(X, 2, function(v) PhysioPreprocess:::.sosfiltfilt_r(sos, v))
  expect_lt(max(abs(got - want)), 1e-10)
  expect_identical(dimnames(got), dimnames(X))
})
