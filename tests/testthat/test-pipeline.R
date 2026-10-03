library(testthat)
library(PhysioPreprocess)

test_that("applyPipeline runs a multi-step pipeline end to end", {
  set.seed(3)
  m <- matrix(rnorm(500 * 4), 500, 4)
  pe <- PhysioExperiment(
    assays = S4Vectors::SimpleList(raw = m),
    colData = S4Vectors::DataFrame(name = paste0("C", 1:4), type = "eeg"),
    samplingRate = 250)

  pipeline <- createPipeline(
    detrend = list(fn = "detrendSignals",    method = "linear"),
    filter  = list(fn = "butterworthFilter", low = 1, high = 40, type = "pass"),
    notch   = list(fn = "notchFilter",       freq = 50)
  )
  expect_s3_class(pipeline, "PhysioPipeline")

  # applyPipeline must pass the working object as the FIRST POSITIONAL argument:
  # the ops take the PhysioExperiment as their first formal (`x`), so passing it
  # as a named `pe =` made every step fail with "unused argument (pe = ...)".
  out <- applyPipeline(pe, pipeline, verbose = FALSE)
  expect_s4_class(out, "PhysioExperiment")
  # each op ran and wrote its result assay (the input `raw` is preserved)
  expect_true(all(c("detrended", "filtered") %in%
                    SummarizedExperiment::assayNames(out)))
  expect_equal(dim(SummarizedExperiment::assay(out, "filtered")), dim(m))
})
