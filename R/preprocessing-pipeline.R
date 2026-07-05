#' Preprocessing Pipeline for PhysioExperiment
#'
#' Functions for creating and applying preprocessing pipelines.

#' Create preprocessing pipeline
#'
#' Creates a preprocessing pipeline specification that can be applied to data.
#'
#' @param ... Preprocessing steps as named function calls.
#' @return A preprocessing pipeline object (list).
#' @references Oppenheim AV, Willsky AS, Nawab SH (1997). "Signals and Systems."
#'   2nd ed. Prentice Hall.
#' @seealso [applyPipeline()] for executing a pipeline on data,
#'   [filterSignals()] and [detrendSignals()] for individual preprocessing
#'   steps that can be included in a pipeline.
#' @export
#' @examples
#' pipeline <- createPipeline(
#'   filter = list(fn = "filterSignals", lowcut = 1, highcut = 40),
#'   detrend = list(fn = "detrendSignals", method = "linear")
#' )
createPipeline <- function(...) {
  steps <- list(...)

  if (length(steps) == 0) {
    stop("Pipeline must contain at least one step", call. = FALSE)
  }

  # Validate steps
  for (i in seq_along(steps)) {
    step <- steps[[i]]
    if (!is.list(step) || !"fn" %in% names(step)) {
      stop("Each pipeline step must be a list with 'fn' element", call. = FALSE)
    }
  }

  structure(
    list(
      steps = steps,
      n_steps = length(steps)
    ),
    class = "PhysioPipeline"
  )
}

#' Apply preprocessing pipeline
#'
#' Applies a preprocessing pipeline to a PhysioExperiment object.
#'
#' @param pe A PhysioExperiment object.
#' @param pipeline A pipeline object created with createPipeline().
#' @param verbose Logical. If TRUE, prints progress messages.
#' @return PhysioExperiment with all pipeline steps applied.
#' @references Oppenheim AV, Willsky AS, Nawab SH (1997). "Signals and Systems."
#'   2nd ed. Prentice Hall.
#' @seealso [createPipeline()] for defining the pipeline steps,
#'   [butterworthFilter()] and [resample()] for common preprocessing
#'   operations used within pipelines.
#' @export
#' @examples
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(1000), nrow = 100, ncol = 10)),
#'   colData = S4Vectors::DataFrame(label = paste0("Ch", 1:10)),
#'   samplingRate = 256
#' )
#' \dontrun{
#' pipeline <- createPipeline(
#'   detrend = list(fn = "detrendSignals", method = "linear")
#' )
#' pe_processed <- applyPipeline(pe, pipeline)
#' }
applyPipeline <- function(pe, pipeline, verbose = FALSE) {
  stopifnot(inherits(pe, "PhysioExperiment"))
  stopifnot(inherits(pipeline, "PhysioPipeline"))

  result <- pe

  for (i in seq_along(pipeline$steps)) {
    step <- pipeline$steps[[i]]
    step_name <- names(pipeline$steps)[i]

    if (verbose) {
      message(sprintf("Applying step %d/%d: %s",
                      i, pipeline$n_steps, step_name))
    }

    fn_name <- step$fn
    fn <- match.fun(fn_name)

    # Get arguments (excluding 'fn')
    args <- step[names(step) != "fn"]
    args <- c(list(pe = result), args)

    result <- do.call(fn, args)
  }

  result
}

#' Print pipeline summary
#'
#' @param x A PhysioPipeline object.
#' @param ... Additional arguments (ignored).
#' @return Invisible x.
#' @seealso [createPipeline()] for creating pipelines,
#'   [applyPipeline()] for executing pipelines.
#' @export
print.PhysioPipeline <- function(x, ...) {
  cat("PhysioPipeline with", x$n_steps, "steps:\n")
  for (i in seq_along(x$steps)) {
    step <- x$steps[[i]]
    step_name <- names(x$steps)[i]
    cat(sprintf("  %d. %s (%s)\n", i, step_name, step$fn))
  }
  invisible(x)
}
