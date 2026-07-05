#' ICA-Based Artifact Removal for PhysioExperiment
#'
#' Functions for Independent Component Analysis (ICA) decomposition
#' and artifact removal using the fastICA algorithm.

#' Perform ICA decomposition
#'
#' Decomposes signals into independent components using the fastICA algorithm.
#'
#' @param pe A PhysioExperiment object.
#' @param n_components Number of components to extract. If NULL, uses all channels.
#' @param assay_name Name of the assay to decompose.
#' @param method ICA algorithm: "fastica" (default).
#' @param max_iter Maximum number of iterations.
#' @param ... Additional arguments passed to fastICA::fastICA.
#' @return A list with components:
#'   \item{S}{Source matrix (time x components)}
#'   \item{A}{Mixing matrix (channels x components)}
#'   \item{W}{Unmixing matrix (components x channels)}
#' @references Hyvarinen A, Oja E (2000). "Independent component analysis:
#'   algorithms and applications." \emph{Neural Networks}, 13(4-5), 411-430.
#' @seealso [removeICAComponents()] for reconstructing signals after removing
#'   artifact components, [icaDecompose()] for the built-in ICA implementation,
#'   [detectBadChannels()] for channel-level artifact detection.
#' @export
#' @examples
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(1000), nrow = 100, ncol = 10)),
#'   colData = S4Vectors::DataFrame(label = paste0("Ch", 1:10)),
#'   samplingRate = 256
#' )
#' \dontrun{
#' ica_result <- runICA(pe, n_components = 5)
#' }
runICA <- function(pe, n_components = NULL, assay_name = NULL,
                   method = "fastica", max_iter = 200, ...) {
  stopifnot(inherits(pe, "PhysioExperiment"))

  if (is.null(assay_name)) {
    assay_name <- defaultAssay(pe)
  }

  data <- SummarizedExperiment::assay(pe, assay_name)

  # Ensure 2D
  if (length(dim(data)) > 2) {
    data <- data[, , 1]
  }

  n_channels <- ncol(data)

  if (is.null(n_components)) {
    n_components <- n_channels
  }

  if (n_components > n_channels) {
    stop("n_components cannot exceed number of channels", call. = FALSE)
  }

  if (!requireNamespace("fastICA", quietly = TRUE)) {
    stop("Package 'fastICA' required for ICA. Install with: install.packages('fastICA')",
         call. = FALSE)
  }

  # Run fastICA
  ica_result <- fastICA::fastICA(data, n.comp = n_components,
                                  maxit = max_iter, ...)

  list(
    S = ica_result$S,
    A = ica_result$A,
    W = ica_result$W,
    n_components = n_components
  )
}

#' Remove ICA components
#'
#' Reconstructs signals after removing specified ICA components (e.g., artifacts).
#'
#' @param pe A PhysioExperiment object.
#' @param ica_result Result from runICA().
#' @param remove_components Integer vector of component indices to remove.
#' @param assay_name Name of the assay to reconstruct.
#' @param output_assay Name for the cleaned output assay.
#' @return PhysioExperiment with cleaned data in output_assay.
#' @references Hyvarinen A, Oja E (2000). "Independent component analysis:
#'   algorithms and applications." \emph{Neural Networks}, 13(4-5), 411-430.
#' @seealso [runICA()] for performing the ICA decomposition,
#'   [icaRemove()] for the alternative component removal implementation,
#'   [detectBadChannels()] for channel-level artifact detection.
#' @export
#' @examples
#' \dontrun{
#' pe <- runICA(pe)
#' pe_clean <- removeICAComponents(pe, ica_result, remove_components = c(1, 3))
#' }
removeICAComponents <- function(pe, ica_result, remove_components,
                                 assay_name = NULL, output_assay = "ica_cleaned") {
  stopifnot(inherits(pe, "PhysioExperiment"))

  if (is.null(assay_name)) {
    assay_name <- defaultAssay(pe)
  }

  # Zero out artifact components
  S_clean <- ica_result$S
  S_clean[, remove_components] <- 0

  # Reconstruct signal
  reconstructed <- S_clean %*% ica_result$A

  # Store in new assay
  SummarizedExperiment::assay(pe, output_assay) <- reconstructed

  pe
}
