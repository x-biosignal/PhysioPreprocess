#' Re-referencing Operations for EEG Data
#'
#' Functions for changing the reference electrode in EEG recordings.
#' Re-referencing is a common preprocessing step that affects the spatial
#' distribution of the signal.
#'
#' @references Nunez, P.L. & Srinivasan, R. (2006). "Electric Fields of the
#'   Brain." 2nd ed. Oxford University Press.

#' Re-reference EEG data
#'
#' Changes the reference electrode for EEG recordings. Supports common
#' re-referencing schemes including average reference, linked mastoids,
#' and single electrode reference.
#'
#' @param x A PhysioExperiment object.
#' @param ref_type Type of re-referencing: "average" (common average reference),
#'   "robust" (PREP robust average reference: iteratively detect and exclude
#'   high-variance bad channels from the average; Bigdely-Shamlo et al., 2015),
#'   "median" (channel-wise median reference, robust to outlier channels),
#'   "channel" (single channel), "channels" (average of specified channels),
#'   or "REST" (Reference Electrode Standardization Technique).
#' @param robust_noise_sd For "robust", the threshold (robust SDs above the
#'   median channel variance) for flagging a bad channel (default: 4).
#' @param robust_max_iter For "robust", the maximum number of PREP iterations
#'   (default: 5).
#' @param ref_channels For "channel" or "channels" type, the channel name(s)
#'   or index/indices to use as reference.
#' @param exclude Channels to exclude from average reference calculation
#'   (e.g., non-EEG channels like EOG, EMG).
#' @param input_assay Input assay name. If NULL, uses default assay.
#' @param output_assay Output assay name. Default is "rereferenced".
#' @param keep_ref Logical. If TRUE, keeps the original reference channel(s)
#'   in the output (zeroed). If FALSE, removes them.
#' @return A `PhysioExperiment` object with re-referenced data stored in a new
#'   assay named `output_assay`. The `reference` and `previous_reference`
#'   metadata fields are updated. When `keep_ref = FALSE` and using channel-based
#'   reference, the reference channel(s) are removed and a new object with
#'   reduced channel count is returned.
#' @references Nunez, P.L. & Srinivasan, R. (2006). "Electric Fields of the
#'   Brain." 2nd ed. Oxford University Press.
#' @seealso [getCurrentReference()] for querying the current reference,
#'   [isAverageReferenced()] for checking average reference status,
#'   [butterworthFilter()] for frequency-domain preprocessing.
#' @details
#' Re-referencing transforms the data by subtracting a reference signal from
#' each channel. The choice of reference affects the spatial distribution
#' and interpretation of the signal.
#'
#' **Average reference** ("average"): Subtracts the mean of all channels at
#' each time point. This is commonly used for high-density EEG and provides
#' a reference-independent measure, but requires good spatial sampling.
#'
#' **Single channel reference** ("channel"): Subtracts the signal from a
#' specified electrode. Common choices include Cz, linked mastoids (A1+A2)/2,
#' or nose reference.
#'
#' **Multi-channel reference** ("channels"): Subtracts the average of multiple
#' specified channels. Useful for linked mastoids or other custom references.
#'
#' @export
#' @examples
#' # Create example EEG data
#' set.seed(123)
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(1000), nrow = 100, ncol = 10)),
#'   colData = S4Vectors::DataFrame(
#'     label = c("Fp1", "Fp2", "F3", "F4", "C3", "C4", "P3", "P4", "O1", "O2"),
#'     type = rep("EEG", 10)
#'   ),
#'   samplingRate = 256
#' )
#'
#' # Apply average reference
#' pe_avg <- rereference(pe, ref_type = "average")
#'
#' # Re-reference to a single channel (Cz)
#' pe_cz <- rereference(pe, ref_type = "channel", ref_channels = "C3")
#'
#' # Re-reference to linked mastoids (if available)
#' # pe_linked <- rereference(pe, ref_type = "channels",
#' #                          ref_channels = c("M1", "M2"))
rereference <- function(x, ref_type = c("average", "robust", "median",
                                        "channel", "channels", "REST"),
                        ref_channels = NULL, exclude = NULL,
                        input_assay = NULL, output_assay = "rereferenced",
                        keep_ref = TRUE, robust_noise_sd = 4,
                        robust_max_iter = 5) {
  stopifnot(inherits(x, "PhysioExperiment"))
  ref_type <- match.arg(ref_type)

  # Get input data
  if (is.null(input_assay)) {
    input_assay <- defaultAssay(x)
  }

  if (is.na(input_assay)) {
    stop("No assays available", call. = FALSE)
  }

  data <- SummarizedExperiment::assay(x, input_assay)
  dims <- dim(data)
  ch_names <- channelNames(x)
  n_channels <- length(ch_names)

  # Convert channel names to indices
  if (is.character(exclude)) {
    exclude_idx <- match(exclude, ch_names)
    exclude_idx <- exclude_idx[!is.na(exclude_idx)]
  } else if (is.numeric(exclude)) {
    exclude_idx <- exclude
  } else {
    exclude_idx <- integer(0)
  }

  if (is.character(ref_channels)) {
    ref_idx <- match(ref_channels, ch_names)
    if (any(is.na(ref_idx))) {
      missing <- ref_channels[is.na(ref_idx)]
      stop("Reference channels not found: ", paste(missing, collapse = ", "), call. = FALSE)
    }
  } else if (is.numeric(ref_channels)) {
    ref_idx <- ref_channels
  } else {
    ref_idx <- NULL
  }

  # Validate ref_type and ref_channels
  if (ref_type == "channel") {
    if (is.null(ref_idx) || length(ref_idx) != 1) {
      stop("ref_type 'channel' requires exactly one ref_channels value", call. = FALSE)
    }
  } else if (ref_type == "channels") {
    if (is.null(ref_idx) || length(ref_idx) < 1) {
      stop("ref_type 'channels' requires at least one ref_channels value", call. = FALSE)
    }
  }

  # Apply re-referencing based on data dimensionality
  excluded_idx <- integer(0)
  rr <- function(mat) {
    .rereferenceMatrix(mat, ref_type, ref_idx, exclude_idx,
                       robust_noise_sd = robust_noise_sd,
                       robust_max_iter = robust_max_iter)
  }
  if (length(dims) == 2) {
    result <- rr(data)
    excluded_idx <- attr(result, "excluded")
    attr(result, "excluded") <- NULL
  } else if (length(dims) == 3) {
    result <- array(NA_real_, dim = dims)
    for (s in seq_len(dims[3])) {
      slice <- rr(data[, , s])
      excluded_idx <- union(excluded_idx, attr(slice, "excluded"))
      result[, , s] <- slice
    }
  } else if (length(dims) == 4) {
    result <- array(NA_real_, dim = dims)
    for (ep in seq_len(dims[3])) {
      for (s in seq_len(dims[4])) {
        slice <- rr(data[, , ep, s])
        excluded_idx <- union(excluded_idx, attr(slice, "excluded"))
        result[, , ep, s] <- slice
      }
    }
  } else {
    stop("Data must be 2D, 3D, or 4D", call. = FALSE)
  }
  excluded_channels <- if (length(excluded_idx)) ch_names[excluded_idx] else NULL

  # Handle reference channel removal
  if (!keep_ref && ref_type %in% c("channel", "channels") && !is.null(ref_idx)) {
    # Remove reference channels from data
    keep_idx <- setdiff(seq_len(n_channels), ref_idx)

    if (length(dims) == 2) {
      result <- result[, keep_idx, drop = FALSE]
    } else if (length(dims) == 3) {
      result <- result[, keep_idx, , drop = FALSE]
    } else if (length(dims) == 4) {
      result <- result[, keep_idx, , , drop = FALSE]
    }

    # Subset colData (channel metadata)
    col_data <- SummarizedExperiment::colData(x)[keep_idx, , drop = FALSE]

    # Subset all existing assays to match new channel count
    old_assays <- SummarizedExperiment::assays(x)
    new_assays <- list()
    for (nm in names(old_assays)) {
      a <- old_assays[[nm]]
      adims <- dim(a)
      if (length(adims) == 2) {
        new_assays[[nm]] <- a[, keep_idx, drop = FALSE]
      } else if (length(adims) == 3) {
        new_assays[[nm]] <- a[, keep_idx, , drop = FALSE]
      } else if (length(adims) == 4) {
        new_assays[[nm]] <- a[, keep_idx, , , drop = FALSE]
      }
    }
    # Add the rereferenced result
    new_assays[[output_assay]] <- result

    # Create new PhysioExperiment with reduced channels
    x <- PhysioExperiment(
      assays = new_assays,
      rowData = SummarizedExperiment::rowData(x),
      colData = col_data,
      metadata = S4Vectors::metadata(x),
      samplingRate = samplingRate(x)
    )
  } else {
    # Store result in existing object
    assays <- SummarizedExperiment::assays(x)
    assays[[output_assay]] <- result
    SummarizedExperiment::assays(x) <- assays
  }

  # Update metadata
  meta <- S4Vectors::metadata(x)
  old_ref <- meta$reference
  if (ref_type == "average") {
    new_ref <- "average"
  } else if (ref_type == "channel") {
    new_ref <- ch_names[ref_idx]
  } else if (ref_type == "channels") {
    new_ref <- paste(ch_names[ref_idx], collapse = "+")
  } else {
    new_ref <- ref_type
  }
  meta$reference <- new_ref
  meta$previous_reference <- old_ref
  meta$reference_excluded_channels <- excluded_channels
  S4Vectors::metadata(x) <- meta

  # Log the re-reference step (method + excluded channels) into provenance.
  x <- PhysioCore::withProvenance(x, x, step = "rereference",
                                  params = list(ref_type = ref_type,
                                                excluded_channels = excluded_channels))
  x
}

#' Variance-based bad-channel indices (modality-agnostic PREP-style detector)
#' @noRd
.robustBadChannels <- function(data, noise_sd, exclude_idx) {
  v <- apply(data, 2, stats::var, na.rm = TRUE)
  v[exclude_idx] <- NA
  med <- stats::median(v, na.rm = TRUE)
  spread <- stats::mad(v, na.rm = TRUE)
  if (!is.finite(spread) || spread <= 0) spread <- stats::sd(v, na.rm = TRUE)
  if (!is.finite(spread) || spread <= 0) return(integer(0))
  which(!is.na(v) & v > med + noise_sd * spread)
}

#' Re-reference a 2D matrix
#' @noRd
.rereferenceMatrix <- function(data, ref_type, ref_idx, exclude_idx,
                               robust_noise_sd = 4, robust_max_iter = 5) {
  n_time <- nrow(data)
  n_channels <- ncol(data)

  result <- data

  if (ref_type == "average") {
    # Calculate average excluding specified channels
    include_idx <- setdiff(seq_len(n_channels), exclude_idx)
    ref_signal <- rowMeans(data[, include_idx, drop = FALSE], na.rm = TRUE)

    # Subtract average from all channels
    for (ch in seq_len(n_channels)) {
      result[, ch] <- data[, ch] - ref_signal
    }
  } else if (ref_type == "robust") {
    # PREP robust average reference: iteratively exclude high-variance bad
    # channels from the average (Bigdely-Shamlo et al., 2015).
    bad_set <- integer(0)
    for (iter in seq_len(robust_max_iter)) {
      new_bad <- .robustBadChannels(result, robust_noise_sd, exclude_idx)
      if (setequal(new_bad, bad_set)) break
      bad_set <- new_bad
      good_idx <- setdiff(seq_len(n_channels), union(bad_set, exclude_idx))
      if (length(good_idx) == 0) {
        stop("Robust reference excluded all channels; loosen robust_noise_sd.",
             call. = FALSE)
      }
      ref_signal <- rowMeans(data[, good_idx, drop = FALSE], na.rm = TRUE)
      for (ch in seq_len(n_channels)) result[, ch] <- data[, ch] - ref_signal
    }
    attr(result, "excluded") <- bad_set
  } else if (ref_type == "median") {
    # Channel-wise median reference (robust to outlier channels).
    ref_signal <- apply(data, 1, stats::median, na.rm = TRUE)
    for (ch in seq_len(n_channels)) {
      result[, ch] <- data[, ch] - ref_signal
    }
  } else if (ref_type == "channel") {
    # Single channel reference
    ref_signal <- data[, ref_idx]

    for (ch in seq_len(n_channels)) {
      result[, ch] <- data[, ch] - ref_signal
    }
  } else if (ref_type == "channels") {
    # Average of multiple channels
    ref_signal <- rowMeans(data[, ref_idx, drop = FALSE], na.rm = TRUE)

    for (ch in seq_len(n_channels)) {
      result[, ch] <- data[, ch] - ref_signal
    }
  } else if (ref_type == "REST") {
    # REST (Reference Electrode Standardization Technique)
    # This is a simplified version - full REST requires lead field matrix
    warning("Full REST requires lead field matrix. Using average reference instead.",
            call. = FALSE)
    include_idx <- setdiff(seq_len(n_channels), exclude_idx)
    ref_signal <- rowMeans(data[, include_idx, drop = FALSE], na.rm = TRUE)

    for (ch in seq_len(n_channels)) {
      result[, ch] <- data[, ch] - ref_signal
    }
  }

  result
}

#' Get current reference
#'
#' Returns the current reference electrode information.
#'
#' @param x A PhysioExperiment object.
#' @return A character string describing the current reference (e.g., "average",
#'   "Cz", "M1+M2"), or `NULL` if no reference has been set.
#' @references Nunez, P.L. & Srinivasan, R. (2006). "Electric Fields of the
#'   Brain." 2nd ed. Oxford University Press.
#' @seealso [rereference()] for changing the reference,
#'   [isAverageReferenced()] for checking average reference status.
#' @export
#' @examples
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(100), nrow = 10, ncol = 10)),
#'   samplingRate = 100
#' )
#' pe <- setReference(pe, "Cz")
#' getCurrentReference(pe)  # "Cz"
getCurrentReference <- function(x) {
  stopifnot(inherits(x, "PhysioExperiment"))
  S4Vectors::metadata(x)$reference
}

#' Check if data is average referenced
#'
#' @param x A PhysioExperiment object.
#' @return A logical scalar: `TRUE` if the reference metadata is set to
#'   "average", `FALSE` otherwise.
#' @references Nunez, P.L. & Srinivasan, R. (2006). "Electric Fields of the
#'   Brain." 2nd ed. Oxford University Press.
#' @seealso [rereference()] for changing the reference,
#'   [getCurrentReference()] for querying the current reference.
#' @export
#' @examples
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(100), nrow = 10, ncol = 10)),
#'   samplingRate = 100
#' )
#' pe <- rereference(pe, ref_type = "average")
#' isAverageReferenced(pe)  # TRUE
isAverageReferenced <- function(x) {
  stopifnot(inherits(x, "PhysioExperiment"))
  ref <- S4Vectors::metadata(x)$reference
  !is.null(ref) && ref == "average"
}
