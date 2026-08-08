# Golden-fixture generator for PhysioPreprocess (WSF-08).
#
# CORE RULE: each golden is an INDEPENDENT reference value, computed from a
# DIFFERENT implementation than the package function under test:
#   - filters:  signal::butter + signal::filtfilt (the reference DSP toolbox)
#   - average reference:  base-R  data - rowMeans(data)  subtraction
#   - ICA:  no external reference obtainable here -> characterization golden,
#           documented as such, pinned to a fixed seed (regression-guard only).
#
# All inputs are built under a FIXED SEED so regeneration is bit-stable.
# Run from the repo root, e.g.:
#   Rscript -e 'devtools::load_all("physio-ecosystem/PhysioPreprocess"); \
#               source("physio-ecosystem/PhysioPreprocess/data-raw/golden.R")'
#
# Writes physio-ecosystem/PhysioPreprocess/tests/testthat/_golden/<key>.rds
# plus a sibling <key>.dcf manifest recording source + tolerance.

suppressMessages({
  library(signal)
  library(PhysioCore)
})

# Locate the golden dir without relying on testthat::test_path (this script is
# sourced outside a testthat run). Mirror helper-golden.R's write_golden().
.pp_root <- local({
  # Prefer the loaded package's source; fall back to a fixed relative path.
  cand <- c(
    file.path("physio-ecosystem", "PhysioPreprocess"),
    "."
  )
  hit <- cand[file.exists(file.path(cand, "DESCRIPTION"))]
  if (length(hit) == 0) stop("Cannot locate PhysioPreprocess package root")
  normalizePath(hit[1])
})
.golden_dir_gen <- file.path(.pp_root, "tests", "testthat", "_golden")

write_golden_gen <- function(value, key, source, tol = 1e-8) {
  dir.create(.golden_dir_gen, showWarnings = FALSE, recursive = TRUE)
  saveRDS(value, file.path(.golden_dir_gen, paste0(key, ".rds")))
  write.dcf(
    data.frame(
      Key = key, Source = source, Tolerance = format(tol),
      Captured = format(Sys.time(), tz = "UTC", usetz = TRUE),
      R = as.character(getRversion()),
      stringsAsFactors = FALSE
    ),
    file.path(.golden_dir_gen, paste0(key, ".dcf"))
  )
  invisible(value)
}

# ---------------------------------------------------------------------------
# Shared deterministic input (used by every filter golden below).
# ---------------------------------------------------------------------------
make_filter_input <- function() {
  set.seed(42)
  n <- 1024L
  sr <- 250
  tt <- seq(0, (n - 1) / sr, length.out = n)
  # 10 Hz component (passband), 50 Hz line noise, plus deterministic jitter.
  sig <- sin(2 * pi * 10 * tt) + 0.5 * sin(2 * pi * 50 * tt) + 0.2 * rnorm(n)
  mat <- cbind(sig, 0.8 * sig + 0.1 * rnorm(n))
  list(mat = mat, n = n, sr = sr)
}

inp <- make_filter_input()
mat <- inp$mat
sr <- inp$sr
nyq <- sr / 2

# apply() inherits column names from cbind(); the package assay is a plain
# unnamed matrix. Strip dimnames so goldens compare on numeric content only.
strip_dn <- function(m) {
  dimnames(m) <- NULL
  m
}

# ---------------------------------------------------------------------------
# 1-3. Butterworth filters  (reference: signal::butter + signal::filtfilt)
#      These target butterworthFilter(..., use_sos = FALSE), the ba-form path,
#      which is an exact zero-phase Butterworth via signal::filtfilt. The
#      reference here is that same DSP toolbox called directly and independently.
# ---------------------------------------------------------------------------

# 1. Bandpass 1-40 Hz, order 4
bf_bp <- signal::butter(n = 4, W = c(1 / nyq, 40 / nyq), type = "pass")
ref_bp <- strip_dn(apply(mat, 2, function(v) as.numeric(signal::filtfilt(bf_bp, v))))
write_golden_gen(
  ref_bp, "butter_bandpass_ba",
  source = "signal::filtfilt(signal::butter(4, c(1,40)/125, 'pass'), x) applied per column; independent DSP reference for butterworthFilter(use_sos=FALSE)",
  tol = 1e-8
)

# 2. Lowpass 30 Hz, order 4
bf_lp <- signal::butter(n = 4, W = 30 / nyq, type = "low")
ref_lp <- strip_dn(apply(mat, 2, function(v) as.numeric(signal::filtfilt(bf_lp, v))))
write_golden_gen(
  ref_lp, "butter_lowpass_ba",
  source = "signal::filtfilt(signal::butter(4, 30/125, 'low'), x) per column; independent DSP reference for butterworthFilter(type='low', use_sos=FALSE)",
  tol = 1e-8
)

# 3. Highpass 0.5 Hz, order 4
bf_hp <- signal::butter(n = 4, W = 0.5 / nyq, type = "high")
ref_hp <- strip_dn(apply(mat, 2, function(v) as.numeric(signal::filtfilt(bf_hp, v))))
write_golden_gen(
  ref_hp, "butter_highpass_ba",
  source = "signal::filtfilt(signal::butter(4, 0.5/125, 'high'), x) per column; independent DSP reference for butterworthFilter(type='high', use_sos=FALSE)",
  tol = 1e-8
)

# ---------------------------------------------------------------------------
# 4. Notch filter (reference: signal::butter(type='stop') + signal::filtfilt)
#    notchFilter internally builds an order-4 Butterworth bandstop over
#    [f - bw/2, f + bw/2] and zero-phase filters it. We reconstruct that exact
#    bandstop design independently from the signal package.
# ---------------------------------------------------------------------------
notch_freq <- 50
notch_bw <- 2
lo <- (notch_freq - notch_bw / 2) / nyq
hi <- (notch_freq + notch_bw / 2) / nyq
bf_notch <- signal::butter(n = 4, W = c(lo, hi), type = "stop")
ref_notch <- strip_dn(apply(mat, 2, function(v) as.numeric(signal::filtfilt(bf_notch, v))))
write_golden_gen(
  ref_notch, "notch_50",
  source = "signal::filtfilt(signal::butter(4, c(49,51)/125, 'stop'), x) per column; independent DSP reference for notchFilter(freq=50, bandwidth=2)",
  tol = 1e-8
)

# ---------------------------------------------------------------------------
# 5. Average re-reference (reference: base-R rowMeans subtraction)
#    Common average reference is, by definition, x_ch(t) - mean_ch x_ch(t).
#    Computed independently in base R.
# ---------------------------------------------------------------------------
set.seed(7)
n_re <- 200L
nch_re <- 8L
mat_re <- matrix(rnorm(n_re * nch_re), n_re, nch_re)
ref_avg <- mat_re - matrix(rowMeans(mat_re), n_re, nch_re)
write_golden_gen(
  ref_avg, "rereference_average",
  source = "base-R common average reference: mat - rowMeans(mat); independent reference for rereference(ref_type='average')",
  tol = 1e-8
)

# ---------------------------------------------------------------------------
# 6. ICA decomposition (characterization golden)
#    FastICA has no closed-form external reference obtainable here (rotation /
#    sign / permutation ambiguity, iterative fixed-point). The package's
#    .fastICA draws its random init via base rnorm, so a set.seed() immediately
#    before icaDecompose() makes the result bit-reproducible. This golden is a
#    REGRESSION GUARD ONLY, not an accuracy check against an external tool.
# ---------------------------------------------------------------------------
set.seed(123)
n_ica <- 200L
s1 <- sin(2 * pi * 3 * seq_len(n_ica) / 250)
s2 <- sign(sin(2 * pi * 7 * seq_len(n_ica) / 250))
s3 <- rnorm(n_ica)
S <- cbind(s1, s2, s3)
A_true <- matrix(c(1, 0.5, 0.3,
                   0.4, 1, 0.6,
                   0.2, 0.7, 1), nrow = 3, byrow = TRUE)
X_ica <- S %*% t(A_true)
pe_ica <- PhysioExperiment(
  assays = S4Vectors::SimpleList(raw = array(X_ica, dim = c(n_ica, 3L))),
  rowData = S4Vectors::DataFrame(t = seq_len(n_ica)),
  colData = S4Vectors::DataFrame(label = paste0("Ch", 1:3)),
  samplingRate = 250
)
set.seed(99)
ica_components <- icaDecompose(pe_ica, n_components = 3L)$components
write_golden_gen(
  ica_components, "ica_characterization",
  source = "characterization: icaDecompose FastICA components on fixed-seed 3-source mixture (set.seed(99) init); regression-guard only, no external reference",
  tol = 1e-6
)

message("PhysioPreprocess goldens written to: ", .golden_dir_gen)
