#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// IIR filtering kernels. These replace pure-R per-sample loops (and the
// ts-object overhead of stats::filter) on the package's hottest paths. Each
// kernel mirrors the semantics of the R reference implementation it replaces
// (kept as .sosfilt_r / .lfilter_r / .sosfiltfilt_r for parity tests).
//
// The cascade is applied sample-by-sample with all section states carried in
// a small local array ("fused" schedule). For an LTI cascade this performs
// exactly the same arithmetic, on the same operands, as filtering the whole
// signal through one section at a time, so results are bit-identical to the
// section-sequential references while touching the data only once per pass.

struct SosCoef {
  double b0, b1, b2, a1, a2;
};

static std::vector<SosCoef> normalize_sos(const NumericMatrix& sos) {
  const int S = sos.nrow();
  std::vector<SosCoef> c(S);
  for (int s = 0; s < S; ++s) {
    const double a0 = sos(s, 3);
    c[s].b0 = sos(s, 0) / a0;
    c[s].b1 = sos(s, 1) / a0;
    c[s].b2 = sos(s, 2) / a0;
    c[s].a1 = sos(s, 4) / a0;
    c[s].a2 = sos(s, 5) / a0;
  }
  return c;
}

// One fused DF2T cascade pass over x[0..n) in place with stride +1 or -1
// (dir = -1 walks the signal backward, equivalent to filtering the reversed
// signal). state has 2 entries per section and is updated in place.
static inline void sos_pass_fused(const std::vector<SosCoef>& c, double* x,
                                  const R_xlen_t n, const int dir,
                                  double* state) {
  const int S = static_cast<int>(c.size());
  double* p = (dir > 0) ? x : x + n - 1;
  for (R_xlen_t m = 0; m < n; ++m, p += dir) {
    double v = *p;
    for (int s = 0; s < S; ++s) {
      const SosCoef& k = c[s];
      double& z1 = state[2 * s];
      double& z2 = state[2 * s + 1];
      const double y = k.b0 * v + z1;
      z1 = k.b1 * v - k.a1 * y + z2;
      z2 = k.b2 * v - k.a2 * y;
      v = y;
    }
    *p = v;
  }
}

// Causal SOS cascade (DF2T) with carried state; mirrors sosfilt()'s R loop.
// [[Rcpp::export]]
List cpp_sosfilt(NumericMatrix sos, NumericVector x, NumericMatrix zi) {
  const int S = sos.nrow();
  std::vector<SosCoef> c = normalize_sos(sos);
  std::vector<double> state(2 * S);
  for (int s = 0; s < S; ++s) {
    state[2 * s] = zi(s, 0);
    state[2 * s + 1] = zi(s, 1);
  }
  NumericVector y = clone(x);
  sos_pass_fused(c, REAL(y), y.size(), 1, state.data());
  NumericMatrix zo(S, 2);
  for (int s = 0; s < S; ++s) {
    zo(s, 0) = state[2 * s];
    zo(s, 1) = state[2 * s + 1];
  }
  return List::create(_["y"] = y, _["zi"] = zo);
}

// General-order transposed direct-form II with carried state; mirrors
// lfilter()'s R loop. b and a must be pre-normalized (a[0] == 1) and padded
// to equal length by the R wrapper.
// [[Rcpp::export]]
List cpp_lfilter(NumericVector b, NumericVector a, NumericVector x,
                 NumericVector zi) {
  const int K = static_cast<int>(zi.size());
  NumericVector z = clone(zi);
  const R_xlen_t nx = x.size();
  NumericVector y(nx);
  if (K == 0) {
    const double b0 = b[0];
    for (R_xlen_t m = 0; m < nx; ++m) y[m] = b0 * x[m];
  } else {
    for (R_xlen_t m = 0; m < nx; ++m) {
      const double xm = x[m];
      const double ym = b[0] * xm + z[0];
      y[m] = ym;
      for (int i = 0; i < K - 1; ++i) {
        z[i] = b[i + 1] * xm + z[i + 1] - a[i + 1] * ym;
      }
      z[K - 1] = b[K] * xm - a[K] * ym;
    }
  }
  return List::create(_["y"] = y, _["zi"] = z);
}

// --- ba-form zero-phase filtering, exact port of signal::filtfilt ----------
//
//   y = filter(b, a, c(x, zeros(2*max(length(a), length(b)))))
//   y = rev(filter(b, a, rev(y)))[1:length(x)]
//
// where filter(b, a, .) is the standard difference equation with zero initial
// state (here computed with the DF2T recursion).

static void ba_normalize(const NumericVector& b, const NumericVector& a,
                         std::vector<double>& bn, std::vector<double>& an) {
  if (a.size() == 0 || b.size() == 0) {
    stop("`b` and `a` must be non-empty coefficient vectors.");
  }
  const int n = static_cast<int>(std::max(b.size(), a.size()));
  bn.assign(n, 0.0);
  an.assign(n, 0.0);
  const double a0 = a[0];
  for (int i = 0; i < b.size(); ++i) bn[i] = b[i] / a0;
  for (int i = 0; i < a.size(); ++i) an[i] = a[i] / a0;
}

// DF2T with zero initial state, in place over x[0..n), walking dir = +/-1.
static void lfilter_zero_pass(const std::vector<double>& b,
                              const std::vector<double>& a, double* x,
                              const R_xlen_t n, const int dir) {
  const int K = static_cast<int>(b.size()) - 1;
  if (K <= 0) {
    for (R_xlen_t m = 0; m < n; ++m) x[m] *= b[0];
    return;
  }
  std::vector<double> z(K, 0.0);
  double* p = (dir > 0) ? x : x + n - 1;
  for (R_xlen_t m = 0; m < n; ++m, p += dir) {
    const double xm = *p;
    const double ym = b[0] * xm + z[0];
    for (int i = 0; i < K - 1; ++i) {
      z[i] = b[i + 1] * xm + z[i + 1] - a[i + 1] * ym;
    }
    z[K - 1] = b[K] * xm - a[K] * ym;
    *p = ym;
  }
}

static void filtfilt_ba_vec(const std::vector<double>& bn,
                            const std::vector<double>& an, const int n_pad,
                            const double* x, const R_xlen_t n, double* out,
                            std::vector<double>& ext) {
  const R_xlen_t ne = n + n_pad;
  ext.assign(ne, 0.0);
  std::copy(x, x + n, ext.begin());
  lfilter_zero_pass(bn, an, ext.data(), ne, 1);   // forward
  lfilter_zero_pass(bn, an, ext.data(), ne, -1);  // == filter(rev(y)) reversed
  std::copy(ext.begin(), ext.begin() + n, out);
}

// [[Rcpp::export]]
NumericVector cpp_filtfilt_ba(NumericVector b, NumericVector a,
                              NumericVector x) {
  std::vector<double> bn, an;
  ba_normalize(b, a, bn, an);
  const int n_pad = 2 * static_cast<int>(std::max(a.size(), b.size()));
  const R_xlen_t n = x.size();
  NumericVector out(n);
  std::vector<double> ext;
  filtfilt_ba_vec(bn, an, n_pad, REAL(x), n, REAL(out), ext);
  return out;
}

// [[Rcpp::export]]
NumericMatrix cpp_filtfilt_ba_mat(NumericVector b, NumericVector a,
                                  NumericMatrix x) {
  std::vector<double> bn, an;
  ba_normalize(b, a, bn, an);
  const int n_pad = 2 * static_cast<int>(std::max(a.size(), b.size()));
  const R_xlen_t n = x.nrow();
  const int nc = x.ncol();
  NumericMatrix out(n, nc);
  const double* xp = REAL(x);
  double* op = REAL(out);
#ifdef _OPENMP
#pragma omp parallel
  {
    std::vector<double> ext;
#pragma omp for schedule(static)
    for (int col = 0; col < nc; ++col) {
      filtfilt_ba_vec(bn, an, n_pad, xp + static_cast<R_xlen_t>(col) * n, n,
                      op + static_cast<R_xlen_t>(col) * n, ext);
    }
  }
#else
  {
    std::vector<double> ext;
    for (int col = 0; col < nc; ++col) {
      filtfilt_ba_vec(bn, an, n_pad, xp + static_cast<R_xlen_t>(col) * n, n,
                      op + static_cast<R_xlen_t>(col) * n, ext);
    }
  }
#endif
  return out;
}

// Zero-phase SOS filtering of one channel: odd-reflection padding of nfact
// samples, zero-state forward pass, zero-state backward pass, unpad. Exact
// port of the R .sosfiltfilt (pad once, both passes start from zero state).
static void sosfiltfilt_vec(const std::vector<SosCoef>& c, const int S,
                            const double* x, const R_xlen_t n, const int nfact,
                            double* out, std::vector<double>& ext,
                            std::vector<double>& state) {
  const R_xlen_t ne = n + 2 * static_cast<R_xlen_t>(nfact);
  ext.resize(ne);
  for (int i = 0; i < nfact; ++i) ext[i] = 2.0 * x[0] - x[nfact - i];
  std::copy(x, x + n, ext.begin() + nfact);
  for (int i = 0; i < nfact; ++i) ext[nfact + n + i] = 2.0 * x[n - 1] - x[n - 2 - i];

  std::fill(state.begin(), state.begin() + 2 * S, 0.0);
  sos_pass_fused(c, ext.data(), ne, 1, state.data());   // forward
  std::fill(state.begin(), state.begin() + 2 * S, 0.0);
  sos_pass_fused(c, ext.data(), ne, -1, state.data());  // backward
  std::copy(ext.begin() + nfact, ext.begin() + nfact + n, out);
}

// [[Rcpp::export]]
NumericVector cpp_sosfiltfilt(NumericMatrix sos, NumericVector x, int nfact) {
  const int S = sos.nrow();
  std::vector<SosCoef> c = normalize_sos(sos);
  const R_xlen_t n = x.size();
  if (nfact < 0) nfact = 0;
  if (nfact >= n) nfact = (n > 0) ? static_cast<int>(n - 1) : 0;
  NumericVector out(n);
  std::vector<double> ext;
  std::vector<double> state(2 * S);
  sosfiltfilt_vec(c, S, REAL(x), n, nfact, REAL(out), ext, state);
  return out;
}

// Column-wise zero-phase SOS filtering of a whole time x channels matrix in
// one call. Channels are independent, so they are processed in parallel when
// OpenMP is available (thread count respects OMP_NUM_THREADS / OMP_THREAD_LIMIT).
// [[Rcpp::export]]
NumericMatrix cpp_sosfiltfilt_mat(NumericMatrix sos, NumericMatrix x,
                                  int nfact) {
  const int S = sos.nrow();
  const std::vector<SosCoef> c = normalize_sos(sos);
  const R_xlen_t n = x.nrow();
  const int nc = x.ncol();
  if (nfact < 0) nfact = 0;
  if (nfact >= n) nfact = (n > 0) ? static_cast<int>(n - 1) : 0;
  NumericMatrix out(n, nc);
  const double* xp = REAL(x);
  double* op = REAL(out);
#ifdef _OPENMP
#pragma omp parallel
  {
    std::vector<double> ext;
    std::vector<double> state(2 * S);
#pragma omp for schedule(static)
    for (int col = 0; col < nc; ++col) {
      sosfiltfilt_vec(c, S, xp + static_cast<R_xlen_t>(col) * n, n, nfact,
                      op + static_cast<R_xlen_t>(col) * n, ext, state);
    }
  }
#else
  {
    std::vector<double> ext;
    std::vector<double> state(2 * S);
    for (int col = 0; col < nc; ++col) {
      sosfiltfilt_vec(c, S, xp + static_cast<R_xlen_t>(col) * n, n, nfact,
                      op + static_cast<R_xlen_t>(col) * n, ext, state);
    }
  }
#endif
  return out;
}
