#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>
using namespace Rcpp;

// Exact 1-D grouping by penalized optimal partitioning (PELT), heteroscedastic
// Gaussian cost. Expected O(n) when groups are plentiful; worst case O(n^2)
// when the data is one giant group (pruning cannot help there).
//
// Model: identical per-point shot-noise sigma as fast_1D_iaac_grouping_cpp:
//     sigma_i = mz_i * p0 * max(1, sqrt(Isat / I_i)) / 1e6        (Da)
// Objective (globally minimized, unlike greedy agglomeration):
//     F = sum_over_segments [ sum_i (mz_i - mu_seg)^2 / sigma_i^2 ] + beta * (#segments)
// with mu_seg the inverse-variance mean of its segment. beta plays the role of
// IAAC's k: larger beta -> fewer, wider groups. beta <= 0 selects the
// BIC-style default beta = 2 * log(n).
//
// Numerical note: prefix-sum O(1) segment costs are the textbook approach but
// are catastrophically cancellative here (sum of (mz/sigma)^2 ~ 1e10 per point
// while segment costs are O(#points)). Instead each surviving PELT candidate
// keeps running sums CENTERED at its segment's first m/z, so competitive
// (short, tight) segments are computed from small, well-conditioned numbers;
// long bad segments lose precision only where their cost is huge and pruned
// anyway. The per-step work is O(|candidates|) either way, so this costs
// nothing asymptotically.

namespace {
struct PeltPoint { double mz, w; int original_idx; };   // w = 1 / sigma^2

struct PeltCand {
  int    tau;        // segment starts at sorted point index tau
  double F;          // optimal cost of points [0, tau)
  double W, V, SS;   // running sums over the open segment, centered at xref
  double xref;       // first m/z of the segment
  double lastCost;   // cost of segment [tau, t) at the current t
};
} // anonymous namespace 

//' Exact 1D grouping via PELT (penalized optimal partitioning)
//' @param mz numeric vector of m/z values (finite, > 0)
//' @param intensity numeric vector (same length as mz; values < 1 are clamped to 1)
//' @param beta penalty per additional segment; beta <= 0 -> auto (2 * log(n))
//' @param p0 floor mass accuracy in ppm (high-intensity limit / systematic floor)
//' @param Isat intensity at which the per-point precision plateaus at p0
//' @param max_seg maximum number of centroids allowed in one group; <= 0 = unlimited.
//'   Caps PELT's worst case (one giant group) at O(n * max_seg) instead of O(n^2).
//'   The result stays the exact global optimum over partitions obeying the cap.
//' @return integer vector of group ids (length(mz)), 1-indexed by min original position
//' @export
// [[Rcpp::export]]
IntegerVector fast_1D_pelt_grouping_cpp(NumericVector mz, NumericVector intensity,
                                        double beta = -1.0, double p0 = 5.0,
                                        double Isat = 5e4, int max_seg = 0) {
  const R_xlen_t nx = mz.size();
  if (intensity.size() != nx)
    stop("`mz` and `intensity` must have the same length");
  if (nx > (R_xlen_t)(std::numeric_limits<int>::max() - 2))
    stop("input too long for this implementation");
  if (std::isnan(beta)) stop("`beta` must not be NA/NaN");
  if (!(p0 > 0.0)   || !std::isfinite(p0))   stop("`p0` must be a finite value > 0");
  if (!(Isat > 0.0) || !std::isfinite(Isat)) stop("`Isat` must be a finite value > 0");
  if (max_seg < 0) max_seg = 0;                          // treat negative as unlimited

  const int n = (int)nx;
  if (n == 0) return IntegerVector(0);
  if (n == 1) {
    if (!std::isfinite(mz[0]) || !(mz[0] > 0.0)) stop("`mz` must be finite and > 0");
    return IntegerVector::create(1);
  }
  if (beta <= 0.0) beta = 2.0 * std::log((double)n);

  // 1. per-point weight from the shot-noise model, then sort by m/z
  std::vector<PeltPoint> pts(n);
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(mz[i]) || !(mz[i] > 0.0))
      stop("`mz` must be finite and > 0 (offending element: %d)", i + 1);
    if (std::isnan(intensity[i]))
      stop("`intensity` must not contain NA/NaN (offending element: %d)", i + 1);
    double I     = std::max(intensity[i], 1.0);
    double mult  = std::max(1.0, std::sqrt(Isat / I));
    double sigma = mz[i] * p0 * mult / 1e6;               // sigma in Da
    pts[i] = { mz[i], 1.0 / (sigma * sigma), i };
  }
  std::stable_sort(pts.begin(), pts.end(),
                   [](const PeltPoint& a, const PeltPoint& b){ return a.mz < b.mz; });

  // 2. PELT: F[t] = min_tau F[tau] + C(tau, t) + beta over surviving candidates
  std::vector<double> F(n + 1);
  std::vector<int>    cp(n + 1, 0);                       // argmin tau for backtracking
  F[0] = 0.0;              // beta counted once per segment (matches the objective above);
                           // the classic F[0] = -beta convention differs only by a constant
                           // offset and yields the identical optimal partition.
  std::vector<PeltCand> cands;
  cands.reserve(256);
  cands.push_back({ 0, F[0], 0.0, 0.0, 0.0, pts[0].mz, 0.0 });

  for (int t = 1; t <= n; ++t) {
    const double x = pts[t - 1].mz, w = pts[t - 1].w;
    double best = std::numeric_limits<double>::infinity();
    int    arg  = 0;
    for (PeltCand& c : cands) {                           // extend every open segment by point t-1
      const double d = x - c.xref;
      c.W += w * d; c.V += w; c.SS += w * d * d;
      double cost = c.SS - c.W * c.W / c.V;               // min_mu weighted SSE, shift-invariant
      if (cost < 0.0) cost = 0.0;                         // clamp float dust
      c.lastCost = cost;
      const double val = c.F + cost;
      if (val < best) { best = val; arg = c.tau; }
    }
    F[t] = best + beta; cp[t] = arg;

    // prune: tau can never win again if F[tau] + C(tau,t) > F[t]  (cost is subadditive)
    // plus the optional length cap: a segment that already spans max_seg points
    // cannot be extended further.
    cands.erase(std::remove_if(cands.begin(), cands.end(),
                  [&](const PeltCand& c){
                    return c.F + c.lastCost > F[t] ||
                           (max_seg > 0 && t - c.tau >= max_seg);
                  }),
                cands.end());
    if (t < n)                                            // open a segment starting at t
      cands.push_back({ t, F[t], 0.0, 0.0, 0.0, pts[t].mz, 0.0 });
  }

  // 3. backtrack segments; id = 1 + min original index within the group
  IntegerVector groups(n);
  for (int t = n; t > 0; ) {
    const int s = cp[t];                                  // segment covers sorted points [s, t-1]
    int minID = pts[s].original_idx + 1;
    for (int i = s + 1; i < t; ++i) minID = std::min(minID, pts[i].original_idx + 1);
    for (int i = s; i < t; ++i)     groups[pts[i].original_idx] = minID;
    t = s;
  }
  return groups;
}
