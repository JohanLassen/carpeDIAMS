// #include <Rcpp.h>
// #include <algorithm>
// #include <vector>
// #include <queue>
// #include <stack>
// #include <cmath>
// using namespace Rcpp;

// // Intensity-adaptive agglomerative grouping (IAAC) of 1-D centroids, O(n log n).
// //
// // Adjacency-constrained agglomeration (priority queue + doubly-linked list).
// // Linkage: inverse-variance consensus z-distance between adjacent clusters,
// // merged smallest-first. Cut: a node is a group iff every internal merge happened
// // at z <= k (the `max_z` field). Each cluster keeps two running sums so the
// // consensus position and uncertainty are O(1) to maintain on merge:
// //     W = sum_i mz_i / sigma_i^2      ->  m_hat   = W / V
// //     V = sum_i 1     / sigma_i^2      ->  sigma_hat = 1 / sqrt(V)
// //
// // NOTE (the sigma-hat floor): a cluster's consensus cannot be more accurate than
// // the instrument's systematic floor p0. Without the floor, sigma_hat shrinks
// // below the real m/z scatter of strong peaks and `max_z` fragments them. The
// // floor `sigma_hat >= m_hat * p0 / 1e6` fixes this; it only affects well-populated
// // strong clusters and leaves singleton / weak behaviour untouched.

// struct Node {
//   int start, end;                 // contiguous range in sorted order
//   int left_child = -1, right_child = -1;
//   double W;                       // sum_i  mz_i / sigma_i^2
//   double V;                       // sum_i  1    / sigma_i^2
//   double max_z;                   // max consensus-merge z anywhere in this subtree
// };

// struct MergeCandidate {
//   double z; int left_node_idx;
//   bool operator>(const MergeCandidate& o) const {
//     if (std::abs(z - o.z) < 1e-15) return left_node_idx > o.left_node_idx;
//     return z > o.z;
//   }
// };

// struct Point { double mz, sigma; int original_idx; };

// // consensus SE, floored at the instrument systematic accuracy p0 (ppm)
// static inline double sigma_hat(const Node& nd, double p0) {
//   double mhat = nd.W / nd.V;
//   double s    = 1.0 / std::sqrt(nd.V);     // inverse-variance consensus SE
//   double flr  = mhat * p0 / 1e6;           // cannot beat p0 ppm
//   return s > flr ? s : flr;
// }
// static inline double consensus_z(const Node& a, const Node& b, double p0) {
//   double sa = sigma_hat(a, p0), sb = sigma_hat(b, p0);
//   return std::abs(a.W / a.V - b.W / b.V) / std::sqrt(sa * sa + sb * sb);
// }

// //' Fast 1D intensity-adaptive grouping (IAAC), O(n log n)
// //' @param mz numeric vector of m/z values
// //' @param intensity numeric vector (same length as mz)
// //' @param k z-threshold; a group is kept iff every internal merge is <= k sigma
// //' @param p0 floor mass accuracy in ppm (high-intensity limit / systematic floor)
// //' @param Isat intensity at which the per-point precision plateaus at p0
// //' @return integer vector of group ids (length(mz)), 1-indexed by min original position
// //' @export
// // [[Rcpp::export]]
// IntegerVector fast_1D_iaac_grouping_cpp(NumericVector mz, NumericVector intensity,
//                                         double k, double p0 = 5.0, double Isat = 5e4) {
//   int n = mz.size();
//   if (n == 0) return IntegerVector(0);
//   if (n == 1) return IntegerVector::create(1);

//   // 1. per-point sigma from the shot-noise model, then sort by m/z
//   std::vector<Point> pts(n);
//   for (int i = 0; i < n; ++i) {
//     double I    = std::max(intensity[i], 1.0);
//     double mult = std::max(1.0, std::sqrt(Isat / I));      // floor at p0
//     pts[i] = { mz[i], mz[i] * p0 * mult / 1e6, i };        // sigma in Da
//   }
//   std::stable_sort(pts.begin(), pts.end(),
//                    [](const Point& a, const Point& b){ return a.mz < b.mz; });

//   // 2. adjacency-constrained agglomeration, smallest consensus-z first
//   std::vector<Node> N; N.reserve(2 * n);
//   std::vector<int> prev_vec(2 * n, -1), next_vec(2 * n, -1);
//   std::vector<char> active(2 * n, 1);
//   std::priority_queue<MergeCandidate, std::vector<MergeCandidate>, std::greater<MergeCandidate>> pq;

//   for (int i = 0; i < n; ++i) {
//     double v = 1.0 / (pts[i].sigma * pts[i].sigma);
//     N.push_back({ i, i, -1, -1, pts[i].mz * v, v, 0.0 });
//     if (i > 0)     prev_vec[i] = i - 1;
//     if (i < n - 1) next_vec[i] = i + 1;
//   }
//   for (int i = 0; i < n - 1; ++i) pq.push({ consensus_z(N[i], N[i + 1], p0), i });

//   int root = 0;
//   while (!pq.empty()) {
//     MergeCandidate top = pq.top(); pq.pop();
//     int l = top.left_node_idx, r = next_vec[l];
//     if (r == -1 || !active[l] || !active[r]) continue;
//     double cur = consensus_z(N[l], N[r], p0);
//     if (std::abs(top.z - cur) > 1e-9 * std::max(1.0, cur)) continue;   // stale: a neighbour changed

//     int idx = (int)N.size();
//     double W = N[l].W + N[r].W, V = N[l].V + N[r].V;
//     double mzn = std::max(cur, std::max(N[l].max_z, N[r].max_z));
//     N.push_back({ N[l].start, N[r].end, l, r, W, V, mzn });
//     active[l] = active[r] = 0; root = idx;

//     int p = prev_vec[l], nn = next_vec[r];
//     prev_vec[idx] = p; next_vec[idx] = nn;
//     if (p  != -1) { next_vec[p]  = idx; pq.push({ consensus_z(N[p], N[idx], p0), p   }); }
//     if (nn != -1) { prev_vec[nn] = idx; pq.push({ consensus_z(N[idx], N[nn], p0), idx }); }
//   }

//   // 3. top-down cut: emit a node as a group iff every internal merge was <= k sigma
//   IntegerVector groups(n);
//   std::stack<int> st; st.push(root);
//   while (!st.empty()) {
//     int cur = st.top(); st.pop();
//     const Node& nd = N[cur];
//     if (nd.max_z <= k) {
//       int minID = pts[nd.start].original_idx + 1;
//       for (int i = nd.start + 1; i <= nd.end; ++i) minID = std::min(minID, pts[i].original_idx + 1);
//       for (int i = nd.start; i <= nd.end; ++i) groups[pts[i].original_idx] = minID;
//     } else {
//       if (nd.right_child != -1) st.push(nd.right_child);
//       if (nd.left_child  != -1) st.push(nd.left_child);
//     }
//   }
//   return groups;
// }

#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <queue>
#include <stack>
#include <cmath>
#include <limits>
using namespace Rcpp;

// Intensity-adaptive agglomerative grouping (IAAC) of 1-D centroids, O(n log n).
//
// Adjacency-constrained agglomeration (priority queue + doubly-linked list).
// Linkage: inverse-variance consensus z-distance between adjacent clusters,
// merged smallest-first. Cut: a node is a group iff every internal merge happened
// at z <= k (the `max_z` field). Each cluster keeps two running sums so the
// consensus position and uncertainty are O(1) to maintain on merge:
//     W = sum_i mz_i / sigma_i^2      ->  m_hat     = W / V
//     V = sum_i 1    / sigma_i^2      ->  sigma_hat = 1 / sqrt(V)
//
// NOTE (the sigma-hat floor): a cluster's consensus cannot be more accurate than
// the instrument's systematic floor p0. Without the floor, sigma_hat shrinks
// below the real m/z scatter of strong peaks and `max_z` fragments them. The
// floor `sigma_hat >= m_hat * p0 / 1e6` fixes this; it only affects well-populated
// strong clusters and leaves singleton / weak behaviour untouched.
//
// Correctness notes on the merge loop:
//   * Nodes are immutable once created (merges make NEW nodes), so a queued
//     candidate's z never goes stale in value; a candidate is invalid only if
//     one of its two nodes has since been consumed by another merge. Storing
//     BOTH node indices makes that check exact (`active` + adjacency), with no
//     floating-point tolerance needed.
//   * Every current adjacency always has an exact live entry in the queue
//     (pushed at the moment the adjacency was created), so the queue always
//     yields the true global minimum over current adjacencies, and the loop
//     performs exactly n-1 merges. Total queue traffic is <= 3n entries.
namespace {
struct Node {
  int start, end;                 // contiguous range in sorted order
  int left_child, right_child;    // -1 for leaves
  double W;                       // sum_i  mz_i / sigma_i^2
  double V;                       // sum_i  1    / sigma_i^2
  double max_z;                   // max consensus-merge z anywhere in this subtree
};

struct MergeCandidate {
  double z;
  int left, right;                // node indices of the adjacent pair when pushed
  // Exact strict weak ordering (no epsilon: a fuzzy comparator is not a valid
  // ordering and is UB in std::priority_queue). Ties break on the node indices
  // so results are fully deterministic.
  bool operator>(const MergeCandidate& o) const {
    if (z != o.z)         return z > o.z;
    if (left != o.left)   return left > o.left;
    return right > o.right;
  }
};

struct Point { double mz, sigma; int original_idx; };

// consensus SE, floored at the instrument systematic accuracy p0 (ppm)
static inline double sigma_hat(const Node& nd, double p0) {
  double mhat = nd.W / nd.V;
  double s    = 1.0 / std::sqrt(nd.V);     // inverse-variance consensus SE
  double flr  = mhat * p0 / 1e6;           // cannot beat p0 ppm
  return s > flr ? s : flr;
}
static inline double consensus_z(const Node& a, const Node& b, double p0) {
  double sa = sigma_hat(a, p0), sb = sigma_hat(b, p0);
  return std::abs(a.W / a.V - b.W / b.V) / std::sqrt(sa * sa + sb * sb);
}
} // anonymous namespace

//' Fast 1D intensity-adaptive grouping (IAAC), O(n log n)
//' @param mz numeric vector of m/z values (finite, > 0)
//' @param intensity numeric vector (same length as mz; values < 1 are clamped to 1)
//' @param k z-threshold; a group is kept iff every internal merge is <= k sigma
//'   (k < 0 yields all singletons; +Inf yields one group; NA is an error)
//' @param p0 floor mass accuracy in ppm (high-intensity limit / systematic floor)
//' @param Isat intensity at which the per-point precision plateaus at p0
//' @return integer vector of group ids (length(mz)), 1-indexed by min original position
//' @export
// [[Rcpp::export]]
IntegerVector fast_1D_iaac_grouping_cpp(NumericVector mz, NumericVector intensity,
                                        double k, double p0 = 5.0, double Isat = 5e4) {
  const R_xlen_t nx = mz.size();
  if (intensity.size() != nx)
    stop("`mz` and `intensity` must have the same length");
  if (nx > (R_xlen_t)(std::numeric_limits<int>::max() / 2 - 2))
    stop("input too long for this implementation");
  if (std::isnan(k))   stop("`k` must not be NA/NaN");
  if (!(p0 > 0.0)   || !std::isfinite(p0))   stop("`p0` must be a finite value > 0");
  if (!(Isat > 0.0) || !std::isfinite(Isat)) stop("`Isat` must be a finite value > 0");

  const int n = (int)nx;
  if (n == 0) return IntegerVector(0);
  if (n == 1) {
    if (!std::isfinite(mz[0]) || !(mz[0] > 0.0)) stop("`mz` must be finite and > 0");
    return IntegerVector::create(1);
  }

  // 1. per-point sigma from the shot-noise model, then sort by m/z
  std::vector<Point> pts(n);
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(mz[i]) || !(mz[i] > 0.0))
      stop("`mz` must be finite and > 0 (offending element: %d)", i + 1);
    if (std::isnan(intensity[i]))
      stop("`intensity` must not contain NA/NaN (offending element: %d)", i + 1);
    double I    = std::max(intensity[i], 1.0);
    double mult = std::max(1.0, std::sqrt(Isat / I));   // sigma floors at mz*p0/1e6
    pts[i] = { mz[i], mz[i] * p0 * mult / 1e6, i };     // sigma in Da
  }
  std::stable_sort(pts.begin(), pts.end(),
                   [](const Point& a, const Point& b){ return a.mz < b.mz; });

  // 2. adjacency-constrained agglomeration, smallest consensus-z first
  std::vector<Node> N; N.reserve(2 * n - 1);            // n leaves + (n-1) merges
  std::vector<int> prev_vec(2 * n - 1, -1), next_vec(2 * n - 1, -1);
  std::vector<char> active(2 * n - 1, 1);

  std::vector<MergeCandidate> pq_buf; pq_buf.reserve(3 * (size_t)n);
  std::priority_queue<MergeCandidate, std::vector<MergeCandidate>,
                      std::greater<MergeCandidate>> pq(std::greater<MergeCandidate>(),
                                                       std::move(pq_buf));

  for (int i = 0; i < n; ++i) {
    double v = 1.0 / (pts[i].sigma * pts[i].sigma);
    N.push_back({ i, i, -1, -1, pts[i].mz * v, v, 0.0 });
    if (i > 0)     prev_vec[i] = i - 1;
    if (i < n - 1) next_vec[i] = i + 1;
  }
  for (int i = 0; i < n - 1; ++i)
    pq.push({ consensus_z(N[i], N[i + 1], p0), i, i + 1 });

  int root = 0;
  while (!pq.empty()) {
    const MergeCandidate top = pq.top(); pq.pop();
    const int l = top.left, r = top.right;
    if (!active[l] || !active[r] || next_vec[l] != r) continue;   // pair consumed

    const int idx = (int)N.size();
    const double W = N[l].W + N[r].W, V = N[l].V + N[r].V;
    const double mzn = std::max(top.z, std::max(N[l].max_z, N[r].max_z));
    N.push_back({ N[l].start, N[r].end, l, r, W, V, mzn });
    active[l] = active[r] = 0; root = idx;

    const int p = prev_vec[l], nn = next_vec[r];
    prev_vec[idx] = p; next_vec[idx] = nn;
    if (p  != -1) { next_vec[p]  = idx; pq.push({ consensus_z(N[p], N[idx], p0), p,   idx }); }
    if (nn != -1) { prev_vec[nn] = idx; pq.push({ consensus_z(N[idx], N[nn], p0), idx, nn  }); }
  }

  // 3. top-down cut: emit a node as a group iff every internal merge was <= k sigma
  //    (a leaf is always a group, so every point gets an id even when k < 0)
  IntegerVector groups(n);
  std::stack<int> st; st.push(root);
  while (!st.empty()) {
    const int cur = st.top(); st.pop();
    const Node& nd = N[cur];
    if (nd.max_z <= k || nd.left_child == -1) {
      int minID = pts[nd.start].original_idx + 1;
      for (int i = nd.start + 1; i <= nd.end; ++i)
        minID = std::min(minID, pts[i].original_idx + 1);
      for (int i = nd.start; i <= nd.end; ++i)
        groups[pts[i].original_idx] = minID;
    } else {
      st.push(nd.right_child);
      st.push(nd.left_child);
    }
  }
  return groups;
}