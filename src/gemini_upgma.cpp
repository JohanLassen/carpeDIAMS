#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <queue>
#include <stack>
#include <cmath>

using namespace Rcpp;

struct Node {
  int start, end;     
  int left_child = -1;
  int right_child = -1;
  double mean;
  long long count;
  double span; // We still need span for your maxDiff grouping logic
};

struct MergeCandidate {
  double mean_dist; 
  int left_node_idx;    
  
  bool operator>(const MergeCandidate& other) const {
    if (std::abs(mean_dist - other.mean_dist) < 1e-15) {
      return left_node_idx > other.left_node_idx; 
    }
    return mean_dist > other.mean_dist;
  }
};

struct Point {
  double val;
  int original_idx;
};

//' Fast 1D UPGMA grouping
//' @param mz numeric vector
//' @param maxDiff numeric scalar
//' @return integer vector of group ids
//' @export
// [[Rcpp::export]]
IntegerVector fast_1D_upgma_grouping_cpp(NumericVector x, double maxDiff) {
  int n = x.size();
  if (n == 0) return IntegerVector(0);
  if (n == 1) return IntegerVector::create(1);
  
  // 1. Sort points
  std::vector<Point> pts(n);
  for (int i = 0; i < n; ++i) pts[i] = {x[i], i};
  std::stable_sort(pts.begin(), pts.end(), [](const Point& a, const Point& b) {
    return a.val < b.val;
  });
  
  // 2. Build UPGMA Tree
  std::vector<Node> nodes;
  nodes.reserve(2 * n); 
  std::vector<int> prev_vec(2 * n, -1);
  std::vector<int> next_vec(2 * n, -1);
  std::vector<bool> active(2 * n, true);
  std::priority_queue<MergeCandidate, std::vector<MergeCandidate>, std::greater<MergeCandidate>> pq;
  
  for (int i = 0; i < n; ++i) {
    nodes.push_back({i, i, -1, -1, pts[i].val, 1, 0.0});
    if (i > 0) prev_vec[i] = i - 1;
    if (i < n - 1) {
      next_vec[i] = i + 1;
      pq.push({std::abs(pts[i+1].val - pts[i].val), i});
    }
  }
  
  int root_idx = 0;
  while (!pq.empty()) {
    MergeCandidate top = pq.top();
    pq.pop();
    
    int l = top.left_node_idx;
    int r = next_vec[l];
    
    if (r == -1 || !active[l] || !active[r]) continue;
    
    // In UPGMA, distance is distance between means
    double current_mean_dist = std::abs(nodes[r].mean - nodes[l].mean);
    if (std::abs(top.mean_dist - current_mean_dist) > 1e-14) continue;
    
    int newNodeIdx = nodes.size();
    long long newCount = nodes[l].count + nodes[r].count;
    // Calculate new mean: (n1*m1 + n2*m2) / (n1 + n2)
    double newMean = (nodes[l].mean * nodes[l].count + nodes[r].mean * nodes[r].count) / newCount;
    
    nodes.push_back({
      nodes[l].start, 
      nodes[r].end, 
      l, r, 
      newMean, 
      newCount, 
      pts[nodes[r].end].val - pts[nodes[l].start].val
    });
    
    active[l] = false; active[r] = false;
    root_idx = newNodeIdx;
    
    int p = prev_vec[l];
    int nn = next_vec[r];
    prev_vec[newNodeIdx] = p; next_vec[newNodeIdx] = nn;
    
    if (p != -1) {
      next_vec[p] = newNodeIdx;
      pq.push({std::abs(nodes[newNodeIdx].mean - nodes[p].mean), p});
    }
    if (nn != -1) {
      prev_vec[nn] = newNodeIdx;
      pq.push({std::abs(nodes[nn].mean - nodes[newNodeIdx].mean), newNodeIdx});
    }
  }
  
  // 3. Top-Down Grouping (Matches your original 'cluster_span' logic)
  IntegerVector groups(n);
  std::stack<int> s;
  s.push(root_idx);
  
  while (!s.empty()) {
    int curr = s.top();
    s.pop();
    
    Node& node = nodes[curr];
    // We still use span (max - min) to decide where to cut, 
    // because that matches your 'groupSimilarityMatrixTree1D' function.
    if (node.span <= maxDiff) {
      int minID = pts[node.start].original_idx + 1;
      for (int i = node.start + 1; i <= node.end; ++i) {
        if (pts[i].original_idx + 1 < minID) minID = pts[i].original_idx + 1;
      }
      for (int i = node.start; i <= node.end; ++i) {
        groups[pts[i].original_idx] = minID;
      }
    } else {
      if (node.right_child != -1) s.push(node.right_child);
      if (node.left_child != -1) s.push(node.left_child);
    }
  }
  return groups;
}