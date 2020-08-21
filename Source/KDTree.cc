//
// (c) Matthew B. Kennel, Institute for Nonlinear Science, UCSD (2004)
//
// Licensed under the Academic Free License version 1.1 found in file LICENSE
// with additional provisions in that same file.

#include <algorithm> 
#include <limits>
#include <queue>
#include <iostream>

#include "Garfield/KDTree.hh"

namespace {

double squared(const double x) { return x * x; }

double dis_from_bnd(const double x, const double amin, const double amax) {
  if (x > amax) {
    return(x-amax); 
  } else if (x < amin)
    return (amin-x);
  else
    return 0.0;
}

}

namespace Garfield {

inline bool operator<(const KDTreeResult& e1, const KDTreeResult& e2) {
  return (e1.dis < e2.dis);
}

// Constructor
KDTree::KDTree(KDTreeArray& data_in)
  : the_data(data_in) {

  const size_t n = data_in.size(); 
  if (!data_in.empty()) {
    dim = data_in[0].size();
  } 

  ind.resize(n);
  for (size_t i = 0; i < n; i++) ind[i] = i; 
  // Build the tree.
  root = build_tree_for_range(0, n - 1, 0); 
  data = &the_data;
}

// Destructor
KDTree::~KDTree() {
  delete root;
}

KDTreeNode* KDTree::build_tree_for_range(int l, int u, KDTreeNode* parent) {

  if (u < l) return nullptr;
  KDTreeNode* node = new KDTreeNode(dim);
  if ((u - l) <= bucketsize) {
    // Create a terminal node. 
    // Always compute true bounding box for terminal node. 
    for (int i = 0; i < dim; i++) {
      node->box[i] = spread_in_coordinate(i, l, u);
    }
    node->cut_dim = 0; 
    node->cut_val = 0.0;
    node->l = l;
    node->u = u;
    node->left = node->right = nullptr;
  } else {
    // Compute an APPROXIMATE bounding box for this node.
    // If parent == nullptr, then this is the root node, and 
    // we compute for all dimensions.
    // Otherwise, we copy the bounding box from the parent for
    // all coordinates except for the parent's cut dimension. 
    // That, we recompute ourself.
    int c = -1;
    double maxspread = 0.0;
    for (int i = 0; i < dim; i++) {
      if (!parent || (parent->cut_dim == i)) {
        node->box[i] = spread_in_coordinate(i, l, u);
      } else {
        node->box[i] = parent->box[i];
      }
      double spread = node->box[i][1] - node->box[i][0]; 
      if (spread > maxspread) {
        maxspread = spread;
        c = i; 
      }
    }

    // Now, c is the identity of which coordinate has the greatest spread.
    double sum = 0.0;
    for (int k = l; k <= u; k++) {
      sum += the_data[ind[k]][c];
    }
    double average = sum / static_cast<double>(u - l + 1);
    int m = select_on_coordinate_value(c, average, l, u);

    // Move the indices around to cut on dim 'c'.
    node->cut_dim = c;
    node->l = l;
    node->u = u;
    node->left = build_tree_for_range(l, m, node);
    node->right = build_tree_for_range(m + 1, u, node);

    if (!node->right) {
      for (int i = 0; i < dim; i++) { 
        node->box[i] = node->left->box[i];
      } 
      node->cut_val = node->left->box[c][1];
      node->cut_val_left = node->cut_val_right = node->cut_val;
    } else if (!node->left) {
      for (int i = 0; i < dim; i++) { 
        node->box[i] = node->right->box[i];
      } 
      node->cut_val = node->right->box[c][1];
      node->cut_val_left = node->cut_val_right = node->cut_val;
    } else {
      node->cut_val_right = node->right->box[c][0];
      node->cut_val_left  = node->left->box[c][1];
      node->cut_val = 0.5 * (node->cut_val_left + node->cut_val_right); 
      
      // Now recompute true bounding box as union of subtree boxes.
      // This is now faster having built the tree, being logarithmic in
      // N, not linear as would be from naive method.
      for (int i = 0; i < dim; i++) {
        node->box[i][1] = std::max(node->left->box[i][1],
                                      node->right->box[i][1]);
        
        node->box[i][0] = std::min(node->left->box[i][0],
                                      node->right->box[i][0]);
      }
    }
  }
  return node;
}

std::array<double, 2> KDTree::spread_in_coordinate(int c, int l, int u) {
  // Return the minimum and maximum of the indexed data between l and u.

  double smin = the_data[ind[l]][c];
  double smax = smin;

  // Process two at a time.
  int i; 
  for (i = l + 2; i <= u; i += 2) {
    double lmin = the_data[ind[i - 1]][c];
    double lmax = the_data[ind[i]][c];
    if (lmin > lmax) std::swap(lmin, lmax); 
    if (smin > lmin) smin = lmin;
    if (smax < lmax) smax = lmax;
  }
  // is there one more element? 
  if (i == u + 1) {
    double last = the_data[ind[u]][c];
    if (smin > last) smin = last;
    if (smax < last) smax = last;
  }
  return {smin, smax};
}

int KDTree::select_on_coordinate_value(int c, double alpha, int l, int u) {
  //  Move indices in ind[l..u] so that the elements in [l .. return]
  //  are <= alpha, and hence are less than the [return + 1 .. u]
  //  elements, viewed across dimension 'c'.
  int lb = l, ub = u;
  while (lb < ub) {
    if (the_data[ind[lb]][c] <= alpha) {
      lb++; // good where it is.
    } else {
      std::swap(ind[lb],ind[ub]); 
      ub--;
    }
  }

  // here ub=lb
  if (the_data[ind[lb]][c] <= alpha)
    return(lb);
  else
    return(lb-1);
  
}

// search record substructure
//
// one of these is created for each search.
// this holds useful information  to be used
// during the search

class SearchRecord {

private:
  friend class KDTree;
  friend class KDTreeNode;

  std::vector<double>& qv; 
  int dim;
  unsigned int nn = 0;
  double ballsize;
  int centeridx, correltime;

  std::priority_queue<KDTreeResult>& result;
  const KDTreeArray* data; 
  const std::vector<int>& ind; 

public:
  SearchRecord(std::vector<double>& qv_in, KDTree& tree_in,
               std::priority_queue<KDTreeResult>& result_in) :
    qv(qv_in),
    result(result_in),
    data(tree_in.data),
    ind(tree_in.ind) { 
    dim = tree_in.dim;
    ballsize = std::numeric_limits<double>::max();
  };

};

// search for n nearest to a given query vector 'qv'.
void KDTree::n_nearest(std::vector<double>& qv, int nn, 
                       std::vector<KDTreeResult>& result) {

  std::priority_queue<KDTreeResult> res; 
  SearchRecord sr(qv, *this, res);
  sr.centeridx = -1;
  sr.correltime = 0;
  sr.nn = nn; 
  root->search(sr); 
  result.clear();
  while (!res.empty()) {
    result.push_back(res.top());
    res.pop();
  }
  if (sort_results) sort(result.begin(), result.end());
  
}

void KDTree::n_nearest_around_point(int idxin, int correltime, int nn,
                                    std::vector<KDTreeResult>& result) {

  std::vector<double> qv(dim);
  for (int i = 0; i < dim; i++) {
    qv[i] = the_data[idxin][i]; 
  }
  std::priority_queue<KDTreeResult> res;
  SearchRecord sr(qv, *this, res);
  sr.centeridx = idxin;
  sr.correltime = correltime;
  sr.nn = nn; 
  root->search(sr); 
  result.clear(); 
  while (!res.empty()) {
    result.push_back(res.top());
    res.pop();
  }
  if (sort_results) sort(result.begin(), result.end());
}

// search for all within a ball of a certain radius
void KDTree::r_nearest(std::vector<double>& qv, double r2, 
                       std::vector<KDTreeResult>& result) {
  std::priority_queue<KDTreeResult> res;
  SearchRecord sr(qv, *this, res);
  sr.centeridx = -1;
  sr.correltime = 0;
  sr.nn = 0; 
  sr.ballsize = r2; 
  root->search(sr); 
  result.clear(); 
  while (!res.empty()) {
    result.push_back(res.top());
    res.pop();
  }
  if (sort_results) sort(result.begin(), result.end());
} 

void KDTree::r_nearest_around_point(int idxin, int correltime, double r2,
                                    std::vector<KDTreeResult>& result) {
  std::vector<double> qv(dim);
  for (int i = 0; i < dim; i++) {
    qv[i] = the_data[idxin][i]; 
  }

  std::priority_queue<KDTreeResult> res;
  SearchRecord sr(qv, *this, res);
  sr.centeridx = idxin;
  sr.correltime = correltime;
  sr.ballsize = r2; 
  sr.nn = 0; 
  root->search(sr); 
  result.clear(); 
  while (!res.empty()) {
    result.push_back(res.top());
    res.pop();
  }
  if (sort_results) sort(result.begin(), result.end());
}

// Constructor
KDTreeNode::KDTreeNode(int dim) : box(dim) {} 

// Destructor
KDTreeNode::~KDTreeNode() {
  if (left) delete left; 
  if (right) delete right; 
}

void KDTreeNode::search(SearchRecord& sr) {
  // the core search routine.
  // This uses true distance to bounding box as the
  // criterion to search the secondary node. 

  if (!left && !right) {
    // We are on a terminal node
    if (sr.nn == 0) {
      process_terminal_node_fixedball(sr);
    } else {
      process_terminal_node(sr);
    }
    return;
  }
  KDTreeNode *ncloser = nullptr;
  KDTreeNode *nfarther = nullptr;

  double extra;
  double qval = sr.qv[cut_dim]; 
  // value of the wall boundary on the cut dimension. 
  if (qval < cut_val) {
    ncloser = left;
    nfarther = right;
    extra = cut_val_right - qval;
  } else {
    ncloser = right;
    nfarther = left;
    extra = qval - cut_val_left; 
  }

  if (ncloser) ncloser->search(sr);
  if ((nfarther) && (squared(extra) < sr.ballsize)) {
    // first cut
    if (nfarther->box_in_search_range(sr)) {
      nfarther->search(sr); 
    }      
  }
}

inline bool KDTreeNode::box_in_search_range(SearchRecord& sr) {

  // does the bounding box, represented by minbox[*],maxbox[*]
  // have any point which is within 'sr.ballsize' to 'sr.qv'??
 
  int dim = sr.dim;
  double dis2 = 0.0; 
  double ballsize = sr.ballsize; 
  for (int i = 0; i < dim; i++) {
    dis2 += squared(dis_from_bnd(sr.qv[i], box[i][0], box[i][1]));
    if (dis2 > ballsize) return false;
  }
  return true;
}

void KDTreeNode::process_terminal_node(SearchRecord& sr) {
  int centeridx = sr.centeridx;
  int correltime = sr.correltime;
  unsigned int nn = sr.nn; 
  int dim = sr.dim;
  double ballsize = sr.ballsize;
 
  const KDTreeArray& data = *sr.data;

  for (int i = l; i <= u; i++) {
    int indexofi = sr.ind[i];
    bool early_exit = false;
    double dis = 0.0;
    for (int k = 0; k < dim; k++) {
      dis += squared(data[indexofi][k] - sr.qv[k]);
      if (dis > ballsize) {
        early_exit = true; 
        break;
      }
    }
    if (early_exit) continue; // next iteration of mainloop
    
    if (centeridx > 0) {
      // Skip points within the decorrelation interval.
      if (abs(indexofi - centeridx) < correltime) continue;
    }

    // Add the point to the list.
    if (sr.result.size() < nn) {
      // The list so far is undersized. 
      KDTreeResult e;
      e.idx = indexofi;
      e.dis = dis;
      sr.result.push(e); 
      // Set the ball radius to the largest on the list (maximum priority).
      if (sr.result.size() == nn) ballsize = sr.result.top().dis;
    } else {
      // if we get here then the current node, has a squared 
      // distance smaller
      // than the last on the list, and belongs on the list.
      KDTreeResult e;
      e.idx = indexofi;
      e.dis = dis;
      sr.result.pop();
      sr.result.push(e); 
      ballsize = sr.result.top().dis;
    }
  } // main loop
  sr.ballsize = ballsize;
}

void KDTreeNode::process_terminal_node_fixedball(SearchRecord& sr) {
  int centeridx = sr.centeridx;
  int correltime = sr.correltime;
  int dim = sr.dim;
  double ballsize = sr.ballsize;

  const KDTreeArray& data = *sr.data;

  for (int i = l; i <= u; i++) {
    int indexofi = sr.ind[i]; 
    bool early_exit = false;
    double dis = 0.0;
    for (int k = 0; k < dim; k++) {
      dis += squared(data[indexofi][k] - sr.qv[k]);
      if (dis > ballsize) {
        early_exit= true; 
        break;
      }
    }
    if (early_exit) continue; // next iteration of mainloop
    
    if (centeridx > 0) {
      // Skip points within the decorrelation interval.
      if (abs(indexofi - centeridx) < correltime) continue;
    }

    KDTreeResult e;
    e.idx = indexofi;
    e.dis = dis;
    // sr.result.push_back(e);
    sr.result.push(e);
  }
}

}
