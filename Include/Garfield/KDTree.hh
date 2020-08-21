#ifndef G_KDTREE2_H
#define G_KDTREE2_H

// (c) Matthew B. Kennel, Institute for Nonlinear Science, UCSD (2004)
//
// Licensed under the Academic Free License version 1.1 found in file LICENSE
// with additional provisions in that same file.


// Implement a kd tree for fast searching of points in a fixed data base
// in k-dimensional Euclidean space.

#include <vector>
#include <algorithm>

namespace Garfield {

typedef std::vector<std::vector<double> > KDTreeArray;

typedef struct {
  double lower, upper;
} interval;

class KDTreeNode; 
class SearchRecord;

struct KDTreeResult {
public:
  double dis;  // its square Euclidean distance
  int idx;    // which neighbor was found
}; 

// KDTree
//
// The main data structure, one for each k-d tree, pointing
// to a tree of an indeterminate number of "KDTreeNode"s.

class KDTree {
public: 
  const KDTreeArray& the_data;   
  // "the_data" is a reference to the underlying 
  // data to be included in the tree.
  //
  // NOTE: this structure does *NOT* own the storage underlying this.
  // Hence, it would be a very bad idea to change the underlying data
  // during use of the search facilities of this tree.
  // Also, the user must deallocate the memory underlying it.

  const int N;   // number of data points
  int dim;
  bool sort_results = false;

public:
  // Constructor.
  KDTree(KDTreeArray& data_in);
  // Destructor.
  ~KDTree();

public:
  // Search for n nearest to a given query vector 'qv'.
  void n_nearest(std::vector<double>& qv, int nn, 
                 std::vector<KDTreeResult>& result);

  // Search for 'nn' nearest to point [idxin] of the input data, excluding
  // neighbors within correltime 
  void n_nearest_around_point(int idxin, int correltime, int nn,
                              std::vector<KDTreeResult>& result);
  
  // Search for all neighbors in ball of size (square Euclidean distance)
  // r2.
  void r_nearest(std::vector<double>& qv, double r2,
                 std::vector<KDTreeResult>& result);

  // Like 'r_nearest', but around existing point, with decorrelation
  // interval. 
  void r_nearest_around_point(int idxin, int correltime, double r2,
                              std::vector<KDTreeResult>& result);

  friend class KDTreeNode;
  friend class SearchRecord;
private:
  KDTreeNode* root = nullptr;
  const KDTreeArray* data = nullptr;

  // Index for the tree leaves. Data in a leaf with bounds [l,u] are
  // in 'the_data[ind[l],*] to the_data[ind[u],*]
  std::vector<int> ind; 

  static const int bucketsize = 12; // global constant. 

private:
  KDTreeNode* build_tree_for_range(int l, int u, KDTreeNode* parent);
  int select_on_coordinate_value(int c, double alpha, int l, int u); 
  void spread_in_coordinate(int c, int l, int u, interval& interv);
};

/// A node in the tree.

class KDTreeNode {
public:
  /// Constructor
  KDTreeNode(int dim);
  /// Destructor
  ~KDTreeNode();

private:
  friend class KDTree;

  // Dimension to cut.
  int cut_dim;                                 
  // Cut value.
  double cut_val, cut_val_left, cut_val_right;  
  // Extents in index array for searching
  int l,u;  
  // [min,max] of the box enclosing all points
  std::vector<interval> box; 

  // Pointers to left and right nodes.
  KDTreeNode *left = nullptr;
  KDTreeNode *right = nullptr;  

  // Recursive innermost core routine for searching.
  void search(SearchRecord& sr); 

  // Return true if the bounding box for this node is within the
  // search range given by the searchvector and maximum ballsize in 'sr'. 
  bool box_in_search_range(SearchRecord& sr);

  // For processing final buckets. 
  void process_terminal_node(SearchRecord& sr);
  void process_terminal_node_fixedball(SearchRecord& sr);

};

}

#endif
