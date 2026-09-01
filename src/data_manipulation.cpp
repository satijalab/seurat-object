#include <RcppEigen.h>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <fstream>
#include <string>
#include <Rinternals.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::Triplet<double> T;

template <typename S>
std::vector<size_t> sort_indexes(const std::vector<S> &v) {
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(),
                   [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  return idx;
}

// [[Rcpp::export(rng = false)
// [[Rcpp::export(rng = false)]]
List GraphToNeighborHelper(Eigen::SparseMatrix<double> mat) {
  mat = mat.transpose();
  //determine the number of neighbors
  int n = 0;
  for(Eigen::SparseMatrix<double>::InnerIterator it(mat, 0); it; ++it) {
    n += 1;
  }
  Eigen::MatrixXd nn_idx(mat.rows(), n);
  Eigen::MatrixXd nn_dist(mat.rows(), n);

  for (int k=0; k<mat.outerSize(); ++k){
    int n_k = 0;
    std::vector<double> row_idx;
    std::vector<double> row_dist;
    row_idx.reserve(n);
    row_dist.reserve(n);
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it) {
      if (n_k > (n-1)) {
        Rcpp::stop("Not all cells have an equal number of neighbors.");
      }
      row_idx.push_back(it.row() + 1);
      row_dist.push_back(it.value());
      n_k += 1;
    }
    if (n_k != n) {
      Rcpp::Rcout << n << ":::" << n_k << std::endl;
      Rcpp::stop("Not all cells have an equal number of neighbors.");
    }
    //order the idx based on dist
    std::vector<size_t> idx_order = sort_indexes(row_dist);
    for(int i = 0; i < n; ++i) {
      nn_idx(k, i) = row_idx[idx_order[i]];
      nn_dist(k, i) = row_dist[idx_order[i]];
    }
  }
  List neighbors = List::create(nn_idx, nn_dist);
  return(neighbors);
}

// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> RowMergeMatricesList(
  List mat_list,
  List mat_rownames,
  std::vector<std::string> all_rownames
) {
  // Convert Rcpp lists to c++ vectors
  std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor>> mat_vec;
  mat_vec.reserve(mat_list.size());
  std::vector<std::vector<std::string>> rownames_vec;
  rownames_vec.reserve(mat_rownames.size());
  std::vector<std::unordered_map<std::string, int>> map_vec;
  map_vec.reserve(mat_list.size());
  int num_cols = 0;
  // Accumulate in 64 bits: summing the non-zero counts of several large layers
  // overflows a signed int, and the negative result then converts to a huge
  // size_t in reserve() below, which throws std::length_error ("vector::reserve")
  int64_t num_nZero = 0;
  // offsets keep track of which column to add in to
  std::vector<int> offsets;
  for (unsigned int i = 0; i < mat_list.size(); i++) {
    mat_vec.emplace_back(Rcpp::as<Eigen::SparseMatrix<double, Eigen::RowMajor>>(mat_list.at(i)));
    rownames_vec.push_back(mat_rownames[i]);
    // Set up hash maps for rowname based lookup
    std::unordered_map<std::string, int> mat_map;
    for (unsigned int j = 0; j < rownames_vec[i].size(); j++) {
      mat_map[rownames_vec[i][j]] = j;
    }
    map_vec.emplace_back(mat_map);
    offsets.push_back(num_cols);
    num_cols += mat_vec[i].cols();
    num_nZero += static_cast<int64_t>(mat_vec[i].nonZeros());
  }
  // set up tripletList for new matrix creation
  std::vector<T> tripletList;
  int num_rows = all_rownames.size();
  // This routine returns a dgCMatrix, whose row indices and column pointers are
  // 32-bit signed integers, so it cannot index more non-zeros than this. Other
  // layer classes are stitched by their own StitchMatrix methods and are not
  // subject to this limit, so point at them rather than failing inside the
  // standard library.
  const int64_t max_nZero = 2147483647LL;
  if (num_nZero > max_nZero) {
    Rcpp::stop(
      "Cannot combine these layers: the result would have " +
      std::to_string(num_nZero) + " non-zero values, more than the " +
      std::to_string(max_nZero) + " that a dgCMatrix can index. Convert the "
      "layers to a matrix class that is not bound by that limit, then join "
      "them again: Seurat's as.DelayedMatrix() keeps them in memory, and "
      "BPCells keeps them on disk. Subsetting to fewer cells or features, or "
      "leaving the layers split, also avoids it."
    );
  }
  tripletList.reserve(static_cast<size_t>(num_nZero));
  // loop over all rows and add nonzero entries to tripletList
  for(int i = 0; i < num_rows; i++) {
    std::string key = all_rownames[i];
    for(int j = 0; j < mat_vec.size(); j++) {
      if (map_vec[j].count(key)) {
        for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it1(mat_vec[j], map_vec[j][key]); it1; ++it1){
          tripletList.emplace_back(i, it1.col() + offsets[j], it1.value());
        }
      }
    }
  }
  Eigen::SparseMatrix<double> combined_mat(num_rows, num_cols);
  combined_mat.setFromTriplets(tripletList.begin(), tripletList.end());
  return combined_mat;
}
