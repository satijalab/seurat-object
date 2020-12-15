#include <RcppEigen.h>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include <Rinternals.h>

using namespace Rcpp;

Eigen::SparseMatrix<double> RowMergeMatricesList(
    List mat_list,
    List mat_rownames,
    std::vector<std::string> all_rownames
);

template <typename S>
std::vector<size_t> sort_indexes(const std::vector<S> &v);
List GraphToNeighborHelper(Eigen::SparseMatrix<double> mat);
