#ifndef SparseHungarian_Defs_H
#define SparseHungarian_Defs_H

#include <Eigen/Dense>
#include <vector>
#include <utility>

namespace SparseHungarian {
  using cost_matrix_t = Eigen::MatrixXf;
  using cost_row_t = Eigen::RowVectorXf;
  using cost_col_t = Eigen::VectorXf;
  using idx_t = Eigen::Index;
  using match_t = std::pair<idx_t, idx_t>;
  using match_vec_t = std::vector<match_t>;
}

#endif //> !SparseHungarian_Defs_H
