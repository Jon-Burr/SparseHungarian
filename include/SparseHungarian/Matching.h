#ifndef SparseHungarian_Matching_H
#define SparseHungarian_Matching_H

#include "Defs.h"
#include "SparseGroup.h"
#include <limits>

namespace SparseHungarian {
  /**
   * \brief Perform a matching without the sparse implementation
   * \param costs The cost matrix defining the problem
   * \param maxCost The maximum cost for a match
   * \return A vector containing any matches that were found
   */
  match_vec_t match(
      const cost_matrix_t& costs,
      float maxCost = std::numeric_limits<float>::infinity() );

  /**
   * \brief Perform a matching using the sparse implementation
   * \param costs The cost matrix defining the problem
   * \param maxCost The maximum cost for a match
   * \return A vector containing any matches that were found
   */
  match_vec_t sparseMatch(
      const cost_matrix_t& costs,
      float maxCost);

  /**
   * \brief Build a match from a list of (disjoint) sparse groups
   * \param The input sparse groups
   */
  match_vec_t matchFromGroups(
      const std::vector<SparseGroup>& groups);

};

#endif //> SparseHungarian_Matching_H
