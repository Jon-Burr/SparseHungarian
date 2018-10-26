#ifndef SparseHungarian_SparseGroup_H
#define SparseHungarian_SparseGroup_H

#include <vector>
#include <set>
//#include "SparseHungarian/Defs.h"
#include "Defs.h"

namespace SparseHungarian {
  /**
   * \brief Describes a subset of the original problem
   *
   * This class is used to split an input problem into sparse subsets. The costs
   * between elements of these subsets are greater than the max cost so they
   * cannot affect each other during the matching process therefore it is safe
   * to match them individually, therefore decreasing the runtime of the
   * matching process.
   */
  class SparseGroup {
    public:
      SparseGroup() {}
      /// The subset of indices from set A in the group
      std::vector<idx_t> indicesA;
      /// The subset of indices from set B in the group
      std::vector<idx_t> indicesB;
      /// The costs between the members of this group
      cost_matrix_t costs;
      /// The maximum cost in this problem
      float maxCost;
      /**
       * \brief build the cost matrix for this group
       * \param fullCosts The cost matrix for the full matching problem
       * \param maxCost The maximum cost in the problem
       */
      void buildCosts(
         const cost_matrix_t& costs,
         float maxCost);

  };

  /**
   * \brief Split a problem into SparseGroups
   * \param costs The costs for this matching problem
   * \param maxCost The maximum cost in this matching problem
   */
  std::vector<SparseGroup> splitProblemIntoSparseGroups(
      const cost_matrix_t& costs,
      float maxCost);
}

#endif //> !SparseHungarian_SparseGroup_H
