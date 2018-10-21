#ifndef SparseHungarian_HungarianSolver_H
#define SparseHungarian_HungarianSolver_H

#include "Defs.h"
#include <vector>
#include <map>

namespace SparseHungarian{
  /**
   * \brief Class containing all of the information necessary to solve a
   * matching problem
   *
   * This class is the one that is actually used to solve the Hungarian
   * algorithm.
   */
  class HungarianSolver {
    public:
      /**
       * \brief Create the solver, this also performs the matching as part of
       * the constructor
       * \param costs The problem's cost matrix. Not required to be square but
       * the number of rows must be less than or equal to the number of columns.
       * \param maxCost If relevant, the maximum cost allowed to count as a
       * matching
       * \param initialMatching Any preliminary attempt at a matching. Supplying
       * this can speed up the algorithm
       */
      HungarianSolver(
          const cost_matrix_t& costs,
          float maxCost = std::numeric_limits<float>::infinity(),
          const match_vec_t& initialMatching = match_vec_t() );

      /// The number of vertices from set A
      const idx_t nVtxA;
      /// The number of vertices from set B
      const idx_t nVtxB;
      /// The solution to this problem
      const match_vec_t& solution() const { return m_solution; }
    private:
      /// The cost matrix, modified to be square
      cost_matrix_t m_costs;
      /// The maximum cost
      const float m_maxCost;
      /// The labels for set A
      std::vector<float> m_labelsA;
      /// The labels for set B
      std::vector<float> m_labelsB;
      /// The solution
      match_vec_t m_solution;
      /// Matches from A to B vertices
      std::vector<idx_t> m_matchA;
      /// Matches from B to A vertices
      std::vector<idx_t> m_matchB;
      /// Try to obtain a solution
      void solve();
      /// Get the slack on an edge
      float getSlack(idx_t a, idx_t b) const;
      /// Perform the depth-first search
      void breadthFirstSearch(idx_t root);
      /// Augment along the provided path
      void augmentPath(const std::map<idx_t, idx_t>& path, idx_t end);
  };
}
#endif //> !SparseHungarian_HungarianSolver_H
