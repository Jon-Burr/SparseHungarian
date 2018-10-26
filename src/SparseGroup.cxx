#include "SparseHungarian/SparseGroup.h"
#include <memory>
#include <algorithm>
#include <queue>

namespace SparseHungarian {
  void SparseGroup::buildCosts(
      const cost_matrix_t& fullCosts,
      float maxCost)
  {
    this->maxCost = maxCost;
    costs.resize(indicesA.size(), indicesB.size() );
    for (idx_t ia = 0; ia < indicesA.size(); ++ia)
      for (idx_t ib = 0; ib < indicesB.size(); ++ib)
        costs(ia, ib) = fullCosts(indicesA[ia], indicesB[ib]);
  }

  std::vector<SparseGroup> splitProblemIntoSparseGroups(
      const Eigen::MatrixXf& costs,
      float maxCost)
  {
    idx_t nVtxA = costs.rows();
    idx_t nVtxB = costs.cols();
    std::vector<SparseGroup> groups;
    // This is essentially a graph partioning problem. Use a breadth first
    // search
    // Keep track of which vertices we've visited
    std::vector<bool> visitedA(nVtxA, false);
    std::vector<bool> visitedB(nVtxB, false);
    // Record the queue for the BFS. As it must contain vertices from both A and
    // B store the B nodes as idx + nVtxA
    std::queue<idx_t> vtxQueue;

    // Step through the A vertices to seed the groups
    idx_t nextVtx = 0;
    while (nextVtx != nVtxA) {
      if (visitedA[nextVtx]) {
        ++nextVtx;
        continue;
      }
      visitedA[nextVtx] = true;
      vtxQueue.push(nextVtx);
      groups.emplace_back();
      SparseGroup& group = groups.back();
      group.indicesA.push_back(nextVtx);
      while (vtxQueue.size() != 0) {
        idx_t current = vtxQueue.front();
        vtxQueue.pop();
        if (current >= nVtxA) {
          // this is a 'B' vertex
          current -= nVtxA;
          // Start from nextVtx+1, we've already visited all the 'A' vertices
          // before this
          for (idx_t ia = nextVtx + 1; ia < nVtxA; ++ia) {
            if (visitedA[ia] || costs.coeff(ia, current) > maxCost)
              continue;
            vtxQueue.push(ia);
            group.indicesA.push_back(ia);
            visitedA[ia] = true;
          }
        }
        else {
          // This is an 'A' vertex
          for (idx_t ib = 0; ib < nVtxB; ++ib) {
            if (visitedB[ib] || costs.coeff(current, ib) > maxCost)
              continue;
            vtxQueue.push(ib+nVtxA);
            group.indicesB.push_back(ib);
            visitedB[ib] = true;
          }
        }
      }
      // Reaching here means that we've completed the group
      if (group.indicesB.size() == 0)
        // We only want to return groups that contain some matchings
        groups.pop_back();
      else
        group.buildCosts(costs, maxCost);
      // Walk on the A index
      ++nextVtx;
    }
    return groups;
  }
}
