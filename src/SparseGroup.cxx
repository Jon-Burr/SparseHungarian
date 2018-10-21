//#include "SparseHungarian/SparseGroup.h"
//#include "../SparseHungarian/SparseGroup.h"
#include "SparseGroup.h"
#include <memory>
#include <iostream>
#include <set>
#include <algorithm>
#include <list>

namespace {
  // check for the intersection of two sets
  template <typename T>
    bool set_intersection(const std::set<T>& set1, const std::set<T>& set2)
    {
      auto itr1 = set1.begin();
      auto itr2 = set2.begin();
      while(itr1 != set1.end() && itr2 != set2.end() )
      {
        if (*itr1 < *itr2)
          ++itr1;
        else if (*itr2 < *itr1)
          ++itr2;
        else
          return true;
      }
      return false;
    }
}

namespace SparseHungarian {

  SparseGroup::SparseGroup(
      const std::set<idx_t>& indicesA,
      const std::set<idx_t>& indicesB,
      const cost_matrix_t& fullCosts,
      float maxCost)
    : indicesA(indicesA.begin(), indicesA.end() ),
      indicesB(indicesB.begin(), indicesB.end() ),
      costs(indicesA.size(), indicesB.size() ),
      maxCost(maxCost)
  {
    for (unsigned int ia = 0; ia < this->indicesA.size(); ++ia) {
      for (unsigned int ib = 0; ib < this->indicesB.size(); ++ib) {
        costs(ia, ib) =
          fullCosts(this->indicesA.at(ia), this->indicesB.at(ib) );
      }
    }
  }

  std::vector<SparseGroup> splitProblemIntoSparseGroups(
      const Eigen::MatrixXf& costs,
      float maxCost)
  {
    // b node associated to each a node
    struct GroupIndices {
      std::set<idx_t> indicesA;
      std::set<idx_t> indicesB;
    };
    // There's going to be a lot of iterator movement and deletion. List
    // *should* be more efficient for this and less likely to invalidate the
    // iterators...
    std::list<GroupIndices> grpIndices;
    for (idx_t ia = 0; ia < costs.rows(); ++ia) {
      grpIndices.emplace_back();
      grpIndices.back().indicesA.insert(ia);
      for (idx_t ib = 0; ib < costs.cols(); ++ib)
        if (costs(ia, ib) <= maxCost)
          grpIndices.back().indicesB.insert(ib);
      if (grpIndices.back().indicesB.size() == 0)
        grpIndices.pop_back();
    }
    // Begin with one group for each 'A' nodes.
    // All grpIndices have at least one of each type of node.
    // Any solo nodes (ones that cannot match to anything) have already been
    // excluded
    auto itr1 = grpIndices.begin();
    ++itr1;
    while (itr1 != grpIndices.end() ) {
      auto itr2 = grpIndices.begin();
      while (itr2 != itr1) {
        // If there is any intersection between the two grpIndices 'A' or 'B' nodes
        // they have to be combined. However, by construction we are increasing
        // the 'A' node index monotonically so there will never be such an
        // intersection
        if (set_intersection(itr2->indicesB, itr1->indicesB) ) {
          itr1->indicesB.insert(itr2->indicesB.begin(), itr2->indicesB.end() );
          itr1->indicesA.insert(itr2->indicesA.begin(), itr2->indicesA.end() );
          itr2 = grpIndices.erase(itr2); // advances the iterator
        }
        else 
          ++itr2;
      }
      // Hit the end of the list, advance itr1 and return itr2 to the beginning
      ++itr1;
    }
    std::vector<SparseGroup> groups;
    for (const auto& indices : grpIndices) {
      groups.emplace_back(indices.indicesA, indices.indicesB, costs, maxCost);
    }
    return groups;
  }
}
