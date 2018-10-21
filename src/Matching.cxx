#include "Matching.h"
#include "HungarianSolver.h"
#include <algorithm>
#include <map>

#include <iostream>
#include <exception>

namespace {
  // helper functions for accessing the search record
  using SparseHungarian::idx_t;
  idx_t& getWithDefault(
      std::map<idx_t, idx_t>& map,
      idx_t index,
      idx_t defaultValue = 0)
  {
    auto itr = map.find(index);
    if (itr == map.end() )
      itr = map.insert(std::make_pair(index, defaultValue) ).first;
    return itr->second;
  }

  // Check the value stored under a key in the map.
  // If the value is missing then return false
  bool checkValueMissing(
      std::map<idx_t, idx_t>& map,
      idx_t index,
      idx_t comparison)
  {
    auto itr = map.find(index);
    if (itr == map.end() )
      return false;
    else
      return itr->second == comparison;
  }
}

namespace SparseHungarian {
  match_vec_t match(
      const cost_matrix_t& costs,
      float maxCost)
  {
    std::cout << "Begin matching: " << std::endl << costs << std::endl;
    // Not required to receive a square matrix, however it's much simpler if we
    // can assume that nRows <= nCols. Therefore if this isn't the case, just
    // flip the cost matrix
    if (costs.rows() > costs.cols() ) {
      match_vec_t matches;
      match_vec_t flippedMatch = match(costs.transpose(), maxCost);
      for (const match_t& match : flippedMatch)
        matches.push_back(std::make_pair(match.second, match.first) );
      return matches;
    }

    // First a quick reminder of notation - 'set A' is the smaller set
    // (represented by the rows) and 'set B' is the larger set (represented by
    // the columns).
    // I will sometimes use 'closest' to refer to the element in the other set
    // with the lowest cost

    std::cout << "Try the simple matching" << std::endl;
    // Start by attempting a very simple matching - just match every element in
    // A to the closest element in B
    match_vec_t matches;
    idx_t nMatchA(costs.rows() ); // number of objects being matched from A
    idx_t nMatchB(costs.cols() ); // number of objects being matched from B
    matches.reserve(nMatchA);
    bool valid = true; // Whether or not the simple match is valid
    std::set<idx_t> matchedIndices;
    for (idx_t ia = 0; ia < nMatchA; ++ia) {
      idx_t minIdx;
      if (costs.row(ia).minCoeff(&minIdx) <= maxCost) {
        // 'second' on the set insert operator tells you if the object was
        // actually inserted (i.e. wasn't already there)
        valid &= matchedIndices.insert(minIdx).second;
        if (!valid) // stop trying the instant a conflict is found
          break;
        matches.push_back(std::make_pair(ia, minIdx) );
      }
    }
    if (valid)
      return matches;

    std::cout << "Simple matching wasn't enough :(" << std::endl;
    HungarianSolver solver(costs, maxCost, matches);
    return solver.solution();
  }

  match_vec_t sparseMatch(
      const cost_matrix_t& cost,
      float maxCost)
  {
    auto groups = splitProblemIntoSparseGroups(cost, maxCost);
    return matchFromGroups(groups);
  }

  match_vec_t matchFromGroups(
      const std::vector<SparseGroup>& groups)
  {
    match_vec_t matches;
    for (const SparseGroup& group : groups) {
      match_vec_t groupMatch = match(group.costs, group.maxCost);
      for (const match_t& match : groupMatch) {
        matches.push_back(std::make_pair(
              group.indicesA.at(match.first),
              group.indicesB.at(match.second) ) );
      }
    }
    return matches;
  }
}
