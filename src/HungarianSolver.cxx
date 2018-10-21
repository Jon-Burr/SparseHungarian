#include "HungarianSolver.h"
#include <exception>
#include <set>
#include <queue>

#include <iostream>

namespace SparseHungarian {
  HungarianSolver::HungarianSolver(
      const cost_matrix_t& costs,
      float maxCost,
      const match_vec_t& initialMatching)
    : 
      nVtxA(costs.rows() ),
      nVtxB(costs.cols() ),
      m_costs(-costs),
      m_maxCost(-maxCost),
      m_labelsA(nVtxB, -maxCost),
      m_labelsB(nVtxB, 0.),
      m_matchA(nVtxB, nVtxB),
      m_matchB(nVtxB, nVtxB)
  {
    // Make sure that the input matrix is correct
    if (nVtxA > nVtxB)
      throw std::runtime_error("Invalid matrix supplied to HungarianSolver"
          "The matrix must have nRows <= nCols!");
    std::cout << "Square the matrix" << std::endl;
    // First square the matrix
    m_costs.conservativeResize(nVtxB, nVtxB);
    for (idx_t ia = nVtxA; ia < nVtxB; ++ia)
      for (idx_t ib = 0; ib < nVtxB; ++ib)
        m_costs(ia, ib) = m_maxCost;
    // Now load the initial matching
    for (const match_t& m : initialMatching) {
      m_matchA[m.first] = m.second;
      m_matchB[m.second] = m.first;
    }
    // Initialise the labels to sensible values
    for (idx_t ia = 0; ia < nVtxA; ++ia)
      m_labelsA[ia] = m_costs.row(ia).maxCoeff();
    solve();
    // Now load the solution into the internal vector
    std::cout << "Build solution" << std::endl;
    for (idx_t ia = 0; ia < nVtxA; ++ia) {
      std::cout << ia << std::endl;
      if (costs.coeff(ia, m_matchA[ia]) < maxCost)
        m_solution.push_back(std::make_pair(ia, m_matchA[ia]) );
    }
    std::cout << "Done" << std::endl;
  }

  void HungarianSolver::solve() 
  {
    std::cout << "In solve" << std::endl;
    // Begin by looking for an unmatched 'A' vertex
    // We can stop this search at nVtxA as we don't care about what the extra
    // dummy vertices match to
    while (true) {
      idx_t ia = 0;
      for (; ia < nVtxA; ++ia)
        if (m_matchA[ia] == nVtxB)
          break;
      if (ia == nVtxA)
        // This means that the matching is complete!
        return;
      breadthFirstSearch(ia);
    }
  }

  float HungarianSolver::getSlack(idx_t a, idx_t b) const
  {
    std::cout << "Get slack" << std::endl;
    std::cout << m_labelsA[a] << " + " << m_labelsB[b] << " - " << m_costs.coeff(a, b) << std::endl;
    return m_labelsA[a] + m_labelsB[b] - m_costs.coeff(a, b);
  }

  void HungarianSolver::breadthFirstSearch(idx_t root)
  {
    // This is a search for any unmatched 'B' node

    // This is really a breadth first search on a slightly modified graph.
    // Vertices which are matched to each other act as single vertices as far as
    // the search is concerned. Therefore we only need to keep track of the root
    // node and which 'B' nodes we have visited.
    std::set<idx_t> visited;

    // The path describes how to go *back* through the tree to the root node -
    // this is the only way we will traverse the tree so it's all we need
    std::map<idx_t, idx_t> path;

    // Keep track of the slacks on the edges between 'A' nodes in the equality
    // subgraph and 'B' nodes outside of it. The index of this vector is the 'B'
    // index
    std::vector<float> slacks;
    slacks.reserve(nVtxB);
    for (idx_t ib = 0; ib < nVtxB; ++ib)
      slacks.push_back(getSlack(root, ib) );
    // Also keep track of which 'A' index that corresponds to. This enables
    // skipping a few steps after updating the labels.
    std::vector<idx_t> minSlackIdx(nVtxB, root);
    // Keep track of the vertices queued up to be inspected
    std::queue<idx_t> vtxQueue;
    vtxQueue.push(root);

    while (true) {
      idx_t current = vtxQueue.front();
      std::cout << "Look at vtx: " << current << std::endl;
      vtxQueue.pop();
      // Find an edge on the equality subgraph leaving from this vertex
      for (idx_t ib = 0; ib < nVtxB; ++ib) {
        float slack = getSlack(current, ib);
        if (slack == 0) { // This is on the equality subgraph
          if (visited.count(ib) ) // but we've already visited
            continue;
          // This is an interesting vertex
          path[ib] = current;
          if (m_matchB[ib] == nVtxB) {
            // it's unmatched! That means we have an alternating augmenting path
            return augmentPath(path, ib);
          }
          else {
            // Add it to the queue and move on...
            visited.insert(ib);
            vtxQueue.push(m_matchB[ib]);
          }
        }
        else if (slack < slacks[ib]) {
          // Update the slacks
          slacks[ib] = slack;
          minSlackIdx[ib] = current;
        }
      }
      if (vtxQueue.size() == 0) {
        std::cout << "Update labels" << std::endl;
        // Being here means that we didn't find the alternating path
        // This means that we need better labelling
        // We find the minimum slack on a vertex heading out of the equality
        // subgraph
        float delta = std::numeric_limits<float>::infinity();
        idx_t minIdx = nVtxB;
        for (idx_t ib = 0; ib < nVtxB; ++ib) {
          if (visited.count(ib) )
            // We only visit nodes on the equality subgraph
            continue;
          if (slacks[ib] < delta) {
            delta = slacks[ib];
            minIdx = ib;
          }
        }
        std::cout << "delta: " << delta << " (" << minIdx << ")" << std::endl;
        // Now we update the labelling. Subtract delta from every 'A' vertex in
        // the equality subgraph and add it to every 'B' vertex in the equality
        // subgraph. This therefore has the effect of adding in a new vertex
        // into the subgraph, the one whose slack we just made 0! That vertex's
        // index is given by minIdx;
        path[minIdx] = minSlackIdx[minIdx];
        if (m_matchB[minIdx] == nVtxB) {
          // It's unmatched!
          return augmentPath(path, minIdx);
        }
        else {
          visited.insert(minIdx);
          vtxQueue.push(m_matchB[minIdx]);
          // And so we go on again :)
        }
      }
    }
  }

  void HungarianSolver::augmentPath(
      const std::map<idx_t, idx_t>& path,
      idx_t end)
  {
    std::cout << "Augment path!" << std::endl;
    do {
      std::cout << "GOOO" << std::endl;
      std::cout << end << std::endl;
      idx_t nextA = path.at(end);
      std::cout << nextA << std::endl;
      idx_t nextB = m_matchA[nextA];
      m_matchB[end] = path.at(end);
      m_matchA[nextA] = end;
      end = nextB;
    }
    // m_matchA[root] = nVtxB
    while (end != nVtxB);
  }
}
