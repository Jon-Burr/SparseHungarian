#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <utility>

//#include "SparseHungarian/SparseGroup.h"
#include "SparseGroup.h"
#include "Matching.h"

namespace {
  float pythagoras(float a, float b) {
    return sqrt(a*a + b*b);
  }
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "No input file found!" << std::endl;
    return 1;
  }
  if (argc < 3) {
    std::cerr << "No output file found!" << std::endl;
    return 1;
  }

  std::ifstream ifs(argv[1]);
  if (!ifs.is_open() ) {
    std::cerr << "Failed to open input file: " << argv[1] << std::endl;
    return 1;
  }
  auto start = ifs.tellg();
  std::ofstream ofs(argv[2]);
  if (!ofs.is_open() ) {
    std::cerr << "Failed to open output file: " << argv[2] << std::endl;
    return 1;
  }
  // Copy over input file
  ofs << ifs.rdbuf();
  ifs.seekg(start);
  std::string word;
  int nPoints;
  float maxCost;
  ifs >> word >> nPoints;
  if (word != "NPoints:")
    throw std::runtime_error(word + "received, 'NPoints:' expected!");
  ifs >> word >> maxCost;
  if (word != "MaxDR:")
    throw std::runtime_error(word + "received, 'MaxDR:' expected!");
  ifs >> word;
  if (word != "SetA:")
    throw std::runtime_error(word + "received, 'SetA:' expected!");
  std::vector<std::pair<float, float> > setA;
  std::vector<std::pair<float, float> > setB;
  for (unsigned int ii = 0; ii < nPoints; ++ii) {
    setA.emplace_back();
    ifs >> setA.back().first >> setA.back().second;
  }
  ifs >> word;
  if (word != "SetB:")
    throw std::runtime_error(word + "received, 'SetB:' expected!");
  for (unsigned int ii = 0; ii < nPoints; ++ii) {
    setB.emplace_back();
    ifs >> setB.back().first >> setB.back().second;
  }

  //std::cout << "Set A:" << std::endl;
  //for (const auto& p : setA)
    //std::cout << p.first << ", " << p.second << std::endl;
  //std::cout << "Set B:" << std::endl;
  //for (const auto& p : setB)
    //std::cout << p.first << ", " << p.second << std::endl;

  std::cout << "Building cost matrix" << std::endl;
  // Let's build the cost matrix
  SparseHungarian::cost_matrix_t costs(setA.size(), setB.size() );
  for (unsigned int ia = 0; ia < setA.size(); ++ia) {
    for (unsigned int ib = 0; ib < setB.size(); ++ib) {
      const auto& pa = setA.at(ia);
      const auto& pb = setB.at(ib);
      costs(ia, ib) = pythagoras(pa.first - pb.first, pa.second - pb.second);
      //std::cout << "(" << pa.first << ", " << pa.second << ") - (" << pb.first << ", " << pb.second << ") = ";
      //std::cout << costs(ia, ib) << std::endl;
    }
  }

  std::cout << "Building groups" << std::endl;

  // Now let's build the groups
  std::vector<SparseHungarian::SparseGroup> groups = 
    SparseHungarian::splitProblemIntoSparseGroups(costs, maxCost);

  ofs << "Groups:" << std::endl;
  //std::cout << "Groups: " << std::endl;
  for (unsigned int ii = 0; ii < groups.size(); ++ii) {
    ofs << "\tGroup " << ii << std::endl;
    ofs << "\tSetA" << std::endl;
    for (auto ia : groups.at(ii).indicesA)
      ofs << "\t\t" << ia << std::endl;
    ofs << "\tSetB" << std::endl;
    for (auto ib : groups.at(ii).indicesB)
      ofs << "\t\t" << ib << std::endl;
  }

  SparseHungarian::match_vec_t matches = SparseHungarian::matchFromGroups(groups);
  //std::cout << costs << std::endl;
  ofs << "Matches:" << std::endl;
  for (const SparseHungarian::match_t& m : matches)
    ofs << m.first << " : " << m.second << std::endl;
    //std::cout << m.first << ": " << m.second << " (" << costs(m.first, m.second) << ")" << std::endl;



  return 0;
}
