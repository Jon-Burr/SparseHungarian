#include "json.hpp"
#include "SparseHungarian/SparseGroup.h"
#include "SparseHungarian/Matching.h"
#include "boost/program_options.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>

namespace {
  const float pi = 3.14159265358979323846;
  float pythagoras(float a, float b) {
    return sqrt(a*a + b*b);
  }
}

int main(int argc, char* argv[]) {
  namespace po = boost::program_options;
  using json = nlohmann::json;

  // The input options
  std::string inputFileName;
  std::string outputFileName;
  po::options_description opts("Allowed options");
  opts.add_options()
    ("help,h", "Produce this message and exit.")
    ("input,i", po::value(&inputFileName), "The input file to read from")
    ("output,o", po::value(&outputFileName), "The output file to write to. "
     "If not set, write to the input file");

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(opts).run(), vm);
  po::notify(vm);

  // Check command line inputs
  if (vm.count("help") ) {
    std::cout << opts << std::endl;
    return 0;
  }

  std::cout << vm.count("input") << std::endl;

  if (inputFileName.empty() ) {
    std::cerr << "No input file name provided!" << std::endl;
    return 1;
  }

  if (outputFileName.empty() )
    outputFileName = inputFileName;

  // Read in the input file
  std::ifstream ifs(inputFileName);
  if (!ifs.is_open() ) {
    std::cerr << "Failed to open input file: " << inputFileName << std::endl;
    return 1;
  }
  json j;
  ifs >> j;
  ifs.close();

  // Make sure we can open the output file
  std::ofstream ofs(outputFileName);
  if (!ofs.is_open() ) {
    std::cerr << "Failed to open output file: " << outputFileName << std::endl;
    return 1;
  }

  auto pointsA = j["PointsA"].get<std::vector<std::pair<float, float>>>();
  auto pointsB = j["PointsB"].get<std::vector<std::pair<float, float>>>();

  float maxCost = j["MaxDR"].get<float>();

  // Now, make the cost matrix
  SparseHungarian::cost_matrix_t costs(pointsA.size(), pointsB.size() );
  for (unsigned int ia = 0; ia < pointsA.size(); ++ia) {
    for (unsigned int ib = 0; ib < pointsB.size(); ++ib) {
      const auto& pa = pointsA.at(ia);
      const auto& pb = pointsB.at(ib);
      double phiDiff = std::fmod(fabs(pa.first - pb.first), 2*pi);
      costs(ia, ib) = pythagoras(
          std::min(phiDiff, fabs(2*pi - phiDiff) ),
          pa.second - pb.second
          );
      // fmod wrt pi to account for detector wrapping
    }
  }

  //std::cout << costs << std::endl;

  auto sparseStart = std::chrono::system_clock::now();
  // Then build the groups
  auto groups = SparseHungarian::splitProblemIntoSparseGroups(costs, maxCost);
  
  // Now solve the problem
  auto sparseMatches = SparseHungarian::matchFromGroups(groups);
  auto sparseEnd = std::chrono::system_clock::now();
  std::chrono::duration<double> sparseDuration = sparseEnd - sparseStart;
  std::cout << "Sparse Hungarian took " << sparseDuration.count() << " seconds." << std::endl;

  // Add a little - solve with the original Hungarian algorithm
  auto origStart = std::chrono::system_clock::now();
  auto origMatches = SparseHungarian::match(costs, maxCost);
  auto origEnd = std::chrono::system_clock::now();
  std::chrono::duration<double> origDuration = origEnd - origStart;
  std::cout << "Original Hungarian took " << origDuration.count() << " seconds." << std::endl;

  // Now, write everything out
  std::vector<std::pair<std::vector<long>, std::vector<long>>>
    groupsOut;
  for (const auto& g : groups) {
    groupsOut.push_back(std::make_pair(g.indicesA, g.indicesB) );
    //groupsOut.emplace_back();
    //groupsOut.back().first.insert(g.indicesA.begin(), g.indicesA.end() );
    //groupsOut.back().second.insert(g.indicesB.begin(), g.indicesB.end() );
  }
  std::ostringstream oss;
  oss << costs;
  j["CostsString"] = oss.str();
  j["Groups"] = groupsOut;
  j["SparseMatches"] = sparseMatches;
  j["HungarianMatches"] = origMatches;

  // Write out the output
  ofs << std::setw(4) << j;
  return 0;
}

