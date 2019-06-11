/*
 * preclustermain.cpp
 *
 *  Created on: 11-jun-2019
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include <lemon/arg_parser.h>
#include <fstream>
#include "precluster.h"

int main(int argc, char** argv)
{
  int nrSegments = 32;
  bool verbose = false;
  int clusterStatistic = 1;
  int nrThreads = -1;
  int timeLimit = -1;
  int memoryLimit = -1;
  double precisionBetaBin = -1;
  std::string stateTreeFilename;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("betaBin", "Beta binomial precision parameter (default: -1, binomial model)", precisionBetaBin)
    .refOption("C", "Clustering statistic:\n" \
                    "     0 -- Cancer Cell Fraction (CCF)\n" \
                    "     1 -- Descendant Cell Fraction (DCF, default)", clusterStatistic)
    .refOption("tl", "Time limit in seconds (default: -1, unlimited)", timeLimit)
    .refOption("t", "Number of threads (default: -1, use all available cores)", nrThreads)
    .refOption("N", "Number of segments (default: 32)", nrSegments, false)
    .refOption("v", "Verbose (default: 0)", verbose, false)
    .refOption("ML", "Memory limit", memoryLimit)
    .refOption("S", "State tree file", stateTreeFilename, false)
    .other("filename");
  ap.parse();
  
  if (!(0 <= clusterStatistic && clusterStatistic < 2))
  {
    std::cerr << "Error: invalid clustering statistic specified" << std::endl;
    return 1;
  }
  
  Solver::ClusterStatisticType statType = static_cast<Solver::ClusterStatisticType>(clusterStatistic);
  
  if (ap.files().empty())
  {
    std::cerr << "Error: missing input filename" << std::endl;
    return 1;
  }
  
  std::ifstream inR(ap.files()[0].c_str());
  if (!inR.good())
  {
    std::cerr << "Error: could not open '" << ap.files()[0] << "' for reading" << std::endl;
    return 1;
  }
  
  ReadMatrix R;
  try
  {
    inR >> R;
  }
  catch (std::runtime_error& e)
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  inR.close();
  
  if (!stateTreeFilename.empty())
  {
    std::ifstream inS(stateTreeFilename.c_str());
    if (!inS.good())
    {
      std::cerr << "Error: could not open '" << stateTreeFilename << "' for reading. Will generate state tree file." << std::endl;
      return 1;
    }
    StateGraph::readStateTrees(inS);
  }
  
  std::cerr << "Arguments: ";
  for (int i = 1; i < argc; ++i)
  {
    std::cerr << " " << argv[i];
  }
  std::cerr << std::endl;
  std::cerr << "Input:      instance with n = " << R.getNrCharacters() << " SNVs and m = " << R.getNrSamples() << " samples" << std::endl;
  
  PreCluster pc(R);
  pc.run(nrSegments, statType, precisionBetaBin, nrThreads, timeLimit, verbose, memoryLimit);
  
  IntMatrix preClustering = pc.getPreClustering();
  
  for (const IntVector& precluster : preClustering)
  {
    bool first = false;
    for (int i : precluster)
    {
      if (first)
        first = false;
      else
        std::cout << " ";
      std::cout << i;
    }
    std::cout << std::endl;
  }
  
  return 0;
}
