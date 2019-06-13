/*
 * precluster.cpp
 *
 *  Created on: 6-jun-2019
 *      Author: M. El-Kebir
 */

#include "precluster.h"
#include "preclusterkmeans.h"

PreCluster::PreCluster(const ReadMatrix& R)
  : _R(R)
  , _nrClusters()
  , _preClustering(R.getNrCharacters(), -1)
{
}

void PreCluster::run(int nrSegments,
                     Solver::ClusterStatisticType statType,
                     double precisionBetaBin,
                     int nrThreads,
                     int timeLimit,
                     bool verbose,
                     int memoryLimit)
{
  IntMatrix partition = _R.partition();
  
  std::cerr << "Precluster: " << partition.size() << " groups of SNVs with identical copy numbers" << std::endl;
  int idx = 0;
  for (const IntVector& snvIndices : partition)
  {
    ++idx;
    std::cerr << "Precluster: clustering group " << idx << "/" << partition.size() << " with " << snvIndices.size() << " SNVs" << std::endl;
    run(snvIndices, nrSegments, statType, precisionBetaBin,
        nrThreads, timeLimit, verbose, memoryLimit);
  }
  
  _invPreClustering = IntMatrix(_nrClusters);
  for (int i = 0; i < _R.getNrCharacters(); ++i)
  {
    _invPreClustering[_preClustering[i]].push_back(i);
  }
}

void PreCluster::run(const IntVector& snvIndices,
                     int nrSegments,
                     Solver::ClusterStatisticType statType,
                     double precisionBetaBin,
                     int nrThreads,
                     int timeLimit,
                     bool verbose,
                     int memoryLimit)
{
  ReadMatrix R(_R.downSampleCharacters(snvIndices));
  
  const int nrObservations = R.getNumberOfObservations();
  
  double min_bic = std::numeric_limits<double>::max();
  IntVector solZ;
  int selected_k = -1;
  
  if (snvIndices.size() == 1)
  {
    selected_k = 1;
    solZ.push_back(0);
  }
  else
  {
    const size_t max_k = std::min((size_t)10, snvIndices.size());
    for (int k = 1; k <= max_k; ++k)
    {
      int nrParameters = R.getNrSamples() * (k+1);
      {
        Solver solver(R, k, 0, statType, precisionBetaBin, false);
        solver.init();
        
        for (int i = 0; i < R.getNrCharacters(); ++i)
        {
          nrParameters += k * solver.getScriptT(i).size() * k;
        }
      }
      
//      PreClusterIlpAlg solver(R, k, nrSegments, statType, precisionBetaBin, true);
      PreClusterKMeans solver(R, k, statType, precisionBetaBin);
      solver.init();
      solver.solve(1, verbose);
//      solver.solve(nrThreads, timeLimit, verbose, memoryLimit);
      
      double b = log(nrObservations) * nrParameters - 2 * solver.getLogLikelihood();
      std::cerr << "k = " << k << " -- Log likelihood " << solver.getLogLikelihood() << " -- BIC " << b << std::endl;
      
  //    if (solver.getLogLikelihood() == INFINITY)
  //    {
  //      continue;
  //    }
      
      if (b == INFINITY)
        continue;
      
      if (b < min_bic)
      {
        min_bic = b;
        solZ = solver.getSolZ();
        selected_k = k;
      }
      else
      {
        break;
      }
    }
  }
  
  for (int i = 0; i < snvIndices.size(); ++i)
  {
    int org_i = snvIndices[i];
    _preClustering[org_i] = _nrClusters + solZ[i];
  }
  _nrClusters += selected_k;
}
