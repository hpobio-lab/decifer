/*
 * preclusterkmeans.cpp
 *
 *  Created on: 13-jun-2019
 *      Author: M. El-Kebir
 */

#include "preclusterkmeans.h"
#include "softpreclusterilpcplex.h"
#include "dkm/dkm.hpp"
#include <tuple>

PreClusterKMeans::PreClusterKMeans(const ReadMatrix& R,
                                   int k,
                                   int nrSegments,
                                   ClusterStatisticType statType,
                                   double precisionBetaBin)
  : Solver(R, k, nrSegments, statType, precisionBetaBin, false)
{
}

void PreClusterKMeans::init()
{
  Solver::init();
}

bool PreClusterKMeans::solve(int nrRestarts,
                             bool verbose)
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  if (m > 20)
  {
    throw std::runtime_error("Error: m > 20 in kmeans clustering");
  }
  
  std::vector<std::array<double, 20> > data;
  for (int i = 0; i < n; ++i)
  {
    data.push_back(std::array<double, 20>{{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}});
    for (int p = 0; p < m; ++p)
    {
      data.back()[p] = _R.getVAF(p, i);
    }
  }
  
  for (int restart = 0; restart < nrRestarts; ++restart)
  {
    std::vector<uint32_t> assignment = std::get<1>(dkm::kmeans_lloyd(data, _k, restart));
    
    IntVector z(_R.getNrCharacters());
    for (int i = 0; i < _R.getNrCharacters(); ++i)
    {
      z[i] = assignment[i];
    }
    
    SoftPreClusterIlpCplex solver(_R, _k, log(_nrSegments) / log(2), _statType, _precisionBetaBin, false);
    solver.init();
    solver.initZ(z);
    solver.solve(1, -1, verbose, -1);
    
    if (_logLikelihood < solver.getLogLikelihood())
    {
      _logLikelihood = solver.getLogLikelihood();
      _solZ = solver.getSolZ();
    }
  }
  
  return true;
}
