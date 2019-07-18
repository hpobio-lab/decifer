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
  
  std::uniform_real_distribution<int> dist(1, 1000000);
  for (int restart = 0; restart < nrRestarts; ++restart)
  {
    auto res = dkm::kmeans_lloyd(data, _k, dist(g_rng));
    
    double logLikelihood = 0;
    for (int i = 0; i < n; ++i)
    {
      for (int p = 0; p < m; ++p)
      {
        int j =  std::get<1>(res)[i];
        int var_pi = _R.getVar(p, i);
        int ref_pi = _R.getRef(p, i);
        double f_pi = std::get<0>(res)[j][p];
        logLikelihood += getLogLikelihood(var_pi, ref_pi, f_pi);
      }
    }
    
    
    
//    SoftPreClusterIlpCplex solver(_R, _k, log(_nrSegments) / log(2), _statType, _precisionBetaBin, true);
//    solver.init();
//    solver.initZ(z);
//    solver.solve(1, -1, verbose, -1);
//
    if (_logLikelihood < logLikelihood)
    {
      _logLikelihood = logLikelihood;
      _solZ = IntVector(n, -1);
      for (int i = 0; i < n; ++i)
      {
        _solZ[i] = std::get<1>(res)[i];
      }
      
      // set solD
      _solD = DoubleMatrix(_k, DoubleVector(m, 0));
      for (int j = 0; j < _k; ++j)
      {
        for (int p = 0; p < m; ++p)
        {
          _solD[j][p] = std::get<0>(res)[j][p];
        }
      }
      
      // set solPi
      _solPi = DoubleVector(_k, 0);
      for (int i = 0; i < n; ++i)
      {
        ++_solPi[_solZ[i]];
      }
      for (int j = 0; j < _k; ++j)
      {
        _solPi[j] /= n;
      }
    }
  }
  
  return true;
}
