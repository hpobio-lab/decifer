/*
 * empre.cpp
 *
 *  Created on: 3-jul-2019
 *      Author: M. El-Kebir
 */

#include "empre.h"
#include "stategraph.h"
#include "dkm/dkm.hpp"
#include <tuple>
#include <iomanip>
#include <lemon/time_measure.h>
#include <algorithm>

EMPre::EMPre(const ReadMatrix& R,
             const IntMatrix& preClustering,
             int k,
             int nrSegments,
             ClusterStatisticType statType,
             double precisionBetaBin,
             bool forceTruncal)
  : Solver(R, k, nrSegments, statType, precisionBetaBin, forceTruncal)
  , _preClustering(preClustering)
  , _solY(R.getNrCharacters())
  , _solT()
  , _initY()
{
}

void EMPre::updatePWLA()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  _dOverallLB = DoubleMatrix(_k, DoubleVector(_R.getNrSamples(), 0));
  _dOverallUB = DoubleMatrix(_k, DoubleVector(_R.getNrSamples(), 1));
  for (int j = 0; j < _k; ++j)
  {
    for (int i = 0; i < n; ++i)
    {
      for (int t = 0; t < _scriptT[i].size(); ++t)
      {
        if (isEnabled(i, t, j))
        {
          for (int p = 0; p < m; ++p)
          {
            if (g_tol.less(_dOverallLB[j][p], _dLB[i][t][p]))
            {
              _dOverallLB[j][p] = _dLB[i][t][p];
            }
            if (g_tol.less(_dUB[i][t][p], _dOverallUB[j][p]))
            {
              _dOverallUB[j][p] = _dUB[i][t][p];
            }
            
            if (!g_tol.different(_dOverallLB[j][p], _dOverallUB[j][p]))
            {
              _dOverallLB[j][p] = _dOverallUB[j][p];
            }
            
            _dOverallLB[j][p] = _dOverallLB[j][p] <= g_thre ? 0 : _dOverallLB[j][p];
            _dOverallUB[j][p] = _dOverallUB[j][p] <= g_thre ? 0 : _dOverallUB[j][p];
            
            assert(_dOverallLB[j][p] <= _dOverallUB[j][p]);
          }
        }
      }
    }
  }
}

void EMPre::initPWLA()
{
  Solver::initPWLA();
  
//  updatePWLA();
}

bool EMPre::initializeD(int seed,
                        int nrThreads,
                        int timeLimit,
                        bool verbose,
                        int nrDownSampledSNVs,
                        int memoryLimit)
{
  const int m = _R.getNrSamples();
  const int n = _R.getNrCharacters();
  const int nrPreClusters = _preClustering.size();
  
  IntVector charToPreCluster(n, -1);
  for (int ii = 0; ii < nrPreClusters; ++ii)
  {
    for (int i : _preClustering[ii])
    {
      charToPreCluster[i] = ii;
    }
  }
  
  IntVector snvIndices;
  for (int ii = 0; ii < nrPreClusters; ++ii)
  {
    std::uniform_int_distribution<> unif(0, _preClustering[ii].size()-1);
    int idx = unif(g_rng);
    snvIndices.push_back(idx);
  }
  
//  IntVector snvIndices(n, 0);
//  for (int i = 1; i < n; ++i)
//  {
//    snvIndices[i] = snvIndices[i-1] + 1;
//  }
//  std::shuffle(snvIndices.begin(), snvIndices.end(), g_rng);
//  snvIndices.erase(snvIndices.begin() + nrDownSampledSNVs, snvIndices.end());
  
  ReadMatrix R = _R.downSampleCharacters(snvIndices);
  
  IntMatrix newPreClustering;
  std::map<int, int> old2new;
  std::map<int, int> old2newSNV;
  int idx = 0;
  for (int i : snvIndices)
  {
    int ii = charToPreCluster[i];
    if (old2new.count(ii) == 0)
    {
      old2new[ii] = newPreClustering.size();
      newPreClustering.push_back(IntVector(1, idx));
    }
    else
    {
      newPreClustering[old2new[ii]].push_back(idx);
    }
    old2newSNV[i] = idx;
    ++idx;
  }
  
  {
    {
      std::unique_ptr<SoftClusterIlp> pILP = createSoftClusterIlpSolver(R, newPreClustering);
//      std::unique_ptr<SoftClusterIlp> pILP = createSoftClusterIlpSolver(_R, _preClustering);
      pILP->init();
      
      for (int i : snvIndices)
      {
        int ii = charToPreCluster[i];
        int new_i = old2newSNV[i];
        
        for (int t = 0; t < _scriptT[i].size(); ++t)
        {
          for (int p = 0; p < m; ++p)
          {
            double dLB = 0;
            double dUB = 1;
            for (int i2 : _preClustering[ii])
            {
              dLB = std::max(dLB, _dLB[i2][t][p]);
              dUB = std::min(dUB, _dUB[i2][t][p]);
            }
            
            pILP->setDLB(new_i, t, p, dLB);
            pILP->setDUB(new_i, t, p, dUB);
          }
        }
      }
      
      if (!pILP->solve(nrThreads, timeLimit, verbose, memoryLimit))
      {
        return false;
      }
      else
      {
        for (int j = 0; j < _k; ++j)
        {
          for (int p = 0; p < m; ++p)
          {
            _solD[j][p] = pILP->getD(p, j);
          }
        }
        
//        _solPi = DoubleVector(_k, 0);
//        for (int i : snvIndices)
//        {
//          int ii = charToPreCluster[i];
//          int new_i = old2newSNV[i];
//          
//          for (int i2 : _preClustering[ii])
//          {
//            for (int t = 0; t < _scriptT[i].size(); ++t)
//            {
//              for (int j = 0; j < _k; ++j)
//              {
//                _solY[i2][t][j] = pILP->getSolY()[new_i][t][j];
//                _solPi[j] += _solY[i2][t][j];
//              }
//            }
//          }
//        }
//        
//        for (int j = 0; j < _k; ++j)
//        {
//          _solPi[j] /= n;
//        }
//        
//        roundD();
       
        pILP = NULL;
      }
    }
  }
  
  return true;
}

bool EMPre::solve(int restart,
                  int seed,
                  int maxIterations,
                  int nrThreads,
                  int timeLimit,
                  bool verbose,
                  int nrDownSampledSNVs,
                  int memoryLimit)
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  lemon::Timer timer;
  initializeD(seed, nrThreads, timeLimit, verbose, nrDownSampledSNVs, memoryLimit);
  
  updateLogLikelihood();
  std::cerr << "k == " << _k << " -- Restart = " << restart << " -- Iteration = 0" << " -- log likelihood = " << _logLikelihood << " -- " << timer.realTime() << " s" << std::endl;
  
  _preClustering = IntMatrix(n);
  for (int i = 0; i < n; ++i)
  {
    _preClustering[i].push_back(i);
  }
  
  initPWLA();
  
  int step = 0;
  while (true)
  {
    if (!stepE(nrThreads, verbose))
    {
      std::cerr << "k == " << _k << " -- Restart = " << restart << " -- Iteration = " << ++step << " -- infeasible solution! (E)" << " -- " << timer.realTime() << " s" << std::endl;
      return false;
    }
    
    double delta = updateLogLikelihood();
    std::cerr << "k == " << _k << " -- Restart = " << restart << " -- Iteration = " << ++step
              << " -- log likelihood = " << _logLikelihood << " -- "
//              << " -- gamma log likelihood = " << getLogLikelihoodGamma() << " -- "
              <<  timer.realTime() << " s (E)" << std::endl;
    
    if (!stepM(nrThreads, verbose))
    {
      std::cerr << "k == " << _k << " -- Restart = " << restart << " -- Iteration = " << ++step << " -- infeasible solution! (M)" << " -- " << timer.realTime() << " s" << std::endl;
      return false;
    }
    roundD();
    
    delta = updateLogLikelihood();
    std::cerr << "k == " << _k << " -- Restart = " << restart << " -- Iteration = " << ++step
              << " -- log likelihood = " << _logLikelihood << " -- "
//              << " -- gamma log likelihood = " << getLogLikelihoodGamma() << " -- "
              <<  timer.realTime() << " s (M)" << std::endl;
    if (!g_tol.nonZero(delta))
    {
      break;
    }
    if (step == maxIterations)
    {
      break;
    }
  }
  
  _solT = PosteriorStateTreeMatrix(n);
  _solZ.clear();
  for (int i = 0; i < n; ++i)
  {
    bool check = false;
    
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        const double y_itj = _solY[i][t][j];
        if (g_tol.nonZero(y_itj))
        {
          check = true;
          
          DoubleVector h_i;
          for (int p = 0; p < m; ++p)
          {
            h_i.push_back((_solD[j][p] - _numerator[i][t][p]) / _denominator[i][p]);
          }
          
          _solT[i].emplace_back(convertToStateTreeFromSNVF(_R, _scriptT[i][t],
                                                           h_i, i),
                                y_itj, t, j);
        }
      }
    }
//    std::sort(_solT[i].begin(),
//              _solT[i].end(),
//              [](const StateTreePosterior& a,
//                 const StateTreePosterior& b)
//                {
//                  return a._gamma > b._gamma;
//                });
    
    if (!check)
    {
        std::cout << "No solution found for mutation: " << i << std::endl;
        return false;
    }
    
    if(_solT[i].empty())
    {
        std::cout << "Mutation index: " << i << std::endl;
        std::cout << "script T size: " << _scriptT[i].size() << std::endl;
    }
    
    assert(!_solT[i].empty());
    _solZ.push_back(_solT[i].back()._j);
  }
  
  return true;
}

double EMPre::getLogLikelihoodGamma(int i) const
{
  const int m = _R.getNrSamples();
  const int n = _R.getNrCharacters();
  
  assert(0 <= i && i < n);
  
  double sum_i = 0;
  const StateEdgeSetVector& scriptT_i = _scriptT[i];
  const int size_scriptT_i = scriptT_i.size();
  const int maxCopyNumber = _R.getMaxCopies(i);
  
  for (int t = 0; t < size_scriptT_i; ++t)
  {
    const StateEdgeSet& T_it = scriptT_i[t];
    for (int j = 0; j < _k; ++j)
    {
      if (_solPi[j] == 0)
        continue;
      
      double prod = _solPi[j] / size_scriptT_i;
      double log_prod = log(prod);
      if (g_tol.nonZero(_solY[i][t][j]))
      {
        for (int p = 0; p < m; ++p)
        {
          const ReadMatrix::CopyNumberStateVector& cnStates_pi = _R.getCopyNumberStates(p, i);
          const int var_pi = _R.getVar(p, i);
          const int ref_pi = _R.getRef(p, i);
          
          assert(isFeasible(_solD[j][p],
                            _xyStar[i][t],
                            cnStates_pi,
                            T_it,
                            maxCopyNumber));
          {
            const double h = (_solD[j][p] - _numerator[i][t][p]) / _denominator[i][p];
            assert(_denominator[i][p] != 0);
            if (!((!(h <= 0 || !g_tol.nonZero(h)) || var_pi == 0) && (!(h >= 1 || g_tol.less(1, h)) || ref_pi == 0)))
            {
              prod = 0;
              break;
            }
            else
            {
              const double val = getLogLikelihood(var_pi, ref_pi, h);
              log_prod += val;
            }
          }
        }
        log_prod *= _solY[i][t][j];
        sum_i += log_prod;
      }
    }
  }
  
  return sum_i;
}

double EMPre::getLogLikelihoodGamma() const
{
  const int n = _R.getNrCharacters();
  
  double sum = 0;
  for (int i = 0; i < n; ++i)
  {
    sum += getLogLikelihoodGamma(i);
  }
  
  return sum;
}

void EMPre::writeMutationProperties(std::ostream& out) const
{
  const int m = _R.getNrSamples();
  const int n = _R.getNrCharacters();
  
  out << "SNV";
  for (int j = 0; j < _k; ++j)
  {
    out << "\t" << "cluster" << j;
  }
  out << "\tcluster";
  
  for (int p = 0; p < m; ++p)
  {
    if (_statType == CLUSTER_DCF)
    {
      out << "\tdcf" << p;
    }
    else
    {
      out << "\tccf" << p;
    }
  }
  out << std::endl;
  
  for (int i = 0; i < n; ++i)
  {
    out << i;
    double max_gamma_j = 0;
    int max_j = -1;
    for (int j = 0; j < _k; ++j)
    {
      double gamma = getGamma(i, j);
      if (max_gamma_j < gamma)
      {
        max_gamma_j = gamma;
        max_j = j;
      }
      out << "\t" << std::setprecision(2) << gamma;
    }
    out << "\t" << max_j;
    
    for (int p = 0; p < m; ++p)
    {
      out << "\t" << std::setprecision(2) << _solD[max_j][p];
    }
    out << std::endl;
  }
}

void EMPre::init()
{
  Solver::init();
  
  const int n = _R.getNrCharacters();
  
  for (int i = 0; i < n; ++i)
  {
    _solY[i] = DoubleMatrix(_scriptT[i].size(), DoubleVector(_k, 1));
  }
}

