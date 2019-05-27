/*
 * em.cpp
 *
 *  Created on: 5-nov-2017
 *      Author: M. El-Kebir
 */

#include "em.h"
#include "stategraph.h"
#include "dkm/dkm.hpp"
#include <tuple>
#include <iomanip>
#include <lemon/time_measure.h>

EM::EM(const ReadMatrix& R,
       int k,
       int nrSegmentBits,
       ClusterStatisticType statType,
       double precisionBetaBin,
       bool forceTruncal)
  : Solver(R, k, (1 << nrSegmentBits) + 1, statType, precisionBetaBin, forceTruncal)
  , _gamma(R.getNrCharacters())
  , _solT()
  , _initY()
{
}

void EM::kMeans(int seed)
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  assert(m <= 10);
  
  std::vector<std::array<double, 10> > data;
  for (int i = 0; i < n; ++i)
  {
    if (_R.isCopyNeutral(i))
    {
      data.push_back(std::array<double, 10>{{0,0,0,0,0,0,0,0,0,0}});
      for (int p = 0; p < m; ++p)
      {
        const int var_pi = _R.getVar(p, i);
        const int ref_pi = _R.getRef(p, i);
        data.back()[p] = double(var_pi) / double(var_pi + ref_pi);
      }
    }
  }
  
  std::cerr << "Identified " << data.size() << " copy-neutral SNVs" << std::endl;
  std::cerr << "Initializing EM algorithm using k-means with seed " << seed << std::endl;
  
  std::vector<std::array<double, 10> > means = std::get<0>(dkm::kmeans_lloyd(data, _k, seed));
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < m; ++p)
    {
      _solD[j][p] = std::min(2 * means[j][p], 1.);
    }
  }
}

void EM::initPWLA()
{
  Solver::initPWLA();
  
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
            //            if (_fLB[i][t][p] > _fOverallLB[j][p])
            if (g_tol.less(_dOverallLB[j][p], _dLB[i][t][p]))
            {
              _dOverallLB[j][p] = _dLB[i][t][p];
            }
            //            if (_fUB[i][t][p] < _fOverallUB[j][p])
            if (g_tol.less(_dUB[i][t][p], _dOverallUB[j][p]))
            {
              _dOverallUB[j][p] = _dUB[i][t][p];
            }
            
            if (!g_tol.different(_dOverallLB[j][p], _dOverallUB[j][p]))
            {
              _dOverallLB[j][p] = _dOverallUB[j][p];
            }
            
            _dOverallLB[j][p] = _dOverallLB[j][p] <= g_thre ? g_thre : _dOverallLB[j][p];
            _dOverallUB[j][p] = _dOverallUB[j][p] <= g_thre ? g_thre : _dOverallUB[j][p];
            
            assert(_dOverallLB[j][p] <= _dOverallUB[j][p]);
          }
        }
      }
    }
  }
}

bool EM::initializeD(int seed,
                     int nrThreads,
                     int timeLimit,
                     bool verbose,
                     int nrDownSampledSNVs,
                     int memoryLimit)
{
  ReadMatrix R = _R.downSampleCharacters(nrDownSampledSNVs);
  
  const int m = R.getNrSamples();
  
  if (seed > 0)
  {
    kMeans(seed);
  }
  else if (seed == 0)
  {
    // hard clustering (distance)
//    std::cerr << "Hard clustering (distance)..." << std::endl;
    std::unique_ptr<ClusterIlp> pILP = createClusterIlpSolver(R, 0);
    pILP->init();
//    if (!_initY.empty())
//    {
//      pILP->initHotStart(_initY);
//    }
    
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
//      _solPi = pILP->getPi();
//      _initY = pILP->getSolY();
    }
  }
  else if (seed % 2 == 0)
  {
    {
//      std::cerr << "Hard clustering (probabilistic)..." << std::endl;
      std::unique_ptr<IncrementalSolver> pILP = createIncrementalSolver(R);
      pILP->init();
//      if (!_initY.empty())
//      {
//        pILP->initHotStart(_initY);
//      }
      
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
//        _solPi = pILP->getPi();
//        _initY = pILP->getSolY();
        
        pILP = NULL;
      }
    }
  }
  else
  {
    double alpha = 0.25;
    double limit = 1. / (1 << 20);
    while (alpha >= limit)
    {
//      std::cerr << "Hard clustering (confidence intervals)..." << std::endl;
      std::unique_ptr<ClusterIlp> pILP = createClusterIlpSolver(R, alpha);
      pILP->init();
      
      if (!pILP->solve(nrThreads, timeLimit, verbose, memoryLimit))
      {
        alpha /= 2;
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
        _solPi = pILP->getPi();
        _initY = pILP->getSolY();
        
        break;
      }
    }
    
    if (alpha < limit)
    {
      return false;
    }
  }
  
  return true;
}

bool EM::solve(int restart,
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
  
  int step = 0;
  while (true)
  {
    stepE();
    if (!stepM(nrThreads, verbose))
    {
      std::cerr << "k == " << _k << " -- Restart = " << restart << " -- Iteration = " << ++step << " -- infeasible solution! " << " -- " << timer.realTime() << " s" << std::endl;
      return false;
    }
    
    double delta = updateLogLikelihood();
    std::cerr << "k == " << _k << " -- Restart = " << restart << " -- Iteration = " << ++step
              << " -- log likelihood = " << _logLikelihood << " -- "
//              << " -- gamma log likelihood = " << getLogLikelihoodGamma() << " -- "
              <<  timer.realTime() << " s" << std::endl;
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
        const double gamma_itj = _gamma[i][t][j];
        if (g_tol.nonZero(gamma_itj))
        {
          check = true;
            
          DoubleVector h_i;
          for (int p = 0; p < m; ++p)
          {
            h_i.push_back((_solD[j][p] - _numerator[i][t][p]) / _denominator[i][p]);
          }
          
          _solT[i].emplace_back(convertToStateTreeFromSNVF(_scriptT[i][t],
                                                           h_i, i),
                                gamma_itj, t, j);
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
      
    if(!check)
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

double EM::getLogLikelihoodGamma(int i) const
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
      if (g_tol.nonZero(_gamma[i][t][j]))
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
        log_prod *= _gamma[i][t][j];
        sum_i += log_prod;
      }
    }
  }
  
  return sum_i;
}

double EM::getLogLikelihoodGamma() const
{
  const int n = _R.getNrCharacters();
  
  double sum = 0;
  for (int i = 0; i < n; ++i)
  {
    sum += getLogLikelihoodGamma(i);
  }
  
  return sum;
}

void EM::writeMutationProperties(std::ostream& out) const
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

void EM::stepE()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  for (int i = 0; i < n; ++i)
  {
    const StateEdgeSetVector& scriptT_i = _scriptT[i];
    const int maxCopyNumber = _R.getMaxCopies(i);
    
    double sum = 0;
    const int size_scriptT_i = scriptT_i.size();
    for (int t = 0; t < size_scriptT_i; ++t)
    {
      const StateEdgeSet& T_it = scriptT_i[t];
      for (int j = 0; j < _k; ++j)
      {
        _gamma[i][t][j] = 0;
        
        if (_solPi[j] == 0)
          continue;
        
        double prod = _solPi[j] / size_scriptT_i;
        double log_prod = log(prod);
        for (int p = 0; p < m; ++p)
        {
          const ReadMatrix::CopyNumberStateVector& cnStates_pi = _R.getCopyNumberStates(p, i);
          const int var_pi = _R.getVar(p, i);
          const int ref_pi = _R.getRef(p, i);
          
          if (isFeasible(_solD[j][p],
                         _xyStar[i][t],
                         cnStates_pi,
                         T_it,
                         maxCopyNumber))
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
          else
          {
            prod = 0;
            break;
          }
        }
        
        if (g_tol.nonZero(prod))
        {
          //sum_i += prod;
          prod = exp(log_prod);
          prod = std::max(std::numeric_limits<double>::min(), prod);
          sum += prod;
          _gamma[i][t][j] = prod;
        }
      }
    }
    
    bool ok = false;
    for (int t = 0; t < scriptT_i.size(); ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        _gamma[i][t][j] /= sum;
        if (!g_tol.nonZero(_gamma[i][t][j]))
        {
          _gamma[i][t][j] = 0;
        }
        else
        {
          if (!g_tol.different(1, _gamma[i][t][j]))
          {
            _gamma[i][t][j] = 1;
          }
          ok = true;
        }
      }
    }
//    if (!ok)
//    {
//      std::cerr << "Houston" << std::endl;
//    }
//    assert(ok);
  }
}

void EM::init()
{
  Solver::init();
  
  const int n = _R.getNrCharacters();
  
  for (int i = 0; i < n; ++i)
  {
    _gamma[i] = DoubleMatrix(_scriptT[i].size(),
                             DoubleVector(_k, 1));
  }
}
