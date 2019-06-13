/*
 * incrementalsolver.cpp
 *
 *  Created on: 13-may-2019
 *      Author: M. El-Kebir
 */

#include "incrementalsolver.h"
#include "softclusterlpcplex.h"

IncrementalSolver::IncrementalSolver(const ReadMatrix& R,
                                     int k,
                                     int maxNrSegmentBits,
                                     ClusterStatisticType statType,
                                     double precisionBetaBin,
                                     bool forceTruncal,
                                     const IntMatrix& preClustering)
  : Solver(R, k, 0, statType, precisionBetaBin, forceTruncal)
  , _maxNrSegmentBits(maxNrSegmentBits)
  , _preClustering(preClustering)
  , _solT()
  , _solY()
{
}

void IncrementalSolver::init()
{
  Solver::init();
}

bool IncrementalSolver::solve(int nrThreads,
                              int timeLimit,
                              bool verbose,
                              int memoryLimit,
                              int nrDownSampledSNVs,
                              int nrRestarts)
{
  for (int restart = 1; restart <= nrRestarts; ++restart)
  {
    ReadMatrix R = _R.downSampleCharacters(nrDownSampledSNVs);
    
    DoubleTensor y;
    DoubleMatrix d;
    DoubleVector pi;
    double logLikelihood;
    if (!solve(R, nrThreads, timeLimit, verbose, memoryLimit, false, d, y, pi, logLikelihood))
    {
      return false;
    }
    
    // resolve with DCF matrix set
    SoftClusterLpCplex solver(_R, _k, _maxNrSegmentBits, _statType, _precisionBetaBin, _forceTruncal, false);
    solver.init();
    solver.initConstraintsDCF(d, 1);
    if (!_preClustering.empty())
    {
      solver.initPreClusteringConstraints(_preClustering);
    }
    
    solver.solve(nrThreads, timeLimit, verbose, memoryLimit);
    logLikelihood = solver.getLogLikelihood();
    y = solver.getSolY();
    d = solver.getD();
    pi = solver.getPi();

    std::cerr << "Resolved -- nrBits = "
              << _maxNrSegmentBits
              << " -- log likelihood = "
              << _logLikelihood << std::endl;

    if (_logLikelihood < logLikelihood)
    {
      _logLikelihood = solver.getLogLikelihood();
      _solY = solver.getSolY();
      _solD = solver.getD();
      _solPi = solver.getPi();
    }
  }
  
  return true;
}

bool IncrementalSolver::solve(int nrThreads,
                              int timeLimit,
                              bool verbose,
                              int memoryLimit)
{
  return solve(_R, nrThreads, timeLimit, verbose,
               memoryLimit, false,
               _solD, _solY, _solPi, _logLikelihood);
}

bool IncrementalSolver::solve(const ReadMatrix& R,
                              int nrThreads,
                              int timeLimit,
                              bool verbose,
                              int memoryLimit,
                              bool includePi,
                              DoubleMatrix& d,
                              DoubleTensor& y,
                              DoubleVector& pi,
                              double& logLikelihood)
{
  // initialization
  {
    SoftClusterLpCplex solver(R, _k, 1,
                              _statType, _precisionBetaBin,
                              _forceTruncal, includePi);
    solver.init();
    solver.solve(nrThreads, timeLimit, verbose, memoryLimit);
    logLikelihood = solver.getLogLikelihood();
    y = solver.getSolY();
    d = solver.getD();
    pi = solver.getPi();
  }
  
  std::cerr << "nrBits = " << 1 << " -- log likelihood = " << logLikelihood << std::endl << std::endl;
  
  for (int nrBits = 2; nrBits <= _maxNrSegmentBits; ++nrBits)
  {
    SoftClusterLpCplex solver(R, _k, nrBits, _statType, _precisionBetaBin, _forceTruncal, includePi);
    solver.init();
    solver.initConstraintsDCF(_solD, 3);
    if (!_preClustering.empty())
    {
      solver.initPreClusteringConstraints(_preClustering);
    }
//    solver.initConstraintsY(y);
    
    solver.solve(nrThreads, timeLimit, verbose, memoryLimit);
    logLikelihood = solver.getLogLikelihood();
    y = solver.getSolY();
    d = solver.getD();
    pi = solver.getPi();
    
    std::cerr << "nrBits = " << nrBits << " -- log likelihood = " << logLikelihood << std::endl << std::endl;
  }
  
  return true;
}
