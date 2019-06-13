/*
 * incrementalpreclustersolver.cpp
 *
 *  Created on: 12-jun-2019
 *      Author: M. El-Kebir
 */

#include "incrementalpreclustersolver.h"
#include "softpreclusterilpcplex.h"

IncrementalPreClusterSolver::IncrementalPreClusterSolver(const ReadMatrix& R,
                                                         int k,
                                                         int maxNrSegmentBits,
                                                         ClusterStatisticType statType,
                                                         double precisionBetaBin)
  : Solver(R, k, 0, statType, precisionBetaBin, false)
  , _maxNrSegmentBits(maxNrSegmentBits)
  , _solT()
  , _solY()
{
}

void IncrementalPreClusterSolver::init()
{
  Solver::init();
}

bool IncrementalPreClusterSolver::solve(int nrThreads,
                                        int timeLimit,
                                        bool verbose,
                                        int memoryLimit)
{
  return solve(_R, nrThreads, timeLimit, verbose,
               memoryLimit, false,
               _solD, _solY, _solPi, _logLikelihood);
}

bool IncrementalPreClusterSolver::solve(const ReadMatrix& R,
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
    SoftPreClusterIlpCplex solver(R, _k, 1,
                                  _statType, _precisionBetaBin,
                                  includePi);
    solver.init();
    solver.solve(nrThreads, timeLimit, verbose, memoryLimit);
    logLikelihood = solver.getLogLikelihood();
    y = solver.getSolY();
    d = solver.getD();
    pi = solver.getPi();
  }
  
//  std::cerr << "nrBits = " << 1 << " -- log likelihood = " << logLikelihood << std::endl << std::endl;
  
  for (int nrBits = 2; nrBits <= _maxNrSegmentBits; ++nrBits)
  {
    SoftPreClusterIlpCplex solver(R, _k, nrBits, _statType, _precisionBetaBin, includePi);
    solver.init();
    solver.initConstraintsDCF(_solD, 3);
//    solver.initConstraintsY(y);
    
    solver.solve(nrThreads, timeLimit, verbose, memoryLimit);
    logLikelihood = solver.getLogLikelihood();
    y = solver.getSolY();
    d = solver.getD();
    pi = solver.getPi();
    
//    std::cerr << "nrBits = " << nrBits << " -- log likelihood = " << logLikelihood << std::endl << std::endl;
  }
  
  return true;
}

