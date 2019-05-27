/*
 * incrementalsolver.h
 *
 *  Created on: 13-may-2019
 *      Author: M. El-Kebir
 */

#ifndef INCREMENTALSOLVER_H
#define INCREMENTALSOLVER_H

#include "utils.h"
#include "readmatrix.h"
#include "stategraph.h"
#include "posteriorstatetree.h"
#include "statetree.h"
#include "solver.h"

class IncrementalSolver : public Solver
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param k Number of clusters
  /// @param maxNrSegmentBits Maximum number of bits for segments for piecewise linear approximation
  /// @param statType Summary statistic to use for clustering
  /// @param precisionBetaBin Precision parameter for beta binomial
  /// @param forceTruncal Force the presence of a dominant truncal cluster
  IncrementalSolver(const ReadMatrix& R,
                    int k,
                    int maxNrSegmentBits,
                    ClusterStatisticType statType,
                    double precisionBetaBin,
                    bool forceTruncal);
  
  /// Destructor
  virtual ~IncrementalSolver()
  {
  };
  
  /// Initialize solver
  virtual void init();
  
  /// Solve
  ///
  /// @param nrThreads Number of threads
  /// @param timeLimit Time limit in seconds
  /// @param verbose Verbose
  /// @param memoryLimit Memory limit
  virtual bool solve(int nrThreads,
                     int timeLimit,
                     bool verbose,
                     int memoryLimit);
  
  /// Solve with downsampling
  ///
  /// @param nrThreads Number of threads
  /// @param timeLimit Time limit in seconds
  /// @param verbose Verbose
  /// @param memoryLimit Memory limit
  /// @param nrDownSampledSNVs Fraction of SNVs to consider
  /// @param nrRestarts Number of restarts
  bool solve(int nrThreads,
             int timeLimit,
             bool verbose,
             int memoryLimit,
             int nrDownSampledSNVs,
             int nrRestarts);
  
  /// Return posterior state trees
  const PosteriorStateTreeMatrix& getPosteriorStateTrees() const
  {
    return _solT;
  }
  
protected:
  /// Solve
  ///
  /// @param R Read matrix
  /// @param nrThreads Number of threads
  /// @param timeLimit Time limit in seconds
  /// @param verbose Verbose
  /// @param memoryLimit Memory limit
  /// @param includePi Include pi in inference
  /// @param d DCF matrix
  /// @param y
  /// @param pi
  /// @param logLikelihood
  bool solve(const ReadMatrix& R,
             int nrThreads,
             int timeLimit,
             bool verbose,
             int memoryLimit,
             bool includePi,
             DoubleMatrix& d,
             DoubleTensor& y,
             DoubleVector& pi,
             double& logLikelihood);
  
protected:
  /// Maximum number of segment bits
  const int _maxNrSegmentBits;
  /// Solution state trees
  PosteriorStateTreeMatrix _solT;
  /// Solution y
  DoubleTensor _solY;
};

#endif // INCREMENTALSOLVER_H
