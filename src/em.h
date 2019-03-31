/*
 * em.h
 *
 *  Created on: 5-nov-2017
 *      Author: M. El-Kebir
 */

#ifndef EM_H
#define EM_H

#include "utils.h"
#include "readmatrix.h"
#include "stategraph.h"
#include "solver.h"
#include "statetree.h"
#include "posteriorstatetree.h"
#include "hardclusterilp.h"
#include "clusterilp.h"

class EM : public Solver
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param k Number of clusters
  /// @param nrSegments Number of segments for piecewise linear approximation
  /// @param statType Summary statistic to use for clustering
  EM(const ReadMatrix& R,
     int k,
     int nrSegments,
     ClusterStatisticType statType);
  
  /// Destructor
  virtual ~EM()
  {
  }
  
  /// Solve
  ///
  /// @param seed Random number generator seed for k-Means.
  /// In case seed < 0, ILP will be used to seed the EM algorithm.
  /// @param maxIterations Maximum number of iterations, ignored when set to -1
  /// @param nrThreads Number of threads
  /// @param timeLimit Time limit in seconds (for initial ILP)
  /// @param verbose Verbose
  /// @param nrDownSampledSNVs Fraction of SNVs to consider
  /// @param memoryLimit Memory limit
  bool solve(int seed,
             int maxIterations,
             int nrThreads,
             int timeLimit,
             bool verbose,
             int nrDownSampledSNVs,
             int memoryLimit);
  
  double getGamma(int i, int j) const
  {
    double gamma = 0;
    
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      gamma += _gamma[i][t][j];
    }
    
    return gamma;
  }
  
  /// Write mutation properties
  ///
  /// @param out Output stream
  virtual void writeMutationProperties(std::ostream& out) const;
  
  /// Initialize solver
  void init();
  
  /// Return posterior state trees
  const PosteriorStateTreeMatrix& getPosteriorStateTrees() const
  {
    return _solT;
  }
  
  /// Initialize hot start
  ///
  /// @param y Clustering and state tree assignment
  void initHotStart(const BoolTensor& y)
  {
    _initY = y;
  }
  
  const BoolTensor& getInitY() const
  {
    return _initY;
  }
  
protected:
  /// E step
  virtual void stepE();
  
  /// M step
  virtual bool stepM(int nrThreads,
                     bool verbose) = 0;
  
  /// Initialize piecewise linear approximation
  virtual void initPWLA();
  
  virtual std::unique_ptr<ClusterIlp> createClusterIlpSolver(const ReadMatrix& R,
                                             double alpha) = 0;
  
  virtual std::unique_ptr<HardClusterIlp> createHardClusterIlpSolver(const ReadMatrix& R) = 0;
  
  void kMeans(int seed);
  
  virtual bool isEnabled(int i, int t, int j) const
  {
    assert(0 <= i && i < _R.getNrCharacters());
    assert(0 <= t && t < _scriptT[i].size());
    assert(0 <= j && j < _k);
    return g_tol.nonZero(_gamma[i][t][j]);
  }
  
  bool initializeD(int seed,
                   int nrThreads,
                   int timeLimit,
                   bool verbose,
                   int nrDownSampledSNVs,
                   int memoryLimit);
  
protected:
  /// _gamma[i][t][j]
  DoubleTensor _gamma;
  /// Solution state trees
  PosteriorStateTreeMatrix _solT;
  /// Initial hot start
  BoolTensor _initY;
};

#endif // EM_H
