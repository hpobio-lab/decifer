/*
 * precluster.h
 *
 *  Created on: 6-jun-2019
 *      Author: M. El-Kebir
 */

#ifndef PRECLUSTER_H
#define PRECLUSTER_H

#include "utils.h"
#include "readmatrix.h"
#include "solver.h"

#ifdef CPLEX
  #include "softpreclusterilpcplex.h"
#else
#endif

class PreCluster
{
public:
  
#ifdef CPLEX
  typedef SoftPreClusterIlpCplex PreClusterIlpAlg;
#endif
  
  /// Constructor
  ///
  /// @param R Read matrix
  PreCluster(const ReadMatrix& R);
  
  /// Run preclustering
  ///
  /// @param nrSegments Number of segments for piecewise linear approximation
  /// @param statType Summary statistic to use for clustering
  /// @param precisionBetaBin Precision parameter for beta binomial
  /// @param nrThreads Number of threads
  /// @param timeLimit Time limit in seconds (unlimited if -1)
  /// @param verbose Verbose
  /// @param memoryLimit Memory limit
  void run(int nrSegments,
           Solver::ClusterStatisticType statType,
           double precisionBetaBin,
           int nrThreads,
           int timeLimit,
           bool verbose,
           int memoryLimit);
  
  /// Return pre clustering
  const IntMatrix& getPreClustering() const
  {
    return _invPreClustering;
  }
  
protected:
  /// Run preclustering for specified cluster
  ///
  /// @param snvIndices SNV indices that comprise cluster
  /// @param nrSegments Number of segments for piecewise linear approximation
  /// @param statType Summary statistic to use for clustering
  /// @param precisionBetaBin Precision parameter for beta binomial
  /// @param nrThreads Number of threads
  /// @param timeLimit Time limit in seconds (unlimited if -1)
  /// @param verbose Verbose
  /// @param memoryLimit Memory limit
  void run(const IntVector& snvIndices,
           int nrSegments,
           Solver::ClusterStatisticType statType,
           double precisionBetaBin,
           int nrThreads,
           int timeLimit,
           bool verbose,
           int memoryLimit);
  
protected:
  /// Read matrix
  const ReadMatrix& _R;
  /// Number of clusters
  int _nrClusters;
  /// Pre clustering
  IntVector _preClustering;
  /// Inverse pre clustering
  IntMatrix _invPreClustering;
};

#endif // PRECLUSTER
