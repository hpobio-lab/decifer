/*
 * preclusterkmeans.h
 *
 *  Created on: 13-jun-2019
 *      Author: M. El-Kebir
 */

#ifndef PRECLUSTERKMEANS_H
#define PRECLUSTERKMEANS_H

#include "utils.h"
#include "readmatrix.h"
#include "stategraph.h"
#include "posteriorstatetree.h"
#include "statetree.h"
#include "solver.h"

class PreClusterKMeans : public Solver
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param k Number of clusters
  /// @param statType Summary statistic to use for clustering
  /// @param precisionBetaBin Precision parameter for beta binomial
  PreClusterKMeans(const ReadMatrix& R,
                   int k,
                   ClusterStatisticType statType,
                   double precisionBetaBin);
  
  /// Destructor
  virtual ~PreClusterKMeans()
  {
  };
  
  /// Initialize solver
  virtual void init();
  
  /// Solve
  ///
  /// @param nrRestarts Number of restarts
  /// @param verbose Verbose
  virtual bool solve(int nrRestarts,
                     bool verbose);
};

#endif // PRECLUSTERKMEANS_H
