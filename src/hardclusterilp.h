/*
 * hardclusterilp.h
 *
 *  Created on: 18-oct-2018
 *      Author: M. El-Kebir
 */

#ifndef HARDCLUSTERILP_H
#define HARDCLUSTERILP_H

#include "utils.h"
#include "readmatrix.h"
#include "stategraph.h"
#include "statetree.h"
#include "solver.h"

class HardClusterIlp : public Solver
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param k Number of clusters
  /// @param nrSegments Number of segments for piecewise linear approximation
  /// @param statType Summary statistic to use for clustering
  /// @param forceTruncal Force the presence of a dominant truncal cluster
  HardClusterIlp(const ReadMatrix& R,
                 int k,
                 int nrSegments,
                 ClusterStatisticType statType,
                 bool forceTruncal);
  
  /// Destructor
  virtual ~HardClusterIlp()
  {
  };
  
  /// Solve
  ///
  /// @param nrThreads Number of threads
  /// @param timeLimit Time limit in seconds (for initial ILP)
  /// @param verbose Verbose
  /// @param memoryLimit Memory limit
  virtual bool solve(int nrThreads,
                     int timeLimit,
                     bool verbose,
                     int memoryLimit) = 0;
  
  /// Export ILP model
  ///
  /// @param filename Filename
  virtual void exportModel(const std::string& filename) = 0;
  
  /// Return state trees
  const StateTreeVector& getStateTrees() const
  {
    return _solT;
  }
  
  /// Initialize solver
  virtual void init();
  
  /// Initialize hot start
  ///
  /// @param y Clustering and state tree assignment
  virtual void initHotStart(const BoolTensor& y) = 0;
  
  const BoolTensor& getSolY() const
  {
    return _solY;
  }
  
protected:
  /// Initialize variables
  virtual void initVariables() = 0;

  /// Initialize constraints
  virtual void initConstraints() = 0;
  
  /// Initialize objective
  virtual void initObjective() = 0;
  
  /// Initialize piecewise linear approximation
  virtual void initPWLA();
  
protected:
  /// hatN[l]
  DoubleVector _hatN;
  /// z[l]
  DoubleVector _z;
  /// Solution state trees
  StateTreeVector _solT;
  /// Solution y
  BoolTensor _solY;
};

#endif // HARDCLUSTERILP_H
