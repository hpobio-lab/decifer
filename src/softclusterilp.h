/*
 * softclusterilp.h
 *
 *  Created on: 13-apr-2019
 *      Author: M. El-Kebir
 */

#ifndef SOFTCLUSTERILP_H
#define SOFTCLUSTERILP_H

#include "utils.h"
#include "readmatrix.h"
#include "stategraph.h"
#include "posteriorstatetree.h"
#include "statetree.h"
#include "solver.h"

class SoftClusterIlp : public Solver
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param k Number of clusters
  /// @param nrSegments Number of segments for piecewise linear approximation
  /// @param statType Summary statistic to use for clustering
  /// @param precisionBetaBin Precision parameter for beta binomial
  /// @param forceTruncal Force the presence of a dominant truncal cluster
  SoftClusterIlp(const ReadMatrix& R,
                 int k,
                 int nrSegments,
                 ClusterStatisticType statType,
                 double precisionBetaBin,
                 bool forceTruncal);
  
  /// Destructor
  virtual ~SoftClusterIlp()
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
  
  /// Return posterior state trees
  const PosteriorStateTreeMatrix& getPosteriorStateTrees() const
  {
    return _solT;
  }
  
  /// Initialize solver
  virtual void init();
  
  /// Initialize hot start
  ///
  /// @param y Clustering and state tree assignment
  virtual void initHotStart(const BoolTensor& y) = 0;
  
  const DoubleTensor& getSolY() const
  {
    return _solY;
  }
  
  /// Initialize DCF constraints
  virtual void initConstraintsDCF(const DoubleMatrix& d,
                                  int multiplicity)
  {
  };
  
protected:
  typedef std::vector<Double5Matrix> Double6Matrix;
  
  /// Initialize piecewise linear approximation
  virtual void initPWLA();
  
  /// Initialize variables
  virtual void initVariables() = 0;

  /// Initialize constraints
  virtual void initConstraints() = 0;
  
  /// Initialize objective
  virtual void initObjective() = 0;
  
  double getCoord(int alpha) const
  {
    assert(0 <= alpha && alpha < _nrSegments);
    return _coord[alpha];
  }
   
protected:
  /// coord[alpha]
  DoubleVector _coord;
  /// coordPi[alpha]
  DoubleVector _coordPi;
  /// G[i][t][j][p][alpha]
  Double5Matrix _hatG;
  /// hatPi[alpha]
  DoubleVector _hatPi;
  /// Solution state trees
  PosteriorStateTreeMatrix _solT;
  /// Solution y
  DoubleTensor _solY;
};

#endif // SOFTCLUSTERILP_H
