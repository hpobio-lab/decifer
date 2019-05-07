/*
 * clusterilp.h
 *
 *  Created on: 10-nov-2017
 *      Author: M. El-Kebir
 */

#ifndef CLUSTERILP_H
#define CLUSTERILP_H

#include "utils.h"
#include "readmatrix.h"
#include "stategraph.h"
#include "statetree.h"
#include "solver.h"

class ClusterIlp : public Solver
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param k Number of clusters
  /// @param alpha Confidence interval width
  /// @param statType Summary statistic to use for clustering
  /// @param precisionBetaBin Precision parameter for beta binomial
  /// @param forceTruncal Force the presence of a dominant truncal cluster
  ClusterIlp(const ReadMatrix& R,
             int k,
             double alpha,
             ClusterStatisticType statType,
             double precisionBetaBin,
             bool forceTruncal);
  
  /// Get VAF lower bound
  ///
  /// @param p Sample
  /// @param i character
  double getVAFlb(int p,
                  int i) const
  {
    return _vafLB[i][p];
  }
  
  /// Get VAF upper bound
  ///
  /// @param p Sample
  /// @param i character
  double getVAFub(int p,
                  int i) const
  {
    return _vafUB[i][p];
  }
  
  const StateTreeVector& getScriptSlb(int i) const
  {
    assert(0 <= i && i < _R.getNrCharacters());
    return _scriptTlb[i];
  }
  
  const StateTreeVector& getScriptSub(int i) const
  {
    assert(0 <= i && i < _R.getNrCharacters());
    return _scriptTub[i];
  }
  
  /// Solve
  ///
  /// @param nrThreads Number of threads
  /// @param timeLimit Time limit in seconds (unlimited if -1)
  /// @param verbose Verbose
  /// @param memoryLimit Memory limit
  virtual bool solve(int nrThreads,
                     int timeLimit,
                     bool verbose,
                     int memoryLimit) = 0;
  
  virtual void exportModel(const std::string& filename) = 0;
  
  /// Return state trees
  const StateTreeVector& getStateTrees() const
  {
    return _solT;
  }
  
  /// Return expected state tree
  const StateTree& Texp(int i,
                        int t) const
  {
    assert(0 <= i && i < _R.getNrCharacters());
    assert(0 <= t && t < _scriptT[i].size());
    
    return _scriptTexp[i][t];
  }
  
  /// Initialize
  virtual void init();
  
  /// Initialize hot start
  ///
  /// @param y Clustering and state tree assignment
  virtual void initHotStart(const BoolTensor& y) = 0;
  
  void writeFeasibleSolutionSpace(std::ostream& out) const;
  
  const BoolTensor& getSolY() const
  {
    return _solY;
  }
  
protected:
  typedef StateGraph::StateEdgeSet StateEdgeSet;
  typedef std::vector<StateEdgeSet> StateEdgeSetVector;
  typedef std::vector<StateEdgeSetVector> StateEdgeSetMatrix;
  
  typedef std::vector<StateTreeVector> StateTreeMatrix;
  
  void initStateTrees();
  
  virtual void initVariables() = 0;
  
  virtual void initConstraints() = 0;
  
  void initVAFs();
  
protected:
  /// Confidence interval width
  const double _alpha;
  /// vafLB[i][p]
  DoubleMatrix _vafLB;
  /// vafUB[i][p]
  DoubleMatrix _vafUB;
  /// State trees \mathcal{T}(\mu_i, f_i)
  StateTreeMatrix _scriptTlb;
  /// State trees \mathcal{T}(\mu_i, f_i)
  StateTreeMatrix _scriptTub;
  /// State trees \mathcal{T}(\mu_i, f_i)
  StateTreeMatrix _scriptTexp;
  /// dcfLB[i][t][p]
  DoubleTensor _dcfLB;
  /// dcfUB[i][t][p]
  DoubleTensor _dcfUB;
  ///
  DoubleTensor _dcfExp;
  /// Solution state trees
  StateTreeVector _solT;
  /// Solution y
  BoolTensor _solY;
};

#endif // CLUSTERILP_H
