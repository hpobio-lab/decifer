/*
 * emcplex.h
 *
 *  Created on: 16-oct-2018
 *      Author: M. El-Kebir
 */

#ifndef EMCPLEX_H
#define EMCPLEX_H

#include "em.h"
#include "hardclusterilpcplex.h"
#include "clusterilpcplex.h"
#include "softclusterlpcplex.h"
#include <ilcplex/ilocplex.h>

class EMCplex : public EM
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
  EMCplex(const ReadMatrix& R,
          int k,
          int nrSegments,
          ClusterStatisticType statType,
          double precisionBetaBin,
          bool forceTruncal);
  
  virtual ~EMCplex()
  {
    _env.end();
  }
  
protected:
  void initPWLA();
  
  bool stepM(int nrThreads,
             bool verbose);
  
  virtual std::unique_ptr<IncrementalSolver> createIncrementalSolver(const ReadMatrix& R)
  {
    return std::unique_ptr<IncrementalSolver>(new IncrementalSolver(R,
                                                                    _k,
                                                                    6,
                                                                    _statType,
                                                                    _precisionBetaBin,
                                                                    _forceTruncal, IntMatrix()));
  }
  
  
  virtual std::unique_ptr<SoftClusterIlp> createSoftClusterIlpSolver(const ReadMatrix& R)
  {
    return std::unique_ptr<SoftClusterIlp>(new SoftClusterLpCplex(R,
                                                                  _k,
                                                                  3,
                                                                  _statType,
                                                                  _precisionBetaBin,
                                                                  _forceTruncal,
                                                                  false));
  }
  
  virtual std::unique_ptr<HardClusterIlp> createHardClusterIlpSolver(const ReadMatrix& R)
  {
    return std::unique_ptr<HardClusterIlp>(new HardClusterIlpCplex(R,
                                                                   _k,
                                                                   _nrSegments,
                                                                   _statType,
                                                                   _precisionBetaBin,
                                                                   _forceTruncal,
                                                                   false));
  }
  
  virtual std::unique_ptr<ClusterIlp> createClusterIlpSolver(const ReadMatrix& R,
                                                             double alpha)
  {
    return std::unique_ptr<ClusterIlp>(new ClusterIlpCplex(R,
                                                           _k,
                                                           alpha,
                                                           _statType,
                                                           _precisionBetaBin,
                                                           _forceTruncal));
  }
  
  typedef IloArray<IloNumVarArray> IloNumVarMatrix;
  typedef IloArray<IloNumVarMatrix> IloNumVar3Matrix;
  typedef IloArray<IloNumVar3Matrix> IloNumVar4Matrix;
  typedef IloArray<IloNumVar4Matrix> IloNumVar5Matrix;
  
private:
  /// Environment
  IloEnv _env;
  /// CPlex model
  IloModel _model;
  /// Solver
  IloCplex _cplex;
  /// lambda[j][p][alpha]
  IloNumVar3Matrix _lambda;
  /// d[j][p]
  IloNumVarMatrix _d;
};

#endif // EMCPLEX_H
