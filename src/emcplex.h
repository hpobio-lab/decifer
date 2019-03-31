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
  EMCplex(const ReadMatrix& R,
          int k,
          int nrSegments,
          ClusterStatisticType statType);
  
  virtual ~EMCplex()
  {
    _env.end();
  }
  
protected:
  void initPWLA();
  
  bool stepM(int nrThreads,
             bool verbose);
  
  virtual std::unique_ptr<HardClusterIlp> createHardClusterIlpSolver(const ReadMatrix& R)
  {
    return std::unique_ptr<HardClusterIlp>(new HardClusterIlpCplex(R,
                                                                   _k,
                                                                   _nrSegments,
                                                                   _statType));
  }
  
  virtual std::unique_ptr<ClusterIlp> createClusterIlpSolver(const ReadMatrix& R,
                                                             double alpha)
  {
    return std::unique_ptr<ClusterIlp>(new ClusterIlpCplex(R,
                                                           _k,
                                                           alpha,
                                                           _statType));
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
  /// lambda[i][t][j][p][l]
  IloNumVar5Matrix _lambdaGRB;
  /// fGRB[j][p]
  IloNumVarMatrix _fGRB;
};

#endif // EMCPLEX_H
