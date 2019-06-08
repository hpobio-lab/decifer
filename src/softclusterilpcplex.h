/*
 * softclusterilpcplex.h
 *
 *  Created on: 15-apr-2019
 *      Author: M. El-Kebir
 */

#ifndef SOFTCLUSTERILPCPLEX_H
#define SOFTCLUSTERILPCPLEX_H

#include "softclusterilp.h"
#include <ilcplex/ilocplex.h>

class SoftClusterIlpCplex : public SoftClusterIlp
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
  SoftClusterIlpCplex(const ReadMatrix& R,
                      int k,
                      int nrSegments,
                      ClusterStatisticType statType,
                      double precisionBetaBin,
                      bool forceTruncal);
  
  /// Export ILP
  ///
  /// @param filename Filename
  virtual void exportModel(const std::string& filename)
  {
    _cplex.exportModel(filename.c_str());
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
                     int memoryLimit);
  
  /// Initialize hot start
  ///
  /// @param y Clustering and state tree assignment
  void initHotStart(const BoolTensor& y);
  
  virtual ~SoftClusterIlpCplex()
  {
    _env.end();
  }
  
protected:
  /// Initialize pre clustering constraint
  ///
  /// @param i1 SNV
  /// @param i2 SNV
  virtual void initPreClusteringConstraint(int i1, int i2);
  
  /// Initialize variables
  virtual void initVariables();
  
  /// Initialize constraints
  virtual void initConstraints();
  
  /// Initialize objective function
  virtual void initObjective();
  
private:
  typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;
  typedef IloArray<IloBoolVarMatrix> IloBoolVar3Matrix;
  typedef IloArray<IloNumVarArray> IloNumVarMatrix;
  typedef IloArray<IloNumVarMatrix> IloNumVar3Matrix;
  typedef IloArray<IloNumVar3Matrix> IloNumVar4Matrix;
  typedef IloArray<IloNumVar4Matrix> IloNumVar5Matrix;
  typedef IloArray<IloNumVar5Matrix> IloNumVar6Matrix;
  
  /// Environment
  IloEnv _env;
  /// CPlex model
  IloModel _model;
  /// Solver
  IloCplex _cplex;
  /// y[i][t][j] -- posterior probability of state tree t and cluster j for SNV i
  IloNumVar3Matrix _y;
  /// yy[i][t][j] = 1 iff y[i][t][j] > 0
  IloBoolVar3Matrix _yy;
  /// d[j][p]
  IloNumVarMatrix _d;
  /// yyd[i][t][j][p] = yy[i][t][j] * d[j][p]
  IloNumVar4Matrix _yyd;
  /// pi[j]
  IloNumVarArray _pi;
  /// gamma[i][t][j][alpha][beta]
  IloNumVar5Matrix _gamma;
  /// lambda[i][t][j][p][alpha][beta
  IloNumVar6Matrix _lambda;
};

#endif // SOFTCLUSTERILPCPLEX_H
