/*
 * minclusterilpcplex.h
 *
 *  Created on: 21-oct-2018
 *      Author: M. El-Kebir
 */

#ifndef MINCLUSTERILPCPLEX_H
#define MINCLUSTERILPCPLEX_H

#include "minclusterilp.h"
#include <ilcplex/ilocplex.h>

class MinClusterIlpCplex : public MinClusterIlp
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param kMax Number of clusters
  /// @param statType Summary statistic to use for clustering
  MinClusterIlpCplex(const ReadMatrix& R,
                     int kMax,
                     ClusterStatisticType statType);
  
  /// Destructor
  virtual ~MinClusterIlpCplex()
  {
    _env.end();
  };
  
  /// Solve
  ///
  /// @param nrThreads Number of threads
  /// @param timeLimit Time limit in seconds (for initial ILP)
  virtual bool solve(int nrThreads,
                     int timeLimit);
  
  /// Export ILP model
  ///
  /// @param filename Filename
  virtual void exportModel(const std::string& filename)
  {
    _cplex.exportModel(filename.c_str());
  }
  
protected:
  /// Initialize variables
  virtual void initVariables();
  
  /// Initialize constraints
  virtual void initConstraints();
  
  /// Initialize objective
  virtual void initObjective();
  
private:
  typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;
  typedef IloArray<IloBoolVarMatrix> IloBoolVar3Matrix;
  typedef IloArray<IloNumVarArray> IloNumVarMatrix;
  typedef IloArray<IloNumVarMatrix> IloNumVar3Matrix;
  typedef IloArray<IloNumVar3Matrix> IloNumVar4Matrix;
  typedef IloArray<IloNumVar4Matrix> IloNumVar5Matrix;
  
  /// Environment
  IloEnv _env;
  /// CPlex model
  IloModel _model;
  /// Solver
  IloCplex _cplex;
  /// y[i][t][j] = 1 iff state tree t and cluster j for SNV i
  IloBoolVar3Matrix _y;
  /// d[j][p]
  IloNumVarMatrix _d;
  /// yd[i][t][j][p] = y[i][t][j] * d[j][p]
  IloNumVar4Matrix _yd;
  /// b[j]
  IloBoolVarArray _b;
  
  /// Solution state trees
  StateTreeVector _solT;
};

#endif // MINCLUSTERILPCPLEX_H

