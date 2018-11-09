/*
 * clusterilpcplex.h
 *
 *  Created on: 16-oct-2018
 *      Author: M. El-Kebir
 */

#ifndef CLUSTERILPCPLEX_H
#define CLUSTERILPCPLEX_H

#include "clusterilp.h"
#include "clusterilpcplexcallback.h"
#include <ilcplex/ilocplex.h>

class ClusterIlpCplex : public ClusterIlp
{
public:
  ClusterIlpCplex(const ReadMatrix& R,
                  int k,
                  double alpha);
  
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
  
  virtual ~ClusterIlpCplex()
  {
    _env.end();
  }
  
protected:
  virtual void initVariables();
  virtual void initConstraints();
  
private:
  typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;
  typedef IloArray<IloBoolVarMatrix> IloBoolVar3Matrix;
  typedef IloArray<IloNumVarArray> IloNumVarMatrix;
  typedef IloArray<IloNumVarMatrix> IloNumVar3Matrix;
  typedef IloArray<IloNumVar3Matrix> IloNumVar4Matrix;
  
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
  /// w[i][t][j][p] = |ff[i][t][j][p] - dcfExp[i][t][p] * y[i][t][j]|
  IloNumVar4Matrix _w;
};

#endif // CLUSTERILPCPLEX_H
