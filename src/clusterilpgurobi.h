/*
 * clusterilpgurobi.h
 *
 *  Created on: 16-oct-2018
 *      Author: M. El-Kebir
 */

#ifndef CLUSTERILPGUROBI_H
#define CLUSTERILPGUROBI_H

#include "clusterilp.h"
#include <gurobi_c++.h>

class ClusterIlpGurobi : public ClusterIlp
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param k Number of clusters
  /// @param alpha Confidence interval
  ClusterIlpGurobi(const ReadMatrix& R,
                   int k,
                   double alpha);

  /// Export ILP
  ///
  /// @param filename Filename
  virtual void exportModel(const std::string& filename)
  {
    _model.write(filename);
  }
  
  /// Solve
  ///
  /// @param nrThreads Number of threads
  /// @param timeLimit Time limit in seconds (unlimited if -1)
  /// @param verbose Verbose
  virtual bool solve(int nrThreads,
                     int timeLimit,
                     bool verbose,
                     int memoryLimit);
  
  /// Initialize hot start
  ///
  /// @param y Clustering and state tree assignment
  void initHotStart(const BoolTensor& y);
  
protected:
  /// Initialize variables
  virtual void initVariables();
  
  /// Initialize constraints
  virtual void initConstraints();
  
private:
  /// Gurobi variable array
  typedef std::vector<GRBVar> VarArray;
  /// Gurobi variable matrix
  typedef std::vector<VarArray> VarMatrix;
  /// Gurobi variable 3D matrix
  typedef std::vector<VarMatrix> Var3Matrix;
  /// Gurobi variable 4D matrix
  typedef std::vector<Var3Matrix> Var4Matrix;
  
  /// Gurobi environment
  GRBEnv _env;
  /// Gurobi model
  GRBModel _model;
  /// y[i][t][j] = 1 iff state tree t and cluster j for SNV i
  Var3Matrix _y;
  /// d[j][p]
  VarMatrix _d;
  /// yd[i][t][j][p] = y[i][t][j] * d[j][p]
  Var4Matrix _yd;
  /// w[i][t][j][p] = |yd[i][t][j][p] - dcfExp[i][t][p] * y[i][t][j]|
  Var4Matrix _w;
};

#endif // CLUSTERILPGUROBI_H
