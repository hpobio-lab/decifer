/*
 * hardclusterilpgurobi.h
 *
 *  Created on: 20-oct-2018
 *      Author: M. El-Kebir
 */

#ifndef HARDCLUSTERILPGUROBI_H
#define HARDCLUSTERILPGUROBI_H

#include "hardclusterilp.h"
#include <gurobi_c++.h>

class HardClusterIlpGurobi : public HardClusterIlp
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param k Number of clusters
  /// @param nrSegments Number of segments for piecewise linear approximation
  HardClusterIlpGurobi(const ReadMatrix& R,
                       int k,
                       int nrSegments);

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
  /// @param memoryLimit Memory limit
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
  
  /// Initialize objective function
  virtual void initObjective();
  
private:
  /// Gurobi variable array
  typedef std::vector<GRBVar> VarArray;
  /// Gurobi variable matrix
  typedef std::vector<VarArray> VarMatrix;
  /// Gurobi variable 3D matrix
  typedef std::vector<VarMatrix> Var3Matrix;
  /// Gurobi variable 4D matrix
  typedef std::vector<Var3Matrix> Var4Matrix;
  /// Gurobi variable 5D matrix
  typedef std::vector<Var4Matrix> Var5Matrix;
  
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
  /// gamma[i][t][j][l]
  Var4Matrix _gamma;
  /// lambda[i][t][j][p][l]
  Var5Matrix _lambda;
};

#endif // HARDCLUSTERILPGUROBI_H
