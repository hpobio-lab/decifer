/*
 * minclusterilpgurobi.h
 *
 *  Created on: 21-oct-2018
 *      Author: M. El-Kebir
 */

#ifndef MINCLUSTERILPGUROBI_H
#define MINCLUSTERILPGUROBI_H

#include "minclusterilp.h"
#include <gurobi_c++.h>

class MinClusterIlpGurobi : public MinClusterIlp
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param kMax Number of clusters
  MinClusterIlpGurobi(const ReadMatrix& R,
                      int kMax);
  
  /// Destructor
  virtual ~MinClusterIlpGurobi()
  {
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
    _model.write(filename);
  }
  
protected:
  /// Initialize variables
  virtual void initVariables();
  
  /// Initialize constraints
  virtual void initConstraints();
  
  /// Initialize objective
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
  /// b[j]
  VarArray _b;
  
  /// Solution state trees
  StateTreeVector _solT;
};

#endif // MINCLUSTERILPGUROBI_H
