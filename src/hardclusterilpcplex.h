/*
 * hardclusterilpcplex.h
 *
 *  Created on: 19-oct-2018
 *      Author: M. El-Kebir
 */

#ifndef HARDCLUSTERILPCPLEX_H
#define HARDCLUSTERILPCPLEX_H

#include "hardclusterilp.h"
#include <ilcplex/ilocplex.h>

class HardClusterIlpCplex : public HardClusterIlp
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
  /// @param includePi Include pi in optimization
  HardClusterIlpCplex(const ReadMatrix& R,
                      int k,
                      int nrSegments,
                      ClusterStatisticType statType,
                      double precisionBetaBin,
                      bool forceTruncal,
                      bool includePi);
  
  /// Initialize pre clustering constraints
  ///
  /// @param preClustering Pre clustering
  virtual void initPreClusteringConstraints(const IntMatrix& preClustering);

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
  
  /// Destructor
  virtual ~HardClusterIlpCplex()
  {
    _env.end();
  }
  
protected:
  /// Initialize variables
  virtual void initVariables();
  
  /// Initialize constraints
  virtual void initConstraints();
  
  /// Initialize objective function
  virtual void initObjective();
  
  /// Initialize pre clustering constraint
  ///
  /// @param i1 SNV
  /// @param i2 SNV
  virtual void initPreClusteringConstraint(int i1, int i2);
  
private:
  typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;
  typedef IloArray<IloBoolVarMatrix> IloBoolVar3Matrix;
  typedef IloArray<IloNumVarArray> IloNumVarMatrix;
  typedef IloArray<IloNumVarMatrix> IloNumVar3Matrix;
  typedef IloArray<IloNumVar3Matrix> IloNumVar4Matrix;
  typedef IloArray<IloNumVar4Matrix> IloNumVar5Matrix;
  
  /// Include pi in optmization
  const bool _includePi;
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
  /// pi[j]
  IloNumVarArray _pi;
  /// gamma[j][alpha]
  IloNumVarMatrix _gamma;
  /// lambda[i][t][j][p][alpha]
  IloNumVar5Matrix _lambda;
  /// llambda[j][p][alpha]
  IloNumVar3Matrix _llambda;
};

#endif // HARDCLUSTERILPCPLEX_H
