/*
 * hardpreclusterilpcplex.h
 *
 *  Created on: 7-jun-2019
 *      Author: M. El-Kebir
 */

#ifndef HARDPRECLUSTERILPCPLEX_H
#define HARDPRECLUSTERILPCPLEX_H

#include "hardclusterilp.h"
#include <ilcplex/ilocplex.h>

class HardPreClusterIlpCplex : public HardClusterIlp
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param k Number of clusters
  /// @param nrSegments Number of segments for piecewise linear approximation
  /// @param statType Summary statistic to use for clustering
  /// @param precisionBetaBin Precision parameter for beta binomial
  /// @param includePi Include pi in optimization
  HardPreClusterIlpCplex(const ReadMatrix& R,
                         int k,
                         int nrSegments,
                         ClusterStatisticType statType,
                         double precisionBetaBin,
                         bool includePi);

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
  void initHotStart(const BoolTensor& y)
  {
  }
  
  /// Destructor
  virtual ~HardPreClusterIlpCplex()
  {
    _env.end();
  }
  
  /// Cluster assignment
  const IntVector& getSolZ() const
  {
    return _solZ;
  }
  
  /// Cluster assignment
  ///
  /// @param i SNV
  const int getSolZ(int i) const
  {
    assert(0 <= i && i < _R.getNrCharacters());
    return _solZ[i];
  }
  
protected:
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
  
  /// Include pi in optmization
  const bool _includePi;
  /// Environment
  IloEnv _env;
  /// CPlex model
  IloModel _model;
  /// Solver
  IloCplex _cplex;
  /// y[i][t][j] = 1 iff state tree t and cluster j for SNV i
  IloNumVar3Matrix _y;
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
  /// w[j][t]
  IloBoolVarMatrix _w;
  /// z[i][j]
  IloBoolVarMatrix _z;
  /// Cluster assignment
  IntVector _solZ;
};

#endif // HARDPRECLUSTERILPCPLEX_H
