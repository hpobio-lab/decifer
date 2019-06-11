/*
 * softpreclusterilpcplex.h
 *
 *  Created on: 11-jun-2019
 *      Author: M. El-Kebir
 */

#ifndef SOFTPRECLUSTERILPCPLEX_H
#define SOFTPRECLUSTERILPCPLEX_H

#include "softclusterilp.h"
#include <ilcplex/ilocplex.h>

class SoftPreClusterIlpCplex : public SoftClusterIlp
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
  SoftPreClusterIlpCplex(const ReadMatrix& R,
                         int k,
                         int nrSegmentBits,
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
  virtual ~SoftPreClusterIlpCplex()
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
  /// Initialize pre clustering constraint
  ///
  /// @param i1 SNV
  /// @param i2 SNV
  virtual void initPreClusteringConstraint(int i1, int i2)
  {
  }
  
  /// Initialize variables
  virtual void initVariables();
  
  /// Initialize constraints
  virtual void initConstraints();
  
  /// Initialize objective function
  virtual void initObjective();
  
  /// Convert binary number to gray
  unsigned int binaryToGray(unsigned int num)
  {
    num = _nrSegments - num - 2;
    unsigned int val = num ^ (num >> 1);
    
    // reverse
    unsigned int reverse_val = 0;
    for (int i = 0; i < _nrSegmentBits; ++i)
    {
      unsigned int bit = (val & (1 << (_nrSegmentBits - i - 1)));
      if (bit)
      {
        reverse_val |= (1 << i);
      }
    }
    
    return reverse_val;
  }
  
private:
  typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;
  typedef IloArray<IloBoolVarMatrix> IloBoolVar3Matrix;
  typedef IloArray<IloNumVarArray> IloNumVarMatrix;
  typedef IloArray<IloNumVarMatrix> IloNumVar3Matrix;
  typedef IloArray<IloNumVar3Matrix> IloNumVar4Matrix;
  typedef IloArray<IloNumVar4Matrix> IloNumVar5Matrix;
  typedef IloArray<IloNumVar5Matrix> IloNumVar6Matrix;
  
  /// Include pi in optmization
  const bool _includePi;
  /// Number of bits to model segments for piecewise linear approximation
  const int _nrSegmentBits;
  /// Environment
  IloEnv _env;
  /// CPlex model
  IloModel _model;
  /// Solver
  IloCplex _cplex;
  /// y[i][t][j] -- posterior probability of state tree t and cluster j for SNV i
  IloNumVar3Matrix _y;
  /// d[j][p]
  IloNumVarMatrix _d;
  /// pi[j]
  IloNumVarArray _pi;
  /// gamma[j][alpha]
  IloNumVarMatrix _gamma;
  /// ggamma[j][alpha]
  IloNumVarMatrix _ggamma;
  /// bggamma[j][beta]
  IloBoolVarMatrix _bggamma;
  /// lambda[i][t][j][p][alpha]
  IloNumVar5Matrix _lambda;
  /// llambda[j][p][alpha]
  IloNumVar3Matrix _llambda;
  /// bllambda[j][p][beta]
  IloBoolVar3Matrix _bllambda;
  /// w[j][t] = 1 iff tree t of cluster j is used
  IloNumVarMatrix _w;
  /// z[i][j] = 1 iff SNV i is assigned to cluster j
  IloBoolVarMatrix _z;
};

#endif // SOFTPRECLUSTERILPCPLEX_H