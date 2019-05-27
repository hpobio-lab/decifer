/*
 * softclusterlpcplexext.h
 *
 *  Created on: 17-apr-2019
 *      Author: M. El-Kebir
 */

#ifndef SOFTCLUSTERLPCPLEXEXT_H
#define SOFTCLUSTERLPCPLEXEXT_H

#include "softclusterilp.h"
#include <ilcplex/ilocplex.h>

class SoftClusterLpCplexExt : public SoftClusterIlp
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param k Number of clusters
  /// @param nrSegments Number of bits to model segments for piecewise linear approximation
  /// @param statType Summary statistic to use for clustering
  /// @param precisionBetaBin Precision parameter for beta binomial
  /// @param forceTruncal Force the presence of a dominant truncal cluster
  /// @param pi Include pi in optimization
  SoftClusterLpCplexExt(const ReadMatrix& R,
                        int k,
                        int nrSegmentBits,
                        ClusterStatisticType statType,
                        double precisionBetaBin,
                        bool forceTruncal,
                        bool pi);
  
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
  
  virtual ~SoftClusterLpCplexExt()
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
  typedef IloArray<IloBoolVar3Matrix> IloBoolVar4Matrix;
  typedef IloArray<IloNumVarArray> IloNumVarMatrix;
  typedef IloArray<IloNumVarMatrix> IloNumVar3Matrix;
  typedef IloArray<IloNumVar3Matrix> IloNumVar4Matrix;
  typedef IloArray<IloNumVar4Matrix> IloNumVar5Matrix;
  typedef IloArray<IloNumVar5Matrix> IloNumVar6Matrix;
  
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
  /// bggamma[j][beta]
  IloBoolVarMatrix _bgamma;
  /// lambda[j][p][alpha]
  IloNumVar3Matrix _lambda;
  /// bllambda[j][p][beta]
  IloBoolVar3Matrix _blambda;
  /// sigma[i][t][j][alpha]
  IloNumVar4Matrix _sigma;
  /// bsigma[i][t][j][beta]
  IloBoolVar4Matrix _bsigma;
  /// ssigma[i][t][j][alpha]
  IloNumVar4Matrix _ssigma;
  /// l[i][t][j]
  IloNumVar3Matrix _l;
  /// Include pi in optmization
  const bool _includePi;
};

#endif // SOFTCLUSTERLPCPLEXEXT_H
