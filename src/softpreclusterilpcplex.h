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
  
  /// Initialize cluster assignment
  void initZ(const IntVector& z)
  {
    const int m = _R.getNrSamples();
    const int n = _R.getNrCharacters();
    assert(z.size() == n);
    
    IntMatrix preClustering(_k);
    for (int i = 0; i < n; ++i)
    {
      assert(0 <= z[i] && z[i] < _k);
      _model.add(_z[i][z[i]] == 1);
      preClustering[z[i]].push_back(i);
    }
    
    for (const IntVector& preCluster : preClustering)
    {
      const int size = preCluster.size();
      for (int i = 1; i < size; ++i)
      {
        int i1 = preCluster[i-1];
        int i2 = preCluster[i];
        
        assert(_scriptT[i1].size() == _scriptT[i2].size());
        const int scriptT_size = _scriptT[i1].size();
        for (int t = 0; t < scriptT_size; ++t)
        {
          for (int j = 0; j < _k; ++j)
          {
            _model.add(_y[i1][t][j] == _y[i2][t][j]);
            
            for (int p = 0; p < m; ++p)
            {
              for (int alpha = 0; alpha < _nrSegments; ++alpha)
              {
                _model.add(_lambda[i1][t][j][p][alpha] == _lambda[i2][t][j][p][alpha]);
              }
            }
          }
        }
      }
    }
  }
  
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
