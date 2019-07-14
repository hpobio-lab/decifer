/*
 * emcplexpre.h
 *
 *  Created on: 2-jul-2019
 *      Author: M. El-Kebir
 */

#ifndef EMCPLEXPRE_H
#define EMCPLEXPRE_H

#include "empre.h"
#include "hardclusterilpcplex.h"
#include "clusterilpcplex.h"
#include "softclusterlpcplexpre.h"
#include "softclusterlpcplex.h"
#include <ilcplex/ilocplex.h>

class EMCplexPre : public EMPre
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param preClustering Pre clustering
  /// @param k Number of clusters
  /// @param nrSegments Number of segments for piecewise linear approximation
  /// @param statType Summary statistic to use for clustering
  /// @param precisionBetaBin Precision parameter for beta binomial
  /// @param forceTruncal Force the presence of a dominant truncal cluster
  EMCplexPre(const ReadMatrix& R,
             const IntMatrix& preClustering,
             int k,
             int nrSegments,
             ClusterStatisticType statType,
             double precisionBetaBin,
             bool forceTruncal);
  
  virtual ~EMCplexPre()
  {
    _env.end();
  }
  
protected:
  void initM();
  
  void initE();
  
  bool stepM(int nrThreads,
             bool verbose);
  
  bool stepE(int nrThreads,
             bool verbose);
  
  virtual std::unique_ptr<SoftClusterIlp> createSoftClusterIlpSolver(const ReadMatrix& R,
                                                                     const IntMatrix& preClustering)
  {
    return std::unique_ptr<SoftClusterIlp>(new SoftClusterLpCplexPre(R,
                                                                     preClustering,
                                                                     _k,
                                                                     3, //log(_nrSegments) / log(2),
                                                                     _statType,
                                                                     _precisionBetaBin,
                                                                     _forceTruncal,
                                                                     true));
  }
  
  typedef IloArray<IloNumVarArray> IloNumVarMatrix;
  typedef IloArray<IloNumVarMatrix> IloNumVar3Matrix;
  typedef IloArray<IloNumVar3Matrix> IloNumVar4Matrix;
  typedef IloArray<IloNumVar4Matrix> IloNumVar5Matrix;
  
private:
  /// Environment
  IloEnv _env;
  /// CPlex model E step (identify y and pi given d)
  IloModel _modelE;
  /// Solver for E step (identify y and pi given d)
  IloCplex _cplexE;
  /// CPlex model M step (identify d given y and pi)
  IloModel _modelM;
  /// Solver for M step (identify d given y and pi)
  IloCplex _cplexM;
  /// lambda[j][p][alpha]
  IloNumVar3Matrix _lambda;
  /// gamma[j][alpha]
  IloNumVarMatrix _gamma;
  /// d[j][p]
  IloNumVarMatrix _d;
  /// yy[ii][t][j]
  IloNumVar3Matrix _y;
  /// pi[j]
  IloNumVarArray _pi;
};

#endif // EMCPLEXPRE_H
