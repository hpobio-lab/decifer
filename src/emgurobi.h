/*
 * emgurobi.h
 *
 *  Created on: 16-oct-2018
 *      Author: M. El-Kebir
 */

#ifndef EMGUROBI_H
#define EMGUROBI_H

#include "em.h"
#include "clusterilpgurobi.h"
#include "hardclusterilpgurobi.h"
#include <gurobi_c++.h>

class EMGurobi : public EM
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param k Number of clusters
  /// @param nrSegments Number of segments for piecewise linear approximation
  EMGurobi(const ReadMatrix& R,
           int k,
           int nrSegments);
  
  virtual ~EMGurobi()
  {
    delete _pModel;
  }
  
protected:
  void initPWLA();
  
  bool stepM(int nrThreads);
  
  virtual std::unique_ptr<HardClusterIlp> createHardClusterIlpSolver(const ReadMatrix& R)
  {
    return std::unique_ptr<HardClusterIlp>(new HardClusterIlpGurobi(R, _k, _nrSegments));
  }
  
  virtual std::unique_ptr<ClusterIlp> createClusterIlpSolver(const ReadMatrix& R,
                                                             double alpha)
  {
    return std::unique_ptr<ClusterIlp>(new ClusterIlpGurobi(R, _k, alpha));
  }
    
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

private:
  /// Gurobi environment
  GRBEnv _env;
  /// Gurobi model
  GRBModel* _pModel;
  /// lambda[i][t][j][p][l]
  Var5Matrix _lambdaGRB;
  /// fGRB[j][p]
  VarMatrix _fGRB;
};

#endif // EMGUROBI_H
