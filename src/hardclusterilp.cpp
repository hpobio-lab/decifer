/*
 * hardclusterilp.cpp
 *
 *  Created on: 19-oct-2018
 *      Author: M. El-Kebir
 */

#include "hardclusterilp.h"

HardClusterIlp::HardClusterIlp(const ReadMatrix& R,
                               int k,
                               int nrSegments,
                               ClusterStatisticType statType,
                               bool forceTruncal)
  : Solver(R, k, nrSegments, statType, forceTruncal)
  , _hatN()
  , _z()
  , _solT()
{
}

void HardClusterIlp::init()
{
  Solver::init();
  
  initPWLA();
  initVariables();
  initConstraints();
  initObjective();
}

void HardClusterIlp::initPWLA()
{
  Solver::initPWLA();
  
  const int n = _R.getNrCharacters();
  
  _hatN = DoubleVector(_nrSegments, 0);
  _z = DoubleVector(_nrSegments, 0);
  
  for (int l = 0; l < _nrSegments; ++l)
  {
    _z[l] = (double(n) / (_nrSegments - 1)) * l;
    _hatN[l] = log(_z[l]);
  }
  
  _hatN[0] = log(g_tol.epsilon());
}
