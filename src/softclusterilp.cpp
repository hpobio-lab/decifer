/*
 * softclusterilp.cpp
 *
 *  Created on: 15-apr-2019
 *      Author: M. El-Kebir
 */

#include "softclusterilp.h"

SoftClusterIlp::SoftClusterIlp(const ReadMatrix& R,
                               int k,
                               int nrSegments,
                               ClusterStatisticType statType,
                               double precisionBetaBin,
                               bool forceTruncal)
  : Solver(R, k, nrSegments, statType, precisionBetaBin, forceTruncal)
  , _coord()
  , _coordPi()
  , _hatG()
  , _hatPi()
  , _solT()
  , _solY()
{
}

void SoftClusterIlp::init()
{
  Solver::init();
  
  initPWLA();
  initVariables();
  initConstraints();
  initObjective();
}

void SoftClusterIlp::initPreClusteringConstraints(const IntMatrix& preClustering)
{
  for (const IntVector& preCluster : preClustering)
  {
    const int size = preCluster.size();
    for (int i = 1; i < size; ++i)
    {
      initPreClusteringConstraint(preCluster[i-1], preCluster[i]);
    }
  }
}

void SoftClusterIlp::initPWLA()
{
  assert(_nrSegments >= 2);
  
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  // compute coord
  _coord = DoubleVector(_nrSegments);
  _coord[0] = 0;
  _coord[_nrSegments - 1] = 1;
  
  const double delta = 1. / (_nrSegments - 1);
  for (int alpha = 1; alpha < _nrSegments - 1; ++alpha)
  {
    _coord[alpha] = delta * alpha;
  }
  
  _coordPi = DoubleVector(_nrSegments);
  _coordPi[0] = g_tol.epsilon();// / (m + 1);
  const double delta2 = (n) / (float) (_nrSegments - 1);
  for (int alpha = 1; alpha < _nrSegments; ++alpha)
  {
    _coordPi[alpha] = delta2 * alpha;
  }
  
  // compute hatG
  const double infeasible_log = -1e300;
  _hatG = Double5Matrix(n);
  for (int i = 0; i < n; ++i)
  {
    const int scriptT_i_size = _scriptT[i].size();
    _hatG[i] = Double4Matrix(scriptT_i_size);
    for (int t = 0; t < scriptT_i_size; ++t)
    {
      _hatG[i][t] = DoubleTensor(_k);
      for (int j = 0; j < _k; ++j)
      {
        _hatG[i][t][j] = DoubleMatrix(m);
        for (int p = 0; p < m; ++p)
        {
          _hatG[i][t][j][p] = DoubleVector(_nrSegments, 0);
          
          const int var_ip = _R.getVar(p, i);
          const int ref_ip = _R.getRef(p, i);
          
          for (int alpha = 0; alpha < _nrSegments; ++alpha)
          {
            double h = (_coord[alpha] - _numerator[i][t][p]) / _denominator[i][p];
            if (h <= 0 || !g_tol.nonZero(h))
            {
              _hatG[i][t][j][p][alpha] = infeasible_log;//log(g_tol.epsilon());
            }
            else if (h >= 1 || g_tol.less(1, h))
            {
              _hatG[i][t][j][p][alpha] = infeasible_log;//log(g_tol.epsilon());
            }
            else
            {
              double likelihood = getLogLikelihood(var_ip, ref_ip, h);
              if (likelihood < infeasible_log)//log(g_tol.epsilon()))
              {
                _hatG[i][t][j][p][alpha] = infeasible_log;//log(g_tol.epsilon());
              }
              else
              {
                _hatG[i][t][j][p][alpha] = likelihood;
              }
            }
          }
        }
      }
    }
  }
  
  // compute hatPi
  _hatPi = DoubleVector(_nrSegments, 0);
  for (int alpha = 0; alpha < _nrSegments; ++alpha)
  {
    double likelihood = log(_coordPi[alpha]);
    _hatPi[alpha] = likelihood;
  }
}
