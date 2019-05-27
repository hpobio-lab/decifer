/*
 * emcplex.cpp
 *
 *  Created on: 16-oct-2018
 *      Author: M. El-Kebir
 */

#include "emcplex.h"

EMCplex::EMCplex(const ReadMatrix& R,
                 int k,
                 int nrSegmentBits,
                 ClusterStatisticType statType,
                 double precisionBetaBin,
                 bool forceTruncal)
  : EM(R, k, nrSegmentBits, statType, precisionBetaBin, forceTruncal)
  , _env()
  , _model(_env)
  , _cplex(_model)
  , _lambda()
  , _d()
{
}

void EMCplex::initPWLA()
{
  EM::initPWLA();
  
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
//  _cplex.clearModel();
  _env.end();
  _env = IloEnv();
  _model = IloModel(_env);
  _cplex = IloCplex(_model);
  
  char buf[1024];
  _lambda = IloNumVar3Matrix(_env, _k);
  for (int j = 0; j < _k; ++j)
  {
    _lambda[j] = IloNumVarMatrix(_env, m);
    for (int p = 0; p < m; ++p)
    {
      _lambda[j][p] = IloNumVarArray(_env, _nrSegments);
      for (int l = 0; l < _nrSegments; ++l)
      {
        snprintf(buf, 1024, "lambda;%d;%d;%d", j, p, l);
        _lambda[j][p][l] = IloNumVar(_env, 0, 1, buf);
      }
    }
  }
  
  _d = IloNumVarMatrix(_env, _k);
  for (int j = 0; j < _k; ++j)
  {
    _d[j] = IloNumVarArray(_env, m);
    for (int p = 0; p < m; ++p)
    {
      snprintf(buf, 1024, "f;%d;%d", j, p);
      _d[j][p] = IloNumVar(_env, 0, 1, buf);
//      _model.add(_fGRB[j][p] == _solF[j][p]);
    }
  }
  
  // initialize constraints
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < m; ++p)
    {
      _model.add(_dOverallLB[j][p] <= _d[j][p]);
      _model.add(_d[j][p] <= _dOverallUB[j][p]);
    }
  }
  
  // compute coord
  DoubleVector coord(_nrSegments);
  coord[0] = 0;
  coord[_nrSegments - 1] = 1;
  
  const double delta = 1. / (_nrSegments - 1);
  for (int alpha = 1; alpha < _nrSegments - 1; ++alpha)
  {
    coord[alpha] = delta * alpha;
  }
  
  IloExpr sum(_env), sum2(_env);
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < m; ++p)
    {
//      _model.add(_d[j][p] == _solD[j][p]);
      for (int l = 0; l < _nrSegments; ++l)
      {
        sum += _lambda[j][p][l];
        sum2 += _lambda[j][p][l] * coord[l];
      }
      _model.add(sum == 1);
      _model.add(_d[j][p] == sum2);
      sum.clear();
      sum2.clear();
    }
  }
  
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < m; ++p)
    {
      for (int alpha = 0; alpha < _nrSegments-1; ++alpha)
      {
        if (g_tol.less(coord[alpha], _dOverallLB[j][p]) && g_tol.less(coord[alpha + 1], _dOverallLB[j][p]))
        {
          _model.add(_lambda[j][p][alpha] == 0);
        }
      }
      
      for (int alpha = 1; alpha < _nrSegments; ++alpha)
      {
        if (g_tol.less(_dOverallUB[j][p], coord[alpha-1]) && g_tol.less(_dOverallUB[j][p], coord[alpha]))
        {
          _model.add(_lambda[j][p][alpha] == 0);
        }
      }
    }
  }
  
  IloExpr obj(_env);
  for (int i = 0; i < n; ++i)
  {
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        if (g_tol.nonZero(_gamma[i][t][j]))
        {
          obj += _gamma[i][t][j] * log(_solPi[j]) - _gamma[i][t][j] * log(_scriptT[i].size());
          assert(!isnan(_gamma[i][t][j]));
          for (int p = 0; p < m; ++p)
          {
            for (int l = 0; l < _nrSegments; ++l)
            {
              assert(!isnan(_hatG[i][t][j][p][l]));
              obj += _gamma[i][t][j] * _hatG[i][t][j][p][l] * _lambda[j][p][l];
            }
          }
        }
      }
    }
  }
  
//  if (_forceTruncal)
//  {
//    // Cluster 0 has largest DCF in all samples
//    for (int j = 1; j < _k; ++j)
//    {
//      for (int p = 0; p < m; ++p)
//      {
//        _model.add(_d[j][p] <= _d[0][p]);
//      }
//    }
//  }
  
  _model.add(IloMaximize(_env, obj));
  obj.end();
  sum.end();
  sum2.end();
}

bool EMCplex::stepM(int nrThreads,
                    bool verbose)
{
  const int m = _R.getNrSamples();
  const int n = _R.getNrCharacters();
  
  {
    double denominator = 0;
    for (int i = 0; i < n; ++i)
    {
      for (int t = 0; t < _scriptT[i].size(); ++t)
      {
        for (int j = 0; j < _k; ++j)
        {
          denominator += _gamma[i][t][j];
        }
      }
    }
    
    for (int j = 0; j < _k; ++j)
    {
      _solPi[j] = 0;
      for (int i = 0; i < n; ++i)
      {
        for (int t = 0; t < _scriptT[i].size(); ++t)
        {
          _solPi[j] += _gamma[i][t][j];
        }
      }
      _solPi[j] /= denominator;
    }
  }
  
  initPWLA();
  if (nrThreads > 0)
  {
    _cplex.setParam(IloCplex::Threads, nrThreads);
  }
  
  if (!verbose)
  {
    _cplex.setOut(_env.getNullStream());
    _cplex.setError(_env.getNullStream());
    _cplex.setWarning(_env.getNullStream());
  }
  
//  _cplex.exportModel("/tmp/infeasible.lp");
  _cplex.solve();
//  std::cerr << "Likelihood: " << _cplex.getObjValue() << std::endl;
  if (_cplex.getStatus() == IloAlgorithm::Infeasible)
  {
    return false;
  }
  
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < m; ++p)
    {
      _solD[j][p] = _cplex.getValue(_d[j][p]);
      
      if (_solD[j][p] == g_tol.epsilon())
      {
        _solD[j][p] = 0;
      }
      
//      for (int alpha = 0; alpha < _nrSegments; ++alpha)
//      {
//        if (g_tol.nonZero(_cplex.getValue(_lambda[j][p][alpha])))
//        {
//          std::cout << _lambda[j][p][alpha] << " = " << _cplex.getValue(_lambda[j][p][alpha]) << std::endl;
//        }
//      }
    }
  }
  
//  for (int i = 0; i < n; ++i)
//  {
//    for (int t = 0; t < _scriptT[i].size(); ++t)
//    {
//      for (int j = 0; j < _k; ++j)
//      {
//        if (g_tol.nonZero(_gamma[i][t][j]))
//        {
//          double sum = _gamma[i][t][j] * log(_solPi[j]) - _gamma[i][t][j] * log(_scriptT[i].size());
//          for (int p = 0; p < m; ++p)
//          {
//            for (int l = 0; l < _nrSegments; ++l)
//            {
//              assert(!isnan(_hatG[i][t][j][p][l]));
//              sum += _gamma[i][t][j] * _hatG[i][t][j][p][l] * _cplex.getValue(_lambda[j][p][l]);
//            }
//          }
//
//          std::cout << i << "," << t << "," << j << " : " << _gamma[i][t][j] << " -- " << getLogLikelihoodGamma(i) << " -- " << sum << std::endl;
//        }
//      }
//    }
//  }
  
  return true;
}
