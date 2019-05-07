/*
 * emcplex.cpp
 *
 *  Created on: 16-oct-2018
 *      Author: M. El-Kebir
 */

#include "emcplex.h"

EMCplex::EMCplex(const ReadMatrix& R,
                 int k,
                 int nrSegments,
                 ClusterStatisticType statType,
                 double precisionBetaBin,
                 bool forceTruncal)
  : EM(R, k, nrSegments, statType, precisionBetaBin, forceTruncal)
  , _env()
  , _model(_env)
  , _cplex(_model)
  , _lambdaGRB()
  , _fGRB()
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
  _lambdaGRB = IloNumVar5Matrix(_env, n);
  for (int i = 0; i < n; ++i)
  {
    _lambdaGRB[i] = IloNumVar4Matrix(_env, _scriptT[i].size());
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      _lambdaGRB[i][t] = IloNumVar3Matrix(_env, _k);
      for (int j = 0; j < _k; ++j)
      {
        if (g_tol.nonZero(_gamma[i][t][j]))
        {
          _lambdaGRB[i][t][j] = IloNumVarMatrix(_env, m);
          for (int p = 0; p < m; ++p)
          {
            _lambdaGRB[i][t][j][p] = IloNumVarArray(_env, _nrSegments);
            for (int l = 0; l < _nrSegments; ++l)
            {
              snprintf(buf, 1024, "lambda;%d;%d;%d;%d;%d", i, t, j, p, l);
              _lambdaGRB[i][t][j][p][l] = IloNumVar(_env, 0, 1, buf);
            }
          }
        }
      }
    }
  }
  
  _fGRB = IloNumVarMatrix(_env, _k);
  for (int j = 0; j < _k; ++j)
  {
    _fGRB[j] = IloNumVarArray(_env, m);
    for (int p = 0; p < m; ++p)
    {
      snprintf(buf, 1024, "f;%d;%d", j, p);
      _fGRB[j][p] = IloNumVar(_env, 0, 1, buf);
//      _model.add(_fGRB[j][p] == _solF[j][p]);
    }
  }
  
  // initialize constraints
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < m; ++p)
    {
      _model.add(_dOverallLB[j][p] <= _fGRB[j][p]);
      _model.add(_fGRB[j][p] <= _dOverallUB[j][p]);
    }
  }
  
  IloExpr sum(_env), sum2(_env);
  for (int i = 0; i < n; ++i)
  {
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        if (g_tol.nonZero(_gamma[i][t][j]))
        {
          for (int p = 0; p < m; ++p)
          {
            for (int l = 0; l < _nrSegments; ++l)
            {
              sum += _lambdaGRB[i][t][j][p][l];
              sum2 += _lambdaGRB[i][t][j][p][l] * _x[i][t][j][p][l];
            }
            _model.add(sum == 1);
            _model.add(_fGRB[j][p] == sum2);
            sum.clear();
            sum2.clear();
          }
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
          assert(!isnan(_gamma[i][t][j]));
          for (int p = 0; p < m; ++p)
          {
            for (int l = 0; l < _nrSegments; ++l)
            {
              assert(!isnan(_hatG[i][t][j][p][l]));
              obj += _gamma[i][t][j] * _hatG[i][t][j][p][l] * _lambdaGRB[i][t][j][p][l];
            }
          }
        }
      }
    }
  }
  
  if (_forceTruncal)
  {
    // Cluster 0 has largest DCF in all samples
    for (int j = 1; j < _k; ++j)
    {
      for (int p = 0; p < m; ++p)
      {
        _model.add(_fGRB[j][p] <= _fGRB[0][p]);
      }
    }
  }
  
  _model.add(IloMaximize(_env, obj));
  obj.end();
  sum.end();
  sum2.end();
}

bool EMCplex::stepM(int nrThreads,
                    bool verbose)
{
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
  if (_cplex.getStatus() == IloAlgorithm::Infeasible)
  {
    return false;
  }
  
  const int m = _R.getNrSamples();
  const int n = _R.getNrCharacters();
  
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < m; ++p)
    {
      _solD[j][p] = _cplex.getValue(_fGRB[j][p]);
      
      if (_solD[j][p] == g_tol.epsilon())
      {
        _solD[j][p] = 0;
      }
    }
  }
  
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
  
  return true;
}
