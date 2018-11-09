/*
 * minclusterilpcplex.cpp
 *
 *  Created on: 21-oct-2018
 *      Author: M. El-Kebir
 */

#include "minclusterilpcplex.h"

MinClusterIlpCplex::MinClusterIlpCplex(const ReadMatrix& R,
                                       int kMax)
  : MinClusterIlp(R, kMax)
  , _env()
  , _model(_env)
  , _cplex(_model)
  , _y()
  , _d()
  , _yd()
  , _b()
  , _solT()
{
}

void MinClusterIlpCplex::initVariables()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  char buf[1024];
  
  _y = IloBoolVar3Matrix(_env, n);
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptT[i].size();
    
    _y[i] = IloBoolVarMatrix(_env, size_T_i);
    for (int t = 0; t < size_T_i; ++t)
    {
      _y[i][t] = IloBoolVarArray(_env, _k);
      for (int j = 0; j < _k; ++j)
      {
        snprintf(buf, 1024, "y;%d;%d;%d", i, t, j);
        _y[i][t][j] = IloBoolVar(_env, buf);
      }
    }
  }
  
  _d = IloNumVarMatrix(_env, _k);
  for (int j = 0; j < _k; ++j)
  {
    _d[j] = IloNumVarArray(_env, m);
    for (int p = 0; p < m; ++p)
    {
      snprintf(buf, 1024, "d;%d;%d", j, p);
      _d[j][p] = IloNumVar(_env, 0, 1, buf);
    }
  }
  
  _yd = IloNumVar4Matrix(_env, n);
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptT[i].size();
    
    _yd[i] = IloNumVar3Matrix(_env, size_T_i);
    for (int t = 0; t < size_T_i; ++t)
    {
      _yd[i][t] = IloNumVarMatrix(_env, _k);
      for (int j = 0; j < _k; ++j)
      {
        _yd[i][t][j] = IloNumVarArray(_env, m);
        for (int p = 0; p < m; ++p)
        {
          snprintf(buf, 1024, "yd;%d;%d;%d;%d", i, t, j, p);
          _yd[i][t][j][p] = IloNumVar(_env, 0, 1, buf);
        }
      }
    }
  }
  
  _b = IloBoolVarArray(_env, _k);
  for (int j = 0; j < _k; ++j)
  {
    snprintf(buf, 1024, "b;%d", j);
    _b[j] = IloBoolVar(_env, buf);
  }
}

void MinClusterIlpCplex::initConstraints()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  IloExpr sum(_env), sum2(_env);
  
  // one state tree and cluster per SNV
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptT[i].size();
    for (int t = 0; t < size_T_i; ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        sum += _y[i][t][j];
      }
    }
    
    _model.add(sum == 1);
    sum.clear();
  }
  
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptT[i].size();
    for (int t = 0; t < size_T_i; ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        for (int p = 0; p < m; ++p)
        {
          _model.add(_y[i][t][j] * _dLB[i][t][p] <= _yd[i][t][j][p]);
          _model.add(_yd[i][t][j][p] <= _dUB[i][t][p] * _y[i][t][j]);
        }
      }
    }
  }
  
  for (int i = 0; i < n; ++i)
  {
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        for (int p = 0; p < m; ++p)
        {
          _model.add(_yd[i][t][j][p] <= _y[i][t][j]);
          _model.add(_yd[i][t][j][p] <= _d[j][p]);
          _model.add(_yd[i][t][j][p] >= _y[i][t][j] + _d[j][p] - 1);
        }
      }
    }
  }
  
  for (int i = 0; i < n; ++i)
  {
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      sum += _y[i][t][0];
    }
  }
  
  for (int j = 1; j < _k; ++j)
  {
    for (int i = 0; i < n; ++i)
    {
      for (int t = 0; t < _scriptT[i].size(); ++t)
      {
        sum2 += _y[i][t][j];
      }
    }
    _model.add(sum >= sum2);
    sum = sum2;
    sum2.clear();
  }
  sum.clear();
  
  for (int j = 0; j < _k; ++j)
  {
    for (int i = 0; i < n; ++i)
    {
      for (int t = 0; t < _scriptT[i].size(); ++t)
      {
        _model.add(_b[j] >= _y[i][t][j]);
        sum += _y[i][t][j];
      }
    }
    _model.add(_b[j] <= sum);
    sum.clear();
  }
  
  sum.end();
  sum2.end();
}

void MinClusterIlpCplex::initObjective()
{
  IloExpr sum(_env);
  
  for (int j = 0; j < _k; ++j)
  {
    sum += _b[j];
  }
  
  _model.add(IloMinimize(_env, sum));
  sum.end();
}

bool MinClusterIlpCplex::solve(int nrThreads, int timeLimit)
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  if (nrThreads > 0)
  {
    _cplex.setParam(IloCplex::Threads, nrThreads);
  }
  if (timeLimit > 0)
  {
    _cplex.setParam(IloCplex::TiLim, timeLimit);
  }
  _cplex.setOut(_env.getNullStream());
  _cplex.setError(_env.getNullStream());
  _cplex.setWarning(_env.getNullStream());
  
  _cplex.solve();
  if (_cplex.getStatus() == IloAlgorithm::Infeasible)
  {
    return false;
  }
  
  _solD = DoubleMatrix(_k, DoubleVector(m, 0));
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < m; ++p)
    {
      _solD[j][p] = _cplex.getValue(_d[j][p]);
    }
  }
  
  _solT.clear();
  _solZ.clear();
  _solPi = DoubleVector(_k);
  for (int i = 0; i < n; ++i)
  {
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        if (_cplex.getValue(_y[i][t][j]) >= 0.4)
        {
          _solT.emplace_back(convertToStateTreeFromDCF(_scriptT[i][t],
                                                       _solD[j], i));
          _solZ.push_back(j);
          ++_solPi[j];
          for (int p = 0; p < m; ++p)
          {
            _solT.back().resetMixtureProportions(p, _solD[j][p]);
          }
        }
      }
    }
  }
  
  for (int j = 0; j < _k; ++j)
  {
    _solPi[j] /= n;
  }
  
  _solMinK = 0;
  for (int j = 0; j < _k; ++j)
  {
    if (_cplex.getValue(_b[j]))
    {
      ++_solMinK;
    }
  }
  
  updateLogLikelihood();
  
  return true;
}

