/*
 * minclusterilpgurobi.cpp
 *
 *  Created on: 21-oct-2018
 *      Author: M. El-Kebir
 */

#include "minclusterilpgurobi.h"

MinClusterIlpGurobi::MinClusterIlpGurobi(const ReadMatrix& R,
                                         int kMax)
  : MinClusterIlp(R, kMax)
  , _env()
  , _model(_env)
  , _y()
  , _d()
  , _yd()
  , _b()
  , _solT()
{
}

void MinClusterIlpGurobi::initVariables()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  char buf[1024];
  
  _y = Var3Matrix(n);
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptT[i].size();
    
    _y[i] = VarMatrix(size_T_i);
    for (int t = 0; t < size_T_i; ++t)
    {
      _y[i][t] = VarArray(_k);
      for (int j = 0; j < _k; ++j)
      {
        snprintf(buf, 1024, "y;%d;%d;%d", i, t, j);
        _y[i][t][j] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
      }
    }
  }
  
  _d = VarMatrix(_k);
  for (int j = 0; j < _k; ++j)
  {
    _d[j] = VarArray(m);
    for (int p = 0; p < m; ++p)
    {
      snprintf(buf, 1024, "d;%d;%d", j, p);
      _d[j][p] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
    }
  }
  
  _yd = Var4Matrix(n);
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptT[i].size();
    
    _yd[i] = Var3Matrix(size_T_i);
    for (int t = 0; t < size_T_i; ++t)
    {
      _yd[i][t] = VarMatrix(_k);
      for (int j = 0; j < _k; ++j)
      {
        _yd[i][t][j] = VarArray(m);
        for (int p = 0; p < m; ++p)
        {
          snprintf(buf, 1024, "yd;%d;%d;%d;%d", i, t, j, p);
          _yd[i][t][j][p] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
        }
      }
    }
  }
  
  _b = VarArray(_k);
  for (int j = 0; j < _k; ++j)
  {
    snprintf(buf, 1024, "b;%d", j);
    _b[j] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
  }
  
  _model.update();
}

void MinClusterIlpGurobi::initConstraints()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  GRBLinExpr sum, sum2;
  
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
    
    _model.addConstr(sum == 1);
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
          _model.addConstr(_y[i][t][j] * _dLB[i][t][p] <= _yd[i][t][j][p]);
          _model.addConstr(_yd[i][t][j][p] <= _dUB[i][t][p] * _y[i][t][j]);
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
          _model.addConstr(_yd[i][t][j][p] <= _y[i][t][j]);
          _model.addConstr(_yd[i][t][j][p] <= _d[j][p]);
          _model.addConstr(_yd[i][t][j][p] >= _y[i][t][j] + _d[j][p] - 1);
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
    _model.addConstr(sum >= sum2);
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
        _model.addConstr(_b[j] >= _y[i][t][j]);
        sum += _y[i][t][j];
      }
    }
    _model.addConstr(_b[j] <= sum);
    sum.clear();
  }
  
  _model.update();
}

void MinClusterIlpGurobi::initObjective()
{
  GRBLinExpr sum;
  
  for (int j = 0; j < _k; ++j)
  {
    sum += _b[j];
  }
  
  _model.setObjective(sum, GRB_MINIMIZE);
  _model.update();
}

bool MinClusterIlpGurobi::solve(int nrThreads, int timeLimit)
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  if (nrThreads > 0)
  {
    _model.getEnv().set(GRB_IntParam_Threads, nrThreads);
  }
  if (timeLimit > 0)
  {
    _model.getEnv().set(GRB_DoubleParam_TimeLimit, timeLimit);
  }
  _model.getEnv().set(GRB_IntParam_LogToConsole, 0);
  
  _model.optimize();
  int status = _model.get(GRB_IntAttr_Status);
  if (status == GRB_INFEASIBLE)
  {
    return false;
  }
  
  _solD = DoubleMatrix(_k, DoubleVector(m, 0));
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < m; ++p)
    {
      _solD[j][p] = _d[j][p].get(GRB_DoubleAttr_X);
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
        if (_y[i][t][j].get(GRB_DoubleAttr_X) >= 0.4)
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
    if (_b[j].get(GRB_DoubleAttr_X) >= 0.4)
    {
      ++_solMinK;
    }
  }
  
  updateLogLikelihood();
  
  return true;
}
