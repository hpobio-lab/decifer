/*
 * clusterilpgurobi.cpp
 *
 *  Created on: 16-oct-2018
 *      Author: M. El-Kebir
 */

#include "clusterilpgurobi.h"

ClusterIlpGurobi::ClusterIlpGurobi(const ReadMatrix& R,
                                   int k,
                                   double alpha)
  : ClusterIlp(R, k, alpha)
  , _env()
  , _model(_env)
  , _y()
  , _d()
  , _yd()
  , _w()
{
}

void ClusterIlpGurobi::initHotStart(const BoolTensor& y)
{
  const int n = _R.getNrCharacters();
  if (y.size() != n)
  {
    return;
  }
  
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptT[i].size();
    if (y[i].size() != size_T_i)
    {
      return;
    }
    
    for (int t = 0; t < size_T_i; ++t)
    {
      if (y[i][t].size() > _k)
      {
        return;
      }
      
      for (int j = 0; j < _k; ++j)
      {
        if (j < y[i][t].size())
        {
          _y[i][t][j].set(GRB_DoubleAttr_Start, y[i][t][j] ? 1 : 0);
        }
        else
        {
          _y[i][t][j].set(GRB_DoubleAttr_Start, 0);
        }
      }
    }
  }
}

void ClusterIlpGurobi::initVariables()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  char buf[1024];
  
  _y = Var3Matrix(n);
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptTlb[i].size();
    assert(size_T_i == _scriptTub[i].size());
    
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
      snprintf(buf, 1024, "f;%d;%d", j, p);
      _d[j][p] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
    }
  }
  
  _yd = Var4Matrix(n);
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptTlb[i].size();
    assert(size_T_i == _scriptTub[i].size());
    
    _yd[i] = Var3Matrix(size_T_i);
    for (int t = 0; t < size_T_i; ++t)
    {
      _yd[i][t] = VarMatrix(_k);
      for (int j = 0; j < _k; ++j)
      {
        _yd[i][t][j] = VarArray(m);
        for (int p = 0; p < m; ++p)
        {
          snprintf(buf, 1024, "ff;%d;%d;%d;%d", i, t, j, p);
          _yd[i][t][j][p] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
        }
      }
    }
  }
  
  _w = Var4Matrix(n);
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptTlb[i].size();
    assert(size_T_i == _scriptTub[i].size());
    assert(size_T_i == _scriptTexp[i].size());
    
    _w[i] = Var3Matrix(size_T_i);
    for (int t = 0; t < size_T_i; ++t)
    {
      _w[i][t] = VarMatrix(_k);
      for (int j = 0; j < _k; ++j)
      {
        _w[i][t][j] = VarArray(m);
        for (int p = 0; p < m; ++p)
        {
          snprintf(buf, 1024, "w;%d;%d;%d;%d", i, t, j, p);
          _w[i][t][j][p] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
        }
      }
    }
  }
  
  _model.update();
}

void ClusterIlpGurobi::initConstraints()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  GRBLinExpr sum;
  
  // one state tree and cluster per SNV
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptTlb[i].size();
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
    const int size_T_i = _scriptTlb[i].size();
    for (int t = 0; t < size_T_i; ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        for (int p = 0; p < m; ++p)
        {
          _model.addConstr(_yd[i][t][j][p] <= _y[i][t][j]);
          _model.addConstr(_yd[i][t][j][p] <= _d[j][p]);
          _model.addConstr(_yd[i][t][j][p] >= _y[i][t][j] + _d[j][p] - 1);
          _model.addConstr(_yd[i][t][j][p] >= _dcfLB[i][t][p] * _y[i][t][j]);
          _model.addConstr(_yd[i][t][j][p] <= _dcfUB[i][t][p] * _y[i][t][j]);
        }
      }
    }
  }
  
//  for (int i = 0; i < n; ++i)
//  {
//    for (int t = 0; t < _scriptT[i].size(); ++t)
//    {
//      sum += _y[i][t][0];
//    }
//  }
//
  GRBLinExpr sum2;
//  for (int j = 1; j < _k; ++j)
//  {
//    for (int i = 0; i < n; ++i)
//    {
//      for (int t = 0; t < _scriptT[i].size(); ++t)
//      {
//        sum2 += _y[i][t][j];
//      }
//    }
//    _model.addConstr(sum >= sum2);
//    sum = sum2;
//    sum2.clear();
//  }
//  sum.clear();
  
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptTlb[i].size();
    for (int t = 0; t < size_T_i; ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        for (int p = 0; p < m; ++p)
        {
          _model.addConstr(_w[i][t][j][p] >= _yd[i][t][j][p] - _dcfExp[i][t][p] * _y[i][t][j]);
          _model.addConstr(_w[i][t][j][p] >= _dcfExp[i][t][p] * _y[i][t][j] - _yd[i][t][j][p]);
          sum += _w[i][t][j][p];
//          sum2 += _w[i][t][j][p];
        }
//        _model.addConstr(sum2 <= m * _y[i][t][j]);
//        sum2.clear();
      }
    }
  }
  
  for (int j = 1; j < _k; ++j)
  {
    _model.addConstr(_d[j-1][0] >= _d[j][0]);
  }
  
  _model.setObjective(sum, GRB_MINIMIZE);
  
  _model.update();
}

bool ClusterIlpGurobi::solve(int nrThreads,
                             int timeLimit,
                             bool verbose,
                             int memoryLimit)
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
  if (memoryLimit > 0)
  {
    _model.getEnv().set(GRB_DoubleParam_NodeLimit, memoryLimit);
  }
  if (!verbose)
  {
    _model.getEnv().set(GRB_IntParam_LogToConsole, 0);
  }

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
  _solY = BoolTensor(n);
  for (int i = 0; i < n; ++i)
  {
    _solY[i] = BoolMatrix(_scriptT[i].size());
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      _solY[i][t] = BoolVector(_k);
      for (int j = 0; j < _k; ++j)
      {
        _solY[i][t][j] = (_y[i][t][j].get(GRB_DoubleAttr_X) >= 0.4);
        if (_solY[i][t][j])
        {
          _solT.push_back(_scriptTexp[i][t]);
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
  
  updateLogLikelihood();
  
  return true;
}
