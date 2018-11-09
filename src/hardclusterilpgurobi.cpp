/*
 * hardclusterilpgurobi.cpp
 *
 *  Created on: 20-oct-2018
 *      Author: M. El-Kebir
 */

#include "hardclusterilpgurobi.h"

HardClusterIlpGurobi::HardClusterIlpGurobi(const ReadMatrix& R,
                                           int k,
                                           int nrSegments)
  : HardClusterIlp(R, k, nrSegments)
  , _env()
  , _model(_env)
  , _y()
  , _d()
  , _yd()
  , _gamma()
  , _lambda()
{
}

void HardClusterIlpGurobi::initHotStart(const BoolTensor& y)
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

void HardClusterIlpGurobi::initVariables()
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
  
  _gamma = Var4Matrix(n);
  for (int i = 0; i < n; ++i)
  {
    _gamma[i] = Var3Matrix(_scriptT[i].size());
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      _gamma[i][t] = VarMatrix(_k);
      for (int j = 0; j < _k; ++j)
      {
        _gamma[i][t][j] = VarArray(_nrSegments);
        for (int l = 0; l < _nrSegments; ++l)
        {
          snprintf(buf, 1024, "gamma;%d;%d;%d;%d", i, t, j, l);
          _gamma[i][t][j][l] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
        }
      }
    }
  }

  _lambda = Var5Matrix(n);
  for (int i = 0; i < n; ++i)
  {
    _lambda[i] = Var4Matrix(_scriptT[i].size());
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      _lambda[i][t] = Var3Matrix(_k);
      for (int j = 0; j < _k; ++j)
      {
        _lambda[i][t][j] = VarMatrix(m);
        for (int p = 0; p < m; ++p)
        {
          _lambda[i][t][j][p] = VarArray(_nrSegments);
          for (int l = 0; l < _nrSegments; ++l)
          {
            snprintf(buf, 1024, "lambda;%d;%d;%d;%d;%d", i, t, j, p, l);
            _lambda[i][t][j][p][l] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
          }
        }
      }
    }
  }
  _model.update();
}

void HardClusterIlpGurobi::initConstraints()
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
          for (int l = 0; l < _nrSegments; ++l)
          {
            sum += _lambda[i][t][j][p][l];
            sum2 += _lambda[i][t][j][p][l] * _x[i][t][j][p][l];
          }
          _model.addConstr(sum == _y[i][t][j]);
          _model.addConstr(sum2 == _yd[i][t][j][p]);
          sum.clear();
          sum2.clear();
        }
      }
    }
  }
  
  GRBLinExpr sum3;
  for (int j = 0; j < _k; ++j)
  {
    for (int ii = 0; ii < n; ++ii)
    {
      for (int tt = 0; tt < _scriptT[ii].size(); ++tt)
      {
        sum3 += _y[ii][tt][j];
      }
    }
    
    for (int i = 0; i < n; ++i)
    {
      for (int t = 0; t < _scriptT[i].size(); ++t)
      {
        for (int l = 0; l < _nrSegments; ++l)
        {
          sum += _gamma[i][t][j][l];
          sum2 += _gamma[i][t][j][l] * _z[l];
        }
        _model.addConstr(sum == _y[i][t][j]);
        _model.addConstr(sum2 <= sum3);
        _model.addConstr(sum2 <= _y[i][t][j]);
        sum.clear();
        sum2.clear();
      }
    }
    sum3.clear();
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
  
  _model.update();
}

void HardClusterIlpGurobi::initObjective()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  GRBLinExpr sum;
  
  for (int i = 0; i < n; ++i)
  {
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        for (int l = 0; l < _nrSegments; ++l)
        {
          sum += _hatN[l] * _gamma[i][t][j][l];
        }
        
        for (int p = 0; p < m; ++p)
        {
          for (int l = 0; l < _nrSegments; ++l)
          {
            sum += _hatG[i][t][j][p][l] * _lambda[i][t][j][p][l];
          }
        }
      }
    }
  }
  
  _model.setObjective(sum, GRB_MAXIMIZE);
  _model.update();
}

bool HardClusterIlpGurobi::solve(int nrThreads,
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
  
  updateLogLikelihood();
  
  return true;
}
