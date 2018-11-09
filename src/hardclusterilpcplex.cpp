/*
 * hardclusterilpcplex.cpp
 *
 *  Created on: 19-oct-2018
 *      Author: M. El-Kebir
 */

#include "hardclusterilpcplex.h"

HardClusterIlpCplex::HardClusterIlpCplex(const ReadMatrix& R,
                                         int k,
                                         int nrSegments)
  : HardClusterIlp(R, k, nrSegments)
  , _env()
  , _model(_env)
  , _cplex(_model)
  , _y()
  , _d()
  , _yd()
  , _gamma()
  , _lambda()
{
}

void HardClusterIlpCplex::initHotStart(const BoolTensor& y)
{
  const int n = _R.getNrCharacters();
  if (y.size() != n)
  {
    return;
  }
  
  IloNumVarArray startVar(_env);
  IloNumArray startVal(_env);

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
          startVar.add(_y[i][t][j]);
          startVal.add(y[i][t][j] ? 1 : 0);
        }
        else
        {
          startVar.add(_y[i][t][j]);
          startVal.add(0);
        }
      }
    }
  }
  
  _cplex.addMIPStart(startVar, startVal);
  startVal.end();
  startVar.end();
}

void HardClusterIlpCplex::initVariables()
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
  
  _gamma = IloNumVar4Matrix(_env, n);
  for (int i = 0; i < n; ++i)
  {
    _gamma[i] = IloNumVar3Matrix(_env, _scriptT[i].size());
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      _gamma[i][t] = IloNumVarMatrix(_env, _k);
      for (int j = 0; j < _k; ++j)
      {
        _gamma[i][t][j] = IloNumVarArray(_env, _nrSegments);
        for (int l = 0; l < _nrSegments; ++l)
        {
          snprintf(buf, 1024, "gamma;%d;%d;%d;%d", i, t, j, l);
          _gamma[i][t][j][l] = IloNumVar(_env, 0, 1, buf);
        }
      }
    }
  }

  _lambda = IloNumVar5Matrix(_env, n);
  for (int i = 0; i < n; ++i)
  {
    _lambda[i] = IloNumVar4Matrix(_env, _scriptT[i].size());
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      _lambda[i][t] = IloNumVar3Matrix(_env, _k);
      for (int j = 0; j < _k; ++j)
      {
        _lambda[i][t][j] = IloNumVarMatrix(_env, m);
        for (int p = 0; p < m; ++p)
        {
          _lambda[i][t][j][p] = IloNumVarArray(_env, _nrSegments);
          for (int l = 0; l < _nrSegments; ++l)
          {
            snprintf(buf, 1024, "lambda;%d;%d;%d;%d;%d", i, t, j, p, l);
            _lambda[i][t][j][p][l] = IloNumVar(_env, 0, 1, buf);
          }
        }
      }
    }
  }
}

void HardClusterIlpCplex::initConstraints()
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
          for (int l = 0; l < _nrSegments; ++l)
          {
            sum += _lambda[i][t][j][p][l];
            sum2 += _lambda[i][t][j][p][l] * _x[i][t][j][p][l];
          }
          _model.add(sum == _y[i][t][j]);
          _model.add(sum2 == _yd[i][t][j][p]);
          sum.clear();
          sum2.clear();
        }
      }
    }
  }
  
  IloExpr sum3(_env);
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
        _model.add(sum == _y[i][t][j]);
        _model.add(sum2 <= sum3);
        _model.add(sum2 <= _y[i][t][j]);
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
    _model.add(sum >= sum2);
    sum = sum2;
    sum2.clear();
  }
  sum.clear();
  
  sum.end();
  sum2.end();
  sum3.end();
}

void HardClusterIlpCplex::initObjective()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  IloExpr sum(_env);
  
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
    
  _model.add(IloMaximize(_env, sum));
  sum.end();
}

bool HardClusterIlpCplex::solve(int nrThreads,
                                int timeLimit,
                                bool verbose,
                                int memoryLimit)
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
  if (memoryLimit > 0)
  {
    _cplex.setParam(IloCplex::WorkMem, memoryLimit);
  }
  
  if (!verbose)
  {
    _cplex.setOut(_env.getNullStream());
    _cplex.setError(_env.getNullStream());
    _cplex.setWarning(_env.getNullStream());
  }
  
  _cplex.solve();
  if (_cplex.getStatus() == IloAlgorithm::Infeasible)
  {
    return false;
  }
  
  if (_cplex.getSolnPoolNsolns() == 0)
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
  _solY = BoolTensor(n);
  for (int i = 0; i < n; ++i)
  {
    _solY[i] = BoolMatrix(_scriptT[i].size());
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      _solY[i][t] = BoolVector(_k);
      for (int j = 0; j < _k; ++j)
      {
        _solY[i][t][j] = (_cplex.getValue(_y[i][t][j]) >= 0.4);
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
