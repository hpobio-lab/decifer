/*
 * clusterilpcplex.cpp
 *
 *  Created on: 16-oct-2018
 *      Author: M. El-Kebir
 */

#include "clusterilpcplex.h"
//#include "clusterilpcplexcallback.h"

ClusterIlpCplex::ClusterIlpCplex(const ReadMatrix& R,
                                 int k,
                                 double alpha,
                                 ClusterStatisticType statType)
  : ClusterIlp(R, k, alpha, statType)
  , _env()
  , _model(_env)
  , _cplex(_model)
  , _y()
  , _d()
  , _yd()
  , _w()
{
}

void ClusterIlpCplex::initHotStart(const BoolTensor& y)
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

void ClusterIlpCplex::initVariables()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  char buf[1024];
  
  _y = IloBoolVar3Matrix(_env, n);
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptTlb[i].size();
    assert(size_T_i == _scriptTub[i].size());
    
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
    const int size_T_i = _scriptTlb[i].size();
    assert(size_T_i == _scriptTub[i].size());
    
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
  
  _w = IloNumVar4Matrix(_env, n);
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptTlb[i].size();
    assert(size_T_i == _scriptTub[i].size());
    assert(size_T_i == _scriptTexp[i].size());
    
    _w[i] = IloNumVar3Matrix(_env, size_T_i);
    for (int t = 0; t < size_T_i; ++t)
    {
      _w[i][t] = IloNumVarMatrix(_env, _k);
      for (int j = 0; j < _k; ++j)
      {
        _w[i][t][j] = IloNumVarArray(_env, m);
        for (int p = 0; p < m; ++p)
        {
          snprintf(buf, 1024, "w;%d;%d;%d;%d", i, t, j, p);
          _w[i][t][j][p] = IloNumVar(_env, 0, 1, buf);
        }
      }
    }
  }
}

void ClusterIlpCplex::initConstraints()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  IloExpr sum(_env);
  
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
    
    _model.add(sum == 1);
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
          _model.add(_yd[i][t][j][p] <= _y[i][t][j]);
//          _model.add(_yd[i][t][j][p] <= _d[j][p]);
          _model.add(_yd[i][t][j][p] >= _y[i][t][j] + _d[j][p] - 1);
          _model.add(_yd[i][t][j][p] >= _dcfLB[i][t][p] * _y[i][t][j]);
          _model.add(_yd[i][t][j][p] <= _dcfUB[i][t][p] * _y[i][t][j]);
        }
      }
    }
  }
  
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptTlb[i].size();
    for (int j = 0; j < _k; ++j)
    {
      for (int p = 0; p < m; ++p)
      {
        for (int t = 0; t < size_T_i; ++t)
        {
          sum += _yd[i][t][j][p];
        }
        _model.add(sum <= _d[j][p]);
        sum.clear();
      }
    }
  }

  for (int i = 0; i < n; ++i)
  {
    for (int p = 0; p < m; ++p)
    {
      const int size_T_i = _scriptTlb[i].size();
      for (int t = 0; t < size_T_i; ++t)
      {
        for (int j = 0; j < _k; ++j)
        {
          sum += _yd[i][t][j][p];
        }
      }
      _model.add(sum <= 1);
      sum.clear();
    }
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
          _model.add(_w[i][t][j][p] >= _yd[i][t][j][p] - _dcfExp[i][t][p] * _y[i][t][j]);
          _model.add(_w[i][t][j][p] >= _dcfExp[i][t][p] * _y[i][t][j] - _yd[i][t][j][p]);
          sum += _w[i][t][j][p];
        }
      }
    }
  }
  
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < m; ++p)
    {
      _model.add(0 <= _d[j][p]);
      _model.add(_d[j][p] <= 1);
    }
  }
  
  for (int j = 1; j < _k; ++j)
  {
    _model.add(_d[j-1][0] >= _d[j][0]);
  }
  
  _model.add(IloMinimize(_env, sum));
  sum.end();
}

bool ClusterIlpCplex::solve(int nrThreads,
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
  
//  IloFastMutex mutex;
//  _cplex.use(IloCplex::Callback(new (_env) ClusterIlpCplexCallback<IloCplex::UserCutCallbackI>(_env, m, n, _k, _y, _yd, _w, _d, _dcfLB, _dcfUB, _dcfExp, &mutex)));
//  _cplex.use(IloCplex::Callback(new (_env) ClusterIlpCplexCallback<IloCplex::LazyConstraintCallbackI>(_env, m, n, _k, _y, _yd, _w, _d, _dcfLB, _dcfUB, _dcfExp, &mutex)));

  try
  {
    _cplex.solve();
  }
  catch (IloException& e)
  {
    std::cerr << e.getMessage() << std::endl;
    e.print(std::cerr);
    return false;
  }
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
  _solY.clear();
  _solPi = DoubleVector(_k);
  _solY = BoolTensor(n);
  for (int i = 0; i < n; ++i)
  {
    _solY[i] = BoolMatrix(_scriptTlb[i].size());
    for (int t = 0; t < _scriptTlb[i].size(); ++t)
    {
      _solY[i][t] = BoolVector(_k);
      for (int j = 0; j < _k; ++j)
      {
        _solY[i][t][j] = (_cplex.getValue(_y[i][t][j]) >= 0.4);
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
