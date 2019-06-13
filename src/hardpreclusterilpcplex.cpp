/*
 * hardpreclusterilpcplex.cpp
 *
 *  Created on: 7-jun-2019
 *      Author: M. El-Kebir
 */

#include "hardpreclusterilpcplex.h"

HardPreClusterIlpCplex::HardPreClusterIlpCplex(const ReadMatrix& R,
                                               int k,
                                               int nrSegments,
                                               ClusterStatisticType statType,
                                               double precisionBetaBin,
                                               bool includePi)
  : HardClusterIlp(R, k, nrSegments, statType, precisionBetaBin, false)
  , _includePi(includePi)
  , _env()
  , _model(_env)
  , _cplex(_model)
  , _y()
  , _d()
  , _gamma()
  , _lambda()
  , _llambda()
  , _w()
  , _z()
{
}

void HardPreClusterIlpCplex::initVariables()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  char buf[1024];
  
  _y = IloNumVar3Matrix(_env, n);
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptT[i].size();
    
    _y[i] = IloNumVarMatrix(_env, size_T_i);
    for (int t = 0; t < size_T_i; ++t)
    {
      _y[i][t] = IloNumVarArray(_env, _k);
      for (int j = 0; j < _k; ++j)
      {
        snprintf(buf, 1024, "y;%d;%d;%d", i, t, j);
        _y[i][t][j] = IloNumVar(_env, 0, 1, buf);
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
  
  if (_includePi)
  {
    _pi = IloNumVarArray(_env, _k);
    for (int j = 0; j < _k; ++j)
    {
      snprintf(buf, 1024, "pi;%d", j);
      _pi[j] = IloNumVar(_env, 0, n, buf);
    }
  
    _gamma = IloNumVarMatrix(_env, _k);
    for (int j = 0; j < _k; ++j)
    {
      _gamma[j] = IloNumVarArray(_env, _nrSegments);
      for (int alpha = 0; alpha < _nrSegments; ++alpha)
      {
        snprintf(buf, 1024, "gamma;%d;%d", j, alpha);
        _gamma[j][alpha] = IloNumVar(_env, 0, n, buf);
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
  
  _llambda = IloNumVar3Matrix(_env, _k);
  for (int j = 0; j < _k; ++j)
  {
    _llambda[j] = IloNumVarMatrix(_env, m);
    for (int p = 0; p < m; ++p)
    {
      _llambda[j][p] = IloNumVarArray(_env, _nrSegments);
      for (int alpha = 0; alpha < _nrSegments; ++alpha)
      {
        snprintf(buf, 1024, "llambda;%d;%d;%d", j, p, alpha);
        _llambda[j][p][alpha] = IloNumVar(_env, 0, 1, buf);
      }
    }
  }
  
  _w = IloBoolVarMatrix(_env, _k);
  for (int j = 0; j < _k; ++j)
  {
    const int scriptT_size = _scriptT[0].size();
    _w[j] = IloBoolVarArray(_env, scriptT_size);
    for (int t = 0; t < scriptT_size; ++t)
    {
      snprintf(buf, 1024, "w;%d;%d", j, t);
      _w[j][t] = IloBoolVar(_env, buf);
    }
  }
  
  _z = IloBoolVarMatrix(_env, n);
  for (int i = 0; i < n; ++i)
  {
    _z[i] = IloBoolVarArray(_env, _k);
    for (int j = 0; j < _k; ++j)
    {
      snprintf(buf, 1024, "z;%d;%d", i, j);
      _z[i][j] = IloBoolVar(_env, buf);
    }
  }
}

void HardPreClusterIlpCplex::initConstraints()
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
  
  if (_includePi)
  {
    for (int j = 0; j < _k; ++j)
    {
      for (int i = 0; i < n; ++i)
      {
        const int size_T_i = _scriptT[i].size();
        for (int t = 0; t < size_T_i; ++t)
        {
          sum += _y[i][t][j];
        }
      }
      _model.add(_pi[j] == sum);
      sum.clear();
    }
    
    for (int j = 0; j < _k; ++j)
    {
      for (int alpha = 0; alpha < _nrSegments; ++alpha)
      {
        sum += _gamma[j][alpha] * _coordPi[alpha];
        sum2 += _gamma[j][alpha];
      }
      _model.add(sum == _pi[j]);
      _model.add(sum2 == 1);
      sum.clear();
      sum2.clear();
    }
  }
  
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < m; ++p)
    {
      for (int alpha = 0; alpha < _nrSegments; ++alpha)
      {
        sum += _llambda[j][p][alpha] * _coord[alpha];
        sum2 += _llambda[j][p][alpha];
      }
      _model.add(sum == _d[j][p]);
      _model.add(sum2 == 1);
      sum.clear();
      sum2.clear();
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
          for (int alpha = 0; alpha < _nrSegments; ++alpha)
          {
            sum += _lambda[i][t][j][p][alpha];
            _model.add(_llambda[j][p][alpha] >= _lambda[i][t][j][p][alpha]);
          }
          _model.add(sum == _y[i][t][j]);
          sum.clear();
          
          for (int alpha = 0; alpha < _nrSegments-1; ++alpha)
          {
            if (_coord[alpha] < _dLB[i][t][p] && _coord[alpha + 1] < _dLB[i][t][p])
            {
              _model.add(_lambda[i][t][j][p][alpha] == 0);
            }
          }
          
          for (int alpha = 1; alpha < _nrSegments; ++alpha)
          {
            if (_coord[alpha-1] > _dUB[i][t][p] && _coord[alpha] > _dUB[i][t][p])
            {
              _model.add(_lambda[i][t][j][p][alpha] == 0);
            }
          }
        }
      }
    }
  }
  
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptT[i].size();
    for (int t = 0; t < size_T_i; ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        _model.add(_y[i][t][j] >= _w[j][t] + _z[i][j] - 1);
      }
    }
  }
  
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < _k; ++j)
    {
      sum += _z[i][j];
    }
    _model.add(sum == 1);
    sum.clear();
  }
  
  const int size_T = _scriptT[0].size();
  for (int j = 0; j < _k; ++j)
  {
    for (int t = 0; t < size_T; ++t)
    {
      sum += _w[j][t];
    }
    _model.add(sum == 1);
    sum.clear();
  }
  
  sum.end();
  sum2.end();
}

void HardPreClusterIlpCplex::initObjective()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  IloExpr sum(_env);
  
  if (_includePi)
  {
    const double log_n = log(n);
    for (int j = 0; j < _k; ++j)
    {
      for (int alpha = 0; alpha < _nrSegments; ++alpha)
      {
        sum += (_hatPi[alpha] - log_n) * _gamma[j][alpha];
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
          for (int alpha = 0; alpha < _nrSegments; ++alpha)
          {
            sum += _hatG[i][t][j][p][alpha] * _lambda[i][t][j][p][alpha];
          }
        }
      }
    }
  }
  
  _model.add(IloMaximize(_env, sum));
  sum.end();
}

bool HardPreClusterIlpCplex::solve(int nrThreads,
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
  else
  {
    _cplex.setOut(std::cerr);
    _cplex.setError(std::cerr);
    _cplex.setWarning(std::cerr);
  }
  
  // TODO: remove me
//  _model.add(_d[0][0] >= .5);
//  _model.add(_w[0][6] == 1);
  
  _cplex.solve();
  if (_cplex.getStatus() == IloAlgorithm::Infeasible)
  {
    return false;
  }
  
  if (_cplex.getSolnPoolNsolns() == 0)
  {
    return false;
  }
  
  if (verbose)
  {
    std::cerr << "Obj value: " << _cplex.getObjValue() << std::endl;
    std::cerr << "Best obj value: " << _cplex.getBestObjValue() << std::endl;
  }
  
  _solD = DoubleMatrix(_k, DoubleVector(m, 0));
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < m; ++p)
    {
      _solD[j][p] = _cplex.getValue(_d[j][p]);
      if (_solD[j][p] == g_tol.epsilon())
      {
        _solD[j][p] = 0;
      }
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
//          double sum = 0;
//          for (int l = 0; l < _nrSegments; ++l)
//          {
//            sum += (_hatPi[l] - log(n)) * _cplex.getValue(_gamma[j][l]);
//          }
//
//          for (int p = 0; p < m; ++p)
//          {
//            for (int l = 0; l < _nrSegments; ++l)
//            {
//              sum += _hatG[i][t][j][p][l] * _cplex.getValue(_lambda[i][t][j][p][l]);
//              if (_cplex.getValue(_lambda[i][t][j][p][l]) > 0.01)
//              {
//                std::cout << _lambda[i][t][j][p][l] << " = " << _cplex.getValue(_lambda[i][t][j][p][l])
//                          << " (" << _hatG[i][t][j][p][l] << ")" << std::endl;
//              }
//            }
//          }
//
//          sum += -log(n);
//          sum += -log(_scriptT[i].size());
//
//          std::cout << i << "," << t << "," << j << ": " << sum << std::endl;
          
          _solT.emplace_back(convertToStateTreeFromDCF(_R,
                                                       _scriptT[i][t],
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
  
  _solZ = IntVector(n, -1);
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < _k; ++j)
    {
      if (_cplex.getValue(_z[i][j]) > .4)
      {
        assert(_solZ[i] == -1);
        _solZ[i] = j;
      }
    }
  }
  
//  const int size_T = _scriptT[0].size();
//  for (int j = 0; j < _k; ++j)
//  {
//    for (int t = 0; t < size_T; ++t)
//    {
//      if (_cplex.getValue(_w[j][t]) > 0.4)
//      {
//        std::cout << _w[j][t] << " = " << _cplex.getValue(_w[j][t]) << std::endl;
//      }
//    }
//  }
  
  for (int i = 0; i < n; ++i)
  {
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        if (g_tol.nonZero(_solY[i][t][j]))
        {
          for (int p = 0; p < m; ++p)
          {
            _solD[j][p] = std::max(_dLB[i][t][p], _solD[j][p]);
            _solD[j][p] = std::min(_dUB[i][t][p], _solD[j][p]);
          }
        }
      }
    }
  }
  
  updateLogLikelihood();
  
  return true;
}
