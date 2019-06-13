/*
 * softclusterilpcplex.cpp
 *
 *  Created on: 15-apr-2019
 *      Author: M. El-Kebir
 */

#include "softclusterilpcplex.h"
#include "softclusterilpcplexcallback.h"

SoftClusterIlpCplex::SoftClusterIlpCplex(const ReadMatrix& R,
                                         int k,
                                         int nrSegments,
                                         ClusterStatisticType statType,
                                         double precisionBetaBin,
                                         bool forceTruncal)
  : SoftClusterIlp(R, k, nrSegments, statType, precisionBetaBin, forceTruncal)
  , _env()
  , _model(_env)
  , _cplex(_model)
  , _y()
  , _yy()
  , _d()
  , _yyd()
  , _pi()
  , _gamma()
  , _lambda()
{
}

void SoftClusterIlpCplex::initPreClusteringConstraint(int i1, int i2)
{
  assert(_scriptT[i1].size() == _scriptT[i2].size());
  
  const int scriptT_size = _scriptT[i1].size();
  for (int t = 0; t < scriptT_size; ++t)
  {
    for (int j = 0; j < _k; ++j)
    {
      _model.add(_y[i1][t][j] == _y[i2][t][j]);
    }
  }
}

void SoftClusterIlpCplex::initHotStart(const BoolTensor& y)
{
}

void SoftClusterIlpCplex::initVariables()
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
  
  _yy = IloBoolVar3Matrix(_env, n);
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptT[i].size();
    
    _yy[i] = IloBoolVarMatrix(_env, size_T_i);
    for (int t = 0; t < size_T_i; ++t)
    {
      _yy[i][t] = IloBoolVarArray(_env, _k);
      for (int j = 0; j < _k; ++j)
      {
        snprintf(buf, 1024, "yy;%d;%d;%d", i, t, j);
        _yy[i][t][j] = IloBoolVar(_env, buf);
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
  
  _yyd = IloNumVar4Matrix(_env, n);
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptT[i].size();
    
    _yyd[i] = IloNumVar3Matrix(_env, size_T_i);
    for (int t = 0; t < size_T_i; ++t)
    {
      _yyd[i][t] = IloNumVarMatrix(_env, _k);
      for (int j = 0; j < _k; ++j)
      {
        _yyd[i][t][j] = IloNumVarArray(_env, m);
        for (int p = 0; p < m; ++p)
        {
          snprintf(buf, 1024, "yyd;%d;%d;%d;%d", i, t, j, p);
          _yyd[i][t][j][p] = IloNumVar(_env, 0, 1, buf);
        }
      }
    }
  }
  
  _pi = IloNumVarArray(_env, _k);
  for (int j = 0; j < _k; ++j)
  {
    snprintf(buf, 1024, "pi;%d", j);
    _pi[j] = IloNumVar(_env, 0, 1, buf);
  }
  
  _gamma = IloNumVar5Matrix(_env, n);
  for (int i = 0; i < n; ++i)
  {
    _gamma[i] = IloNumVar4Matrix(_env, _scriptT[i].size());
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      _gamma[i][t] = IloNumVar3Matrix(_env, _k);
      for (int j = 0; j < _k; ++j)
      {
        _gamma[i][t][j] = IloNumVarMatrix(_env, _nrSegments);
        for (int alpha = 0; alpha < _nrSegments; ++alpha)
        {
          _gamma[i][t][j][alpha] = IloNumVarArray(_env, _nrSegments);
          for (int beta = 0; beta < 2; ++beta)
          {
            snprintf(buf, 1024, "gamma;%d;%d;%d;%d;%d", i, t, j, alpha, beta);
            _gamma[i][t][j][alpha][beta] = IloNumVar(_env, 0, 1, buf);
          }
        }
      }
    }
  }

  _lambda = IloNumVar6Matrix(_env, n);
  for (int i = 0; i < n; ++i)
  {
    _lambda[i] = IloNumVar5Matrix(_env, _scriptT[i].size());
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      _lambda[i][t] = IloNumVar4Matrix(_env, _k);
      for (int j = 0; j < _k; ++j)
      {
        _lambda[i][t][j] = IloNumVar3Matrix(_env, m);
        for (int p = 0; p < m; ++p)
        {
          _lambda[i][t][j][p] = IloNumVarMatrix(_env, _nrSegments);
          for (int alpha = 0; alpha < _nrSegments; ++alpha)
          {
            _lambda[i][t][j][p][alpha] = IloNumVarArray(_env, _nrSegments);
            for (int beta = 0; beta < 2; ++beta)
            {
              snprintf(buf, 1024, "lambda;%d;%d;%d;%d;%d;%d", i, t, j, p, alpha, beta);
              _lambda[i][t][j][p][alpha][beta] = IloNumVar(_env, 0, 1, buf);
            }
          }
        }
      }
    }
  }
}

void SoftClusterIlpCplex::initConstraints()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  IloExpr sum(_env), sum2(_env), sum3(_env);
  
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
  
  for (int j = 0; j < _k; ++j)
  {
    sum += _pi[j];
  }
  _model.add(sum == 1);
  sum.clear();
  
  for (int i = 0; i < n; ++i)
  {
    const int size_T_i = _scriptT[i].size();
    for (int t = 0; t < size_T_i; ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        _model.add(_yy[i][t][j] >= 0);
//        _model.add(_yy[i][t][j] >= _y[i][t][j]);
//        for (int p = 0; p < m; ++p)
//        {
//          _model.add(_yy[i][t][j] * _dLB[i][t][p] <= _yyd[i][t][j][p]);
//          _model.add(_yyd[i][t][j][p] <= _dUB[i][t][p] * _yy[i][t][j]);
//        }
      }
    }
  }
//
  for (int i = 0; i < n; ++i)
  {
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        for (int p = 0; p < m; ++p)
        {
          _model.add(_yyd[i][t][j][p] >= 0);
//          _model.add(_yyd[i][t][j][p] <= _yy[i][t][j]);
//          _model.add(_yyd[i][t][j][p] <= _d[j][p]);
//          _model.add(_yyd[i][t][j][p] >= _yy[i][t][j] + _d[j][p] - 1);
          for (int alpha = 0; alpha < _nrSegments; ++alpha)
          {
            for (int beta = 0; beta < 2; ++beta)
            {
              sum += _lambda[i][t][j][p][alpha][beta];
              sum2 += _lambda[i][t][j][p][alpha][beta] * _coord[alpha];
              if (beta == 1)
              {
                sum3 += _lambda[i][t][j][p][alpha][beta];
              }
            }
          }
          _model.add(sum == 1);
          _model.add(sum2 == _d[j][p]);
          _model.add(sum3 == _y[i][t][j]);
          sum.clear();
          sum2.clear();
          sum3.clear();
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
        for (int alpha = 0; alpha < _nrSegments; ++alpha)
        {
          for (int beta = 0; beta < _nrSegments; ++beta)
          {
            sum += _gamma[i][t][j][alpha][beta];
            sum2 += _gamma[i][t][j][alpha][beta] * _coord[alpha];
            if (beta == 1)
            {
              sum3 += _gamma[i][t][j][alpha][beta];
            }
          }
        }
        _model.add(sum == 1);
        _model.add(sum2 == _pi[j]);
        _model.add(sum3 == _y[i][t][j]);
        sum.clear();
        sum2.clear();
        sum3.clear();
      }
    }
  }
  
  // Cluster sizes are nonincreasing
//  for (int i = 0; i < n; ++i)
//  {
//    for (int t = 0; t < _scriptT[i].size(); ++t)
//    {
//      sum += _y[i][t][0];
//    }
//  }
//  for (int j = 1; j < _k; ++j)
//  {
//    for (int i = 0; i < n; ++i)
//    {
//      for (int t = 0; t < _scriptT[i].size(); ++t)
//      {
//        sum2 += _y[i][t][j];
//      }
//    }
//    _model.add(sum >= sum2);
//    sum = sum2;
//    sum2.clear();
//  }
//  sum.clear();
  for (int j = 1; j < _k; ++j)
  {
    _model.add(_d[j][0] <= _d[j-1][0]);
  }
  
  
  if (_forceTruncal)
  {
    // Cluster 0 has largest DCF in all samples
    for (int j = 1; j < _k; ++j)
    {
      for (int p = 0; p < m; ++p)
      {
        _model.add(_d[j][p] <= _d[0][p]);
      }
    }
  }
  
  sum.end();
  sum2.end();
  sum3.end();
}

void SoftClusterIlpCplex::initObjective()
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
        for (int alpha = 0; alpha < _nrSegments; ++alpha)
        {
            sum += _hatPi[alpha] * _gamma[i][t][j][alpha][1];
        }
        
        for (int p = 0; p < m; ++p)
        {
          for (int alpha = 0; alpha < _nrSegments; ++alpha)
          {
            sum += _hatG[i][t][j][p][alpha] * _lambda[i][t][j][p][alpha][1];
          }
        }
      }
    }
  }
  
  _model.add(IloMaximize(_env, sum));
  sum.end();
}

bool SoftClusterIlpCplex::solve(int nrThreads,
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
  
  IloFastMutex mutex;
//  _cplex.use(IloCplex::Callback(new (_env) SoftClusterIlpCplexCallback<IloCplex::UserCutCallbackI>(_env, m, n, _k, _y, _yy, _yyd, _d, _dLB, _dUB, &mutex)));
  _cplex.use(IloCplex::Callback(new (_env) SoftClusterIlpCplexCallback<IloCplex::LazyConstraintCallbackI>(_env, m, n, _k, _y, _yy, _yyd, _d, _dLB, _dUB, &mutex)));
  
  try {
    _cplex.solve();
  } catch (IloException& e) {
    std::cerr << e.getMessage() << std::endl;
    return false;
  }
//  _cplex.solve();
  if (_cplex.getStatus() == IloAlgorithm::Infeasible)
  {
    return false;
  }
  
  if (_cplex.getSolnPoolNsolns() == 0)
  {
    return false;
  }
  
//  std::cout << "Obj value: " << _cplex.getObjValue() << std::endl;
//  std::cout << "Best obj value: " << _cplex.getBestObjValue() << std::endl;
  
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
  
  _solT = PosteriorStateTreeMatrix(n);
  _solZ.clear();
  _solPi = DoubleVector(_k);
  _solY = DoubleTensor(n);
  for (int i = 0; i < n; ++i)
  {
    _solY[i] = DoubleMatrix(_scriptT[i].size());
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      _solY[i][t] = DoubleVector(_k);
      for (int j = 0; j < _k; ++j)
      {
        _solY[i][t][j] = _cplex.getValue(_y[i][t][j]);
//        if (_solY[i][t][j])
        {
          double sum = 0;
          for (int alpha = 0; alpha < _nrSegments; ++alpha)
          {
            sum += _hatPi[alpha] * _cplex.getValue(_gamma[i][t][j][alpha][1]);
          }
          
          for (int p = 0; p < m; ++p)
          {
            for (int alpha = 0; alpha < _nrSegments; ++alpha)
            {
              sum += _hatG[i][t][j][p][alpha] * _cplex.getValue(_lambda[i][t][j][p][alpha][1]);
            }
          }
          
          DoubleVector h_i;
          for (int p = 0; p < m; ++p)
          {
            h_i.push_back((_solD[j][p] - _numerator[i][t][p]) / _denominator[i][p]);
          }
          
          _solT[i].emplace_back(convertToStateTreeFromSNVF(_R, _scriptT[i][t],
                                                           h_i, i),
                                _solY[i][t][j], t, j);
          
          std::cout << i << "," << t << "," << j << ": " << _solY[i][t][j] << " -- " << sum << std::endl;
        }
      }
    }
    
//    std::sort(_solT[i].begin(),
//              _solT[i].end(),
//              [](const PosteriorStateTree& a,
//                 const PosteriorStateTree& b)
//                {
//                  return a._gamma > b._gamma;
//                });

    // TODO: shouldn't we sort?
    assert(!_solT[i].empty());
    _solZ.push_back(_solT[i].back()._j);
  }
  
  for (int j = 0; j < _k; ++j)
  {
    _solPi[j] = _cplex.getValue(_pi[j]);
  }
  
  updateLogLikelihood();
  
  return true;
}
