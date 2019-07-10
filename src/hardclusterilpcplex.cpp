/*
 * hardclusterilpcplex.cpp
 *
 *  Created on: 19-oct-2018
 *      Author: M. El-Kebir
 */

#include "hardclusterilpcplex.h"

HardClusterIlpCplex::HardClusterIlpCplex(const ReadMatrix& R,
                                         int k,
                                         int nrSegments,
                                         ClusterStatisticType statType,
                                         double precisionBetaBin,
                                         bool forceTruncal,
                                         bool includePi)
  : HardClusterIlp(R, k, nrSegments, statType, precisionBetaBin, forceTruncal)
  , _includePi(includePi)
  , _env()
  , _model(_env)
  , _cplex(_model)
  , _y()
  , _d()
  , _gamma()
  , _lambda()
{
}

void HardClusterIlpCplex::initPreClusteringConstraints(const IntMatrix& preClustering)
{
  for (const IntVector& preCluster : preClustering)
  {
    const int size = preCluster.size();
    for (int i = 1; i < size; ++i)
    {
      initPreClusteringConstraint(preCluster[i-1], preCluster[i]);
    }
  }
}

void HardClusterIlpCplex::initPreClusteringConstraint(int i1, int i2)
{
  assert(_scriptT[i1].size() == _scriptT[i2].size());
  
  const int m = _R.getNrSamples();
  const int scriptT_size = _scriptT[i1].size();
  for (int t = 0; t < scriptT_size; ++t)
  {
    for (int j = 0; j < _k; ++j)
    {
      _model.add(_y[i1][t][j] == _y[i2][t][j]);
      
      for (int p = 0; p < m; ++p)
      {
        for (int alpha = 0; alpha < _nrSegments; ++alpha)
        {
          _model.add(_lambda[i1][t][j][p][alpha] == _lambda[i2][t][j][p][alpha]);
        }
      }
    }
  }
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
      bool ok = true;
      for (int p = 0; p < m; ++p)
      {
        if (_dUB[i][t][p] == 0 && _R.getVar(p, i) > 0)
        {
          ok = false;
        }
      }
      
      if (!ok)
      {
        for (int j = 0; j < _k; ++j)
        {
          _model.add(_y[i][t][j] == 0);
        }
      }
    }
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
//  if (_includePi)
//  {
//    for (int j = 1; j < _k; ++j)
//    {
//      _model.add(_pi[j] <= _pi[j-1]);
//    }
//  }
   
  sum.end();
  sum2.end();
}

void HardClusterIlpCplex::initObjective()
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
//            sum += _hatN[l] * _cplex.getValue(_gamma[i][t][j][l]);
//          }
          
//          for (int p = 0; p < m; ++p)
//          {
//            for (int l = 0; l < _nrSegments; ++l)
//            {
//              sum += _hatG[i][t][j][p][l] * _cplex.getValue(_lambda[i][t][j][p][l]);
//            }
//          }
//
//          sum += -log(n);
//          sum += -log(_scriptT[i].size());
          
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
  
  updateLogLikelihood();
  
  return true;
}
