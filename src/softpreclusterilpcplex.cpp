/*
 * softpreclusterilpcplex.cpp
 *
 *  Created on: 11-jun-2019
 *      Author: M. El-Kebir
 */

#include "softpreclusterilpcplex.h"

SoftPreClusterIlpCplex::SoftPreClusterIlpCplex(const ReadMatrix& R,
                                               int k,
                                               int nrSegmentBits,
                                               ClusterStatisticType statType,
                                               double precisionBetaBin,
                                               bool includePi)
  : SoftClusterIlp(R, k, (1 << nrSegmentBits) + 1, statType, precisionBetaBin, false)
  , _includePi(includePi)
  , _nrSegmentBits(nrSegmentBits)
  , _env()
  , _model(_env)
  , _cplex(_model)
  , _y()
  , _d()
  , _pi()
  , _gamma()
  , _ggamma()
  , _bggamma()
  , _lambda()
  , _llambda()
  , _bllambda()
  , _w()
  , _z()
{
}

void SoftPreClusterIlpCplex::initVariables()
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
  }
  
  if (_includePi)
  {
    _ggamma = IloNumVarMatrix(_env, _k);
    _gamma = IloNumVarMatrix(_env, _k);
    for (int j = 0; j < _k; ++j)
    {
      _ggamma[j] = IloNumVarArray(_env, _nrSegments);
      _gamma[j] = IloNumVarArray(_env, _nrSegments);
      for (int alpha = 0; alpha < _nrSegments; ++alpha)
      {
        snprintf(buf, 1024, "ggamma;%d;%d", j, alpha);
        _ggamma[j][alpha] = IloNumVar(_env, 0, 1, buf);
        
        snprintf(buf, 1024, "gamma;%d;%d", j, alpha);
        _gamma[j][alpha] = IloNumVar(_env, 0, n, buf);
      }
    }
    
    _bggamma = IloBoolVarMatrix(_env, _k);
    for (int j = 0; j < _k; ++j)
    {
      _bggamma[j] = IloBoolVarArray(_env, _nrSegmentBits);
      for (int beta = 0; beta < _nrSegmentBits; ++beta)
      {
        snprintf(buf, 1024, "bggamma;%d;%d", j, beta);
        _bggamma[j][beta] = IloBoolVar(_env, buf);
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
          for (int alpha = 0; alpha < _nrSegments; ++alpha)
          {
            snprintf(buf, 1024, "lambda;%d;%d;%d;%d;%d", i, t, j, p, alpha);
            _lambda[i][t][j][p][alpha] = IloNumVar(_env, 0, 1, buf);
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
  
  _bllambda = IloBoolVar3Matrix(_env, _k);
  for (int j = 0; j < _k; ++j)
  {
    _bllambda[j] = IloBoolVarMatrix(_env, m);
    for (int p = 0; p < m; ++p)
    {
      _bllambda[j][p] = IloBoolVarArray(_env, _nrSegments);
      for (int beta = 0; beta < _nrSegmentBits; ++beta)
      {
        snprintf(buf, 1024, "bllambda;%d;%d;%d", j, p, beta);
        _bllambda[j][p][beta] = IloBoolVar(_env, buf);
      }
    }
  }
  
  _w = IloNumVarMatrix(_env, _k);
  for (int j = 0; j < _k; ++j)
  {
    const int scriptT_size = _scriptT[0].size();
    _w[j] = IloNumVarArray(_env, scriptT_size);
    for (int t = 0; t < scriptT_size; ++t)
    {
      snprintf(buf, 1024, "w;%d;%d", j, t);
      _w[j][t] = IloBoolVar(_env, 0., 1., buf);
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

void SoftPreClusterIlpCplex::initConstraints()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  IloExpr sum(_env), sum2(_env), sum3(_env), sum4(_env), sum5(_env);
  
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
  }
  
  IloNumArray vals(_env, _nrSegments);
  for (int alpha = 0; alpha < _nrSegments; ++alpha)
  {
    vals[alpha] = _coord[alpha];
  }
  
  if (_includePi)
  {
    for (int j = 0; j < _k; ++j)
    {
      //      IloSOS2 sos2(_env, _ggamma[j], vals);
      //      _model.add(sos2);
      
      // adjacency condition
      for (int beta = 0; beta < _nrSegmentBits; ++beta)
      {
        for (int alpha = 0; alpha < _nrSegments; ++alpha)
        {
          if (alpha == 0)
          {
            unsigned int grayCode = binaryToGray(0);
            if (grayCode & (1 << beta))
            {
              sum += _ggamma[j][alpha];
            }
            else
            {
              sum2 += _ggamma[j][alpha];
            }
          }
          else if (alpha == _nrSegments)
          {
            unsigned int grayCode = binaryToGray(_nrSegments - 1);
            if (grayCode & (1 << beta))
            {
              sum += _ggamma[j][alpha];
            }
            else
            {
              sum2 += _ggamma[j][alpha];
            }
          }
          else
          {
            unsigned int grayCode1 = binaryToGray(alpha - 1);
            unsigned int grayCode2 = binaryToGray(alpha);
            if ((grayCode1 & (1 << beta)) && (grayCode2 & (1 << beta)))
            {
              sum += _ggamma[j][alpha];
            }
            else if (!(grayCode1 & (1 << beta)) && !(grayCode2 & (1 << beta)))
            {
              sum2 += _ggamma[j][alpha];
            }
          }
        }
        _model.add(sum <= _bggamma[j][beta]);
        _model.add(sum2 <= 1 - _bggamma[j][beta]);
        sum.clear();
        sum2.clear();
      }
      
      
      for (int alpha = 0; alpha < _nrSegments; ++alpha)
      {
        sum += _ggamma[j][alpha] * _coordPi[alpha];
        sum2 += _ggamma[j][alpha];
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
      //      IloSOS2 sos2(_env, _llambda[j][p], vals);
      //      _model.add(sos2);
      
      // adjacency condition
      for (int beta = 0; beta < _nrSegmentBits; ++beta)
      {
        for (int alpha = 0; alpha < _nrSegments; ++alpha)
        {
          if (alpha == 0)
          {
            unsigned int grayCode = binaryToGray(0);
            if (grayCode & (1 << beta))
            {
              sum += _llambda[j][p][alpha];
            }
            else
            {
              sum2 += _llambda[j][p][alpha];
            }
          }
          else if (alpha == _nrSegments)
          {
            unsigned int grayCode = binaryToGray(_nrSegments - 1);
            if (grayCode & (1 << beta))
            {
              sum += _llambda[j][p][alpha];
            }
            else
            {
              sum2 += _llambda[j][p][alpha];
            }
          }
          else
          {
            unsigned int grayCode1 = binaryToGray(alpha - 1);
            unsigned int grayCode2 = binaryToGray(alpha);
            if ((grayCode1 & (1 << beta)) && (grayCode2 & (1 << beta)))
            {
              sum += _llambda[j][p][alpha];
            }
            else if (!(grayCode1 & (1 << beta)) && !(grayCode2 & (1 << beta)))
            {
              sum2 += _llambda[j][p][alpha];
            }
          }
        }
        _model.add(sum <= _bllambda[j][p][beta]);
        _model.add(sum2 <= 1 - _bllambda[j][p][beta]);
        
        sum.clear();
        sum2.clear();
      }
      
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
        }
      }
    }
  }
  
  if (_includePi)
  {
    for (int j = 0; j < _k; ++j)
    {
      for (int alpha = 0; alpha < _nrSegments; ++alpha)
      {
        sum3 += _gamma[j][alpha];
        _model.add(n * _ggamma[j][alpha] >= _gamma[j][alpha]);
      }
      _model.add(sum3 == _pi[j]);
      sum3.clear();
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
          for (int alpha = 0; alpha < _nrSegments-1; ++alpha)
          {
            if (g_tol.less(_coord[alpha], _dLB[i][t][p]) && g_tol.less(_coord[alpha + 1], _dLB[i][t][p]))
              //            if (_coord[alpha] < _dLB[i][t][p] && _coord[alpha + 1] < _dLB[i][t][p])
            {
              _model.add(_lambda[i][t][j][p][alpha] == 0);
            }
          }
          
          for (int alpha = 1; alpha < _nrSegments; ++alpha)
          {
            if (g_tol.less(_dUB[i][t][p], _coord[alpha-1]) && g_tol.less(_dUB[i][t][p], _coord[alpha]))
              //            if (_coord[alpha-1] > _dUB[i][t][p] && _coord[alpha] > _dUB[i][t][p])
            {
              _model.add(_lambda[i][t][j][p][alpha] == 0);
            }
          }
        }
      }
    }
  }
  
  //  if (_includePi)
  //  {
  //    for (int j = 1; j < _k; ++j)
  //    {
  //      _model.add(_pi[j-1] >= _pi[j]);
  //    }
  //  }
  
  
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
  //  else
  //  {
  //    for (int j = 1; j < _k; ++j)
  //    {
  //      _model.add(_pi[j] >= _pi[j-1]);
  //    }
  //  }
  
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
  sum3.end();
  sum4.end();
}

void SoftPreClusterIlpCplex::initObjective()
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

bool SoftPreClusterIlpCplex::solve(int nrThreads,
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
  
  std::cout << "Obj value: " << _cplex.getObjValue() << std::endl;
  std::cout << "Best obj value: " << _cplex.getBestObjValue() << std::endl;
  
  _solD = DoubleMatrix(_k, DoubleVector(m, 0));
  _solPi = DoubleVector(_k);
  
  // get d
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

  // get y
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
        if (!g_tol.nonZero(_solY[i][t][j]))
        {
          _solY[i][t][j] = 0;
        }
      }
    }
  }
  
  // get pi
  if (_includePi)
  {
    for (int j = 0; j < _k; ++j)
    {
      _solPi[j] = _cplex.getValue(_pi[j]) / n;
    }
  }
  else
  {
    // compute maximum likelihood pi
    for (int j = 0; j < _k; ++j)
    {
      for (int i = 0; i < n; ++i)
      {
        for (int t = 0; t < _scriptT[i].size(); ++t)
        {
          _solPi[j] += _solY[i][t][j];
        }
      }
      
      _solPi[j] /= n;
    }
  }
  
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

  _solT = PosteriorStateTreeMatrix(n);
  _solZ.clear();
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
        if (g_tol.nonZero(_solY[i][t][j]))
        {
          DoubleVector h_i;
          for (int p = 0; p < m; ++p)
          {
            h_i.push_back((_solD[j][p] - _numerator[i][t][p]) / _denominator[i][p]);
          }
          
          _solT[i].emplace_back(convertToStateTreeFromSNVF(_scriptT[i][t],
                                                           h_i, i),
                                _solY[i][t][j], t, j);
        }
      }
    }

    // TODO: shouldn't we sort?
    assert(!_solT[i].empty());
    _solZ.push_back(_solT[i].back()._j);
  }
  
  updateLogLikelihood();
  return true;
}
