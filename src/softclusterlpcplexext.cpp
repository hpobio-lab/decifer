/*
 * softclusterlpcplexext.cpp
 *
 *  Created on: 17-apr-2019
 *      Author: M. El-Kebir
 */

#include "softclusterlpcplexext.h"
#include "ilconcert/ilolinear.h"
#include "softclusterlpcplexcallback.h"

SoftClusterLpCplexExt::SoftClusterLpCplexExt(const ReadMatrix& R,
                                             int k,
                                             int nrSegmentBits,
                                             ClusterStatisticType statType,
                                             double precisionBetaBin,
                                             bool forceTruncal,
                                             bool pi)
  : SoftClusterIlp(R, k, (1 << nrSegmentBits) + 1, statType, precisionBetaBin, forceTruncal)
  , _nrSegmentBits(nrSegmentBits)
  , _env()
  , _model(_env)
  , _cplex(_model)
  , _y()
  , _d()
  , _pi()
  , _gamma()
  , _bgamma()
  , _lambda()
  , _blambda()
  , _sigma()
  , _bsigma()
  , _ssigma()
  , _l()
  , _includePi(pi)
{
}

void SoftClusterLpCplexExt::initPreClusteringConstraint(int i1, int i2)
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

void SoftClusterLpCplexExt::initHotStart(const BoolTensor& y)
{
}

void SoftClusterLpCplexExt::initVariables()
{
  for (unsigned int i = 0; i < _nrSegments - 1; ++i)
  {
    std::cout << i << ":";
    unsigned int gray = binaryToGray(i);

    for (int bit = 0; bit < 32; ++bit)
    {
      std::cout << " ";
      if (gray & (1 << bit))
        std::cout << 1;
      else
        std::cout << 0;
    }
    std::cout << std::endl;
  }
  
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
    _gamma = IloNumVarMatrix(_env, _k);
    for (int j = 0; j < _k; ++j)
    {
      _gamma[j] = IloNumVarArray(_env, _nrSegments);
      for (int alpha = 0; alpha < _nrSegments; ++alpha)
      {
        snprintf(buf, 1024, "gamma;%d;%d", j, alpha);
        _gamma[j][alpha] = IloNumVar(_env, 0, 1, buf);
      }
    }
    
    _bgamma = IloBoolVarMatrix(_env, _k);
    for (int j = 0; j < _k; ++j)
    {
      _bgamma[j] = IloBoolVarArray(_env, _nrSegmentBits);
      for (int beta = 0; beta < _nrSegmentBits; ++beta)
      {
        snprintf(buf, 1024, "bgamma;%d;%d", j, beta);
        _bgamma[j][beta] = IloBoolVar(_env, buf);
      }
    }
  }
  
  _lambda = IloNumVar3Matrix(_env, _k);
  for (int j = 0; j < _k; ++j)
  {
    _lambda[j] = IloNumVarMatrix(_env, m);
    for (int p = 0; p < m; ++p)
    {
      _lambda[j][p] = IloNumVarArray(_env, _nrSegments);
      for (int alpha = 0; alpha < _nrSegments; ++alpha)
      {
        snprintf(buf, 1024, "lambda;%d;%d;%d", j, p, alpha);
        _lambda[j][p][alpha] = IloNumVar(_env, 0, 1, buf);
      }
    }
  }
  
  _blambda = IloBoolVar3Matrix(_env, _k);
  for (int j = 0; j < _k; ++j)
  {
    _blambda[j] = IloBoolVarMatrix(_env, m);
    for (int p = 0; p < m; ++p)
    {
      _blambda[j][p] = IloBoolVarArray(_env, _nrSegments);
      for (int beta = 0; beta < _nrSegmentBits; ++beta)
      {
        snprintf(buf, 1024, "blambda;%d;%d;%d", j, p, beta);
        _blambda[j][p][beta] = IloBoolVar(_env, buf);
      }
    }
  }
  
  _sigma = IloNumVar4Matrix(_env, n);
  _bsigma = IloBoolVar4Matrix(_env, n);
  _ssigma = IloNumVar4Matrix(_env, n);
  _l = IloNumVar3Matrix(_env, n);
  for (int i = 0; i < n; ++i)
  {
    _sigma[i] = IloNumVar3Matrix(_env, _scriptT[i].size());
    _bsigma[i] = IloBoolVar3Matrix(_env, _scriptT[i].size());
    _ssigma[i] = IloNumVar3Matrix(_env, _scriptT[i].size());
    _l[i] = IloNumVarMatrix(_env, _scriptT[i].size());
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      _sigma[i][t] = IloNumVarMatrix(_env, _k);
      _bsigma[i][t] = IloBoolVarMatrix(_env, _k);
      _ssigma[i][t] = IloNumVarMatrix(_env, _k);
      _l[i][t] = IloNumVarArray(_env, _k);
      for (int j = 0; j < _k; ++j)
      {
        snprintf(buf, 1024, "l;%d;%d;%d", i, t, j);
        _l[i][t][j] = IloNumVar(_env, -1e20, 0, buf);
        
        _sigma[i][t][j] = IloNumVarArray(_env, _nrSegments);
        _ssigma[i][t][j] = IloNumVarArray(_env, _nrSegments);
        for (int alpha = 0; alpha < _nrSegments; ++alpha)
        {
          snprintf(buf, 1024, "sigma;%d;%d;%d;%d", i, t, j, alpha);
          _sigma[i][t][j][alpha] = IloNumVar(_env, 0, 1, buf);
          snprintf(buf, 1024, "ssigma;%d;%d;%d;%d", i, t, j, alpha);
          _ssigma[i][t][j][alpha] = IloNumVar(_env, 0, 1, buf);
        }
        
        _bsigma[i][t][j] = IloBoolVarArray(_env, _nrSegmentBits);
        for (int beta = 0; beta < _nrSegmentBits; ++beta)
        {
          snprintf(buf, 1024, "bsigma;%d;%d;%d;%d", i, t, j, beta);
          _bsigma[i][t][j][beta] = IloBoolVar(_env, 0, 1, buf);
        }
      }
    }
  }
}

void SoftClusterLpCplexExt::initConstraints()
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
              sum += _gamma[j][alpha];
            }
            else
            {
              sum2 += _gamma[j][alpha];
            }
          }
          else if (alpha == _nrSegments)
          {
            unsigned int grayCode = binaryToGray(_nrSegments - 1);
            if (grayCode & (1 << beta))
            {
              sum += _gamma[j][alpha];
            }
            else
            {
              sum2 += _gamma[j][alpha];
            }
          }
          else
          {
            unsigned int grayCode1 = binaryToGray(alpha - 1);
            unsigned int grayCode2 = binaryToGray(alpha);
            if ((grayCode1 & (1 << beta)) && (grayCode2 & (1 << beta)))
            {
              sum += _gamma[j][alpha];
            }
            else if (!(grayCode1 & (1 << beta)) && !(grayCode2 & (1 << beta)))
            {
              sum2 += _gamma[j][alpha];
            }
          }
        }
        _model.add(sum <= _bgamma[j][beta]);
        _model.add(sum2 <= 1 - _bgamma[j][beta]);
        sum.clear();
        sum2.clear();
      }
      
      
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
              sum += _lambda[j][p][alpha];
            }
            else
            {
              sum2 += _lambda[j][p][alpha];
            }
          }
          else if (alpha == _nrSegments)
          {
            unsigned int grayCode = binaryToGray(_nrSegments - 1);
            if (grayCode & (1 << beta))
            {
              sum += _lambda[j][p][alpha];
            }
            else
            {
              sum2 += _lambda[j][p][alpha];
            }
          }
          else
          {
            unsigned int grayCode1 = binaryToGray(alpha - 1);
            unsigned int grayCode2 = binaryToGray(alpha);
            if ((grayCode1 & (1 << beta)) && (grayCode2 & (1 << beta)))
            {
              sum += _lambda[j][p][alpha];
            }
            else if (!(grayCode1 & (1 << beta)) && !(grayCode2 & (1 << beta)))
            {
              sum2 += _lambda[j][p][alpha];
            }
          }
        }
        _model.add(sum <= _blambda[j][p][beta]);
        _model.add(sum2 <= 1 - _blambda[j][p][beta]);
        
        sum.clear();
        sum2.clear();
      }
      
      for (int alpha = 0; alpha < _nrSegments; ++alpha)
      {
        sum += _lambda[j][p][alpha] * _coord[alpha];
        sum2 += _lambda[j][p][alpha];
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
                sum += _sigma[i][t][j][alpha];
              }
              else
              {
                sum2 += _sigma[i][t][j][alpha];;
              }
            }
            else if (alpha == _nrSegments)
            {
              unsigned int grayCode = binaryToGray(_nrSegments - 1);
              if (grayCode & (1 << beta))
              {
                sum += _sigma[i][t][j][alpha];
              }
              else
              {
                sum2 += _sigma[i][t][j][alpha];
              }
            }
            else
            {
              unsigned int grayCode1 = binaryToGray(alpha - 1);
              unsigned int grayCode2 = binaryToGray(alpha);
              if ((grayCode1 & (1 << beta)) && (grayCode2 & (1 << beta)))
              {
                sum += _sigma[i][t][j][alpha];
              }
              else if (!(grayCode1 & (1 << beta)) && !(grayCode2 & (1 << beta)))
              {
                sum2 += _sigma[i][t][j][alpha];
              }
            }
          }
          _model.add(sum <= _bsigma[i][t][j][beta]);
          _model.add(sum2 <= 1 - _bsigma[i][t][j][beta]);
          
          sum.clear();
          sum2.clear();
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
        if (_includePi)
        {
          for (int alpha = 0; alpha < _nrSegments; ++alpha)
          {
            sum += _gamma[j][alpha] * log(_coordPi[alpha] / n);
          }
        }
        
        for (int p = 0; p < m; ++p)
        {
          for (int alpha = 0; alpha < _nrSegments; ++alpha)
          {
            sum += _lambda[j][p][alpha] * _hatG[i][t][j][p][alpha];
          }
        }
        
        _model.add(_l[i][t][j] == sum);
        sum.clear();
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
          if (alpha == 0)
          {
            sum += _sigma[i][t][j][alpha] * (m+1) * log(g_tol.epsilon());
          }
          else
          {
            sum += _sigma[i][t][j][alpha] * log(_coordPi[alpha] / n);
          }
          sum2 += _sigma[i][t][j][alpha];
          sum3 += _ssigma[i][t][j][alpha];
          _model.add(_sigma[i][t][j][alpha] >= _ssigma[i][t][j][alpha]);
        }
        _model.add(sum == _l[i][t][j]);
        _model.add(sum2 == 1);
        _model.add(sum3 == _y[i][t][j]);
        sum.clear();
        sum2.clear();
        sum3.clear();
      }
    }
  }
  
//  // Cluster sizes are nonincreasing
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
  
  /*
   3 #clusters
   2 #samples
   0.140431 0.373801
   0.524305 0.811727
   0.865575 0.96332
   0.46
   0.22
   0.32
  */
//  _model.add(_d[0][0] == 0.140431);
//  _model.add(_d[0][1] == 0.373801);
//  _model.add(_d[1][0] == 0.524305);
//  _model.add(_d[1][1] == 0.811727);
//  _model.add(_d[2][0] == 0.865575);
//  _model.add(_d[2][1] == 0.96332);
//
//  _model.add(_pi[0] == 0.46);
//  _model.add(_pi[1] == 0.22);
//  _model.add(_pi[2] == 0.32);
  
  // PWLA
  IloConstraintArray cuts(_env);



//  try
//  {
//    _cplex.addLazyConstraints(cuts);
//  }
//  catch (IloCplex::InvalidCutException& e)
//  {
//    std::cerr << e.getCut() << std::endl;
//    std::cerr << e.getMessage() << std::endl;
//    exit(1);
//  }
  
  sum.end();
  sum2.end();
  sum3.end();
  sum4.end();
}

void SoftClusterLpCplexExt::initObjective()
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
          sum += log(_coordPi[alpha] / n) * _ssigma[i][t][j][alpha];
        }
      }
    }
  }
  
  _model.add(IloMaximize(_env, sum));
  sum.clear();
  sum.end();
}

bool SoftClusterLpCplexExt::solve(int nrThreads,
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
  
//  _cplex.setParam(IloCplex::PreInd, 0);
  
//  IloFastMutex mutex;
//  _cplex.use(IloCplex::Callback(new (_env) SoftClusterLpCplexCallback<IloCplex::LazyConstraintCallbackI>(_env, m, n, _k, _y, _lambda, _llambda, _dLB, _dUB, &mutex)));
//  _cplex.use(IloCplex::Callback(new (_env) SoftClusterLpCplexCallback<IloCplex::UserCutCallbackI>(_env, m, n, _k, _y, _lambda, _llambda, _dLB, _dUB, &mutex)));
  
  if (!verbose)
  {
    _cplex.setOut(_env.getNullStream());
    _cplex.setError(_env.getNullStream());
    _cplex.setWarning(_env.getNullStream());
  }
  
  try {
    _cplex.solve();
  } catch (IloException& e) {
    std::cerr << e.getMessage() << std::endl;
    return false;
  }
//  _cplex.solve();
  if (_cplex.getStatus() == IloAlgorithm::Infeasible)
  {
    std::cerr << "Error: Infeasible. More clusters or more segments are needed." << std::endl;
    exportModel("/tmp/infeasible.lp");
    return false;
  }
  
  if (_cplex.getSolnPoolNsolns() == 0)
  {
    return false;
  }
  
//  std::cout << "Obj value: " << _cplex.getObjValue() << std::endl;
//  std::cout << "Best obj value: " << _cplex.getBestObjValue() << std::endl;
  
  _solD = DoubleMatrix(_k, DoubleVector(m, 0));
  _solPi = DoubleVector(_k);
  
  std::cout << "d:" << std::endl;
  for (int j = 0; j < _k; ++j)
  {
//    _solPi[j] = _cplex.getValue(_pi[j]);
    for (int p = 0; p < m; ++p)
    {
      _solD[j][p] = _cplex.getValue(_d[j][p]);
      if (_solD[j][p] == g_tol.epsilon())
      {
        _solD[j][p] = 0;
      }
      std::cout << _solD[j][p] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  
//  for (int alpha = 0; alpha < _nrSegments; ++alpha)
//  {
//    std::cout << _llambda[2][0][alpha] << " = " << _cplex.getValue(_llambda[2][0][alpha]) << std::endl;
//  }
  
  std::cout << "pi :";
  for (int j = 0; j < _k; ++j)
  {
    std::cout << " " << _solPi[j];
  }
  std::cout << std::endl;
  
  std::cout << "bpD :" << std::endl;
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < m; ++p)
    {
      for (int alpha = 0; alpha < _nrSegments; ++alpha)
      {
        if (g_tol.nonZero(_cplex.getValue(_lambda[j][p][alpha])))
        {
          std::cout << _lambda[j][p][alpha] << " == " << _cplex.getValue(_lambda[j][p][alpha]) << std::endl;
        }
      }
    }
  }
  std::cout << std::endl;

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
          std::cout << i << "," << t << "," << j << ": " << _solY[i][t][j] << " -- " << _cplex.getValue(_l[i][t][j]) << std::endl;
          
          for (int alpha = 0; alpha < _nrSegments; ++alpha)
          {
            if (g_tol.nonZero(_cplex.getValue(_sigma[i][t][j][alpha])))
            {
              std::cout << _sigma[i][t][j][alpha] << " = " << _cplex.getValue(_sigma[i][t][j][alpha]) << std::endl;
            }
            if (g_tol.nonZero(_cplex.getValue(_ssigma[i][t][j][alpha])))
            {
              std::cout << _ssigma[i][t][j][alpha] << " = " << _cplex.getValue(_sigma[i][t][j][alpha]) << std::endl;
            }
          }
          
          std::cout << std::endl;
          
          DoubleVector h_i;
          for (int p = 0; p < m; ++p)
          {
            h_i.push_back((_solD[j][p] - _numerator[i][t][p]) / _denominator[i][p]);
          }
          
          _solT[i].emplace_back(convertToStateTreeFromSNVF(_R, _scriptT[i][t],
                                                           h_i, i),
                                _solY[i][t][j], t, j);
          
          for (int p = 0; p < m; ++p)
          {
            std::cout << "[" << _dLB[i][t][p] << "," << _dUB[i][t][p] << "] ";
          }
          std::cout << std::endl;
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
  

  if (_includePi)
  {
    for (int j = 0; j < _k; ++j)
    {
      _solPi[j] = _cplex.getValue(_pi[j]);
    }
  }
  
  updateLogLikelihood();
  
  return true;
}
