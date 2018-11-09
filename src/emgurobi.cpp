/*
 * emgurobi.cpp
 *
 *  Created on: 16-oct-2018
 *      Author: M. El-Kebir
 */

#include "emgurobi.h"

EMGurobi::EMGurobi(const ReadMatrix& R,
                   int k,
                   int nrSegments)
  : EM(R, k, nrSegments)
  , _env()
  , _pModel(NULL)
  , _lambdaGRB()
  , _fGRB()
{
}

void EMGurobi::initPWLA()
{
  EM::initPWLA();
  
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  assert(!_pModel);
  _pModel = new GRBModel(_env);
  
  char buf[1024];
  
  _lambdaGRB = Var5Matrix(n);
  for (int i = 0; i < n; ++i)
  {
    _lambdaGRB[i] = Var4Matrix(_scriptT[i].size());
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      _lambdaGRB[i][t] = Var3Matrix(_k);
      for (int j = 0; j < _k; ++j)
      {
        if (g_tol.nonZero(_gamma[i][t][j]))
        {
          _lambdaGRB[i][t][j] = VarMatrix(m);
          for (int p = 0; p < m; ++p)
          {
            _lambdaGRB[i][t][j][p] = VarArray(_nrSegments);
            for (int l = 0; l < _nrSegments; ++l)
            {
              snprintf(buf, 1024, "lambda;%d;%d;%d;%d;%d", i, t, j, p, l);
              _lambdaGRB[i][t][j][p][l] = _pModel->addVar(0, 1, 0, GRB_CONTINUOUS, buf);
            }
          }
        }
      }
    }
  }
  
  _fGRB = VarMatrix(_k);
  for (int j = 0; j < _k; ++j)
  {
    _fGRB[j] = VarArray(m);
    for (int p = 0; p < m; ++p)
    {
      snprintf(buf, 1024, "f;%d;%d", j, p);
      _fGRB[j][p] = _pModel->addVar(0, 1, 0, GRB_CONTINUOUS, buf);
    }
  }
  
  _pModel->update();
  
  // initialize constraints
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < m; ++p)
    {
      _pModel->addConstr(_dOverallLB[j][p] <= _fGRB[j][p]);
      _pModel->addConstr(_fGRB[j][p] <= _dOverallUB[j][p]);
    }
  }
  
  GRBLinExpr sum, sum2;
  for (int i = 0; i < n; ++i)
  {
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        if (g_tol.nonZero(_gamma[i][t][j]))
        {
          for (int p = 0; p < m; ++p)
          {
            for (int l = 0; l < _nrSegments; ++l)
            {
              sum += _lambdaGRB[i][t][j][p][l];
              sum2 += _lambdaGRB[i][t][j][p][l] * _x[i][t][j][p][l];
            }
            _pModel->addConstr(sum == 1);
            _pModel->addConstr(_fGRB[j][p] == sum2);
            sum.clear();
            sum2.clear();
          }
        }
      }
    }
  }
  
  GRBLinExpr obj;
  for (int i = 0; i < n; ++i)
  {
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        if (g_tol.nonZero(_gamma[i][t][j]))
        {
          assert(!isnan(_gamma[i][t][j]));
          for (int p = 0; p < m; ++p)
          {
            for (int l = 0; l < _nrSegments; ++l)
            {
              assert(!isnan(_hatG[i][t][j][p][l]));
              obj += _gamma[i][t][j] * _hatG[i][t][j][p][l] * _lambdaGRB[i][t][j][p][l];
            }
          }
        }
      }
    }
  }
  
  try
  {
    _pModel->setObjective(obj, GRB_MAXIMIZE);
  }
  catch (GRBException& e)
  {
    std::cerr << e.getErrorCode() << " " << e.getMessage() << std::endl;
    abort();
  }
  _pModel->update();
}

bool EMGurobi::stepM(int nrThreads)
{
  initPWLA();
  if (nrThreads > 0)
  {
    _pModel->getEnv().set(GRB_IntParam_Threads, nrThreads);
  }
  
  _pModel->getEnv().set(GRB_IntParam_LogToConsole, 0);
  _pModel->optimize();
  int status = _pModel->get(GRB_IntAttr_Status);
  if (status == GRB_INFEASIBLE)
  {
    //    try {
    //        _pModel->computeIIS();
    //    } catch (GRBException e) {
    //        std::cout << e.getMessage() << std::endl;
    //    }
    //    _pModel->write("/tmp/error.ilp");
    return false;
  }
  
  const int m = _R.getNrSamples();
  const int n = _R.getNrCharacters();
  
  try
  {
    for (int j = 0; j < _k; ++j)
    {
      for (int p = 0; p < m; ++p)
      {
        _solD[j][p] = _fGRB[j][p].get(GRB_DoubleAttr_X);
      }
    }
  }
  catch (GRBException& e)
  {
    std::cerr << e.getErrorCode() << " " << e.getMessage() << std::endl;
    abort();
  }
  
  {
    double denominator = 0;
    for (int i = 0; i < n; ++i)
    {
      for (int t = 0; t < _scriptT[i].size(); ++t)
      {
        for (int j = 0; j < _k; ++j)
        {
          denominator += _gamma[i][t][j];
        }
      }
    }
    
    for (int j = 0; j < _k; ++j)
    {
      _solPi[j] = 0;
      for (int i = 0; i < n; ++i)
      {
        for (int t = 0; t < _scriptT[i].size(); ++t)
        {
          _solPi[j] += _gamma[i][t][j];
        }
      }
      _solPi[j] /= denominator;
    }
  }
  
  delete _pModel;
  _pModel = NULL;
  
  return true;
}
