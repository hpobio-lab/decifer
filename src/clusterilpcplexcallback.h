/*
 * clusterilpcplexcallback.h
 *
 *  Created on: 24-oct-2018
 *      Author: M. El-Kebir
 */

#ifndef CLUSTERILPCPLEXCALLBACK_H
#define CLUSTERILPCPLEXCALLBACK_H

#include <ilcplex/ilocplex.h>
#include <ilconcert/ilothread.h>
#include "utils.h"

template<class T>
class ClusterIlpCplexCallback : public T
{
private:
  typedef std::vector<int> IntVector;
  typedef std::vector<IntVector> IntMatrix;
  typedef std::vector<IntMatrix> Int3Matrix;
  typedef std::vector<Int3Matrix> Int4Matrix;
  
public:
  typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;
  typedef IloArray<IloBoolVarMatrix> IloBoolVar3Matrix;
  typedef IloArray<IloNumVarArray> IloNumVarMatrix;
  typedef IloArray<IloNumVarMatrix> IloNumVar3Matrix;
  typedef IloArray<IloNumVar3Matrix> IloNumVar4Matrix;
  
  ClusterIlpCplexCallback(IloEnv env,
                          int m,
                          int n,
                          int k,
                          const IloBoolVar3Matrix& y,
                          const IloNumVar4Matrix& yd,
                          const IloNumVar4Matrix& w,
                          const IloNumVarMatrix& d,
                          const DoubleTensor& dLB,
                          const DoubleTensor& dUB,
                          const DoubleTensor& dExp,
                          IloFastMutex* pMutex)
    : T(env)
    , _m(m)
    , _n(n)
    , _k(k)
    , _dLB(dLB)
    , _dUB(dUB)
    , _dExp(dExp)
    , _y()
    , _yd()
    , _w()
    , _d()
    , _idx_y()
    , _idx_yd()
    , _idx_w()
    , _idx_d()
    , _maxIterations(100)
    , _currentIterations(0)
    , _nodeId()
    , _pMutex(pMutex)
  {
    int count_y = 0;
    int count_yd = 0;
    int count_w = 0;
    int count_d = 0;
    _idx_y = Int3Matrix(n);
    _idx_yd = Int4Matrix(n);
    _idx_w = Int4Matrix(n);
    for (int i = 0; i < n; ++i)
    {
      _idx_y[i] = IntMatrix(y[i].getSize());
      _idx_yd[i] = Int3Matrix(y[i].getSize());
      _idx_w[i] = Int3Matrix(y[i].getSize());
      
      for (int t = 0; t < y[i].getSize(); ++t)
      {
        _idx_y[i][t] = IntVector(_k);
        _idx_yd[i][t] = IntMatrix(_k);
        _idx_w[i][t] = IntMatrix(_k);
        for (int j = 0; j < _k; ++j)
        {
          _idx_y[i][t][j] = count_y;
          ++count_y;
          
          _idx_yd[i][t][j] = IntVector(_m);
          _idx_w[i][t][j] = IntVector(_m);
          for (int p = 0; p < _m; ++p)
          {
            _idx_yd[i][t][j][p] = count_yd;
            _idx_w[i][t][j][p] = count_w;
            ++count_yd;
            ++count_w;
          }
        }
      }
    }
    
    _idx_d = IntMatrix(_k);
    for (int j = 0; j < _k; ++j)
    {
      _idx_d[j] = IntVector(_m);
      for (int p = 0; p < m; ++p)
      {
        _idx_d[j][p] = count_d;
        ++count_d;
      }
    }
    
    _y = IloBoolVarArray(env, count_y);
    _yd = IloNumVarArray(env, count_y * m);
    _w = IloNumVarArray(env, count_y * m);
    _d = IloNumVarArray(env, count_d);
    
    for (int i = 0; i < n; ++i)
    {
      for (int t = 0; t < y[i].getSize(); ++t)
      {
        for (int j = 0; j < _k; ++j)
        {
          _y[_idx_y[i][t][j]] = y[i][t][j];
          for (int p = 0; p < _m; ++p)
          {
            _yd[_idx_yd[i][t][j][p]] = yd[i][t][j][p];
            _w[_idx_w[i][t][j][p]] = w[i][t][j][p];
          }
        }
      }
    }
    
    for (int j = 0; j < _k; ++j)
    {
      for (int p = 0; p < _m; ++p)
      {
        _d[_idx_d[j][p]] = d[j][p];
      }
    }
  }
  
  IloCplex::CallbackI *duplicateCallback() const
  {
    return (new (IloCplex::UserCutCallbackI::getEnv()) ClusterIlpCplexCallback(*this));
  }
  
  void main();
  
  void separate();
  
private:
  const int _m;
  const int _n;
  const int _k;
  const DoubleTensor& _dLB;
  const DoubleTensor& _dUB;
  const DoubleTensor& _dExp;
  IloBoolVarArray _y;
  IloNumVarArray _yd;
  IloNumVarArray _w;
  IloNumVarArray _d;
  Int3Matrix _idx_y;
  Int4Matrix _idx_yd;
  Int4Matrix _idx_w;
  IntMatrix _idx_d;

  int _maxIterations;
  int _currentIterations;
  IloCplex::MIPCallbackI::NodeId _nodeId;
  IloFastMutex* _pMutex;
};

template<class T>
inline void ClusterIlpCplexCallback<T>::separate()
{
  IloNumArray y_vals = IloNumArray(T::getEnv(), _y.getSize());
  IloNumArray yd_vals = IloNumArray(T::getEnv(), _yd.getSize());
  IloNumArray w_vals = IloNumArray(T::getEnv(), _w.getSize());
  IloNumArray d_vals = IloNumArray(T::getEnv(), _d.getSize());
  
  _pMutex->lock();
  T::getValues(y_vals, _y);
  T::getValues(yd_vals, _yd);
  T::getValues(w_vals, _w);
  T::getValues(d_vals, _d);
  _pMutex->unlock();
  
  for (int i = 0; i < _n; ++i)
  {
    const int size_T_i = _idx_y[i].size();
    for (int t = 0; t < size_T_i; ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        for (int p = 0; p < _m; ++p)
        {
          int idx_yd_itjp = _idx_yd[i][t][j][p];
          int idx_y_itj = _idx_y[i][t][j];
          int idx_d_jp = _idx_d[j][p];
          if (g_tol.different(yd_vals[idx_yd_itjp], y_vals[idx_y_itj] * d_vals[idx_d_jp]))
          {
            T::add(_yd[idx_yd_itjp] <= _y[idx_y_itj], IloCplex::UseCutPurge).end();
            T::add(_yd[idx_yd_itjp] <= _d[idx_d_jp], IloCplex::UseCutPurge).end();
            T::add(_yd[idx_yd_itjp] >= _y[idx_y_itj] + _d[idx_d_jp] - 1, IloCplex::UseCutPurge).end();
          }
//          if (g_tol.less(yd_vals[idx_yd_itjp], _dLB[i][t][p] * y_vals[idx_y_itj]))
//          {
//            T::add(_yd[idx_yd_itjp] >= _dLB[i][t][p] * _y[idx_y_itj], IloCplex::UseCutPurge).end();
//          }
//          if (g_tol.less(_dUB[i][t][p] * y_vals[idx_y_itj], yd_vals[idx_yd_itjp]))
//          {
//            T::add(_yd[idx_yd_itjp] <= _dUB[i][t][p] * _y[idx_y_itj], IloCplex::UseCutPurge).end();
//          }
        }
      }
    }
  }
}

template<class T>
inline void ClusterIlpCplexCallback<T>::main()
{
  if (_currentIterations == _maxIterations)
  {
    return;
  }
  
  IloCplex::MIPCallbackI::NodeId newId = T::getNodeId();
  if (_currentIterations == 0 || newId != _nodeId)
  {
    _nodeId = newId;
  }
  _currentIterations++;
  
  separate();
}

#endif // CLUSTERILPCPLEXCALLBACK_H
