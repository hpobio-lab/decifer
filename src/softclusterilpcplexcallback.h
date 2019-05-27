/*
 * softclusterilpcplexcallback.h
 *
 *  Created on: 24-oct-2018
 *      Author: M. El-Kebir
 */

#ifndef SOFTCLUSTERILPCPLEXCALLBACK_H
#define SOFTCLUSTERILPCPLEXCALLBACK_H

#include <ilcplex/ilocplex.h>
#include <ilconcert/ilothread.h>
#include "utils.h"

template<class T>
class SoftClusterIlpCplexCallback : public T
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
  typedef IloArray<IloNumVar4Matrix> IloNumVar5Matrix;
  
  SoftClusterIlpCplexCallback(IloEnv env,
                              int m,
                              int n,
                              int k,
                              const IloNumVar3Matrix& y,
                              const IloBoolVar3Matrix& yy,
                              const IloNumVar4Matrix& yyd,
                              const IloNumVarMatrix& d,
                              const DoubleTensor& dLB,
                              const DoubleTensor& dUB,
                              IloFastMutex* pMutex)
  : T(env)
  , _m(m)
  , _n(n)
  , _k(k)
  , _dLB(dLB)
  , _dUB(dUB)
  , _y()
  , _yy()
  , _yyd()
  , _d()
  , _idx_y()
  , _idx_yyd()
  , _idx_d()
  , _maxIterations(100)
  , _currentIterations(0)
  , _nodeId()
  , _pMutex(pMutex)
  {
    int count_y = 0;
    int count_yd = 0;
    _idx_y = Int3Matrix(n);
    _idx_yyd = Int4Matrix(n);
    for (int i = 0; i < n; ++i)
    {
      _idx_y[i] = IntMatrix(y[i].getSize());
      _idx_yyd[i] = Int3Matrix(y[i].getSize());
      
      for (int t = 0; t < y[i].getSize(); ++t)
      {
        _idx_y[i][t] = IntVector(_k);
        _idx_yyd[i][t] = IntMatrix(_k);
        for (int j = 0; j < _k; ++j)
        {
          _idx_y[i][t][j] = count_y;
          ++count_y;
          
          _idx_yyd[i][t][j] = IntVector(_m);
          for (int p = 0; p < _m; ++p)
          {
            _idx_yyd[i][t][j][p] = count_yd;
            ++count_yd;
          }
        }
      }
    }
    
    int count_d = 0;
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
    
    _y = IloNumVarArray(env, count_y);
    _yy = IloBoolVarArray(env, count_y);
    _yyd = IloNumVarArray(env, count_y * m);
    _d = IloNumVarArray(env, count_d);
    
    for (int i = 0; i < n; ++i)
    {
      for (int t = 0; t < y[i].getSize(); ++t)
      {
        for (int j = 0; j < _k; ++j)
        {
          _y[_idx_y[i][t][j]] = y[i][t][j];
          _yy[_idx_y[i][t][j]] = yy[i][t][j];
          for (int p = 0; p < _m; ++p)
          {
            _yyd[_idx_yyd[i][t][j][p]] = yyd[i][t][j][p];
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
    return (new (IloCplex::UserCutCallbackI::getEnv()) SoftClusterIlpCplexCallback(*this));
  }
  
  void main();
  
  void separate();
  
private:
  const int _m;
  const int _n;
  const int _k;
  const DoubleTensor& _dLB;
  const DoubleTensor& _dUB;
  IloNumVarArray _y;
  IloBoolVarArray _yy;
  IloNumVarArray _yyd;
  IloNumVarArray _d;
  Int3Matrix _idx_y;
  Int4Matrix _idx_yyd;
  IntMatrix _idx_d;
  
  int _maxIterations;
  int _currentIterations;
  IloCplex::MIPCallbackI::NodeId _nodeId;
  IloFastMutex* _pMutex;
};

template<class T>
inline void SoftClusterIlpCplexCallback<T>::separate()
{
  assert(_y.getSize() == _yy.getSize());
  
  IloNumArray y_vals = IloNumArray(T::getEnv(), _y.getSize());
  IloNumArray yy_vals = IloNumArray(T::getEnv(), _y.getSize());
  
  _pMutex->lock();
  T::getValues(y_vals, _y);
  T::getValues(yy_vals, _yy);
  
  int count = 0;
  for (int i = 0; i < _n; ++i)
  {
    const int size_T_i = _idx_y[i].size();
    for (int t = 0; t < size_T_i; ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        const int idx_y_itj = _idx_y[i][t][j];
        const double yy = yy_vals[idx_y_itj];
        const double y = y_vals[idx_y_itj];
        if (yy < y)
        {
          ++count;
//          T::add(_yy[idx_y_itj] >= _y[idx_y_itj], IloCplex::UseCutPurge).end();
          T::add(_yy[idx_y_itj] >= _y[idx_y_itj]);
          for (int p = 0; p < _m; ++p)
          {
            const int idx_yyd_itjp = _idx_yyd[i][t][j][p];
            const int idx_d_jp = _idx_d[j][p];
            T::add(_yy[idx_y_itj] * _dLB[i][t][p] <= _yyd[idx_yyd_itjp]);
            T::add(_yyd[idx_yyd_itjp] <= _dUB[i][t][p] * _yy[idx_y_itj]);
            
            T::add(_yyd[idx_yyd_itjp] <= _yy[idx_y_itj]);
            T::add(_yyd[idx_yyd_itjp] <= _d[idx_d_jp]);
            T::add(_yyd[idx_yyd_itjp] >= _yy[idx_y_itj] + _d[idx_d_jp] - 1);

//            T::add(_yy[idx_y_itj] * _dLB[i][t][p] <= _yyd[idx_yyd_itjp], IloCplex::UseCutPurge).end();
//            T::add(_yyd[idx_yyd_itjp] <= _dUB[i][t][p] * _yy[idx_y_itj], IloCplex::UseCutPurge).end();
//
//            T::add(_yyd[idx_yyd_itjp] <= _yy[idx_y_itj], IloCplex::UseCutPurge).end();
//            T::add(_yyd[idx_yyd_itjp] <= _d[idx_d_jp], IloCplex::UseCutPurge).end();
//            T::add(_yyd[idx_yyd_itjp] >= _yy[idx_y_itj] + _d[idx_d_jp] - 1, IloCplex::UseCutPurge).end();
          }
        }
      }
    }
  }
  _pMutex->unlock();
}

template<class T>
inline void SoftClusterIlpCplexCallback<T>::main()
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

#endif // SOFTCLUSTERILPCPLEXCALLBACK_H

