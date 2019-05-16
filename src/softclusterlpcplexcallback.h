/*
 * softclusterlpcplexcallback.h
 *
 *  Created on: 24-oct-2018
 *      Author: M. El-Kebir
 */

#ifndef SOFTCLUSTERLPCPLEXCALLBACK_H
#define SOFTCLUSTERLPCPLEXCALLBACK_H

#include <ilcplex/ilocplex.h>
#include <ilconcert/ilothread.h>
#include "utils.h"

template<class T>
class SoftClusterLpCplexCallback : public T
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
  
  SoftClusterLpCplexCallback(IloEnv env,
                             int m,
                             int n,
                             int k,
                             const IloNumVar3Matrix& y,
                             const IloNumVar5Matrix& lambda,
                             const IloNumVar3Matrix& llambda,
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
  , _lambda(lambda)
  , _llambda(llambda)
  , _idx_y()
  , _maxIterations(100)
  , _currentIterations(0)
  , _nodeId()
  , _pMutex(pMutex)
  {
    int count_y = 0;
    _idx_y = Int3Matrix(n);
    for (int i = 0; i < n; ++i)
    {
      _idx_y[i] = IntMatrix(y[i].getSize());
      
      for (int t = 0; t < y[i].getSize(); ++t)
      {
        _idx_y[i][t] = IntVector(_k);
        for (int j = 0; j < _k; ++j)
        {
          _idx_y[i][t][j] = count_y;
          ++count_y;
        }
      }
    }
    
    _valsY = IloNumArray(env, count_y);
    _y = IloNumVarArray(env, count_y);
    for (int i = 0; i < n; ++i)
    {
      for (int t = 0; t < y[i].getSize(); ++t)
      {
        for (int j = 0; j < _k; ++j)
        {
          const int idx = _idx_y[i][t][j];
          _y[idx] = y[i][t][j];
          _valsY[idx] = 0;
        }
      }
    }
  }
  
  IloCplex::CallbackI *duplicateCallback() const
  {
    return (new (IloCplex::UserCutCallbackI::getEnv()) SoftClusterLpCplexCallback(*this));
  }
  
  void main();
  
  void separate();
  
private:
  const int _m;
  const int _n;
  const int _k;
  const DoubleTensor& _dLB;
  const DoubleTensor& _dUB;
  /// y[i][t][j]
  IloNumVarArray _y;
  ///
  IloNumArray _valsY;
  /// lambda[i][t][j][p][alpha]
  IloNumVar5Matrix _lambda;
  /// llambda[j][p][alpha]
  IloNumVar3Matrix _llambda;

  Int3Matrix _idx_y;
  
  int _maxIterations;
  int _currentIterations;
  IloCplex::MIPCallbackI::NodeId _nodeId;
  IloFastMutex* _pMutex;
};

template<class T>
inline void SoftClusterLpCplexCallback<T>::separate()
{
  IloNumArray y_vals = IloNumArray(T::getEnv(), _y.getSize());
  
  _pMutex->lock();
  T::getValues(y_vals, _y);
  _pMutex->unlock();
  
  const int nrSegments = _llambda[0][0].getSize();
  IloExpr sum(T::getEnv());
  
//  static int count = 0;
  for (int i = 0; i < _n; ++i)
  {
    const int size_T_i = _idx_y[i].size();
    for (int t = 0; t < size_T_i; ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        const int idx_y_itj = _idx_y[i][t][j];
        const double curr_val = y_vals[idx_y_itj];
//        const double prev_val = _valsY[idx_y_itj];
//        if (!g_tol.nonZero(prev_val) && g_tol.nonZero(curr_val))
        if (g_tol.nonZero(curr_val))
        {
//          _valsY[idx_y_itj] = curr_val;
//          ++count;
          for (int p = 0; p < _m; ++p)
          {
            for (int alpha = 0; alpha < nrSegments; ++alpha)
            {
//              sum += _lambda[i][t][j][p][alpha];
              T::add(_llambda[j][p][alpha] >= _lambda[i][t][j][p][alpha], IloCplex::UseCutPurge);
            }
//            T::add(sum == _y[idx_y_itj], IloCplex::UseCutPurge);
            sum.clear();
          }
        }
      }
    }
  }
//  std::cout << count << "/" << _y.getSize() << std::endl;
}

template<class T>
inline void SoftClusterLpCplexCallback<T>::main()
{
//  if (_currentIterations == _maxIterations)
//  {
//    return;
//  }
  
  IloCplex::MIPCallbackI::NodeId newId = T::getNodeId();
  if (_currentIterations == 0 || newId != _nodeId)
  {
    _nodeId = newId;
  }
  _currentIterations++;
  
  separate();
}

#endif // SOFTCLUSTERILPCPLEXCALLBACK_H

