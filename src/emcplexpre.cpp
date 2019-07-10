/*
 * emcplexpre.cpp
 *
 *  Created on: 3-jul-2019
 *      Author: M. El-Kebir
 */

#include "emcplexpre.h"

EMCplexPre::EMCplexPre(const ReadMatrix& R,
                       const IntMatrix& preClustering,
                       int k,
                       int nrSegments,
                       ClusterStatisticType statType,
                       double precisionBetaBin,
                       bool forceTruncal)
  : EMPre(R, preClustering, k, nrSegments, statType, precisionBetaBin, forceTruncal)
  , _env()
  , _modelE(_env)
  , _cplexE(_modelE)
  , _modelM(_env)
  , _cplexM(_modelM)
  , _lambda()
  , _d()
  , _y()
  , _pi()
{
}

void EMCplexPre::initE()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  const int nrPreClusters = _preClustering.size();

  // set character to pre clusters map
  IntVector charToPreCluster(n, -1);
  for (int ii = 0; ii < nrPreClusters; ++ii)
  {
    for (int i : _preClustering[ii])
    {
      charToPreCluster[i] = ii;
    }
  }
  
  _env.end();
  _env = IloEnv();
  _modelE = IloModel(_env);
  _cplexE = IloCplex(_modelE);
  
  IloExpr sum(_env);
  IloExpr sum2(_env);
  
  char buf[1024];
  _y = IloNumVar3Matrix(_env, nrPreClusters);
  for (int ii = 0; ii < nrPreClusters; ++ii)
  {
    assert(!_preClustering[ii].empty());
    const int size_T_ii = _scriptT[_preClustering[ii].front()].size();
    
    _y[ii] = IloNumVarMatrix(_env, size_T_ii);
    for (int t = 0; t < size_T_ii; ++t)
    {
      _y[ii][t] = IloNumVarArray(_env, _k);
      for (int j = 0; j < _k; ++j)
      {
        snprintf(buf, 1024, "y;%d;%d;%d", ii, t, j);
        _y[ii][t][j] = IloNumVar(_env, 0, 1, buf);
      }
    }
  }
  
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
      _gamma[j][alpha] = IloNumVar(_env, 0, 1, buf);
    }
  }
  
  // one state tree and cluster per SNV
  for (int ii = 0; ii < nrPreClusters; ++ii)
  {
    assert(!_preClustering[ii].empty());
    const int size_T_ii = _scriptT[_preClustering[ii].front()].size();
    
    for (int t = 0; t < size_T_ii; ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        sum += _y[ii][t][j];
      }
    }
    
    _modelE.add(sum == 1);
    sum.clear();
  }
  
  for (int j = 0; j < _k; ++j)
  {
    for (int i = 0; i < n; ++i)
    {
      const int ii = charToPreCluster[i];
      const int size_T_i = _scriptT[i].size();
      for (int t = 0; t < size_T_i; ++t)
      {
        sum += _y[ii][t][j];
      }
    }
    _modelE.add(_pi[j] == sum);
    sum.clear();
  }
  
  DoubleVector coordPi(_nrSegments);
  coordPi[0] = g_tol.epsilon();
  const double delta2 = (n) / (float) (_nrSegments - 1);
  for (int alpha = 1; alpha < _nrSegments; ++alpha)
  {
    coordPi[alpha] = delta2 * alpha;
  }
  
  // compute hatPi
  DoubleVector hatPi(_nrSegments, 0);
  for (int alpha = 0; alpha < _nrSegments; ++alpha)
  {
    double likelihood = log(coordPi[alpha]);
    hatPi[alpha] = likelihood;
  }
  
  for (int j = 0; j < _k; ++j)
  {
    for (int alpha = 0; alpha < _nrSegments; ++alpha)
    {
      sum += _gamma[j][alpha] * coordPi[alpha];
      sum2 += _gamma[j][alpha];
    }
    _modelE.add(sum == _pi[j]);
    _modelE.add(sum2 == 1);
    sum.clear();
    sum2.clear();
  }
  
  IloExpr obj(_env);
  const double log_n = log(n);
  for (int j = 0; j < _k; ++j)
  {
    for (int alpha = 0; alpha < _nrSegments; ++alpha)
    {
      obj += (hatPi[alpha] - log_n) * _gamma[j][alpha];
    }
  }
  
  const double infeasible_log = -1e300;
  BoolTensor feasible(nrPreClusters);
  for (int ii = 0; ii < nrPreClusters; ++ii)
  {
    const int size_T_ii = _scriptT[_preClustering[ii].front()].size();
    feasible[ii] = BoolMatrix(size_T_ii, BoolVector(_k, true));
    for (int i : _preClustering[ii])
    {
      for (int t = 0; t < _scriptT[i].size(); ++t)
      {
        for (int j = 0; j < _k; ++j)
        {
          for (int p = 0; p < m; ++p)
          {
            const double h = (_solD[j][p] - _numerator[i][t][p]) / _denominator[i][p];
            const int var_ip = _R.getVar(p, i);
            const int ref_ip = _R.getRef(p, i);
            
            const double likelihood = getLogLikelihood(var_ip, ref_ip, h);
            if (likelihood < infeasible_log || isnan(likelihood) || _solD[j][p] < _dLB[i][t][p] || _solD[j][p] > _dUB[i][t][p])
            {
              if (feasible[ii][t][j])
              {
                _modelE.add(_y[ii][t][j] == 0);
              }
              feasible[ii][t][j] = false;
            }
          }
        }
      }
    }
  }
  
  for (int i = 0; i < n; ++i)
  {
    const int ii = charToPreCluster[i];
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        if (feasible[ii][t][j])
        {
          obj += - _y[ii][t][j] * log(_scriptT[i].size());
          for (int p = 0; p < m; ++p)
          {
            const double h = (_solD[j][p] - _numerator[i][t][p]) / _denominator[i][p];
            const int var_ip = _R.getVar(p, i);
            const int ref_ip = _R.getRef(p, i);
            
            const double likelihood = getLogLikelihood(var_ip, ref_ip, h);
            assert(likelihood >= infeasible_log);
            obj += _y[ii][t][j] * likelihood;
          }
        }
      }
    }
  }
  
  for (int ii = 0; ii < nrPreClusters; ++ii)
  {
    assert(!_preClustering[ii].empty());
    const int size_T_ii = _scriptT[_preClustering[ii].front()].size();
    
    IntSet dominatingStateTrees;
    
    for (int i : _preClustering[ii])
    {
      Digraph G;
      IntNodeMap node2idx(G);
      DoubleVectorNodeMap node2dcfExp(G);
      for (int t = 0; t < size_T_ii; ++t)
      {
        Node v = G.addNode();
        node2idx[v] = t;
        node2dcfExp[v] = DoubleVector(m, 0);
        
        // compute dcfExp
        
        DoubleVector f_i;
        for (int p = 0; p < m; ++p)
        {
          f_i.push_back(_R.getVAF(p, i));
        }
        
        StateTree TT = Solver::convertToStateTreeFromSNVF(_R, _scriptT[i][t], f_i, i);
        for (int p = 0; p < m; ++p)
        {
          node2dcfExp[v][p] = TT.dcf(p);
        }
      }
      
      for (NodeIt v(G); v != lemon::INVALID; ++v)
      {
        bool dominating = true;
        for (NodeIt w(G); w != lemon::INVALID; ++w)
        {
          if (v == w)
            continue;
          
          bool edge = true;
          for (int p = 0; p < m; ++p)
          {
            edge &= g_tol.less(node2dcfExp[v][p], node2dcfExp[w][p]) && feasible[ii][node2idx[w]][p];
          }
          
          if (edge)
          {
            G.addArc(v, w);
            dominating = false;
          }
        }
        
        if (dominating)
        {
          dominatingStateTrees.insert(node2idx[v]);
        }
      }
    }
    
//    std::cerr << "SNV precluster: " << ii << " -- number of state trees: "
//              << dominatingStateTrees.size() << "/" << size_T_ii << std::endl;
    if (!dominatingStateTrees.empty())
    {
      for (int t = 0; t < size_T_ii; ++t)
      {
        if (dominatingStateTrees.count(t) == 0)
        {
          for (int j = 0; j < _k; ++j)
          {
            _modelE.add(_y[ii][t][j] == 0);
          }
        }
      }
    }
  }
  
  _modelE.add(IloMaximize(_env, obj));
  obj.end();
  sum.end();
  sum2.end();
}

void EMCplexPre::initM()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  updatePWLA();
  
  _env.end();
  _env = IloEnv();
  _modelM = IloModel(_env);
  _cplexM = IloCplex(_modelM);
  
  char buf[1024];
  _lambda = IloNumVar3Matrix(_env, _k);
  for (int j = 0; j < _k; ++j)
  {
    _lambda[j] = IloNumVarMatrix(_env, m);
    for (int p = 0; p < m; ++p)
    {
      _lambda[j][p] = IloNumVarArray(_env, _nrSegments);
      for (int l = 0; l < _nrSegments; ++l)
      {
        snprintf(buf, 1024, "lambda;%d;%d;%d", j, p, l);
        _lambda[j][p][l] = IloNumVar(_env, 0, 1, buf);
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
  
  // initialize constraints
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < m; ++p)
    {
      _modelM.add(_dOverallLB[j][p] <= _d[j][p]);
      _modelM.add(_d[j][p] <= _dOverallUB[j][p]);
    }
  }
  
  // compute coord
  DoubleVector coord(_nrSegments);
  coord[0] = 0;
  coord[_nrSegments - 1] = 1;
  
  const double delta = 1. / (_nrSegments - 1);
  for (int alpha = 1; alpha < _nrSegments - 1; ++alpha)
  {
    coord[alpha] = delta * alpha;
  }
  
  IloExpr sum(_env), sum2(_env);
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < m; ++p)
    {
      for (int l = 0; l < _nrSegments; ++l)
      {
        sum += _lambda[j][p][l];
        sum2 += _lambda[j][p][l] * coord[l];
      }
      _modelM.add(sum == 1);
      _modelM.add(_d[j][p] == sum2);
      sum.clear();
      sum2.clear();
    }
  }
  
//  _modelM.add(_d[0][0] == 0.548387);
//  _modelM.add(_d[0][1] == 0.516129);
//  _modelM.add(_d[0][2] == 0.67);
//  _modelM.add(_d[1][0] == 0.322581);
//  _modelM.add(_d[1][1] == 0.0322581);
//  _modelM.add(_d[1][2] == 0.612903);
//  _modelM.add(_d[2][0] == 0.0322581);
//  _modelM.add(_d[2][1] == 0.258065);
//  _modelM.add(_d[2][2] == 0.0322581);
//  _modelM.add(_d[3][0] == 0.0967742);
//  _modelM.add(_d[3][1] == 0.0322581);
//  _modelM.add(_d[3][2] == 0.0967742);
//  _modelM.add(_d[4][0] == 0.258065);
//  _modelM.add(_d[4][1] == 0.387097);
//  _modelM.add(_d[4][2] == 0.387097);
  

  IloExpr obj(_env);
  for (int i = 0; i < n; ++i)
  {
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        if (g_tol.nonZero(_solY[i][t][j]))
        {
          obj += _solY[i][t][j] * log(_solPi[j]) - _solY[i][t][j] * log(_scriptT[i].size());
          assert(!isnan(_solY[i][t][j]));
          for (int p = 0; p < m; ++p)
          {
            for (int l = 0; l < _nrSegments; ++l)
            {
              assert(!isnan(_hatG[i][t][j][p][l]));
              obj += _solY[i][t][j] * _hatG[i][t][j][p][l] * _lambda[j][p][l];
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
        _modelM.add(_d[j][p] <= _d[0][p]);
      }
    }
  }
  
  _modelM.add(IloMaximize(_env, obj));
  obj.end();
  sum.end();
  sum2.end();
}

bool EMCplexPre::stepM(int nrThreads,
                       bool verbose)
{
  const int m = _R.getNrSamples();
  
  initM();
  if (nrThreads > 0)
  {
    _cplexM.setParam(IloCplex::Threads, nrThreads);
  }
  
  if (!verbose)
  {
    _cplexM.setOut(_env.getNullStream());
    _cplexM.setError(_env.getNullStream());
    _cplexM.setWarning(_env.getNullStream());
  }
  
  _cplexM.solve();
//  std::cerr << "Likelihood: " << _cplex.getObjValue() << std::endl;
  if (_cplexM.getStatus() == IloAlgorithm::Infeasible)
  {
    _cplexM.exportModel("/tmp/infeasible.lp");
    _cplexM.end();
    return false;
  }
  
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < m; ++p)
    {
      _solD[j][p] = _cplexM.getValue(_d[j][p]);
      
      if (_solD[j][p] == g_tol.epsilon())
      {
        _solD[j][p] = 0;
      }
    }
  }
  
  _cplexM.end();
  
  return true;
}

bool EMCplexPre::stepE(int nrThreads,
                       bool verbose)
{
  const int m = _R.getNrSamples();
  const int n = _R.getNrCharacters();
  const int nrPreClusters = _preClustering.size();
  
  initE();
  if (nrThreads > 0)
  {
    _cplexE.setParam(IloCplex::Threads, nrThreads);
  }
  
  if (!verbose)
  {
    _cplexE.setOut(_env.getNullStream());
    _cplexE.setError(_env.getNullStream());
    _cplexE.setWarning(_env.getNullStream());
  }
  
  _cplexE.solve();
  //  std::cerr << "Likelihood: " << _cplex.getObjValue() << std::endl;
  if (_cplexE.getStatus() == IloAlgorithm::Infeasible)
  {
    _cplexE.exportModel("/tmp/infeasible.lp");
    _cplexE.end();
    return false;
  }
  
  // TODO: round down to zero?
  for (int j = 0; j < _k; ++j)
  {
    _solPi[j] = _cplexE.getValue(_pi[j]);
    _solPi[j] /= n;
  }
  
  _solT = PosteriorStateTreeMatrix(n);
  for (int ii = 0; ii < nrPreClusters; ++ii)
  {
    const int scriptT_size = _scriptT[_preClustering[ii].front()].size();
    for (int t = 0; t < scriptT_size; ++t)
    {
      for (int j = 0; j < _k; ++j)
      {
        double y_iitj = _cplexE.getValue(_y[ii][t][j]);
        if (!g_tol.nonZero(y_iitj))
          y_iitj = 0;
        
        for (int i : _preClustering[ii])
        {
          _solY[i][t][j] = y_iitj;
          
          DoubleVector f_i;
          for (int p = 0; p < m; ++p)
          {
            f_i.push_back((_solD[j][p] - _numerator[i][t][p]) / _denominator[i][p]);
          }
          
          _solT[i].emplace_back(convertToStateTreeFromSNVF(_R, _scriptT[i][t],
                                                           f_i, i),
                                y_iitj, t, j);
        }
      }
    }
  }
  
  // sort and update _solZ
  for (int i = 0; i < n; ++i)
  {
    std::sort(_solT[i].begin(),
              _solT[i].end(),
              [](const PosteriorStateTree& a,
                 const PosteriorStateTree& b)
                {
                  return a._gamma > b._gamma;
                });
    _solZ[i] = _solT[i][0]._j;
  }
  
  _cplexE.end();
  
  return true;
}
