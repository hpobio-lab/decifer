/*
 * clusterilp.cpp
 *
 *  Created on: 10-nov-2017
 *      Author: M. El-Kebir
 */

#include "clusterilp.h"
#include <boost/math/distributions/beta.hpp>

ClusterIlp::ClusterIlp(const ReadMatrix& R,
                       int k,
                       double alpha,
                       ClusterStatisticType statType)
  : Solver(R, k, 0, statType)
  , _alpha(alpha)
  , _vafLB()
  , _vafUB()
  , _scriptTlb()
  , _scriptTub()
  , _scriptTexp()
  , _dcfLB()
  , _dcfUB()
  , _dcfExp()
  , _solT()
  , _solY()
{
}

void ClusterIlp::init()
{
  Solver::init();
  
  initVAFs();
  initStateTrees();
  initVariables();
  initConstraints();
}

void ClusterIlp::writeFeasibleSolutionSpace(std::ostream& out) const
{
  const int m = _R.getNrSamples();
  
  out << "SNV\tstate_tree";
  for (int p = 0; p < m; ++p)
  {
    out << "\tvaf" << p;
  }
  
  for (int p = 0; p < m; ++p)
  {
    out << "\tcf" << p;
  }
  
  for (int p = 0; p < m; ++p)
  {
    out << "\tdcf" << p;
  }
  out << std::endl;
  
  const int n = _R.getNrCharacters();
  for (int i = 0; i < n; ++i)
  {
    for (int t = 0; t < _scriptTexp[i].size(); ++t)
    {
      const StateTree& T_it = _scriptTexp[i][t];
      
      out << i << "\t" << t;
      for (int p = 0; p < m; ++p)
      {
        out << "\t" << T_it.vaf(p);
      }
      for (int p = 0; p < m; ++p)
      {
        out << "\t" << T_it.cf(p);
      }
      for (int p = 0; p < m; ++p)
      {
        out << "\t" << T_it.dcf(p);
      }
      out << std::endl;
    }
  }
}

void ClusterIlp::initVAFs()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  _vafLB = DoubleMatrix(n, DoubleVector(m, 0.));
  _vafUB = DoubleMatrix(n, DoubleVector(m, 1.));
  if (_alpha > 0)
  {
    for (int i = 0; i < n; ++i)
    {
      for (int p = 0; p < m; ++p)
      {
        boost::math::beta_distribution<> b_ip(1 + _R.getVar(p, i), 1 + _R.getRef(p, i));
        double h_ip_lb = boost::math::quantile(b_ip, _alpha / 2);
        double h_ip_ub = boost::math::quantile(b_ip, 1. - _alpha / 2);
        
        _vafLB[i][p] = h_ip_lb;
        _vafUB[i][p] = h_ip_ub;
      }
    }
  }
  else if (_alpha < 0)
  {
    double alpha = -_alpha;
    for (int i = 0; i < n; ++i)
    {
      for (int p = 0; p < m; ++p)
      {
        double vaf = _R.getVAF(p, i);
        double vaf_lb = std::max(0., vaf - alpha);
        double vaf_ub = std::min(1., vaf + alpha);
        
        _vafLB[i][p] = vaf_lb;
        _vafUB[i][p] = vaf_ub;
      }
    }
  }
}

void ClusterIlp::initStateTrees()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  _scriptTlb = _scriptTub = _scriptTexp = StateTreeMatrix(n);
  _dcfLB = DoubleTensor(n);
  _dcfUB = DoubleTensor(n);
  _dcfExp = DoubleTensor(n);
  for (int i = 0; i < n; ++i)
  {
    // 0. obtain L
    const ReadMatrix::CopyNumberStateVector& cnStates = _R.getCopyNumberStates(0, i);
    // + 1 is done on purpose, allowing the following state tree for instance
    // (1,1,0) -> (2,2,0) -> (2,2,1) -> (2,0,1)
    const int maxCopies = _R.getMaxCopies(i) + 1;
    
    IntPairSet L;
    for (const ReadMatrix::CopyNumberState& cnState : cnStates)
    {
      L.insert(IntPair(cnState._x, cnState._y));
    }
    
    // 1. enumerate all state trees
    const auto& res = StateGraph::getStateTrees(L, maxCopies);
    
    // 2. check if ok
    for (const StateEdgeSet& T_it : res)
    {
      IntVector pi;
      std::map<StateGraph::CnaTriple, int> vertices;
      
      StateGraph::CnaTripleSet verticesPreMut, verticesPostMut;
      StateGraph::CnaTriple mutationVertex;
      StateGraph::partition(T_it, verticesPreMut, verticesPostMut, mutationVertex);
      
      assert(mutationVertex._x != -1);
      assert(mutationVertex._y != -1);
      
      vertices[StateGraph::CnaTriple(1, 1, 0)] = 0;
      pi.push_back(-1);
      for (const StateGraph::StateEdge& edge : T_it)
      {
        if (vertices.count(edge.first) == 0)
        {
          vertices[edge.first] = pi.size();
          pi.push_back(-1);
        }
        if (vertices.count(edge.second) == 0)
        {
          vertices[edge.second] = pi.size();
          pi.push_back(vertices[edge.first]);
        }
        pi[vertices[edge.second]] = vertices[edge.first];
      }
      
      StateTree S_it_lb(pi);
      StateTree S_it_ub(pi);
      StateTree S_it_exp(pi);
      
      bool ok = true;
      for (int p = 0; p < m && ok; ++p)
      {
        double muTotal = 0;
        for (const ReadMatrix::CopyNumberState& cnState : _R.getCopyNumberStates(p, i))
        {
          muTotal += (cnState._x + cnState._y) * cnState._mu;
        }
        
        double muMut = 0;
        for (const auto& kv : vertices)
        {
          if (kv.first._x != mutationVertex._x || kv.first._y != mutationVertex._y)
          {
            muMut += kv.first._z * _R.getMu(p, i, kv.first._x, kv.first._y);
          }
        }
        
        double muStar = _R.getMu(p, i, mutationVertex._x, mutationVertex._y);
        
        const double h_ip = _R.getVAF(p, i);
        const double h_ip_lb = std::max(_vafLB[i][p], muMut / muTotal);
        const double h_ip_ub = std::min(_vafUB[i][p], (muMut + muStar) / muTotal);
        
        ok &= h_ip_lb <= h_ip_ub;
        if (ok)
        {
          for (const auto& kv : vertices)
          {
            if (p == 0)
            {
              S_it_lb.setCnState(kv.second, kv.first._x, kv.first._y, kv.first._z);
              S_it_ub.setCnState(kv.second, kv.first._x, kv.first._y, kv.first._z);
              S_it_exp.setCnState(kv.second, kv.first._x, kv.first._y, kv.first._z);
            }
            
            if (kv.first._x != mutationVertex._x || kv.first._y != mutationVertex._y)
            {
              S_it_lb.setMixtureProportion(kv.second, _R.getMu(p, i, kv.first._x, kv.first._y));
              S_it_ub.setMixtureProportion(kv.second, _R.getMu(p, i, kv.first._x, kv.first._y));
              S_it_exp.setMixtureProportion(kv.second, _R.getMu(p, i, kv.first._x, kv.first._y));
            }
            else if (kv.first._x == mutationVertex._x && kv.first._y == mutationVertex._y && kv.first._z == 1)
            {
              double s_mut_lb = h_ip_lb * muTotal - muMut;
              double s_mut_ub = h_ip_ub * muTotal - muMut;
              double s_mut_exp = h_ip * muTotal - muMut;
              
              S_it_lb.setMixtureProportion(kv.second, s_mut_lb);
              S_it_ub.setMixtureProportion(kv.second, s_mut_ub);
              S_it_exp.setMixtureProportion(kv.second, s_mut_exp);
            }
            else
            {
              assert(kv.first._x == mutationVertex._x && kv.first._y == mutationVertex._y && kv.first._z == 0);
              
              double s_nonMut_lb = muStar - (h_ip_ub * muTotal - muMut);
              double s_nonMut_ub = muStar - (h_ip_lb * muTotal - muMut);
              double s_nonMut_exp = muStar - (h_ip * muTotal - muMut);
              
              S_it_lb.setMixtureProportion(kv.second, s_nonMut_lb);
              S_it_ub.setMixtureProportion(kv.second, s_nonMut_ub);
              S_it_exp.setMixtureProportion(kv.second, s_nonMut_exp);
            }
          }
        }
      }
      
      if (ok)
      {
        _scriptTlb[i].push_back(S_it_lb);
        _scriptTub[i].push_back(S_it_ub);
        _scriptTexp[i].push_back(S_it_exp);
        
        _dcfLB[i].push_back(DoubleVector(m, 0));
        _dcfUB[i].push_back(DoubleVector(m, 0));
        _dcfExp[i].push_back(DoubleVector(m, 0));
        for (int p = 0; p < m; ++p)
        {
          const int h_ip = _R.getVAF(p, i);
          
          _dcfLB[i].back()[p] = S_it_lb.dcf(p);
          _dcfUB[i].back()[p] = S_it_ub.dcf(p);
          if (_vafLB[i][p] <= h_ip && h_ip <= _vafUB[i][p])
          {
            if (_statType == CLUSTER_CCF)
            {
              _dcfExp[i].back()[p] = S_it_exp.cf(p);
            }
            else
            {
              _dcfExp[i].back()[p] = S_it_exp.dcf(p);
            }
          }
          else if (h_ip < _vafLB[i][p])
          {
            if (_statType == CLUSTER_CCF)
            {
              _dcfExp[i].back()[p] = S_it_lb.cf(p);
            }
            else
            {
              _dcfExp[i].back()[p] = S_it_lb.dcf(p);
            }
          }
          else
          {
            if (_statType == CLUSTER_CCF)
            {
              _dcfExp[i].back()[p] = S_it_ub.cf(p);
            }
            else
            {
              _dcfExp[i].back()[p] = S_it_ub.dcf(p);
            }
          }
        }
      }
    }
  }
}
