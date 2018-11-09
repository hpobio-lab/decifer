/*
 * generativemodel.cpp
 *
 *  Created on: 23-oct-2017
 *      Author: M. El-Kebir
 */

#include "generativemodel.h"
#include "statetreessampler.h"
#include <iomanip>

GenerativeModel::GenerativeModel(int m,
                                 int n,
                                 int k,
                                 int maxXY,
                                 int maxCN,
                                 int coverage)
  : _m(m)
  , _n(n)
  , _k(k)
  , _maxXY(maxXY)
  , _maxCN(maxCN)
  , _coverage(coverage)
  , _purityVector(m)
  , _f(m, DoubleVector(k, 0))
  , _z(n, -1)
  , _R(m, n)
  , _S()
  , _clsT()
  , _rootClsT(lemon::INVALID)
  , _nodeToClusterClsT(_clsT, -1)
  , _clusterToNodeClsT(_k, lemon::INVALID)
  , _pPhyloT(NULL)
{
}

GenerativeModel::~GenerativeModel()
{
  delete _pPhyloT;
}

void GenerativeModel::writeSamplePurities(std::ostream& out) const
{
  out << _m << " #samples" << std::endl;
  for (int p = 0; p < _m; ++p)
  {
    out << _purityVector[p] << std::endl;
  }
}

void GenerativeModel::generateRootedTree(const int k,
                                         Digraph& T,
                                         Node& root)
{
  if (k == 0)
  {
    return;
  }
  else if (k == 1)
  {
    root = T.addNode();
  }
  else if (k == 2)
  {
    root = T.addNode();
    Node child = T.addNode();
    T.addArc(root, child);
  }
  else
  {
    IntVector seq(k - 2, 0);
    std::uniform_int_distribution<> unif_0k(0, k - 1);
    
    for (int i = 0; i < k - 2; ++i)
    {
      seq[i] = unif_0k(g_rng);
    }
    
    Graph unrootedT;
    converPruferToTree(seq, unrootedT);
    
    int rootIdx = unif_0k(g_rng);
    Graph::Node rootUnrootedT;
    for (Graph::NodeIt v(unrootedT); v != lemon::INVALID; ++v)
    {
      if (--rootIdx < 0)
      {
        rootUnrootedT = v;
        break;
      }
    }
    
    Graph::NodeMap<bool> visited(unrootedT, false);
    visited[rootUnrootedT] = true;
    T.clear();
    root = T.addNode();
    rootTree(unrootedT, rootUnrootedT, visited, T, root);
  }
}

void GenerativeModel::rootTree(const Graph& unrootedT,
                               const Graph::Node vertexUnrootedT,
                               Graph::NodeMap<bool>& visited,
                               Digraph& rootedT,
                               Node vertexRootedT)
{
  for (Graph::IncEdgeIt a(unrootedT, vertexUnrootedT); a != lemon::INVALID; ++a)
  {
    Graph::Node v = unrootedT.oppositeNode(vertexUnrootedT, a);
    if (visited[v])
      continue;
    
    visited[v] = true;
    Node vv = rootedT.addNode();
    rootedT.addArc(vertexRootedT, vv);
    
    rootTree(unrootedT, v, visited, rootedT, vv);
  }
}

void GenerativeModel::converPruferToTree(const IntVector& pruferSequence,
                                         Graph& T)
{
  const int k = pruferSequence.size() + 2;
  
  for (int a_i : pruferSequence)
  {
    assert(0 <= a_i && a_i < k);
  }
  
  typedef std::vector<Graph::Node> NodeVector;
  NodeVector nodes;
  
  // The tree will have k nodes.
  for (int i = 0; i < k; ++i)
  {
    nodes.push_back(T.addNode());
  }
  
  // For each node set its degree to the number of times it
  // appears in the sequence plus 1
  typedef Graph::NodeMap<int> IntNodeMap;
  IntNodeMap deg(T, 1);
  for (int i : pruferSequence)
  {
    ++deg[nodes[i]];
  }
  
  // Next, for each number i in the sequence, find the first
  // (lowest-numbered) node, v_j, with degree equal to 1,
  // add the edge (v_i, v_j) to the tree,
  // and decrement the degrees of v_i and v_j
  for (int i : pruferSequence)
  {
    Graph::Node v_i = nodes[i];
    for (Graph::Node v_j : nodes)
    {
      if (deg[v_j] == 1)
      {
        T.addEdge(v_i, v_j);
        --deg[v_i];
        --deg[v_j];
        break;
      }
    }
  }
  
  // Two nodes with degree 1 remain (call them u, v).
  // Add the edge (u,v) to the tree
  Graph::Node v_i = lemon::INVALID;
  Graph::Node v_j = lemon::INVALID;
  for (Graph::Node v_l : nodes)
  {
    if (deg[v_l] == 1)
    {
      if (v_i == lemon::INVALID)
      {
        v_i = v_l;
      }
      else
      {
        assert(v_j == lemon::INVALID);
        v_j = v_l;
      }
    }
  }
  T.addEdge(v_i, v_j);
}

void GenerativeModel::generateDCF(Node v_j,
                                  const DoubleMatrix& usages)
{
  const int j = _nodeToClusterClsT[v_j];
  for (OutArcIt a_jl(_clsT, v_j); a_jl != lemon::INVALID; ++a_jl)
  {
    Node v_l = _clsT.target(a_jl);
    generateDCF(v_l, usages);
  }
  
  for (int p = 0; p < _m; ++p)
  {
    double sum = 0;
    for (OutArcIt a_jl(_clsT, v_j); a_jl != lemon::INVALID; ++a_jl)
    {
      Node v_l = _clsT.target(a_jl);
      int l = _nodeToClusterClsT[v_l];
      sum += _f[p][l];
    }
    _f[p][j] = sum + usages[p][j];
  }
}

void GenerativeModel::writeDOT(std::ostream& out) const
{
  out << "digraph T {" << std::endl;
  
  for (NodeIt v_j(_clsT); v_j != lemon::INVALID; ++v_j)
  {
    const int j = _nodeToClusterClsT[v_j];
    out << "\t" << _clsT.id(v_j) << " [label=\"" << _nodeToClusterClsT[v_j];
    for (int p = 0; p < _m; ++p)
    {
       out << "\\n" << std::setprecision(2) << _f[p][j];
    }
    out << "\"]" << std::endl;
  }
  
  for (ArcIt a_jl(_clsT); a_jl != lemon::INVALID; ++a_jl)
  {
    Node v_j = _clsT.source(a_jl);
    Node v_l = _clsT.target(a_jl);
    
    out << "\t" << _clsT.id(v_j) << " -> " << _clsT.id(v_l) << std::endl;
  }
  
  out << "}" << std::endl;
}

void GenerativeModel::writeMutationProperties(std::ostream& out) const
{
  out << "SNV\tcluster";
  for (int p = 0; p < _m; ++p)
  {
    out << "\tvaf" << p << "\tcf" << p << "\tdcf" << p;
  }
  out << std::endl;

  for (int i = 0; i < _n; ++i)
  {
    const StateTree& S_i = _S[i];
    out << i << "\t" << _z[i];

    for (int p = 0; p < _m; ++p)
    {
      out << "\t" << S_i.vaf(p)
          << "\t" << S_i.cf(p)
          << "\t" << S_i.dcf(p);
    }
    out << std::endl;
  }
}

void GenerativeModel::generate()
{
  std::uniform_real_distribution<> unif_01(0, 1);
  std::uniform_int_distribution<> unif_0k(0, _k - 1);
  std::poisson_distribution<> poisson(_coverage);
  sftrabbit::beta_distribution<> beta(36, 4);
  std::gamma_distribution<> gamma(1, 1);
  
  // 0a. generate cluster tree
  generateRootedTree(_k, _clsT, _rootClsT);
  int j = 0;
  for (NodeIt v_j(_clsT); v_j != lemon::INVALID; ++v_j, ++j)
  {
    _clusterToNodeClsT[j] = v_j;
    _nodeToClusterClsT[v_j] = j;
  }
  
  // 0b. generate purities
  for (int p = 0; p < _m; ++p)
  {
    _purityVector[p] = beta(g_rng);
  }
  
  // 0c. generate usage matrix
  DoubleMatrix clusterUsages(_m, DoubleVector(_k, 0));
  for (int p = 0; p < _m; ++p)
  {
    double sum = 0;
    for (int j = 0; j < _k; ++j)
    {
      clusterUsages[p][j] = gamma(g_rng);
      sum += clusterUsages[p][j];
    }
    
    for (int j = 0; j < _k; ++j)
    {
      clusterUsages[p][j] /= sum;
      clusterUsages[p][j] *= _purityVector[p];
    }
  }
  
  // 0d. update DCFs
  generateDCF(_rootClsT, clusterUsages);
  
  // 1. generate cluster assignments
  for (int i = 0; i < _n; ++i)
  {
    _z[i] = unif_0k(g_rng);
  }
  
  // 2. generate state trees
  StateTreesSampler sampler(_maxXY, _n, _maxCN);
  _S = sampler.sample();
  
  // 3. generate phylogenetic tree
  _pPhyloT = new PhylogeneticTree(_S, _z,
                                  _clsT, _rootClsT,
                                  _nodeToClusterClsT,
                                  _clusterToNodeClsT);
  
  // 4. generate samples from phylogenetic tree
  _pPhyloT->sample(_purityVector, clusterUsages, _S);
  
  // 5. compute DCFs and generate read counts
  for (int i = 0; i < _n; ++i)
  {
    const StateTree& S_i = _S[i];
    for (int p = 0; p < _m; ++p)
    {
//      _f[p][_z[i]] = S_i.dcf(p);
      double h_pi = S_i.vaf(p);
      int d_pi = poisson(g_rng);
      std::binomial_distribution<> binom(d_pi, h_pi);
      int a_pi = binom(g_rng);
      _R.setVar(p, i, a_pi);
      _R.setRef(p, i, d_pi - a_pi);
      
      std::map<IntPair, double> muMap;
      for (NodeIt v_j(S_i.S()); v_j != lemon::INVALID; ++v_j)
      {
        int j = S_i.state(v_j);
        if (S_i.isPresent(j))
        {
          IntPair xy(S_i.x(j), S_i.y(j));
          if (muMap.count(xy) == 0)
          {
            muMap[xy] = S_i.s(j, p);
          }
          else
          {
            muMap[xy] += S_i.s(j, p);
          }
        }
      }
      
      for (const auto& kv : muMap)
      {
        _R.setCopyNumberState(p, i,
                              kv.first.first,
                              kv.first.second,
                              kv.second);
      }
    }
  }
}

void GenerativeModel::writeClustering(std::ostream& out) const
{
  for (int c = 0; c < _n; ++c)
  {
    out << c << "\t" << _R.indexToCharacter(c) << "\t" << _z[c] << std::endl;
  }
}

void GenerativeModel::writeFandPi(std::ostream& out) const
{
  out << _k << " #clusters" << std::endl;
  out << _m << " #samples" << std::endl;
  
  for (int j = 0; j < _k; ++j)
  {
    for (int p = 0; p < _m; ++p)
    {
      if (p != 0)
      {
        out << " ";
      }
      out << _f[p][j];
    }
    out << std::endl;
  }
  
  DoubleVector pi(_k, 0);
  for (int i = 0; i < _n; ++i)
  {
    ++pi[_z[i]];
  }
  for (int j = 0; j < _k; ++j)
  {
    pi[j] /= _n;
  }
  
  for (int j = 0; j < _k; ++j)
  {
    out << pi[j] << std::endl;
  }
}
