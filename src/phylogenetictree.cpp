/*
 * phylogenetictree.cpp
 *
 *  Created on: 20-jan-2018
 *      Author: M. El-Kebir
 */

#include "phylogenetictree.h"
#include <lemon/bfs.h>
#include <lemon/lgf_reader.h>
#include <lemon/lgf_writer.h>
#include <iomanip>

PhylogeneticTree::PhylogeneticTree()
  : _T()
  , _rootT(lemon::INVALID)
  , _nodeToClusterT(_T)
  , _clusterToNodeT()
  , _nodeToCharStateT(_T)
  , _charStateToNodeT()
  , _isMutationNodeT(_T)
  , _clusterCoverT(_T)
  , _U(_T)
{
}

PhylogeneticTree::PhylogeneticTree(const StateTreeVector& S,
                                   const IntVector& z,
                                   const Digraph& mutT,
                                   const Node rootMutT,
                                   const IntNodeMap& nodeToClusterMutT,
                                   const NodeVector& clusterToNodeMutT)
  : _T()
  , _rootT(_T.addNode())
  , _nodeToClusterT(_T, -1)
  , _clusterToNodeT()
  , _nodeToCharStateT(_T)
  , _charStateToNodeT()
  , _isMutationNodeT(_T, false)
  , _clusterCoverT(_T, -1)
  , _U(_T)
{
  const int nrClusters = *(std::max_element(z.begin(), z.end())) + 1;
  _clusterToNodeT = NodeVector(nrClusters, lemon::INVALID);
  
  init(S, z,
       mutT,
       nodeToClusterMutT,
       clusterToNodeMutT,
       rootMutT);
  updateClusterCover(_rootT, -1);
}

void PhylogeneticTree::updateClusterCover(Node v, int cluster)
{
  _clusterCoverT[v] = cluster;
  
  for (OutArcIt a(_T, v); a != lemon::INVALID; ++a)
  {
    Node w = _T.target(a);
    if (_isMutationNodeT[w])
    {
      assert(0 <= _nodeToClusterT[w] && _nodeToClusterT[w] < _clusterToNodeT.size());
      updateClusterCover(w, _nodeToClusterT[w]);
    }
    else
    {
      updateClusterCover(w, cluster);
    }
  }
}

void PhylogeneticTree::sample(const DoubleVector& purityVector,
                              const DoubleMatrix& clusterUsages,
                              StateTreeVector& S)
{
  assert(!purityVector.empty());
  const int nrSamples = purityVector.size();
  const int nrClusters = _clusterToNodeT.size();
  
  std::gamma_distribution<> gamma(1, 1);
  
  // 0. initialize usages of all vertices
  for (NodeIt v(_T); v != lemon::INVALID; ++v)
  {
    _U[v] = DoubleVector(nrSamples, 0);
  }
  
  // 1. set purity of root vertex
  for (int p = 0; p < nrSamples; ++p)
  {
    _U[_rootT][p] = 1. - purityVector[p];
  }
  
  // 2. sample remaining vertices
  for (int clusterIdx = 0; clusterIdx < nrClusters; ++clusterIdx)
  {
    DoubleVector sum(nrSamples, 0.);
    for (NodeIt v(_T); v != lemon::INVALID; ++v)
    {
      if (_clusterCoverT[v] == clusterIdx)
      {
        for (int p = 0; p < nrSamples; ++p)
        {
          _U[v][p] = gamma(g_rng);
          sum[p] += _U[v][p];
        }
      }
    }
    
    for (NodeIt v(_T); v != lemon::INVALID; ++v)
    {
      if (_clusterCoverT[v] == clusterIdx)
      {
        for (int p = 0; p < nrSamples; ++p)
        {
          _U[v][p] /= sum[p];
          _U[v][p] *= clusterUsages[p][clusterIdx];
        }
      }
    }
  }
  
  // 3. update state trees
  const int nrCharacters = getNrCharacters();
  for (int c = 0; c < nrCharacters; ++c)
  {
    S[c].resetMixtureProportions(nrSamples);
    updateStateProportions(_rootT,
                           c, S[c].rootState(),
                           S[c]);
  }
  
  // TODO: 4. prune unused state tree vertices!
}

void PhylogeneticTree::updateStateProportions(Node v,
                                              int c, int i,
                                              StateTree& S_c)
{
  const int nrSamples = _U[v].size();
  for (int p = 0; p < nrSamples; ++p)
  {
    S_c.incrementMixtureProportion(p, i, _U[v][p]);
  }
  
  for (OutArcIt a(_T, v); a != lemon::INVALID; ++a)
  {
    Node w = _T.target(a);
    int j = i;
    for (IntPair dj : _nodeToCharStateT[w])
    {
      if (dj.first == c)
      {
        j = dj.second;
        break;
      }
    }
    updateStateProportions(w, c, j, S_c);
  }
}

Node PhylogeneticTree::pickDescendant(const Digraph& T,
                                      Node v) const
{
  BoolNodeMap reachedNodeMap(_T, false);
  lemon::bfs(_T).reachedMap(reachedNodeMap).run(v);
  
  NodeVector visitedNodes;
  for (NodeIt u(_T); u != lemon::INVALID; ++u)
  {
    if (reachedNodeMap[u])
    {
      visitedNodes.push_back(u);
    }
  }
  
  std::uniform_int_distribution<> unif(0, visitedNodes.size() - 1);
  int i = unif(g_rng);
  
  return visitedNodes[i];
}

void PhylogeneticTree::init(const StateTreeVector& S,
                            const IntVector& z,
                            const Digraph& mutT,
                            const IntNodeMap& nodeToClusterMutT,
                            const NodeVector& clusterToNodeMutT,
                            const Node nodeMutT)
{
  // 1. Collect mutations
  int clusterIdx = nodeToClusterMutT[nodeMutT];
  IntVector characters;
  
  const int n = z.size();
  for (int i = 0; i < n; ++i)
  {
    if (z[i] == clusterIdx)
    {
      characters.push_back(i);
    }
  }
  
  // 2. Identify insertion node (in _T) for mutations in clusterIdx
  Node insertionNode;
  if (InArcIt(mutT, nodeMutT) == lemon::INVALID)
  {
    insertionNode = _rootT;
  }
  else
  {
    int parentClusterIdx = nodeToClusterMutT[mutT.source(InArcIt(mutT, nodeMutT))];
    Node ancestorInsertionNode = _clusterToNodeT[parentClusterIdx];
    assert(ancestorInsertionNode != lemon::INVALID);
    
    insertionNode = pickDescendant(_T, ancestorInsertionNode);
  }
  
  // 3. Extend _T
  Node mutationNode = extend(S, clusterIdx,
                             characters, insertionNode);
  _clusterToNodeT[clusterIdx] = mutationNode;
  _nodeToClusterT[mutationNode] = clusterIdx;
  
  // 4. Recurse on children of the mutation tree
  for (OutArcIt a(mutT, nodeMutT); a != lemon::INVALID; ++a)
  {
    Node childNodeMutT = mutT.target(a);
    init(S,
         z,
         mutT,
         nodeToClusterMutT,
         clusterToNodeMutT,
         childNodeMutT);
  }
}

Node PhylogeneticTree::pickParent(const StateTree& S_c,
                                  const IntPair& ci) const
{
  const int parentState = S_c.parent(ci.second);
  
  assert(_charStateToNodeT.count(IntPair(ci.first, parentState)) == 1);
  Node ancestorNodeT = _charStateToNodeT.find(IntPair(ci.first, parentState))->second;
  
  // 1. Identify black list
  NodeSet blackList;
  for (OutArcIt a(S_c.S(), S_c.node(parentState));
       a != lemon::INVALID; ++a)
  {
    int j = S_c.state(S_c.S().target(a));
    IntPair cj(ci.first, j);
    if (j != ci.second && _charStateToNodeT.count(cj) == 1)
    {
      Node blackListNode = _charStateToNodeT.find(cj)->second;
      BoolNodeMap visited(_T, false);
      lemon::bfs(_T).reachedMap(visited).run(blackListNode);
      for (NodeIt v(_T); v != lemon::INVALID; ++v)
      {
        if (visited[v])
        {
          blackList.insert(v);
        }
      }
    }
  }
  
  // 2. Identify white list
  NodeVector whiteList;
  
  BoolNodeMap reachedNodeMap(_T, false);
  lemon::bfs(_T).reachedMap(reachedNodeMap).run(ancestorNodeT);
  
  for (NodeIt u(_T); u != lemon::INVALID; ++u)
  {
    if (reachedNodeMap[u] && blackList.count(u) == 0)
    {
      whiteList.push_back(u);
    }
  }
  
  assert(!whiteList.empty());
  
  // 3. Pick white list node
  std::uniform_int_distribution<> unif(0, whiteList.size() - 1);
  int i = unif(g_rng);
  
  return whiteList[i];
}

Node PhylogeneticTree::extend(const StateTreeVector& S,
                              const int clusterIdx,
                              const IntVector& characters,
                              Node insertionNode)
{
  const int n = S.size();
  
  IntPairVector frontier;
  
  // 1. Determine mutation states for each state tree
  IntVector mutationStateVector(n, -1);
  for (const int c : characters)
  {
    assert(0 <= c && c < n);
    const StateTree& S_c = S[c];
    
    mutationStateVector[c] = S_c.mutationState();
  }
  
  // 2. Initialize pre-mutation path and frontier
  {
    std::vector<IntList> ancestorList(n);
    
    for (const int c : characters)
    {
      const StateTree& S_c = S[c];
      
      int rootState_c = S_c.rootState();
      _charStateToNodeT[IntPair(c, rootState_c)] = _rootT;
      _nodeToCharStateT[_rootT].insert(IntPair(c, rootState_c));
      
      int i = mutationStateVector[c];
      while ((i = S_c.parent(i)) != -1)
      {
        // skip the root state
        if (i != rootState_c)
        {
          ancestorList[c].push_front(i);
          _charStateToNodeT[IntPair(c, i)] = lemon::INVALID;
        }
      }
      
      for (OutArcIt a(S_c.S(), S_c.node(rootState_c)); a != lemon::INVALID; ++a)
      {
        int childState_c = S_c.state(S_c.S().target(a));
        if (childState_c != mutationStateVector[c]
            && _charStateToNodeT.count(IntPair(c, childState_c)) == 0)
        {
          frontier.push_back(std::make_pair(c, childState_c));
        }
      }
    }
    
    bool done = false;
    IntVector newCharacters = characters;
    while (!done)
    {
      done = true;
      
      std::shuffle(newCharacters.begin(),
                   newCharacters.end(),
                   g_rng);
      
      for (const int c : newCharacters)
      {
        const StateTree& S_c = S[c];
        int rootState_c = S_c.rootState();
        
        if (!ancestorList[c].empty())
        {
          done = false;
          int i = ancestorList[c].front();
          ancestorList[c].pop_front();
          
          Node newInsertionNode = _T.addNode();
          _isMutationNodeT[newInsertionNode] = false;
          _nodeToClusterT[newInsertionNode] = -1;
          _T.addArc(insertionNode, newInsertionNode);
          _nodeToCharStateT[newInsertionNode].insert(IntPair(c, i));
          _charStateToNodeT[IntPair(c, i)] = newInsertionNode;
          insertionNode = newInsertionNode;
          
          for (OutArcIt a(S_c.S(), S_c.node(rootState_c)); a != lemon::INVALID; ++a)
          {
            int childState_c = S_c.state(S_c.S().target(a));
            if (childState_c != mutationStateVector[c]
                && _charStateToNodeT.count(IntPair(c, childState_c)) == 0)
            {
              frontier.push_back(std::make_pair(c, childState_c));
            }
          }
        }
      }
    }
  }
  
  // 3. Add mutation vertex
  Node mutationNode = _T.addNode();
  _isMutationNodeT[mutationNode] = true;
  _nodeToClusterT[mutationNode] = clusterIdx;
  _clusterToNodeT[clusterIdx] = mutationNode;
  _T.addArc(insertionNode, mutationNode);
  
  for (const int c : characters)
  {
    int i = mutationStateVector[c];
    _nodeToCharStateT[mutationNode].insert(IntPair(c, i));
    _charStateToNodeT[IntPair(c, i)] = mutationNode;
    
    const StateTree& S_c = S[c];
    for (OutArcIt a(S_c.S(), S_c.node(i)); a != lemon::INVALID; ++a)
    {
      int childState_c = S_c.state(S_c.S().target(a));
      assert(_charStateToNodeT.count(IntPair(c, childState_c)) == 0);
      frontier.push_back(std::make_pair(c, childState_c));
    }
  }
  
  // 4. Add remaining vertices
  while (!frontier.empty())
  {
    int idx = std::uniform_int_distribution<>(0, frontier.size() - 1)(g_rng);
    IntPair ci = frontier[idx];
    frontier.erase(frontier.begin() + idx);
    
    Node parentNode = pickParent(S[ci.first], ci);
    Node newNode = _T.addNode();
    _isMutationNodeT[newNode] = false;
    _nodeToClusterT[newNode] = -1;
    _T.addArc(parentNode, newNode);
    
    _nodeToCharStateT[newNode].insert(ci);
    _charStateToNodeT[ci] = newNode;

    const StateTree& S_c = S[ci.first];
    for (OutArcIt a(S_c.S(), S_c.node(ci.second));
         a != lemon::INVALID; ++a)
    {
      int childState_c = S_c.state(S_c.S().target(a));
      assert(_charStateToNodeT.count(IntPair(ci.first, childState_c)) == 0);
      frontier.push_back(IntPair(ci.first, childState_c));
    }
  }
  
  return mutationNode;
}

void PhylogeneticTree::writeDOT(std::ostream& out) const
{
  out << "digraph T {" << std::endl;
  
  for (NodeIt v(_T); v != lemon::INVALID; ++v)
  {
    out << "\t" << _T.id(v) << " [label=\"";
    for (double u : _U[v])
    {
      out << std::setprecision(2) << u << "\\n";
    }
    out << "\"";
    if (_isMutationNodeT[v])
    {
      out << ",color=red";
    }
    out << "]" << std::endl;
  }
  
  for (ArcIt a(_T); a != lemon::INVALID; ++a)
  {
    Node v = _T.source(a);
    Node w = _T.target(a);
    
    out << "\t" << _T.id(v) << " -> " << _T.id(w) << " [label=\"";
    
    bool first = true;
    for (const IntPair& ci : _nodeToCharStateT[w])
    {
      if (first)
      {
        first = false;
      }
      else
      {
        out << "\\n";
      }
      
      out << ci.first << ", " << ci.second;
    }
    
    out << "\"]";
    out << std::endl;
  }
  
  out << "}" << std::endl;
}

void PhylogeneticTree::writeDOT(const StateTreeVector& S,
                                std::ostream& out) const
{
  out << "digraph T {" << std::endl;
  
  for (NodeIt v(_T); v != lemon::INVALID; ++v)
  {
    out << "\t" << _T.id(v) << " [label=\"";
    for (double u : _U[v])
    {
      out << std::setprecision(2) << u << "\\n";
    }
    out << "\"";
    if (_isMutationNodeT[v])
    {
      out << ",color=red";
    }
    out << "]" << std::endl;
  }
  
  for (ArcIt a(_T); a != lemon::INVALID; ++a)
  {
    Node v = _T.source(a);
    Node w = _T.target(a);
    
    out << "\t" << _T.id(v) << " -> " << _T.id(w) << " [label=\"";
    
    bool first = true;
    for (const IntPair& ci : _nodeToCharStateT[w])
    {
      if (first)
      {
        first = false;
      }
      else
      {
        out << "\\n";
      }
      
      int x = -1, y = -1, z = -1;
      S[ci.first].getCnState(ci.second, x, y, z);
      out << ci.first << ", (" << x << "," << y << "," << z << ")";
    }
    
    out << "\"]";
    out << std::endl;
  }
  
  out << "}" << std::endl;
}

std::ostream& operator<<(std::ostream& out,
                         const PhylogeneticTree& phyloT)
{
  lemon::digraphWriter(phyloT._T, out)
    .node("rootT", phyloT._rootT)
    .nodeMap("nodeToClusterT", phyloT._nodeToClusterT)
    .nodeMap("nodeToCharStateT", phyloT._nodeToCharStateT)
    .nodeMap("isMutationNodeT", phyloT._isMutationNodeT)
    .nodeMap("clusterCoverT", phyloT._clusterCoverT)
    .nodeMap("U", phyloT._U)
    .run();
  
  return out;
}

std::istream& operator>>(std::istream& in,
                         PhylogeneticTree& phyloT)
{
  lemon::digraphReader(phyloT._T, in)
    .node("rootT", phyloT._rootT)
    .nodeMap("nodeToClusterT", phyloT._nodeToClusterT)
    .nodeMap("nodeToCharStateT", phyloT._nodeToCharStateT)
    .nodeMap("isMutationNodeT", phyloT._isMutationNodeT)
    .nodeMap("clusterCoverT", phyloT._clusterCoverT)
    .nodeMap("U", phyloT._U)
    .run();
  
  // Update _clusterToNodeT and _charStateToNodeT
  const int nrClusters = lemon::mapMaxValue(phyloT._T, phyloT._nodeToClusterT);
  phyloT._clusterToNodeT = NodeVector(nrClusters, lemon::INVALID);
  for (NodeIt v(phyloT._T); v != lemon::INVALID; ++v)
  {
    const int j = phyloT._nodeToClusterT[v];
    if (j != -1)
    {
      phyloT._clusterToNodeT[j] = v;
    }
    
    for (const IntPair& ci : phyloT._nodeToCharStateT[v])
    {
      phyloT._charStateToNodeT[ci] = v;
    }
  }
  
  return in;
}

std::ostream& operator<<(std::ostream& out, const IntPairSet& s)
{
  bool first = true;
  for (const IntPair& ci : s)
  {
    if (first)
    {
      first = false;
    }
    else
    {
      out << " ";
    }
    out << ci.first << " " << ci.second;
  }
  
  return out;
}

std::istream& operator>>(std::istream& in, IntPairSet& s)
{
  while (in.good())
  {
    int c = -1, i = -1;
    in >> c >> i;
    
    if (c == -1 || i == -1)
    {
      throw std::runtime_error("Parsing error while reading input file");
    }
    
    s.insert(IntPair(c, i));
  }
  
  return in;
}

std::ostream& operator<<(std::ostream& out, const DoubleVector& s)
{
  bool first = true;
  for (double v : s)
  {
    if (first)
    {
      first = false;
    }
    else
    {
      out << " ";
    }

    out << v;
  }
  
  return out;
}

std::istream& operator>>(std::istream& in,
                         DoubleVector& s)
{
  while (in.good())
  {
    double v = NAN;
    in >> v;
    
    if (v == NAN)
    {
      throw std::runtime_error("Parsing error while reading input file");
    }
    
    s.push_back(v);
  }
  
  return in;
}
