/*
 * stategraph.cpp
 *
 *  Created on: 11-jan-2016
 *      Author: M. El-Kebir
 */

#include <lemon/bfs.h>
#include "stategraph.h"

StateGraph::Dictionary StateGraph::_dict = StateGraph::Dictionary();

const StateGraph::StateEdgeSetSet& StateGraph::getStateTrees(const IntPairSet& L,
                                                             int maxCopyNumber)
{
  StateGraph G(maxCopyNumber);
  
  if (_dict.find(L) == _dict.end())
  {
    G.enumerate(L);
    _dict[L] = G._result;
  }
  
  return _dict[L];
}
  
StateGraph::StateGraph(int maxCopyNumber)
  : _maxCopyNumber(maxCopyNumber)
  , _G()
  , _x(_G, 0)
  , _y(_G, 0)
  , _xbar(_G, 0)
  , _ybar(_G, 0)
  , _type(_G)
  , _toNode(_maxCopyNumber + 1, Node3Matrix(_maxCopyNumber + 1,
                                            NodeMatrix(_maxCopyNumber + 1,
                                                       NodeVector(_maxCopyNumber + 1,
                                                                  lemon::INVALID))))
  , _root(lemon::INVALID)
  , _result()
{
  init();
}
  
void StateGraph::init()
{
  // vertices
  for (int x = 0; x <= _maxCopyNumber; ++x)
  {
    for (int y = 0; y <= x; ++y)
    {
      for (int z = 0; z <= x; ++z)
      {
        Node v_xyz0 = _G.addNode();
        _toNode[x][y][z][0] = v_xyz0;
        _x[v_xyz0] = x;
        _y[v_xyz0] = y;
        _xbar[v_xyz0] = z;
        _ybar[v_xyz0] = 0;
        
        if (0 < z && z <= y)
        {
          Node v_xy0z = _G.addNode();
          _toNode[x][y][0][z] = v_xy0z;
          _x[v_xy0z] = x;
          _y[v_xy0z] = y;
          _xbar[v_xy0z] = 0;
          _ybar[v_xy0z] = z;
        }
      }
    }
  }
  
  _root = _toNode[1][1][0][0];
  
  // edges
  for (NodeIt v(_G); v != lemon::INVALID; ++v)
  {
    const int x = _x[v];
    const int y = _y[v];
    const int xbar = _xbar[v];
    const int ybar = _ybar[v];
    
    // mutation edges
    if (xbar == 0 && ybar == 0 && x > 0)
    {
      addEdge(v, MUTATION, x, y, xbar + 1, ybar);
    }
    if (xbar == 0 && ybar == 0 && y > 0)
    {
      addEdge(v, MUTATION, x, y, xbar, ybar + 1);
    }
    // amplification edges
    if (0 < x && x < _maxCopyNumber && xbar < x)
    {
      addEdge(v, AMPLIFICATION, x+1, y, xbar, ybar);
    }
    if (0 < x && x < _maxCopyNumber && xbar > 0)
    {
      addEdge(v, AMPLIFICATION, x+1, y, xbar+1, ybar);
    }
    if (0 < y && y < x && ybar < y)
    {
      addEdge(v, AMPLIFICATION, x, y+1, xbar, ybar);
    }
    if (0 < y && y < x && ybar > 0)
    {
      addEdge(v, AMPLIFICATION, x, y+1, xbar, ybar+1);
    }
    if (x == y && xbar != ybar && 0 < y && ybar < y && y < _maxCopyNumber)
    {
      addEdge(v, AMPLIFICATION, y+1, x, ybar, xbar);
    }
    if (x == y && xbar != ybar && 0 < y && ybar > 0 && y < _maxCopyNumber)
    {
      addEdge(v, AMPLIFICATION, y+1, x, ybar+1, xbar);
    }
    // deletion edges
    if (x > xbar && x > y)
    {
      addEdge(v, DELETION, x-1, y, xbar, ybar);
    }
    if (xbar > 0 && x > y)
    {
      addEdge(v, DELETION, x-1, y, xbar-1, ybar);
    }
    if (y > ybar)
    {
      addEdge(v, DELETION, x, y-1, xbar, ybar);
    }
    if (ybar > 0)
    {
      addEdge(v, DELETION, x, y-1, xbar, ybar-1);
    }
    if (x == y && xbar != ybar && x > xbar)
    {
      addEdge(v, DELETION, y, x-1, ybar, xbar);
    }
    if (x == y && xbar != ybar && xbar > 0)
    {
      addEdge(v, DELETION, y, x-1, ybar, xbar-1);
    }
  }
}
  
void StateGraph::writeDOT(std::ostream& out) const
{
  IntNodeMap level(_G, 0);
  lemon::bfs(_G).distMap(level).run(_toNode[1][1][0][0]);
  
  out << "digraph G {" << std::endl;
  
  for (NodeIt v(_G); v != lemon::INVALID; ++v)
  {
    out << "\t" << _G.id(v) << " [label=\"("
        << _x[v] << "," << _y[v] << ","
        << _xbar[v] << "," << _ybar[v] << ")\\n" << level[v] << "\"]" << std::endl;
  }
  
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    Node v = _G.source(a);
    Node w = _G.target(a);
    
    out << "\t" << _G.id(v) << " -> " << _G.id(w) << " [color=";
    switch (_type[a])
    {
      case MUTATION:
        out << "black";
        break;
      case AMPLIFICATION:
        out << "green";
        break;
      case DELETION:
        out << "red";
        break;
    }
    out << "]" << std::endl;
  }
  out << "}" << std::endl;
}
  
void StateGraph::init(SubDigraph& subG,
                      SubDigraph& T,
                      ArcList& F)
{
  T.enable(_root);
  F.clear();
  
  for (OutArcIt a(_G, _root); a != lemon::INVALID; ++a)
  {
    F.push_back(a);
  }
}
  
int StateGraph::enumerate(const IntPairSet& L)
{
  BoolNodeMap filterNodesT(_G, false);
  BoolArcMap filterArcsT(_G, false);
  SubDigraph T(_G, filterNodesT, filterArcsT);
  
  BoolNodeMap filterNodesG(_G, true);
  BoolArcMap filterArcsG(_G, true);
  SubDigraph subG(_G, filterNodesG, filterArcsG);
  
  ArcList F;
  init(subG, T, F);
  
  Arc a = lemon::INVALID;
  IntMatrix V(_maxCopyNumber + 1,
              IntVector(_maxCopyNumber + 1, 0));
  V[1][1] = 1;
  
  BoolMatrix LL(_maxCopyNumber + 1,
                BoolVector(_maxCopyNumber + 1, false));
  for (const IntPair& xy : L)
  {
    LL[xy.first][xy.second] = true;
  }
  
  _result.clear();
  grow(L, LL, subG, T, F, V, a);
  
  return _result.size();
}
  
bool StateGraph::isValid(const IntPairSet& L,
                         const BoolMatrix& LL,
                         const SubDigraph& T,
                         const IntMatrix& V) const
{
  for (const IntPair& xy : L)
  {
    if (V[xy.first][xy.second] == 0)
      return false;
  }
  
  // check #2: all leaves are in L
  for (SubNodeIt v(T); v != lemon::INVALID; ++v)
  {
    bool leaf = SubOutArcIt(T, v) == lemon::INVALID;
    if (leaf && !LL[_x[v]][_y[v]])
    {
      return false;
    }
  }
    
  return true;
}
  
bool StateGraph::isValid(const SubDigraph& T) const
{
  if (!T.status(_root))
    return false;
  
  // at most one mutation edge
  int mutCount = 0;
  Arc mutationEdge = lemon::INVALID;
  for (SubArcIt a_cidj(T); a_cidj != lemon::INVALID; ++a_cidj)
  {
    if (_type[a_cidj] == MUTATION)
    {
      ++mutCount;
      mutationEdge = a_cidj;
    }
  }
  
  if (mutCount > 1)
    return false;
  
  int mut_x = _x[T.source(mutationEdge)];
  int mut_y = _y[T.source(mutationEdge)];
  
  // inf sites on copy states
  // this boils down to checking that |V_(x,y)| = 1
  IntMatrix V(_maxCopyNumber + 1, IntVector(_maxCopyNumber + 1, 0));
  for (SubNodeIt v(T); v != lemon::INVALID; ++v)
  {
    int x = _x[v];
    int y = _y[v];
    
    if (x == mut_x && y == mut_y)
      continue;
    
    if (V[x][y] == 0)
    {
      V[x][y] = 1;
    }
    else
    {
      return false;
    }
  }
  
  return true;
}
  
bool StateGraph::isValid(const SubDigraph& T,
                         Arc a)
{
  Node w = T.target(a);
  
  assert(T.status(T.source(a)));
  assert(!T.status(a));
  assert(!T.status(w));
  
  T.enable(a);
  T.enable(w);
  
  bool res = isValid(T);

  T.disable(a);
  T.disable(w);
  
  return res;
}
  
void StateGraph::writeDOT(const SubDigraph& T, std::ostream& out) const
{
  out << "digraph S {" << std::endl;

  for (SubNodeIt v(T); v != lemon::INVALID; ++v)
  {
    out << "\t" << T.id(v) << " [label=\"("
        << _x[v] << "," << _y[v] << ","
        << _xbar[v] << "," << _ybar[v] << ")\"]" << std::endl;
  }
  
  for (SubArcIt a(T); a != lemon::INVALID; ++a)
  {
    out << "\t" << T.id(T.source(a)) << " -> " << T.id(T.target(a)) << std::endl;
  }
  
  out << "}" << std::endl;
}
  
void StateGraph::writeDOT(const SubDigraph& T, const ArcList& F, std::ostream& out) const
{
  out << "digraph S {" << std::endl;
  
  for (SubNodeIt v(T); v != lemon::INVALID; ++v)
  {
    out << "\t" << T.id(v) << " [label=\"("
        << _x[v] << "," << _y[v] << ","
        << _xbar[v] << "," << _ybar[v] << ")\"]" << std::endl;
  }
  
  for (ArcListIt it = F.begin(); it != F.end(); ++it)
  {
    Arc st = *it;
    Node t = _G.target(st);
    
    out << "\t" << T.id(t) << " [label=\"("
        << _x[t] << "," << _y[t] << ","
        << _xbar[t] << "," << _ybar[t] << ")\",color=red]" << std::endl;
  }
  
  for (SubArcIt a(T); a != lemon::INVALID; ++a)
  {
    out << "\t" << T.id(T.source(a)) << " -> " << T.id(T.target(a)) << std::endl;
  }
  
  for (ArcListIt it = F.begin(); it != F.end(); ++it)
  {
    Arc st = *it;
    out << "\t" << _G.id(_G.source(st)) << " -> " << _G.id(_G.target(st)) << " [color=red]" << std::endl;
  }
  
  out << "}" << std::endl;
}
  
bool StateGraph::finalize(const IntPairSet& L,
                          const BoolMatrix& LL,
                          const SubDigraph& T,
                          const IntMatrix& V,
                          const Arc& mutationEdge)
{
//  static int idx = 0;
//  ++idx;
  
  // convert to triples
  CnaTripleNodeMap toTriple(_G, CnaTriple(0, 0, 0));
  
  for (SubNodeIt v(T); v != lemon::INVALID; ++v)
  {
    int x_v = _x[v];
    int y_v = _y[v];
    int xbar_v = _xbar[v];
    int ybar_v = _ybar[v];
    
    toTriple[v]._x = x_v;
    toTriple[v]._y = y_v;
    toTriple[v]._z = std::max(xbar_v, ybar_v);
  }
  
  // mutation edge must be in LL
  if (mutationEdge == lemon::INVALID ||
      !LL[_x[T.source(mutationEdge)]][_y[T.source(mutationEdge)]])
  {
    return false;
  }
  
  StateEdgeSet S;
  
  // shortcircuit inner nodes not in LL
  for (SubNodeIt v(T); v != lemon::INVALID; ++v)
  {
    int x_v = _x[v];
    int y_v = _y[v];
    
    if (LL[x_v][y_v])
    {
      // find parent
      SubInArcIt a(T, v);
      while (a != lemon::INVALID)
      {
        Node u = T.source(a);
        
        int x_u = _x[u];
        int y_u = _y[u];
        
        if (LL[x_u][y_u])
        {
          S.insert(std::make_pair(toTriple[u], toTriple[v]));
          break;
        }
        else
        {
          a = SubInArcIt(T, u);
        }
      }
    }
  }
  
  bool ok = false;
  for (const auto st : S)
  {
    ok |= st.second._z > 0;
  }
  
//  if (idx == 19)
//  {
//    writeDOT(T, std::cout);
//  }
  
  if (ok)
  {
    _result.insert(S);
    return false;
  }
  
  return false;
}
  
void StateGraph::grow(const IntPairSet& L,
                      const BoolMatrix& LL,
                      SubDigraph& G,
                      SubDigraph& T,
                      ArcList& F,
                      IntMatrix& V,
                      Arc& mutationEdge)
{
  if (isValid(L, LL, T, V) && mutationEdge != lemon::INVALID)
  {
    // report
    finalize(L, LL, T, V, mutationEdge);
//    writeDOT(T, std::cout);
  }
  else
  {    
    ArcList FF;
    while (!F.empty())
    {
      assert(!F.empty());
      
      Arc uv = F.back();
      F.pop_back();
      
      Node u = G.source(uv);
      Node v = G.target(uv);
      
      int x_v = _x[v];
      int y_v = _y[v];
      
      // mutation edge => V[x_v][y_v] == 0
      assert(mutationEdge == lemon::INVALID || V[x_v][y_v] == 0);
      assert(V[x_v][y_v] <= 1);

      // update V
      ++V[x_v][y_v];
      if (_type[uv] == MUTATION)
      {
        mutationEdge = uv;
      }
      
      assert(T.status(u));
      assert(!T.status(v));
      assert(!T.status(uv));
      
      // add uv to T
      T.enable(v);
      T.enable(uv);

      ArcList newF = F;
      
      // remove arcs from F
      for (ArcListNonConstIt it = newF.begin(); it != newF.end();)
      {
        Arc st = *it;
        Node s = T.source(st);
        Node t = T.target(st);
        
        int x_t = _x[t];
        int y_t = _y[t];
        
        if (v == t || (v != s && x_v == x_t && y_v == y_t) || (mutationEdge == uv && V[x_t][y_t] == 1))
        {
          assert(T.status(s));
          it = newF.erase(it);
        }
        else
        {
          ++it;
        }
      }
      
      // push each arc vw where w not in V(T) onto F
      for (SubOutArcIt vw(G, v); vw != lemon::INVALID; ++vw)
      {
        Node w = G.target(vw);
        if (T.status(w))
          continue;
        
        int x_w = _x[w];
        int y_w = _y[w];
        
        if ((mutationEdge == lemon::INVALID && x_v == x_w && y_v == y_w)
            || V[x_w][y_w] == 0)
        {
          newF.push_back(vw);
        }
      }
      
//      writeDOT(T, newF, std::cout);
      grow(L, LL, G, T, newF, V, mutationEdge);
      
      G.disable(uv);
      if (uv == mutationEdge)
      {
        mutationEdge = lemon::INVALID;
      }
      
      T.disable(uv);
      T.disable(v);
      --V[x_v][y_v];
      
      FF.push_back(uv);
    }

    for (ArcListRevIt it = FF.rbegin(); it != FF.rend(); ++it)
    {
      Arc a = *it;
      assert(!G.status(a));
      
      F.push_back(*it);
      G.enable(a);
    }
  }
}

void StateGraph::partition(const StateEdgeSet& S,
                           const CnaTriple& vertex,
                           bool& preMutationFlag,
                           CnaTriple& mutation,
                           CnaTripleSet& preMutation,
                           CnaTripleSet& postMutation)
{
  if (preMutationFlag)
  {
    preMutation.insert(vertex);
  }
  else
  {
    postMutation.insert(vertex);
  }
  
  for (const StateEdge& st : S)
  {
    if (st.first == vertex)
    {
      if (preMutationFlag && st.second._z == 1)
      {
        bool newPreMutationFlag = false;
        partition(S, st.second, newPreMutationFlag, mutation, preMutation, postMutation);
        mutation = st.second;
      }
      else
      {
        partition(S, st.second, preMutationFlag, mutation, preMutation, postMutation);
      }
    }
  }
}

void StateGraph::partition(const StateEdgeSet& S,
                           CnaTripleSet& preMutation,
                           CnaTripleSet& postMutation,
                           CnaTriple& mutation)
{
  CnaTriple root(1, 1, 0);
  bool preMutationFlag = true;
  mutation = CnaTriple(-1, -1, -1);
  partition(S, root, preMutationFlag, mutation, preMutation, postMutation);
}
  
bool operator<(const StateGraph::CnaTriple& lhs, const StateGraph::CnaTriple& rhs)
{
  return lhs._x < rhs._x || (lhs._x == rhs._x && lhs._y < rhs._y) || (lhs._x == rhs._x && lhs._y == rhs._y && lhs._z < rhs._z);
}
