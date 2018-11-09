/*
 * coordinateascent.cpp
 *
 *  Created on: 25-oct-2017
 *      Author: M. El-Kebir
 */

#include "coordinateascent.h"

CoordinateAscent::CoordinateAscent(const ReadMatrix& R,
                                   int k,
                                   int nrGridPoints)
  : _R(R)
  , _k(k)
  , _logBinomCoeff(R.getNrSamples(), DoubleVector(R.getNrCharacters(), 0))
  , _gridF()
  , _stateTrees(R.getNrCharacters())
  , _f(R.getNrSamples(), DoubleVector(k, 0))
  , _z(R.getNrCharacters(), -1)
  , _obj(0)
{
  const double delta = 1. / (nrGridPoints - 1);
  for (int i = 0; i < nrGridPoints; ++i)
  {
    _gridF.push_back(delta * i);
  }
  
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  // for each sample, dcf and mutation enumerate all state trees
  for (int p = 0; p < m; ++p)
  {
    for (double f : _gridF)
    {
      for (int i = 0; i < n; ++i)
      {
        enumerateStateTrees(_R.getCopyNumberStates(p, i),
                            f,
                            _R.getMaxCopies(i),
                            _stateTrees[i][std::make_pair(p, f)]);
      }
    }
  }
  
  // init binomial coefficient
  DoubleVector logFactorial(2, 0);
  const int maxDepth = _R.getMaxReadDepth();
  for (int i = 1; i <= maxDepth; ++i)
  {
    logFactorial.push_back(logFactorial.back() + log(i));
  }
  
  for (int p = 0; p < m; ++p)
  {
    for (int i = 0; i < n; ++i)
    {
      int var_pi = _R.getVar(p, i);
      int ref_pi = _R.getRef(p, i);
      int tot_pi = var_pi + ref_pi;
      _logBinomCoeff[p][i] = logFactorial[tot_pi] - (logFactorial[var_pi] + logFactorial[ref_pi]);
    }
  }
}

void CoordinateAscent::solve(IntVector z0)
{
  const int n = _R.getNrCharacters();
  assert(z0.size() == n);
  _z = z0;
  
  _obj = 0;
  for (int iteration = 0; iteration < 100; ++iteration)
  {
    double prevObj = _obj;
    solveF();
    std::cout << "Iteration: " << iteration << " F-step; objective: " << _obj << std::endl;
    solveZ();
    std::cout << "Iteration: " << iteration << " Z-step; objective: " << _obj << std::endl;
    if (prevObj == _obj)
      break;
  }
  
  const int m = _R.getNrSamples();
  
  std::cout << "z =";
  for (int i = 0; i < n; ++i)
    std::cout << " " << _z[i];
  std::cout << std::endl << std::endl;
  
  std::cout << "F =" << std::endl;
  for (int p = 0; p < m; ++p)
  {
    for (int j = 0; j < _k; ++j)
    {
      if (j != 0)
      {
        std::cout << " ";
      }
      std::cout << _f[p][j];
    }
    std::cout << std::endl;
  }
}

void CoordinateAscent::solve(int nrIterations, int seed)
{
  // Random number generator
  std::mt19937 rng(seed);

  std::uniform_int_distribution<> unif_0k(0, _k - 1);

  // 1. sample z
  for (int c = 0; c < _R.getNrCharacters(); ++c)
  {
    _z[c] = unif_0k(rng);
  }
  
  _obj = 0;
  for (int iteration = 0; iteration < nrIterations; ++iteration)
  {
    double prevObj = _obj;
    solveF();
    std::cout << "Iteration: " << iteration << " F-step; objective: " << _obj << std::endl;
    solveZ();
    std::cout << "Iteration: " << iteration << " Z-step; objective: " << _obj << std::endl;
    if (prevObj == _obj)
      break;
  }
  
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();

  std::cout << "z =";
  for (int i = 0; i < n; ++i)
    std::cout << " " << _z[i];
  std::cout << std::endl << std::endl;
  
  std::cout << "F =" << std::endl;
  for (int p = 0; p < m; ++p)
  {
    for (int j = 0; j < _k; ++j)
    {
      if (j != 0)
      {
        std::cout << " ";
      }
      std::cout << _f[p][j];
    }
    std::cout << std::endl;
  }
}

void CoordinateAscent::enumerateStateTrees(const ReadMatrix::CopyNumberStateVector& cnStates,
                                           const double dcf,
                                           const int maxCopyNumber,
                                           FrequencyMap& mapToFreqs) const
{
  // 0. obtain L
  IntPairSet L;
  for (const ReadMatrix::CopyNumberState& cnState : cnStates)
  {
    L.insert(IntPair(cnState._x, cnState._y));
  }
  
  // 1. enumerate all state trees
  const StateGraph::StateEdgeSetSet& setS = StateGraph::getStateTrees(L, maxCopyNumber);
  
  // 2. now check for each S in setS whether S is feasible
  for (const StateGraph::StateEdgeSet& S : setS)
  {
    DoubleTensor f(maxCopyNumber + 1,
                   DoubleMatrix(maxCopyNumber + 1,
                                DoubleVector(maxCopyNumber + 1, 0)));
    
    // 2a. get vertices of S
    StateGraph::CnaTripleSet verticesS;
    verticesS.insert(StateGraph::CnaTriple(1,1,0));
    
    int mut_x = -1, mut_y = -1;
    StateGraph::StateEdge mutEdge;
    for (const StateGraph::StateEdge& st : S)
    {
      verticesS.insert(st.first);
      verticesS.insert(st.second);
      
      // check whether st is a mutation edge
      if (st.first._x == st.second._x && st.first._y == st.second._y)
      {
        mut_x = st.first._x;
        mut_y = st.first._y;
        mutEdge = st;
      }
    }
    
    if (mut_x == -1 && mut_y == -1) continue;
    
    StateGraph::CnaTripleSet verticesPreMutS, verticesPostMutS;
    StateGraph::CnaTriple mutationVertex;
    StateGraph::partition(S, verticesPreMutS, verticesPostMutS, mutationVertex);
    
    double massPreMut = 0, massPostMut = 0, massMut = 0;
    
    for (const ReadMatrix::CopyNumberState& cnState : cnStates)
    {
      if (cnState._x == mut_x && cnState._y == mut_y)
      {
        massMut = cnState._mu;
      }
    }
    
    for (const StateGraph::CnaTriple& xyz : verticesS)
    {
      if (xyz != mutEdge.first && xyz != mutEdge.second)
      {
        for (const ReadMatrix::CopyNumberState& cnState : cnStates)
        {
          if (cnState._x == xyz._x && cnState._y == xyz._y)
          {
            if (verticesPostMutS.count(xyz) == 1)
            {
              massPostMut += cnState._mu;
              f[xyz._x][xyz._y][xyz._z] = massPostMut;
            }
            else
            {
              assert(verticesPreMutS.count(xyz) == 1);
              massPreMut += cnState._mu;
              f[xyz._x][xyz._y][xyz._z] = massPreMut;
            }
          }
        }
      }
    }
    
    if (massPostMut > dcf)
    {
      continue;
    }
    
    double slack = dcf - massPostMut;
    if (slack > massMut)
    {
      continue;
    }
    
    f[mutEdge.second._x][mutEdge.second._y][mutEdge.second._z] = slack;
    f[mutEdge.first._x][mutEdge.first._y][mutEdge.first._z] = massMut - slack;
    
    mapToFreqs[S] = f;
  }
}

bool CoordinateAscent::next(const DoubleMatrix& allowedDCFs, IntVector& indices) const
{
  const int m = _R.getNrSamples();
  for (int p = 0; p < m; ++p)
  {
    if (indices[p] < allowedDCFs[p].size() - 1)
    {
      ++indices[p];
      
      for (int q = 0; q < p; ++q)
      {
        indices[q] = 0;
      }
      
      return true;
    }
  }
  return false;
}

void CoordinateAscent::identifyDCFs(const IntVector& Z_j,
                                    DoubleMatrix& dcf) const
{
  const int m = _R.getNrSamples();

  dcf = DoubleMatrix(m);
  for (int p = 0; p < m; ++p)
  {
    for (double f : _gridF)
    {
      bool f_ok = true;
      for (int i : Z_j)
      {
        f_ok &= !_stateTrees[i].find(IntDoublePair(p, f))->second.empty();
      }
      if (f_ok)
      {
        dcf[p].push_back(f);
      }
    }
  }
}

void CoordinateAscent::intersect(const DoubleVector& DCFs,
                                 const int i,
                                 StateGraph::StateEdgeSetSet& commonStateTrees) const
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  assert(0 <= i && i < n);
  
  for (const auto& kv : _stateTrees[i].find(IntDoublePair(0, DCFs[0]))->second)
  {
    commonStateTrees.insert(kv.first);
  }
  
  for (int p = 1; p < m; ++p)
  {
    const FrequencyMap& freqMap = _stateTrees[i].find(IntDoublePair(p, DCFs[p]))->second;
    for (StateGraph::StateEdgeSetSetNonConstIt it = commonStateTrees.begin(); it != commonStateTrees.end();)
    {
      if (freqMap.count(*it) == 0)
      {
        it = commonStateTrees.erase(it);
      }
      else
      {
        ++it;
      }
    }
  }
}

void CoordinateAscent::solveF()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();
  
  _obj = 1;
  for (int j = 0; j < _k; ++j)
  {
    // collection all mutations in cluster j
    IntVector Z_j;
    for (int i = 0; i < n; ++i)
    {
      if (_z[i] == j)
      {
        Z_j.push_back(i);
      }
    }
    
    // pick a set of DCFs for each sample that result in a non-empty set of state trees
    DoubleMatrix allowedDCFs;
    identifyDCFs(Z_j, allowedDCFs);
    
    double maxProb = 0;
    DoubleVector maxDCF;
    
    IntVector f_j_indices(m, 0);
    do
    {
      DoubleVector DCF(m, 0);
      for (int p = 0; p < m; ++p)
      {
        DCF[p] = allowedDCFs[p][f_j_indices[p]];
      }
      
      StateGraph::StateEdgeSetSet TT;
      double prob = 1;
      for (int i : Z_j)
      {
        StateGraph::StateEdgeSetSet commonStateTrees;
        intersect(DCF, i, commonStateTrees);
        
        double prob_i = 0;
        for (const StateGraph::StateEdgeSet& T_i : commonStateTrees)
        {
          double log_prob_i = 0;
          for (int p = 0; p < m; ++p)
          {
            const int var_pi = _R.getVar(p, i);
            const int ref_pi = _R.getRef(p, i);
            
            const DoubleTensor& S_pi = _stateTrees[i].find(IntDoublePair(p, DCF[p]))->second.find(T_i)->second;
            const double h_pi = getVAF(S_pi);
                       
            log_prob_i += _logBinomCoeff[p][i] + var_pi * log(h_pi) + ref_pi * log(1 - h_pi);
          }
          prob_i += exp(log_prob_i);
        }
        prob *= prob_i;
      }
      
      if (prob > maxProb)
      {
        maxProb = prob;
        maxDCF = DCF;
      }
    } while (next(allowedDCFs, f_j_indices));
    
    for (int p = 0; p < m; ++p)
    {
      _f[p][j] = maxDCF[p];
    }
    
    _obj *= maxProb;
  }
}

void CoordinateAscent::solveZ()
{
  const int n = _R.getNrCharacters();
  const int m = _R.getNrSamples();

  double obj = 1;
  for (int i = 0; i < n; ++i)
  {
    double max_prob_i = 0;
    for (int j = 0; j < _k; ++j)
    {
      DoubleVector DCF(m, 0);
      for (int p = 0; p < m; ++p)
      {
        DCF[p] = _f[p][j];
      }
      
      StateGraph::StateEdgeSetSet commonStateTrees;
      intersect(DCF, i, commonStateTrees);
      
      double prob_i = 0;
      for (const StateGraph::StateEdgeSet& T_i : commonStateTrees)
      {
        double log_prob_i = 0;
        for (int p = 0; p < m; ++p)
        {
          const int var_pi = _R.getVar(p, i);
          const int ref_pi = _R.getRef(p, i);
          
          const DoubleTensor& S_pi = _stateTrees[i].find(IntDoublePair(p, DCF[p]))->second.find(T_i)->second;
          const double h_pi = getVAF(S_pi);
          
          log_prob_i += _logBinomCoeff[p][i] + var_pi * log(h_pi) + ref_pi * log(1 - h_pi);
        }
        prob_i += exp(log_prob_i);
      }
      
      if (prob_i > max_prob_i)
      {
        max_prob_i = prob_i;
        _z[i] = j;
      }
    }
    obj *= max_prob_i;
  }
  _obj = obj;
}

