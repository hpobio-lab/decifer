/*
 *  python.cpp
 *
 *   Created on: 10-nov-2017
 *       Author: M. El-Kebir
 */

#include <boost/scoped_array.hpp>
#include <boost/python.hpp>
#include <iostream>
#include "utils.h"
#include "stategraph.h"
#include "statetree.h"
#include "statetreessampler.h"
#include "readmatrix.h"
#include <fstream>
#include "solver.h"
#include "posteriorstatetree.h"
#include "clusterilp.h"
#include "phylogenetictree.h"

#ifdef CPLEX
  #include "clusterilpcplex.h"
  typedef ClusterIlpCplex ClusterIlpAlg;
#else
  #include "clusterilpgurobi.h"
  typedef ClusterIlpGurobi ClusterIlpAlg;
#endif

namespace p = boost::python;

p::list computeSangerCCFs(const p::str& inputFilenameR,
                          const p::str& inputFilenamePurities)
{
  p::list res;
  
  std::string cppFilename = p::extract<std::string>(inputFilenameR);
  std::ifstream inR(cppFilename.c_str());
  if (!inR.good())
  {
    throw std::runtime_error("Error: could not open '" + cppFilename + "' for reading");
  }
  
  ReadMatrix R;
  inR >> R;
  inR.close();
  
  cppFilename = p::extract<std::string>(inputFilenamePurities);
  std::ifstream inPurities(cppFilename.c_str());
  if (!inPurities.good())
  {
    throw std::runtime_error("Error: could not open '" + cppFilename + "' for reading");
  }
  
  DoubleVector purities;
  std::string line;
  getline(inPurities, line);
  for (int p = 0; p < R.getNrSamples(); ++p)
  {
    getline(inPurities, line);
    std::stringstream ss(line);
    double puritiy = 0;
    ss >> puritiy;
    purities.push_back(puritiy);
  }
  
  const int n = R.getNrCharacters();
  const int m = R.getNrSamples();

  for (int i = 0; i < n; ++i)
  {
    p::list ccfs_i;
    for (int p = 0; p < m; ++p)
    {
      IntDoublePair nChr_CCF =
        R.computeSangerCCF(R.getRef(p, i),
                           R.getVar(p, i),
                           purities[p],
                           R.getCopyNumberStates(p, i));
      
      ccfs_i.append(nChr_CCF.second);
    }
    res.append(ccfs_i);
  }
  
  return res;
}

p::tuple computeSangerCCF(int ref,
                          int var,
                          double purity,
                          const p::list& l)
{
  // 1. transform l into CopyNumberStateVector
  ReadMatrix::CopyNumberStateVector cnStates;
  
  for (int i = 0; i < p::len(l); ++i)
  {
    ReadMatrix::CopyNumberState cnState;
    p::tuple tpl = p::extract<p::tuple>(l[i]);
    
    cnState._x = p::extract<int>(tpl[0]);
    cnState._y = p::extract<int>(tpl[1]);
    cnState._mu = p::extract<double>(tpl[2]);
    
    cnStates.push_back(cnState);
  }
  
  // 2. compute CCF and n_chr
  IntDoublePair nChr_CCF = ReadMatrix::computeSangerCCF(ref, var,
                                                        purity,
                                                        cnStates);
  // 3. return solution
  return p::make_tuple(nChr_CCF.first, nChr_CCF.second);
}

p::list enumerateStateTrees(const p::list& l)
{
  IntPairSet L;
  int max_xy = 0;
  for (int i = 0; i < p::len(l); ++i)
  {
    p::tuple tpl = p::extract<p::tuple>(l[i]);
    int x = p::extract<int>(tpl[0]);
    int y = p::extract<int>(tpl[1]);
    max_xy = std::max(max_xy, std::max(x, y));
    L.insert(IntPair(x, y));
    std::cout << "Hello (" << x << "," << y << ")" << std::endl;
  }
  
  p::list res;
  const StateGraph::StateEdgeSetSet& S = StateGraph::getStateTrees(L, max_xy);
  for (const StateGraph::StateEdgeSet& T_i : S)
  {
    StateGraph::CnaTripleSet vertices;
    p::list edges_i;
    for (const StateGraph::StateEdge edge : T_i)
    {
      auto s = p::make_tuple(edge.first._x, edge.first._y, edge.first._z);
      auto t = p::make_tuple(edge.second._x, edge.second._y, edge.second._z);
      edges_i.append(p::make_tuple(s, t));
    }
    res.append(edges_i);
  }
  return res;
}

std::string visualizeStateTree(const p::str& inputFilename,
                               const int idx,
                               const int m)
{
  std::string cppFilename = p::extract<std::string>(inputFilename);
  std::ifstream inT(cppFilename.c_str());
  if (!inT.good())
  {
    throw std::runtime_error("Error: could not open '" + cppFilename + "' for reading");
  }
  
  StateTreeVector T;
  inT >> T;

  if (!(0 <= idx && idx < T.size()))
  {
    throw std::runtime_error("Error: invalid index");
  }
  
  std::stringstream ss;
  T[idx].writeDOT(m, ss);
  return ss.str();
}

void visualizeStateTrees(const p::str& inputFilename,
                         const p::str& outputPrefix,
                         const int m)
{
  std::string cppFilename = p::extract<std::string>(inputFilename);
 
  std::ifstream inT(cppFilename.c_str());
  if (!inT.good())
  {
    throw std::runtime_error("Error: could not open '" + cppFilename + "' for reading");
  }
  
  StateTreeVector T;
  inT >> T;

  std::string cppOutputPrefix = p::extract<std::string>(outputPrefix);
  char buf[1024];
  for (int idx = 0; idx < T.size(); ++idx)
  {
    snprintf(buf, 1024, "%s.S%d.dot", cppOutputPrefix.c_str(), idx);
    std::ofstream outT(buf);
    if (!outT.good())
    {
      throw std::runtime_error("Error: could not open '" + std::string(buf) + "' for reading");
    }
    T[idx].writeDOT(m, outT);
    outT.close();
  }
}

double computeSnvLogLikelihood(const p::str& readMatrixFilename,
                               const p::str& solutionFilename,
                               int i)
{
  std::string cppReadMatrixFilename = p::extract<std::string>(readMatrixFilename);
  std::string cppSolutionFilename = p::extract<std::string>(solutionFilename);
  
  std::ifstream inR(cppReadMatrixFilename.c_str());
  if (!inR.good())
  {
    throw std::runtime_error("Error: could not open '" + cppReadMatrixFilename + "' for reading");
  }
  
  ReadMatrix R;
  inR >> R;
 
  assert(0 <= i && i < R.getNrCharacters());
  
  std::ifstream inF(cppSolutionFilename.c_str());
  if (!inF.good())
  {
    throw std::runtime_error("Error: could not open '" + cppSolutionFilename + "' for reading");
  }
  
  DoubleMatrix f;
  DoubleVector pi;
  Solver::readSolution(inF, f, pi);
  Solver solver(R, pi.size(), 0);
  solver.init(i);
  solver.setSolution(f, pi);
  
  return solver.getLogLikelihood(i);
}

double computeLogLikelihood(const p::str& readMatrixFilename,
                            const p::str& solutionFilename)
{
  std::string cppReadMatrixFilename = p::extract<std::string>(readMatrixFilename);
  std::string cppSolutionFilename = p::extract<std::string>(solutionFilename);
  
  std::ifstream inR(cppReadMatrixFilename.c_str());
  if (!inR.good())
  {
    throw std::runtime_error("Error: could not open '" + cppReadMatrixFilename + "' for reading");
  }
  
  ReadMatrix R;
  inR >> R;
  
  std::ifstream inF(cppSolutionFilename.c_str());
  if (!inF.good())
  {
    throw std::runtime_error("Error: could not open '" + cppSolutionFilename + "' for reading");
  }
  
  DoubleMatrix f;
  DoubleVector pi;
  Solver::readSolution(inF, f, pi);
  Solver solver(R, pi.size(), 0);
  solver.init();
  solver.setSolution(f, pi);
  
  return solver.getLogLikelihood();
}

std::string sampleStateTree(const p::list& l,
                            double dcf,
                            double purity)
{
  IntPairSet L;
  int max_xy = 0;
  for (int i = 0; i < p::len(l); ++i)
  {
    p::tuple tpl = p::extract<p::tuple>(l[i]);
    int x = p::extract<int>(tpl[0]);
    int y = p::extract<int>(tpl[1]);
    max_xy = std::max(max_xy, std::max(x, y));
    L.insert(IntPair(x, y));
  }
  
  StateTreesSampler sampler(max_xy, 1, L.size());
  
  StateTreeVector S = sampler.sample();
  
  std::uniform_int_distribution<> unif(0, S.size() - 1);
  int idx = unif(g_rng);
  
  S[idx].sampleMixtureProportions(dcf, purity);
  
  std::stringstream ss;
  S[idx].writeDOT(1, ss);
  return ss.str();
}

p::tuple compareStateTrees(const p::str& trueFilename,
                           const p::str& inferredFilename)
{
  std::string cppTrueFilename = p::extract<std::string>(trueFilename);
  std::string cppInferredFilename = p::extract<std::string>(inferredFilename);
  
  std::ifstream inTrueT(cppTrueFilename.c_str());
  if (!inTrueT.good())
  {
    throw std::runtime_error("Error: could not open '" + cppTrueFilename + "' for reading");
  }
  
  StateTreeVector trueT;
  inTrueT >> trueT;
  
  std::ifstream inInfT(cppInferredFilename.c_str());
  if (!inInfT.good())
  {
    throw std::runtime_error("Error: could not open '" + cppInferredFilename + "' for reading");
  }
  
  PosteriorStateTreeMatrix infT;
  inInfT >> infT;
  
  if (trueT.size() != infT.size())
  {
    throw std::runtime_error("Error: invalid sizes");
  }
  
  const int n = trueT.size();
  typedef std::tuple<int, int, int> CnTriple;
  typedef std::pair<CnTriple, CnTriple> CnTriplePair;
  typedef std::set<CnTriplePair> CnTripleSet;
  
  double equal = 0;
  p::list res;
  for (int i = 0; i < n; ++i)
  {
    CnTripleSet trueEdges, inferredEdges;
    
    // extract true edges
    const StateTree& trueT_i = trueT[i];
    for (ArcIt a(trueT_i.S()); a != lemon::INVALID; ++a)
    {
      Node u = trueT_i.S().source(a);
      Node v = trueT_i.S().target(a);
      
      int state_u = trueT_i.state(u);
      int state_v = trueT_i.state(v);
      
      CnTriplePair edge;
      trueT_i.getCnState(state_u,
                         std::get<0>(edge.first),
                         std::get<1>(edge.first),
                         std::get<2>(edge.first));
      trueT_i.getCnState(state_v,
                         std::get<0>(edge.second),
                         std::get<1>(edge.second),
                         std::get<2>(edge.second));
      trueEdges.insert(edge);
    }
    
    // identify most likely inferred state tree
    const PosteriorStateTreeVector& posterior_T_i = infT[i];
    double max_gamma = 0;
    int max_t = -1, t = 0;
    for (const PosteriorStateTree& T : posterior_T_i)
    {
      if (T._gamma > max_gamma)
      {
        max_gamma = T._gamma;
        max_t = t;
      }
      ++t;
    }
    
    // extract true edges
    const StateTree& infT_i = posterior_T_i[max_t]._T;
    for (ArcIt a(infT_i.S()); a != lemon::INVALID; ++a)
    {
      Node u = infT_i.S().source(a);
      Node v = infT_i.S().target(a);
      
      int state_u = infT_i.state(u);
      int state_v = infT_i.state(v);
      
      CnTriplePair edge;
      infT_i.getCnState(state_u,
                        std::get<0>(edge.first),
                        std::get<1>(edge.first),
                        std::get<2>(edge.first));
      infT_i.getCnState(state_v,
                        std::get<0>(edge.second),
                        std::get<1>(edge.second),
                        std::get<2>(edge.second));
      inferredEdges.insert(edge);
    }
    
    res.append(inferredEdges == trueEdges);
    if (inferredEdges == trueEdges)
    {
      ++equal;
    }
  }
  
  double recall = equal / n;
  
  return p::make_tuple(recall, res);
}

std::string visualizePosteriorStateTrees(const p::str& filename,
                                         int index)
{
  std::string cppFilename = p::extract<std::string>(filename);
  
  std::ifstream inT(cppFilename.c_str());
  if (!inT.good())
  {
    throw std::runtime_error("Error: could not open '" + cppFilename + "' for reading");
  }
  
  PosteriorStateTreeMatrix T;
  
  inT >> T;
  
  std::stringstream ss;
  PosteriorStateTree::writeDOT(T[index], ss);  
  return ss.str();
}

std::string visualize(const p::str& inputFilename,
                      int k,
                      int m,
                      double alpha,
                      int i,
                      int t)
{
  std::string cppInputFilename = p::extract<std::string>(inputFilename);
  
  std::ifstream inR(cppInputFilename.c_str());
  if (!inR.good())
  {
    throw std::runtime_error("Error: could not open '" + cppInputFilename + "' for reading");
  }
  
  ReadMatrix R;
  inR >> R;
  inR.close();
  
  ClusterIlpAlg clusterIlp(R, k, alpha);
  clusterIlp.init();
  
  const StateTree& T_it = clusterIlp.Texp(i, t);

  std::stringstream ss;
  T_it.writeDOT(m, ss);
  return ss.str();
}

void exportFeasibleSolutionSpace(const p::str& inputFilename,
                                 int k,
                                 double alpha,
                                 const p::str& outputFilename)
{
  std::string cppInputFilename = p::extract<std::string>(inputFilename);
  
  std::ifstream inR(cppInputFilename.c_str());
  if (!inR.good())
  {
    throw std::runtime_error("Error: could not open '" + cppInputFilename + "' for reading");
  }
  
  ReadMatrix R;
  inR >> R;
  inR.close();
  
  ClusterIlpAlg clusterIlp(R, k, alpha);
  clusterIlp.init();
  
  std::string cppOutputFilename = p::extract<std::string>(outputFilename);
  std::ofstream outT(cppOutputFilename.c_str());
  clusterIlp.writeFeasibleSolutionSpace(outT);
  outT.close();
}

std::string visualizePhylogeny(const std::string& phylogenyFilename,
                               const std::string& stateTreesFilename)
{
  std::ifstream inPhyloT(phylogenyFilename.c_str());
  if (!inPhyloT.good())
  {
    throw std::runtime_error("Error: could not open '" + phylogenyFilename + "' for reading");
  }
  
  PhylogeneticTree phyloT;
  inPhyloT >> phyloT;
  inPhyloT.close();
  
  std::ifstream inS(stateTreesFilename.c_str());
  if (!inS.good())
  {
    throw std::runtime_error("Error: could not open '" + stateTreesFilename + "' for reading");
  }
  
  StateTreeVector S;
  inS >> S;
  inS.close();
  
  std::stringstream ss;
  phyloT.writeDOT(S, ss);
  return ss.str();
}

BOOST_PYTHON_MODULE(decifermod)
{
  p::def("enumerateStateTrees", enumerateStateTrees);
  p::def("visualizePhylogeny", visualizePhylogeny);
  p::def("visualizeStateTree", visualizeStateTree);
  p::def("visualizeStateTrees", visualizeStateTrees);
  p::def("computeLogLikelihood", computeLogLikelihood);
  p::def("computeSnvLogLikelihood", computeSnvLogLikelihood);
  p::def("sampleStateTree", sampleStateTree);
  p::def("computeSangerCCF", computeSangerCCF);
  p::def("computeSangerCCFs", computeSangerCCFs);
  p::def("visualizePosteriorStateTrees", visualizePosteriorStateTrees);
  p::def("compareStateTrees", compareStateTrees);
  p::def("exportFeasibleSolutionSpace", exportFeasibleSolutionSpace);
  p::def("visualize", visualize);
}
