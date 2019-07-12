/*
 *  posteriorstatetree.cpp
 *
 *   Created on: 19-nov-2017
 *       Author: M. El-Kebir
 */

#include "posteriorstatetree.h"
#include <iomanip>

void PosteriorStateTree::writeDOT(const PosteriorStateTreeVector& T,
                                  std::ostream& out)
{
  out << "digraph T {" << std::endl;
  
  for (int t = 0; t < T.size(); ++t)
  {
    T[t]._T.writeClusterDOT(T[t]._gamma, t, out);
  }
  
  out << "}" << std::endl;
}

std::ostream& operator<<(std::ostream& out,
                         const PosteriorStateTreeMatrix& matT)
{
  const int n = matT.size();
  out << n << " #characters" << std::endl;
  for (int i = 0; i < n; ++i)
  {
    out << matT[i].size() << " #state trees" << std::endl;
    for (const PosteriorStateTree& T : matT[i])
    {
      out << T._gamma << std::endl;
      out << T._j << std::endl;
      out << T._t << std::endl;
      out << T._T << std::endl;
    }
  }
  
  return out;
}

std::istream& operator>>(std::istream& in,
                         PosteriorStateTreeMatrix& matT)
{
  g_lineNumber = 0;
  
  std::string line;
  getline(in, line);
  std::stringstream ss(line);
  
  int n = -1;
  ss >> n;
  
  if (n < 1)
  {
    throw std::runtime_error(getLineNumber() + "Error: incorrect number of characters");
  }
  
  matT = PosteriorStateTreeMatrix(n);
  for (int i = 0; i < n; ++i)
  {
    int nrStateTrees = -1;
    getline(in, line);
    ss.clear();
    ss.str(line);
    ss >> nrStateTrees;
    
    if (nrStateTrees < 1)
    {
      throw std::runtime_error(getLineNumber() + "Error: incorrect number of state trees");
    }
    
    for (int l = 0; l < nrStateTrees; ++l)
    {
      double gamma = -1;
      getline(in, line);
      ss.clear();
      ss.str(line);
      ss >> gamma;
      if (g_tol.less(gamma, 0.) || g_tol.less(1., gamma))
      {
        throw std::runtime_error(getLineNumber() + "Error: invalid posterior probability");
      }
      
      int j = -1;
      getline(in, line);
      ss.clear();
      ss.str(line);
      ss >> j;
      if (j < 0)
      {
        throw std::runtime_error(getLineNumber() + "Error: invalid cluster assignment");
      }
      
      int t = -1;
      getline(in, line);
      ss.clear();
      ss.str(line);
      ss >> t;
      if (t < 0)
      {
        throw std::runtime_error(getLineNumber() + "Error: invalid tree index");
      }
      
      StateTree T;
      in >> T;
      
      matT[i].emplace_back(T, gamma, t, j);
      getline(in, line);
    }
  }
  
  return in;
}

void PosteriorStateTree::writeSummary(const ReadMatrix& R,
                                      const PosteriorStateTreeMatrix& sol,
                                      std::ostream& out)
{
  const int n = R.getNrCharacters();
  const int m = R.getNrSamples();
  
  out << "SNV_index\t"
      << "SNV_label\t"
      << "sample_index\t"
      << "sample_label\t"
      << "ref\t"
      << "alt\t"
      << "VAF\t"
      << "sanger_n_mut\t"
      << "sanger_CCF\t"
      << "decifer_CCF\t"
      << "decifer_DCF\t"
      << "decifer_cluster_CCF\t"
      << "decifer_cluster_DCF\t"
      << "decifer_cluster_index\t"
      << "decifer_state_tree_index\t"
      << "nr_nonzero_z"
      << std::endl;

  for (int i = 0; i < n; ++i)
  {
    for (int p = 0; p < m; ++p)
    {
      int ref = R.getRef(p, i);
      int var = R.getVar(p, i);
      double vaf = R.getVAF(p, i);
      double purity = 1;
      auto sanger = R.computeSangerCCF(ref, var,
                                       purity, R.getCopyNumberStates(p, i));
      int sanger_n_mut = sanger.first;
      double sangerCCF = sanger.second;
      
      double max_gamma = -1;
      int max_t = -1;
      for (int t = 0; t < sol[i].size(); ++t)
      {
        if (sol[i][t]._gamma > max_gamma)
        {
          max_t = t;
          max_gamma = sol[i][t]._gamma;
        }
      }
      
      IntSet zSet;
      IntPairSet mutTumorSet;
      IntPairSet nonMutTumorSet;
      const StateTree& T_it = sol[i][max_t]._T;
      for (int j = 0; j != T_it.numVertices(); ++j)
      {
        if (T_it.isPresent(j))
        {
          int x_j = T_it.x(j);
          int y_j = T_it.y(j);
          int z_j = T_it.z(j);
          if (z_j > 0)
          {
            zSet.insert(z_j);
          }
          if (x_j != 1 || y_j != 1)
          {
            if (z_j > 0)
            {
              mutTumorSet.insert(IntPair(x_j, y_j));
            }
            else
            {
              nonMutTumorSet.insert(IntPair(x_j, y_j));
            }
          }
        }
      }
      
      double deciferCCF = sol[i][max_t]._T.maxLikelihoodCCF(p, vaf);
      double deciferDCF = sol[i][max_t]._T.maxLikelihoodDCF(p, vaf);
      double deciferClusterDCF = sol[i][max_t]._T.dcf(p);
      double deciferClusterCCF = sol[i][max_t]._T.cf(p);
      
      out << i << "\t" << R.indexToCharacter(i) << "\t"
          << p << "\t" << R.indexToSample(p) << "\t"
          << ref << "\t" << var << "\t"
          << R.getVAF(p, i) << "\t"
          << sanger_n_mut << "\t" << sangerCCF << "\t"
          << deciferCCF << "\t" << deciferDCF << "\t"
          << deciferClusterCCF << "\t" << deciferClusterDCF << "\t"
          << sol[i][max_t]._j << "\t" << max_t << "\t"
          << zSet.size() //<< "\t"
//          << nonMutTumorSet.size() << "\t"
//          << mutTumorSet.size()
          << std::endl;
    }
  }
}
