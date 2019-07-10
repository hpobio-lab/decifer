/*
 * analyzemain.cpp
 *
 *  Created on: 12-jun-2019
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "stategraph.h"
#include "readmatrix.h"
#include <lemon/arg_parser.h>
#include <fstream>
#include "solver.h"

void writeStats(const ReadMatrix& R,
                std::ostream& out)
{
  Solver solver(R, 1, 0, Solver::CLUSTER_DCF, -1, false);
  solver.init();
  
  out << "SNV_index\tSNV_label\tCNA_loss\tnr_state_trees\tnr_distinct_ccfs\tnr_distinct_dcfs\tidentical_ccf_dcf";
  for (int p = 0; p < R.getNrSamples(); ++p)
  {
    out << "\tmax_dcf_sample" << p;
  }
  out << std::endl;
  for (int i = 0; i < R.getNrCharacters(); ++i)
  {
    bool loss = false;
    for (const ReadMatrix::CopyNumberState& cnState : R.getCopyNumberStates(0, i))
    {
      if (cnState._y == 0)
      {
        loss = true;
      }
    }
    
    const Solver::StateEdgeSetVector& scriptT_i = solver.getStateTrees(i);
    
    const int nrSamples = R.getNrSamples();
    
    DoubleVector f_i;
    for (int p = 0; p < nrSamples; ++p)
    {
      f_i.push_back(R.getVAF(p, i));
    }
    
    std::map<IntVector, int> ccf_i;
    std::map<IntVector, int> dcf_i;
    bool identical = true;
    IntVector max_dcf(R.getNrSamples(), 0);
    for (const Solver::StateEdgeSet& T : scriptT_i)
    {
      bool ok = true;
      IntVector ccf, dcf;
      StateTree TT = Solver::convertToStateTreeFromSNVF(R, T, f_i, i);
      
      bool update = true;
      for (int p = 0; p < nrSamples; ++p)
      {
        ccf.push_back(TT.cf(p) * 100);
        dcf.push_back(TT.dcf(p) * 100);
        
        if (ccf.back() < 0 || dcf.back() < 0 || ccf.back() > 100 || dcf.back() > 100)
        {
          ok = false;
        }
        
        if (dcf.back() < max_dcf[p])
        {
          update = false;
        }
      }
      
      if (update)
      {
        max_dcf = dcf;
      }
      
      if (ok)
      {
        if (ccf != dcf)
        {
          identical = false;
        }
        
        ++ccf_i[ccf];
        ++dcf_i[dcf];
      }
    }
    
    out << i << "\t" << R.indexToCharacter(i)
        << "\t" << (loss ? "1" : "0")
        << "\t" << scriptT_i.size()
        << "\t" << ccf_i.size()
        << "\t" << dcf_i.size()
        << "\t" << (identical ? "1" : "0");
    
    for (int p = 0; p < R.getNrSamples(); ++p)
    {
      out << "\t" << max_dcf[p] / 100.;
    }
    
    out << std::endl;
  }
}

int main(int argc, char** argv)
{
  std::string inputStateTreeFilename;
  bool summary = false;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("S", "Input state tree file", inputStateTreeFilename, false)
    .refOption("s", "Summarize", summary)
    .other("Read_matrix");
  ap.parse();
  
  if (!inputStateTreeFilename.empty())
  {
    std::ifstream inS(inputStateTreeFilename.c_str());
    if (!inS.good())
    {
      std::cerr << "Error: could not open '" << inputStateTreeFilename << "' for reading. Will generate state tree file." << std::endl;
    }
    StateGraph::readStateTrees(inS);
  }
  
  if (ap.files().empty())
  {
    std::cerr << "Error: missing input filename" << std::endl;
    return 1;
  }
  
  std::ifstream inR(ap.files()[0].c_str());
  if (!inR.good())
  {
    std::cerr << "Error: could not open '" << ap.files()[0] << "' for reading" << std::endl;
    return 1;
  }
  
  ReadMatrix R;
  try
  {
    inR >> R;
  }
  catch (std::runtime_error& e)
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  inR.close();
  
  if (summary)
  {
    
  }
  else
  {
    writeStats(R, std::cout);
  }
  
  return 0;
}
