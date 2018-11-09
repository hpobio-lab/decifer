/*
 * summarizemain.cpp
 *
 *  Created on: 23-oct-2018
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include <lemon/arg_parser.h>
#include <fstream>
#include "readmatrix.h"
#include "posteriorstatetree.h"

int main(int argc, char** argv)
{
  bool hard = false;
  lemon::ArgParser ap(argc, argv);
  ap.other("read_matrix")
    .other("decifer_state_trees")
    .refOption("H", "Hard clustering result", hard);
  ap.parse();
  
  if (ap.files().size() != 2)
  {
    std::cerr << "Error: expected two input files" << std::endl;
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
  
  std::ifstream inSol(ap.files()[1].c_str());
  if (!inSol.good())
  {
    std::cerr << "Error: could not open '" << ap.files()[1] << "' for reading" << std::endl;
    return 1;
  }
  
  if (!hard)
  {
    PosteriorStateTreeMatrix sol;
    try
    {
      inSol >> sol;
    }
    catch (std::runtime_error& e)
    {
      std::cerr << e.what() << std::endl;
      return 1;
    }
    inSol.close();
    
    PosteriorStateTree::writeSummary(R, sol, std::cout);
  }
  else
  {
    StateTreeVector sol;
    try
    {
      inSol >> sol;
    }
    catch (std::runtime_error& e)
    {
      std::cerr << e.what() << std::endl;
      return 1;
    }
    inSol.close();
    
    StateTree::writeSummary(R, sol, std::cout);
  }
  return 0;
}
