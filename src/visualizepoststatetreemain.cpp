/*
 * visualizestatetree.cpp
 *
 *  Created on: 10-nov-2017
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "posteriorstatetree.h"
#include <lemon/arg_parser.h>
#include <fstream>

int main(int argc, char** argv)
{
  int m = -1;
  int i = -1;
  int t = -1;
  lemon::ArgParser ap(argc, argv);
  ap.refOption("i", "SNV index", i, true)
    .refOption("t", "State tree index", t, true)
    .refOption("m", "Number of samples", m)
    .other("decifer_posterior_state_trees");
  ap.parse();
  
  if (ap.files().size() != 1)
  {
    std::cerr << "Error: expected input file" << std::endl;
    return 1;
  }
  
  std::ifstream inSol(ap.files()[0].c_str());
  if (!inSol.good())
  {
    std::cerr << "Error: could not open '" << ap.files()[0] << "' for reading" << std::endl;
    return 1;
  }
  
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
  
  if (!(0 <= i && i < sol.size()))
  {
    std::cerr << "Incorrect SNV index" << std::endl;
    return 1;
  }
  
  if (!(0 <= t && t < sol[i].size()))
  {
    std::cerr << "Incorrect state tree index" << std::endl;
    return 1;
  }
  
  if (m == -1)
  {
    sol[i][t]._T.writeDOT(std::cout);
  }
  else
  {
    sol[i][t]._T.writeDOT(m, std::cout);
  }
  
  return 0;
}
