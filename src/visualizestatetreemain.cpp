/*
 * visualizestatetree.cpp
 *
 *  Created on: 10-nov-2017
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "statetree.h"
#include <lemon/arg_parser.h>
#include <fstream>

int main(int argc, char** argv)
{
  int idx = -1;
  int m = -1;
  lemon::ArgParser ap(argc, argv);
  ap.refOption("idx", "State tree index", idx, true)
    .refOption("m", "Number of samples", m)
    .other("filename");
  ap.parse();
  
  std::ifstream inT(ap.files()[0].c_str());
  if (!inT.good())
  {
    std::cerr << "Error: could not open '" << ap.files()[0] << "' for reading" << std::endl;
    return 1;
  }
  
  StateTreeVector T;
  try
  {
    inT >> T;
  }
  catch (std::runtime_error& e)
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  
  if (!(0 <= idx && idx < T.size()))
  {
    std::cerr << "Error: invalid index" << std::endl;
    return 1;
  }
  
  if (m == -1)
  {
    T[idx].writeDOT(std::cout);
  }
  else
  {
    T[idx].writeDOT(m, std::cout);
  }
  
  return 0;
}
