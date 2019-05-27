/*
 * generatestatetreesmain.cpp
 *
 *  Created on: 28-mar-2019
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "stategraph.h"
#include "readmatrix.h"
#include <lemon/arg_parser.h>
#include <fstream>
#include "solver.h"

bool next(const int maxXY,
          const int maxCN,
          const int nrCnStates,
          const IntVector& minEnumState,
          const IntVector& maxEnumState,
          IntVector& enumState)
{
  assert(enumState.size() == maxCN);
  
  IntVector newEnumState = enumState;
  
  int i = maxCN - 1;
  for (; i >= 0; --i)
  {
    if (newEnumState[i] < maxEnumState[i] - 1)
    {
      ++newEnumState[i];
      break;
    }
  }
  
  if (i < 0)
  {
    return false;
  }
  else
  {
    for (int j = i + 1; j < maxCN; ++j)
    {
      newEnumState[j] = newEnumState[j-1] + 1;
    }
  }
  
  enumState = newEnumState;
  
  return true;
}

int main(int argc, char** argv)
{
  int maxXY = 3;
  int maxCN = 2;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("maxXY", "Maximum number of maternal/paternal copies (default: 3)", maxXY, false)
    .refOption("maxCN", "Maximum number of copy number events (default: 2)", maxCN, false);
  ap.parse();
  
  if (!ap.files().empty())
  {
    for (const std::string& filename : ap.files())
    {
      std::ifstream inR(filename.c_str());
      if (!inR.good())
      {
        std::cerr << "Error: could not open '" << filename << "' for reading" << std::endl;
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
      
      Solver solver(R, 1, 1, Solver::CLUSTER_DCF, -1, false);
      solver.init();
    }
  }
  else
  {
    std::vector<IntPair> cnStates;
    for (int x = 0; x <= maxXY; ++x)
    {
      for (int y = x; y <= maxXY; ++y)
      {
        cnStates.push_back(IntPair(x, y));
        if (x != y)
        {
          cnStates.push_back(IntPair(y, x));
        }
      }
    }
    
    IntVector enumState(maxCN);
    for (int i = 0; i < maxCN; ++i)
    {
      enumState[i] = i;
    }
    
    IntVector minEnumState = enumState;
    
    IntVector maxEnumState(maxCN);
    for (int i = 0; i < maxCN; ++i)
    {
      maxEnumState[i] = cnStates.size() - (maxCN - i - 1);
    }
    
    do
    {
      IntPairSet L;
      int maxCopy = 0;
      for (int i = 0; i < maxCN; ++i)
      {
        L.insert(cnStates[enumState[i]]);
        maxCopy = std::max(maxCopy,
                           std::max(cnStates[enumState[i]].first,
                                    cnStates[enumState[i]].second));
      }
      
      StateGraph::getStateTrees(L, maxCopy + 1);
      std::cerr << enumState[0] << " " << enumState[1] << std::endl;
    } while (next(maxXY, maxCN, cnStates.size(), minEnumState, maxEnumState, enumState));
  }
  
  StateGraph::writeStateTrees(std::cout);
  
  return 0;
}
