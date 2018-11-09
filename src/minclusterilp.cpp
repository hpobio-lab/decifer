/*
 * minclusterilp.cpp
 *
 *  Created on: 21-oct-2018
 *      Author: M. El-Kebir
 */

#include "minclusterilp.h"

MinClusterIlp::MinClusterIlp(const ReadMatrix& R,
                             int kMax)
  : Solver(R, kMax, 0)
  , _solMinK(0)
{
}

void MinClusterIlp::init()
{
  Solver::init();
  
  initVariables();
  initConstraints();
  initObjective();
}
