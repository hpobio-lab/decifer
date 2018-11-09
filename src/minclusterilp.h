/*
 * minclusterilp.h
 *
 *  Created on: 21-oct-2018
 *      Author: M. El-Kebir
 */

#ifndef MINCLUSTERILP_H
#define MINCLUSTERILP_H

#include "solver.h"

class MinClusterIlp : public Solver
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param kMax Number of clusters
  MinClusterIlp(const ReadMatrix& R,
                int kMax);
  
  /// Destructor
  virtual ~MinClusterIlp()
  {
  };
  
  /// Solve
  ///
  /// @param nrThreads Number of threads
  /// @param timeLimit Time limit in seconds (for initial ILP)
  virtual bool solve(int nrThreads,
                     int timeLimit) = 0;
  
  /// Export ILP model
  ///
  /// @param filename Filename
  virtual void exportModel(const std::string& filename) = 0;
  
  /// Initialize solver
  virtual void init();
  
  /// Get minimum cluster size
  int getMinK() const
  {
    return _solMinK;
  }
  
protected:
  /// Initialize variables
  virtual void initVariables() = 0;
  
  /// Initialize constraints
  virtual void initConstraints() = 0;
  
  /// Initialize objective
  virtual void initObjective() = 0;
  
protected:
  /// Minimum cluster size
  int _solMinK;
};

#endif // MINCLUSTERILP_H
