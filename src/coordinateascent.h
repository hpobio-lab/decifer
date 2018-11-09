/*
 * coordinateascent.h
 *
 *  Created on: 25-oct-2017
 *      Author: M. El-Kebir
 */

#ifndef COORDINATEASCENT_H
#define COORDINATEASCENT_H

#include "utils.h"
#include "readmatrix.h"
#include "stategraph.h"

/// This class models a coordinate ascent algorithm
class CoordinateAscent
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param k Number of clusters
  /// @param nrGridPoints Number of grid points
  CoordinateAscent(const ReadMatrix& R,
                   int k,
                   int nrGridPoints);
  
  /// Solve
  ///
  /// @param nrIterations Maximum number of iterations
  /// @param seed Random number generator seed
  void solve(int nrIterations, int seed);
  
  /// Solve using intial cluster assignment
  ///
  /// @param z0 Initial cluster assignment
  void solve(IntVector z0);
  
private:
  /// Solve for F
  void solveF();
  
  /// Solve for Z
  void solveZ();
  
  typedef std::map<StateGraph::StateEdgeSet, DoubleTensor> FrequencyMap;
  typedef std::vector<FrequencyMap> FrequencyMapVector;
  typedef std::pair<int, double> IntDoublePair;
  typedef std::map<IntDoublePair, FrequencyMap> FrequencyMapType;
  typedef std::vector<FrequencyMapType> FrequencyMapTypeVector;
  
  void identifyDCFs(const IntVector& Z_j,
                    DoubleMatrix& dcf) const;
  
  void intersect(const DoubleVector& DCFs,
                 const int i,
                 StateGraph::StateEdgeSetSet& commonStateTrees) const;
  
  double getVAF(const DoubleTensor& S_pi) const
  {
    double numerator = 0;
    double denominator = 0;

    for (int x = 0; x < S_pi.size(); ++x)
    {
      for (int y = 0; y < S_pi[x].size(); ++y)
      {
        for (int z = 0; z < S_pi[x][y].size(); ++z)
        {
          numerator += z * S_pi[x][y][z];
          denominator += (x + y) * S_pi[x][y][z];
        }
      }
    }
    
    return numerator / denominator;
  }
  
  bool next(const DoubleMatrix& allowedDCFs, IntVector& indices) const;
  
  void enumerateStateTrees(const ReadMatrix::CopyNumberStateVector& cnStates,
                           const double dcf,
                           const int maxCopyNumber,
                           FrequencyMap& mapToFreqs) const;
  
private:
  /// Read matrix
  const ReadMatrix& _R;
  /// Number of clusters
  const int _k;
  /// \binom(a_pi, d_pi)
  DoubleMatrix _logBinomCoeff;
  /// Grid
  DoubleVector _gridF;
  /// Feasible state trees
  FrequencyMapTypeVector _stateTrees;
  /// Descendant cell fractions
  DoubleMatrix _f;
  /// Cluster assignment
  IntVector _z;
  /// Objective
  double _obj;
};

#endif // COORDINATEASCENT_H
