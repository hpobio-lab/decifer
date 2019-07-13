/*
 * empre.h
 *
 *  Created on: 2-jul-2019
 *      Author: M. El-Kebir
 */

#ifndef EMPRE_H
#define EMPRE_H

#include "utils.h"
#include "readmatrix.h"
#include "stategraph.h"
#include "solver.h"
#include "statetree.h"
#include "posteriorstatetree.h"
#include "hardclusterilp.h"
#include "softclusterilp.h"
#include "clusterilp.h"
#include "incrementalsolver.h"

class EMPre : public Solver
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param preClustering Pre clustering
  /// @param k Number of clusters
  /// @param nrSegments Number of segments for piecewise linear approximation
  /// @param statType Summary statistic to use for clustering
  /// @param precisionBetaBin Precision parameter for beta binomial
  /// @param forceTruncal Force the presence of a dominant truncal cluster
  EMPre(const ReadMatrix& R,
        const IntMatrix& preClustering,
        int k,
        int nrSegments,
        ClusterStatisticType statType,
        double precisionBetaBin,
        bool forceTruncal);
  
  /// Destructor
  virtual ~EMPre()
  {
  }
  
  /// Solve
  ///
  /// @param restart Restart number
  /// @param seed Random number generator seed for k-Means.
  /// In case seed < 0, ILP will be used to seed the EM algorithm.
  /// @param maxIterations Maximum number of iterations, ignored when set to -1
  /// @param nrThreads Number of threads
  /// @param timeLimit Time limit in seconds (for initial ILP)
  /// @param verbose Verbose
  /// @param nrDownSampledSNVs Fraction of SNVs to consider for initialization
  /// @param memoryLimit Memory limit
  bool solve(int restart,
             int seed,
             int maxIterations,
             int nrThreads,
             int timeLimit,
             bool verbose,
             int nrDownSampledSNVs,
             int memoryLimit);
  
  /// Get sum of _y[i][*][y]
  ///
  /// @param i SNV
  /// @param j Cluster
  double getGamma(int i, int j) const
  {
    double gamma = 0;
    
    for (int t = 0; t < _scriptT[i].size(); ++t)
    {
      gamma += _solY[i][t][j];
    }
    
    return gamma;
  }
  
  /// Write mutation properties
  ///
  /// @param out Output stream
  virtual void writeMutationProperties(std::ostream& out) const;
  
  /// Initialize solver
  void init();
  
  /// Return posterior state trees
  const PosteriorStateTreeMatrix& getPosteriorStateTrees() const
  {
    return _solT;
  }
  
  /// Initialize hot start
  ///
  /// @param y Clustering and state tree assignment
  void initHotStart(const DoubleTensor& y)
  {
    _initY = y;
  }
  
  const DoubleTensor& getInitY() const
  {
    return _initY;
  }
  
  double getLogLikelihoodGamma(int i) const;
  
  double getLogLikelihoodGamma() const;
  
  double roundD()
  {
    const int n = _R.getNrCharacters();
    const int m = _R.getNrSamples();
    
    for (int i = 0; i < n; ++i)
    {
      for (int t = 0; t < _scriptT[i].size(); ++t)
      {
        for (int j = 0; j < _k; ++j)
        {
          if (g_tol.nonZero(_solY[i][t][j]))
          {
            for (int p = 0; p < m; ++p)
            {
              _solD[j][p] = std::max(_dLB[i][t][p], _solD[j][p]);
              _solD[j][p] = std::min(_dUB[i][t][p], _solD[j][p]);
            }
          }
        }
      }
    }
  }
  
protected:
  /// E step, estimate y and pi
  ///
  /// @param nrThreads Number of threads to use
  /// @param verbose Verbose
  virtual bool stepE(int nrThreads,
                     bool verbose) = 0;
  
  /// M step, estimate d
  ///
  /// @param nrThreads Number of threads to use
  /// @param verbose Verbose
  virtual bool stepM(int nrThreads,
                     bool verbose) = 0;
  
  /// Initialize piecewise linear approximation
  virtual void initPWLA();
  
  virtual std::unique_ptr<SoftClusterIlp> createSoftClusterIlpSolver(const ReadMatrix& R,
                                                                     const IntMatrix& preClustering) = 0;
  
  
  virtual bool isEnabled(int i, int t, int j) const
  {
    assert(0 <= i && i < _R.getNrCharacters());
    assert(0 <= t && t < _scriptT[i].size());
    assert(0 <= j && j < _k);
    return g_tol.nonZero(_solY[i][t][j]);
  }
  
  void updatePWLA();
  
  bool initializeD(int seed,
                   int nrThreads,
                   int timeLimit,
                   bool verbose,
                   int nrDownSampledSNVs,
                   int memoryLimit);
  
protected:
  /// Pre clustering
  IntMatrix _preClustering;
  /// _solY[i][t][j]
  DoubleTensor _solY;
  /// Solution state trees
  PosteriorStateTreeMatrix _solT;
  /// Initial hot start
  DoubleTensor _initY;
};

#endif // EM_H

