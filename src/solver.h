/*
 * solver.h
 *
 *  Created on: 12-nov-2017
 *      Author: M. El-Kebir
 */

#ifndef SOLVER_H
#define SOLVER_H

#include "utils.h"
#include "readmatrix.h"
#include "stategraph.h"
#include "statetree.h"

class Solver
{
public:
  /// Constructor
  ///
  /// @param R Read matrix
  /// @param k Number of clusters
  /// @param nrSegments Number of segments for piecewise linear approximation
  Solver(const ReadMatrix& R,
         int k,
         int nrSegments);
  
  /// Destructor
  virtual ~Solver()
  {
  }
  
  /// Return DCF
  ///
  /// @param p Sample
  /// @param j Cluster
  double getD(int p, int j) const
  {
    assert(0 <= p && p < _R.getNrSamples());
    assert(0 <= j && j < _k);
    
    return _solD[j][p];
  }
  
  /// Return DCF
  const DoubleMatrix& getF() const
  {
    return _solD;
  }
  
  /// Return pi
  ///
  /// @param j Cluster
  int getPi(int j) const
  {
    assert(0 <= j && j < _k);
    
    return _solPi[j];
  }
  
  /// Return pi
  const DoubleVector& getPi() const
  {
    return _solPi;
  }
  
  /// Return log likelihood
  double getLogLikelihood() const
  {
    return _logLikelihood;
  }
  
  /// Write mutation clustering
  ///
  /// @param out Output stream
  virtual void writeClustering(std::ostream& out) const;
  
  /// Write mutation properties
  ///
  /// @param out Output stream
  virtual void writeMutationProperties(std::ostream& out) const;
  
  /// Initialize solver
  virtual void init();
  
  /// Initialize solver
  ///
  /// @param i SNV
  virtual void init(int i);
  
  /// Set solution
  void setSolution(const DoubleMatrix& d, const DoubleVector& pi)
  {
    assert(d.size() == _k);
    assert(d[0].size() == _R.getNrSamples());
    assert(pi.size() == _k);
    
    _solD = d;
    _solPi = pi;
    
    updateLogLikelihood();
  }
  
  /// Write solution
  ///
  /// @param out Output stream
  void writeSolution(std::ostream& out) const;
  
  /// Read solution
  ///
  /// @param in Input stream
  bool readSolution(std::istream& in);

  /// Read solution
  ///
  /// @param in Input stream
  /// @param outF Output matrix f
  /// @param outPi Output vector pi
  static void readSolution(std::istream& in, DoubleMatrix& outF, DoubleVector& outPi);
  
  /// Compute log likelihood
  ///
  /// @param i SNV
  double getLogLikelihood(int i) const;
  
protected:
  typedef ReadMatrix::CopyNumberStateVector CopyNumberStateVector;
  typedef std::vector<IntPair> IntPairVector;
  typedef std::vector<IntPairVector> IntPairMatrix;
  typedef StateGraph::StateEdgeSet StateEdgeSet;
  typedef std::vector<StateEdgeSet> StateEdgeSetVector;
  typedef std::vector<StateEdgeSetVector> StateEdgeSetMatrix;
  typedef std::vector<DoubleTensor> Double4Matrix;
  typedef std::vector<Double4Matrix> Double5Matrix;
  
  /// Update log likelihood
  double updateLogLikelihood();
  
  /// Return whether provided DCF is feasible
  ///
  /// @param d_jp DCF
  /// @param xyStar Copy number state where SNV was introduced
  /// @param cnStates Copy number states
  /// @param T_it State tree
  /// @param maxCopyNumber Maximum copy number
  bool isFeasible(const double d_jp,
                  const IntPair& xyStar,
                  const CopyNumberStateVector& cnStates,
                  const StateEdgeSet& T_it,
                  const int maxCopyNumber) const;
  
  /// Construct a state tree from SNVF vector
  ///
  /// @param T_it State tree
  /// @param f_i SNVF
  /// @param i SNV
  StateTree convertToStateTreeFromSNVF(const StateEdgeSet& T_it,
                                       const DoubleVector& f_i,
                                       const int i) const;
  
  /// Construct a state tree from DCF vector
  ///
  /// @param T_it State tree
  /// @param d_i DCF
  /// @param i SNV
  StateTree convertToStateTreeFromDCF(const StateEdgeSet& T_it,
                                      const DoubleVector& d_i,
                                      const int i) const;
  
  /// Return whether state tree T_(i,t) is enabled for SNV i
  ///
  /// @param i SNV
  /// @param t State tree
  /// @param j Cluster
  virtual bool isEnabled(int i, int t, int j) const
  {
    return false;
  }
  
  /// Initialize piecewise linear approximation
  virtual void initPWLA();
  
protected:
  /// Read matrix
  const ReadMatrix& _R;
  /// Number of clusters
  const int _k;
  /// Number of segments for piecewise linear approximation
  const int _nrSegments;
  /// Log factorial look-up table
  DoubleVector _logFactorial;
  /// State trees \f$\mathcal{T}(\mu_i)\f$
  StateEdgeSetMatrix _scriptT;
  
  /// x[i][t][j][p][l]
  Double5Matrix _x;
  /// G[i][t][j][p][l]
  Double5Matrix _hatG;
  /// _dOverallLB[j][p]
  DoubleMatrix _dOverallLB;
  /// _dOverallUB[j][p]
  DoubleMatrix _dOverallUB;
  
  /// _dLB[i][t][p]
  DoubleTensor _dLB;
  /// _dUB[i][t][p]
  DoubleTensor _dUB;
  /// denominator[i][p]
  DoubleMatrix _denominator;
  /// numerator[i][t][p]
  DoubleTensor _numerator;
  /// _xyStar[i][t];
  IntPairMatrix _xyStar;
  /// \f$\log \binom{a_{p,i}}{d_{p,i}}\f$
  DoubleMatrix _logBinomCoeff;
  /// _solD[j][p]
  DoubleMatrix _solD;
  /// _solPi[j]
  DoubleVector _solPi;
  /// _solZ[i]
  IntVector _solZ;
  /// Log likelihood
  double _logLikelihood;
  
public:
  /// State tree set
  ///
  /// @param i SNV
  const StateEdgeSetVector& getScriptT(int i) const
  {
    assert(0 <= i && i < _R.getNrCharacters());
    return _scriptT[i];
  }
};

#endif // SOLVER_H
