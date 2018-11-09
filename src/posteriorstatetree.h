/*
 *  posteriorstatetree.h
 *
 *   Created on: 19-nov-2017
 *       Author: M. El-Kebir
 */

#ifndef POSTERIORSTATETREE_H
#define POSTERIORSTATETREE_H

#include "utils.h"
#include "statetree.h"
#include "readmatrix.h"

class PosteriorStateTree;

typedef std::vector<PosteriorStateTree> PosteriorStateTreeVector;

typedef std::vector<PosteriorStateTreeVector> PosteriorStateTreeMatrix;

class PosteriorStateTree
{
public:
  PosteriorStateTree(const StateTree& T,
                     double gamma,
                     int t,
                     int j)
    : _T(T)
    , _gamma(gamma)
    , _t(t)
    , _j(j)
  {
  }

  const StateTree _T;
  const double _gamma;
  const int _t;
  const int _j;
  
  static void writeDOT(const PosteriorStateTreeVector& T,
                       std::ostream& out);
  
  static void writeSummary(const ReadMatrix& R,
                           const PosteriorStateTreeMatrix& sol,
                           std::ostream& out);
};

/// Output state tree posterior matrix
///
/// @param out Output stream
/// @param matT State tree posterior matrix
std::ostream& operator<<(std::ostream& out,
                         const PosteriorStateTreeMatrix& matT);

/// Input state tree posterior matrix
///
/// @param in Input stream
/// @param vecS State tree posterior matrix
std::istream& operator>>(std::istream& in,
                         PosteriorStateTreeMatrix& matT);

#endif // POSTERIORSTATETREE_H
