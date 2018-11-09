/*
 * statetreessampler.h
 *
 *  Created on: 6-feb-2016
 *      Author: M. El-Kebir
 */

#ifndef STATETREESSAMPLER_H
#define STATETREESSAMPLER_H

#include "utils.h"
#include "stategraph.h"
#include "statetree.h"
#include <random>

/// This class samples n state trees
class StateTreesSampler
{
public:
  /// Constructor
  ///
  /// @param maxXY Maximum number of maternal and/or paternal copies
  /// @param n Number of mutations
  /// @param maxCopyEvents Maximum number of copy number events
  StateTreesSampler(int maxXY,
                    int n,
                    int maxCopyEvents);
  
  typedef std::vector<StateTree> StateTreeVector;
  typedef std::vector<IntPair> IntPairVector;
  typedef std::map<StateGraph::CnaTriple, int> CnaTripleMap;
  typedef CnaTripleMap::const_iterator CnaTripleMapIt;
  typedef StateGraph::CnaTripleVector CnaTripleVector;
  typedef std::vector<std::string> StringVector;
  typedef std::vector<StateGraph::StateEdgeSet> StateEdgeSetVector;

  StateTreeVector sample();

private:
  void initSampleSpace();
  IntPair sampleCopyState();

private:
  /// Maximum number of maternal and/or paternal copies
  const int _maxXY;
  /// Number of mutations
  const int _n;
  /// Maximum number of copy number events
  const int _maxCopyEvents;
  /// Vector of copy number states (biased towards (1,1))
  IntPairVector _sampleSpace;
};

#endif // STATETREESSAMPLER_H
