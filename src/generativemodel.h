/*
 * generativemodel.h
 *
 *  Created on: 23-oct-2017
 *      Author: M. El-Kebir
 */

#ifndef GENERATIVEMODEL_H
#define GENERATIVEMODEL_H

#include "generativemodel.h"
#include "statetree.h"
#include "readmatrix.h"
#include "beta_distribution.hpp"
#include "phylogenetictree.h"

class GenerativeModel
{
public:
  /// Constructor
  ///
  /// @param m Number of samples
  /// @param n Number of mutations
  /// @param k Number of clusters
  /// @param maxXY Maximum number of maternal/paternal copies
  /// @param maxCN Maximum number of copy number events
  /// @param coverage Read depth
  GenerativeModel(int m,
                  int n,
                  int k,
                  int maxXY,
                  int maxCN,
                  int coverage);
  
  /// Destructor
  virtual ~GenerativeModel();
  
  /// Generate instance
  void generate();
  
  /// Return read matrix
  const ReadMatrix& R() const
  {
    return _R;
  }
  
  /// Return state trees
  const StateTreeVector& S() const
  {
    return _S;
  }
  
  /// Return phylogenetic tree
  const PhylogeneticTree& phyloT() const
  {
    assert(_pPhyloT);
    return *_pPhyloT;
  }
  
  /// Write sample purities
  ///
  /// @param out Output stream
  void writeSamplePurities(std::ostream& out) const;
  
  /// Write mutation clustering
  ///
  /// @param out Output stream
  void writeClustering(std::ostream& out) const;
  
  /// Write mutation properties
  ///
  /// @param out Output stream
  void writeMutationProperties(std::ostream& out) const;
  
  /// Write mutation cluster tree
  ///
  /// @param out Output stream
  void writeDOT(std::ostream& out) const;
  
  /// Write F and Pi
  ///
  /// @param out Output stream
  void writeFandPi(std::ostream& out) const;
  
protected:
  /// Number of samples
  const int _m;
  /// Number of mutations
  const int _n;
  /// Number of clusters
  const int _k;
  /// Maximum number of maternal/paternal copies
  const int _maxXY;
  /// Maximum number of copy number events
  const int _maxCN;
  /// Read depth
  const int _coverage;
  /// Sample purities
  DoubleVector _purityVector;
  /// Descendant cell fractions
  DoubleMatrix _f;
  /// Cluster assignment
  IntVector _z;
  /// Copy number mixture proportions and read counts
  ReadMatrix _R;
  /// State tree vector
  StateTreeVector _S;
  
  /// Cluster tree
  Digraph _clsT;
  /// Root of phylogenetic tree
  Node _rootClsT;
  /// Phylogenetic tree node to cluster index
  IntNodeMap _nodeToClusterClsT;
  /// Cluster index to phylogenetic tree node
  NodeVector _clusterToNodeClsT;
  
  /// Phylogenetic tree
  PhylogeneticTree* _pPhyloT;

  typedef lemon::ListGraph Graph;
  
  void generateDCF(Node v_j,
                   const DoubleMatrix& usages);
  
  /// Convert Prufer sequence to an unrooted tree
  ///
  /// @param seq Prufer sequence
  /// @param T Resulting tree
  static void converPruferToTree(const IntVector& seq, Graph& T);
  
  /// Root unrooted tree
  static void rootTree(const Graph& unrootedT,
                       const Graph::Node vertexUnrootedT,
                       Graph::NodeMap<bool>& visited,
                       Digraph& rootedT,
                       Node vertexRootedT);
  
  /// Generate rooted tree with k vertices uniformly at random
  ///
  /// @param k Number of nodes
  /// @param T Resulting tree
  /// @param root Resulting root node
  static void generateRootedTree(const int k,
                                 Digraph& T,
                                 Node& root);
};

#endif // GENERATIVEMODEL_H
