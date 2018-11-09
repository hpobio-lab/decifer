/*
 * phylogenetictree.h
 *
 *  Created on: 20-jan-2018
 *      Author: M. El-Kebir
 */

#ifndef PHYLOGENETICTREE_H
#define PHYLOGENETICTREE_H

#include "utils.h"
#include "statetree.h"

class PhylogeneticTree
{
public:
  /// Default constructor, only used for deserialization
  PhylogeneticTree();
  
  /// Constructor
  ///
  /// @param S State tree vector
  /// @param z Cluster assignment
  /// @param mutTree Tree relating mutations
  /// @param rootMutTree Mutation tree node to use to extend phylogenetic tree
  /// @param nodeToClusterMutTree Mutation tree node to cluster index
  /// @param clusterToNodeMutTree Cluster index to phylogenetic tree node
  PhylogeneticTree(const StateTreeVector& S,
                   const IntVector& z,
                   const Digraph& mutTree,
                   const Node rootMutTree,
                   const IntNodeMap& nodeToClusterMutTree,
                   const NodeVector& clusterToNodeMutTree);
  
  /// Write phylogenetic tree in DOT format
  ///
  /// @param S State tree vector
  /// @param out Output stream
  void writeDOT(const StateTreeVector& S,
                std::ostream& out) const;
  
  /// Write phylogenetic tree in DOT format
  ///
  /// @param out Output stream
  void writeDOT(std::ostream& out) const;
  
  /// Write mutation clustering
  ///
  /// @param out Output stream
  void writeClustering(std::ostream& out) const;
  
  /// Write mutation tree
  ///
  /// @param out Output stream
  void writeMutationTreeDOT(std::ostream& out) const;
  
  /// Sample vertices
  ///
  /// @param purityVector Vector of sample purities
  /// @param clusterUsages Matrix of mutation cluster usages (samples by clusters)
  /// @param S State trees to be updated
  void sample(const DoubleVector& purityVector,
              const DoubleMatrix& clusterUsages,
              StateTreeVector& S);
  
private:
  typedef std::vector<IntPair> IntPairVector;
  typedef Digraph::NodeMap<IntPairVector> StateVectorNodeMap;
  typedef Digraph::NodeMap<IntPairSet> StateSetNodeMap;
  
  /// Randomly initialize phylogenetic tree
  ///
  /// @param S State tree vector
  /// @param z Cluster assignment
  /// @param mutT Mutation tree
  /// @param nodeToClusterMutT Mutation tree node to cluster index
  /// @param clusterToNodeMutT Cluster index to phylogenetic tree node
  /// @param nodeMutT Mutation tree node to use to extend phylogenetic tree
  void init(const StateTreeVector& S,
            const IntVector& z,
            const Digraph& mutT,
            const IntNodeMap& nodeToClusterMutT,
            const NodeVector& clusterToNodeMutT,
            const Node nodeMutT);
  
  /// Uniformly at random pick a descendant of specified node and tree
  ///
  /// @param T Tree
  /// @param v Ancestral node
  Node pickDescendant(const Digraph& T,
                      Node v) const;
  
  /// Uniformly at random pick a parental node in T
  ///
  /// @param S_c State tree
  /// @param ci Character-state pair of new child
  Node pickParent(const StateTree& S_c,
                  const IntPair& ci) const;
  
  /// Extend the phylogenetic tree _T with the specified characters.
  /// Return the node on whose incoming edge the mutations occur.
  ///
  /// @param S State tree vector
  /// @param clusterIndex Cluster Index
  /// @param characters Set of characters (mutations) that occur on the same phylogenetic branch
  /// @param insertionNode Root of the subtree
  Node extend(const StateTreeVector& S,
              const int clusterIdx,
              const IntVector& characters,
              Node insertionNode);
  
  /// Update cluster cover recursively.
  /// Initial call updateCluster(_rootT, -1)
  ///
  /// @param v Node
  /// @param cluster Cluster index
  void updateClusterCover(Node v, int cluster);
  
  /// Update state proportions recursively.
  /// Initial call updateStateProportions(_rootT)
  ///
  /// @param v Node
  void updateStateProportions(Node v,
                              int c,
                              int i,
                              StateTree& S_c);
  
  int getNrCharacters() const
  {
    IntSet set;
    for (const auto& kv : _charStateToNodeT)
    {
      set.insert(kv.first.first);
    }
    return set.size();
  }
  
private:
  /// Phylogenetic tree
  Digraph _T;
  /// Root node
  Node _rootT;
  /// Tree node to cluster index
  IntNodeMap _nodeToClusterT;
  /// Cluster index to tree node
  NodeVector _clusterToNodeT;
  /// Arc labeling
  StateSetNodeMap _nodeToCharStateT;
  /// Character-state pair to node
  std::map<IntPair, Node> _charStateToNodeT;
  /// Mutation node indicator
  BoolNodeMap _isMutationNodeT;
  /// Indicates the most recently introduced mutation cluster for each node in T
  IntNodeMap _clusterCoverT;
  /// Usages
  DoubleVectorNodeMap _U;
  
  friend std::ostream& operator<<(std::ostream& out, const PhylogeneticTree& phyloT);
  friend std::istream& operator>>(std::istream& in, PhylogeneticTree& phyloT);
};

/// Output phylogenetic tree
///
/// @param out Output stream
/// @param phyloT Phylogenetic tree
std::ostream& operator<<(std::ostream& out, const PhylogeneticTree& phyloT);

/// Input phylogenetic tree
///
/// @param in Input stream
/// @param phyloT Phylogenetic tree
std::istream& operator>>(std::istream& in, PhylogeneticTree& phyloT);

std::ostream& operator<<(std::ostream& out, const IntPairSet& s);

std::istream& operator>>(std::istream& in, IntPairSet& s);

std::ostream& operator<<(std::ostream& out, const DoubleVector& s);

std::istream& operator>>(std::istream& in, DoubleVector& s);

#endif // PHYLOGENETICTREE_H
