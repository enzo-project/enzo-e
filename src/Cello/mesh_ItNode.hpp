// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_ItNode.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-11
/// @brief    [\ref Mesh] Declaration of the ItNode iterator class
///

#ifndef MESH_IT_NODE_HPP
#define MESH_IT_NODE_HPP

class Tree;
class NodeTrace;
class Node;

class ItNode {

  /// @class    ItNode
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Iterator over Nodes in a Patch

public: // interface

  /// Constructor to iterate over the entire tree
  ItNode(Tree * tree) throw();

  /// Constructor to iterate over one level in the tree
  ItNode(Tree * tree, int level) throw();

  /// Constructor to iterate between two levels in the tree
  ItNode(Tree * tree, int level_lower, int level_upper) throw();

  /// Destructor
  virtual ~ItNode() throw();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    p | lower_level_;
    p | upper_level_;
    WARNING("ItNode::pup","skipping tree_ (aliased)");
    // p | *tree_;
    WARNING("ItNode::pup","skipping node_trace_");
    p | node_trace_;
    p | reset_;
  }

  /// Iterate through all leaf Nodes in the Tree
  Node * next_leaf () throw();

  /// Iterate through all Nodes in the Tree
  Node * next_node () throw();

  /// Return whether the iteration is complete
  bool done() const throw();

  /// Return the nodeTrace corresponding to the current Node
  const NodeTrace * node_trace() const
  { return & node_trace_; };

private: // functions

  /// Go to next common ancestor that with unvisited siblings
  void seek_next_fork_()  //  e.g. 01001111 -> 0100
  {
    while ((node_trace_.index() + 1 == tree_->num_children()) ) {
      node_trace_.pop();
    }
  };

  /// Go to next unvisited sibling of the node trace
  void seek_next_sibling_()  //  e.g. 0100 -> 0101000
  {
    int index = node_trace_.index();
    node_trace_.pop();
    node_trace_.push(index + 1);
  };

  /// Go to next unvisited leaf in the node trace
  void seek_next_leaf_()  //  e.g. 0100 -> 0101000
  {
    while ( ! node_trace_.node()->is_leaf()  &&
	    node_trace_.level() < upper_level_) {
      node_trace_.push(0);
    }
  };


private: // attributes

  /// Do not iterate over nodes in levels lower than lower_level_ (root = 0)
  int lower_level_;

  /// do not iterate over nodes in levels greater than upper_level_ (root = 0)
  int upper_level_;

  /// The Tree being iterated over
  Tree * tree_;

  /// The trace of the current Node
  NodeTrace node_trace_;

  /// Signal we've reached the last Node
  bool reset_;

};

#endif /* MESH_IT_NODE_HPP */

