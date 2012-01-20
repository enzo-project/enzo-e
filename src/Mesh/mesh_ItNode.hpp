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

class ItNode : public It<Node> {

  /// @class    ItNode
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Iterator over Nodes in a Patch

public: // interface

  /// Constructor
  ItNode(Tree * tree) throw();

  /// Destructor
  virtual ~ItNode() throw();

  /// Iterate through all Nodes in the Tree
  Node * operator++ () throw();

  /// Return whether the iteration is complete
  bool done() const throw();

  /// Return the nodeTrace corresponding to the current Node
  const NodeTrace * node_trace() const
  { return & node_trace_; };

private: // functions

  /// Go to first leaf in the Tree
  void seek_first_leaf_()
  {
    node_trace_.reset();
    while ( ! node_trace_.node()->is_leaf() ) {
      node_trace_.push(0);
    }
  };

  /// Go to next common ancestor with unvisited children
  void seek_next_fork_()  //  e.g. 01001111 -> 0100
  {
    while ((node_trace_.level() > 0) && 
	   (node_trace_.index() + 1 == tree_->num_children()) ) {
      node_trace_.pop();
    }
  };

  /// Go to next unvisited leaf in 
  void seek_next_leaf_()  //  e.g. 0100 -> 0101000
  {
    int index = node_trace_.index();
    node_trace_.pop();
    node_trace_.push(index + 1);
    while ( ! node_trace_.node()->is_leaf() ) {
      node_trace_.push(0);
    }
  };


private: // attributes

  /// The Tree being iterated over
  Tree * tree_;

  /// The trace of the current Node
  NodeTrace node_trace_;

  /// Signal we've reached the last Node
  bool reset_;

};

#endif /* MESH_IT_NODE_HPP */

