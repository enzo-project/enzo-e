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

private: // attributes

  /// The Tree being iterated over
  Tree * tree_;

  /// The trace of the current Node
  NodeTrace node_trace_;

  /// Signal we've reached the last Node
  bool reset_;

};

#endif /* MESH_IT_NODE_HPP */

