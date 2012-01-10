// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Node.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    [\ref Mesh] Declaration of the Node class
///

#ifndef MESH_NODE_HPP
#define MESH_NODE_HPP

class Node {

  /// @class    Node
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  Node() throw();

  /// Destructor
  ~Node() throw();

  /// Set the data payload for the Node
  void set_data(void * data)
  { data_ = data; }

  /// Return the data for the given Node
  void * data() const
  { return data_; }

  /// Refine the node into multiple child nodes
  void refine (int num_children);

  /// Remove all descendent nodes, if any
  void coarsen (int num_children);

  /// return the kth child Node
  Node * child (int k) const;

  /// return whether the node is a leaf node
  bool is_leaf () const
  { return (child_ == 0); }

private: // functions


private: // attributes

  /// Pointer to the data payload for the Node
  void * data_;

  /// Array of child nodes
  Node * child_;

};

#endif /* MESH_NODE_HPP */

