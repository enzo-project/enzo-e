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

  void pup (PUP::er &p);

  /// Set the data payload for the Node
  void set_data(void * data)
  { data_ = data;
    have_data_ = data != NULL;
  }

  /// Return the data for the given Node
  void * data() const
  { return data_; }

  /// Refine the node into multiple child nodes, and return the number
  /// of new nodes created
  int refine (int num_children);

  /// Remove all descendent nodes, and return the number of nodes deleted
  int coarsen (int num_children);

  /// return the kth child Node
  Node * child (int k);

  /// return whether the node is a leaf node
  bool is_leaf () const
  { return (child_array_.size() == 0); }

private: // functions


private: // attributes

  /// Whether data_ is null
  bool have_data_;

  /// Pointer to the data payload for the Node
  void * data_;

  /// Array of child nodes
  std::vector<Node> child_array_;

};

#endif /* MESH_NODE_HPP */

