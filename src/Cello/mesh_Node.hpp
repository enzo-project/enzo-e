// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Node.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    [\ref Mesh] Declaration of the Node class
///

#ifndef MESH_NODE_HPP
#define MESH_NODE_HPP

#ifdef REMOVE_PATCH
#else /* REMOVE_PATCH */
class Patch;
#endif

class Node {

  /// @class    Node
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  Node() throw();

  /// Destructor
  ~Node() throw();

#ifdef CONFIG_USE_CHARM
  void pup (PUP::er &p);
#endif

  /// Set the data payload for the Node
  void set_data(void * data)
  { data_ = data;
#ifdef CONFIG_USE_CHARM
    have_data_ = data != NULL;
#endif
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
  Node * child (int k) const;

  /// return whether the node is a leaf node
  bool is_leaf () const
  { return (child_array_ == 0); }

private: // functions


private: // attributes

#ifdef CONFIG_USE_CHARM
  /// Whether data_ is null
  bool have_data_;
#endif

  /// Pointer to the data payload for the Node
  void * data_;

#ifdef CONFIG_USE_CHARM
  /// Number of children: used for CHARM++ pup()
  char size_; 
#endif

  /// Array of child nodes
  Node * child_array_;

};

#endif /* MESH_NODE_HPP */

