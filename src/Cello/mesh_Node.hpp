// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Node.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    [\ref Mesh] Declaration of the Node class
///

#ifndef MESH_NODE_HPP
#define MESH_NODE_HPP

class Patch;

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
  void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    TRACE1("data_ = %p",data_);
    // int * test = new int[1];
    // *test = 2;
    // int *& test2 = test;
    // TRACE1 ("test2 = %d",*test2);
    // Patch *& data_alias = (Patch *) data_;
    p | *((CProxy_Patch *)data_);
    TRACE0;
    p | size_;
    if (p.isUnpacking()) child_array_ = new Node[size_];
    PUParray(p,child_array_,size_);
  };
#endif

  /// Set the data payload for the Node
  void set_data(void * data)
  { data_ = data; }

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

