// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Tree.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    [\ref Mesh] Declaration of the Tree class
///

#ifndef MESH_TREE_HPP
#define MESH_TREE_HPP

class Tree {

  /// @class    Tree
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  Tree(int d, int r) throw();

  /// Destructor
  ~Tree() throw();

  /// Copy constructor
  Tree(const Tree & tree) throw();

  /// Assignment operator
  Tree & operator= (const Tree & tree) throw();

  /// Return the dimensionality of the tree
  int dimension()
  { return d_; }

  /// Return the refinement level of the tree
  int refinement() 
  { return r_; }

  /// The number of allocated Nodes in the tree
  int num_nodes ()
  { return num_nodes_; }

  /// Return the number of levels in the Tree
  void num_levels ();

  /// Return the root Node of the tree
  Node * root_node()
  { return root_; }

  /// Find the node containing the given point
  Node * find_node (double x, double y, double z);

  /// Refine the given node
  void refine_node (Node * node);

  /// Delete the given node and all its children from the tree
  void delete_node (Node * node);

  /// Balance the Node with respect to its neighbors
  void balance_node (Node * node);

  /// Return the neighboring Node in the specified direction
  Node * node_neighbor (Node * node);

  /// Return the parent node
  Node * node_parent (Node * node);

  /// Return the specified child Node
  Node * node_child (Node * node);

private: // functions

  
private: // attributes

  /// Dimensionality of the Tree
  int d_;

  /// Refinement level of the Tree
  int r_;

  /// Pointer to the root node
  Node * root_;

  /// Number of nodes in the tree
  int num_nodes_;
  

};

#endif /* MESH_TREE_HPP */

