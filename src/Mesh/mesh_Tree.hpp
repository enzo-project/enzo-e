// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Tree.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-06
/// @brief    [\ref Mesh] Declaration of the Tree class for r^d-trees
///

#ifndef MESH_TREE_HPP
#define MESH_TREE_HPP

class Tree {

  /// @class    Tree
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Create a r-way tree in d-dimensions, i.e. Nodes have r^d children
  Tree(int d, int r) throw();

  /// Destructor
  ~Tree() throw();

  /// Return the dimensionality of the tree
  int dimension() const
  { return d_; }

  /// Return the refinement level of the tree
  int refinement() const
  { return r_; }

  /// The number of allocated Nodes in the tree
  int num_nodes () const
  { return num_nodes_; }

  /// Return the maximum Node level in the Tree, with root == 0
  int max_level () const
  { return max_level_; }

  /// Return the root Node of the tree
  Node * root_node() const
  { return root_; }

  /// Refine the given node
  void refine_node (Node * node);

  /// Delete the given node and all its children from the tree
  void delete_node (Node * node);

  /// Balance the Node with respect to its neighbors
  void balance_node (Node * node);

  /// Return the neighboring Node in the specified direction
  Node * node_neighbor (Node * node) const;

  /// Return the parent node
  Node * node_parent (Node * node) const;

  /// Return the specified child Node
  Node * node_child (Node * node, int k) const;

private: // functions

  
private: // attributes

  /// Dimensionality of the Tree
  int d_;

  /// Refinement ratio of levels in the Tree
  int r_;

  /// Pointer to the root node
  Node * root_;

  /// Number of nodes in the tree
  int num_nodes_;
  
  /// Maximum levels in the Tree, starting at 0
  int max_level_;
  

};

#endif /* MESH_TREE_HPP */

