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

  /// Return the number of child nodes in each node (d**r)
  int num_children() const
  { return c_; }

  /// The number of allocated Nodes in the tree
  int num_nodes () const
  { return num_nodes_; }

  /// Return the root Node of the tree
  Node * root_node() const
  { return root_; }

  /// Refine the given node
  void refine_node (NodeTrace *);

  /// Delete the given node and all its children from the tree
  void delete_node (NodeTrace *);

  /// Balance the Node with respect to its neighbors
  void balance_node (NodeTrace *);

  /// Return the neighboring Node in the specified direction
  NodeTrace * node_neighbor (NodeTrace *) const;

  /// Return the parent node
  NodeTrace * node_parent (NodeTrace *) const;

  /// Return the specified child Node
  NodeTrace * node_child (NodeTrace *, int k) const;

private: // functions

  
private: // attributes

  /// Dimensionality of the Tree
  int d_;

  /// Refinement ratio of levels in the Tree
  int r_;

  /// Number of children per node
  int c_;

  /// Pointer to the root node
  Node * root_;

  /// Number of nodes in the tree
  int num_nodes_;

};

#endif /* MESH_TREE_HPP */

