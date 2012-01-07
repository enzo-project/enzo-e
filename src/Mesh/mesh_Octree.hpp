// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Octree.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    [\ref Mesh] Declaration of the Octree class
///

#ifndef MESH_OCTREE_HPP
#define MESH_OCTREE_HPP

class Octree {

  /// @class    Octree
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  Octree() throw();

  /// Destructor
  ~Octree() throw();

  /// Copy constructor
  Octree(const Octree & octree) throw();

  /// Assignment operator
  Octree & operator= (const Octree & octree) throw();

  /// Return the dimensionality of the tree
  int dimension()
  { return 2; }

  /// Return the refinement level of the tree
  int refinement() 
  { return 2; }

  /// The number of allocated Nodes in the tree
  int num_nodes ()
  { return num_nodes_; }

  /// Return the number of levels in the Octree
  void num_levels ();

  /// Return the root Node of the tree
  Node * root_node()
  { return root_; }

  /// Find the node containing the given point
  // Node * find_node (double x, double y, double z)
  // { return root_->find_node(x,y,z); }

  /// Create a duplicate copy of the tree starting at the given node.
  // void tree_copy ();


private: // functions

  
private: // attributes

  /// Pointer to the root node
  Node * root_;

  /// Number of nodes in the tree
  int num_nodes_;
  

};

#endif /* MESH_OCTREE_HPP */

