// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_TREE3K_HPP
#define MESH_TREE3K_HPP

/// @file     mesh_Tree3K.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-29
/// @todo     Decide on either k or r as refinement factor variable name
/// @brief    [\ref Mesh] Interface for the Tree3K class

class Tree3K : public TreeK {

  /// @class    Tree3K
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Generalized k^3-tree

public: // interface

  /// Initialize a Tree3K with given refinement factor k
  Tree3K(int k);

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~Tree3K();

  /// Copy constructor
  Tree3K(const Tree3K & tree3k) throw();

  /// Assignment operator
  Tree3K & operator= (const Tree3K & tree3k) throw();

  /// Return the number of nodes in the tree
  int num_nodes()
  { return root_->num_nodes(); };

  /// Refine down to array
  void refine
    (const int * level_array, 
     int ndx, int ndy, int ndz,
     int max_level, 
     bool full_nodes = true
     );

  /// Refine nodes to remove level jumps
  void balance(bool full_nodes = true);

  /// Replace uniformly-refined patch with single node
  void optimize();

  /// Create an image of levels
  float * create_image (int n, int line_width, int axis);

  /// Create a geomview file of levels
  void geomview (std::string filename);

private: // attributes

  /// Root of the tree
  Node3K * root_;
};

#endif /* TREE_K_HPP */
