// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_TREE4_HPP
#define MESH_TREE4_HPP

/// @file     mesh_tree4.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-09-12
/// @bug      Implementations of Tree4 and Tree16 are essentially identical
/// @brief    Declaration of the Tree4 class

#include "error.hpp"

/// @def      MAX_LEVELS
/// @brief    Maximum number of tree levels
#define MAX_LEVELS 64

class Tree4 {

  /// @class    Tree4
  /// @ingroup  Mesh
  /// @brief    Generalized 2^2-tree (quadtree)

public: // interface

  /// Initialize the Tree4 object
  Tree4();

  /// Deallocate the Tree4 object
  ~Tree4() { delete root_; };

  /// Copy constructor
  Tree4(const Tree4 & tree4) throw()
  { INCOMPLETE_MESSAGE("Tree4::Tree4",""); };

  /// Assignment operator
  Tree4 & operator= (const Tree4 & tree4) throw()
  { INCOMPLETE_MESSAGE("Tree4::operator =",""); 
    return *this; };

  /// Return the number of nodes in the tree
  int num_nodes()
  { return root_->num_nodes(); };

  /// Refine down to array
  void refine
    (const int * level_array, 
     int nd0, int nd1, 
     int max_level, 
     bool full_nodes = true
     );

  /// Refine nodes to remove level jumps
  void balance(bool full_nodes = true);

  /// Replace uniformly-refined patch with single node
  void optimize();

  /// Create an image of levels
 float * create_image (int n, int line_width);

  /// Return the number of levels
  int levels() { return levels_; }

private: // attributes

  /// Number of levels in the tree
  int levels_;

  /// Root node of the tree
  Node4 * root_;

};

#endif /* MESH_TREE4_HPP */
