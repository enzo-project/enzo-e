// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_TREE16_HPP
#define MESH_TREE16_HPP

/// @file     mesh_tree16.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-09-18
/// @bug      Implementations of Tree4 and Tree16 are essentially identical
/// @brief    [\ref Mesh] Declaration of the Tree16 class 

/// @def      MAX_LEVELS
/// @brief    Maximum number of tree levels
#define MAX_LEVELS 64

class Tree16 {

  /// @class    Tree16
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Generalized 4^2-tree

public: // interface

  /// Initialize the  Tree16 object
  Tree16();

  /// Destructor
  ~Tree16() 
  { delete root_; };

  /// Copy constructor
  Tree16(const Tree16 & tree16) throw()
  { INCOMPLETE("Tree16::Tree16"); };

  /// Assignment operator
  Tree16 & operator= (const Tree16 & tree16) throw()
  { INCOMPLETE("Tree16::operator ="); 
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

 private:

  /// Number of levels in the tree
  int levels_;

  /// Root node of the tree
  Node16 * root_;

};

#endif /* MESH_TREE16_HPP */
