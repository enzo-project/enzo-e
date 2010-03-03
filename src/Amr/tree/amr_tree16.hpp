// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef AMR_TREE16_HPP
#define AMR_TREE16_HPP

/// @file     amr_tree16.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-09-18
/// @brief    Declaration of the Tree16 class 

/// @def      MAX_LEVELS
/// @brief    Maximum number of tree levels
#define MAX_LEVELS 64

class Tree16 {

  /// @class    Tree16
  /// @ingroup  Amr
  /// @brief    Generalized 4^2-tree

public: // interface

  /// Initialize the  Tree16 object
  Tree16();

  /// Deallocate the Tree16 object
  ~Tree16() { delete root_; };

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

#endif
