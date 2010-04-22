// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef AMR_TREE3K_HPP
#define AMR_TREE3K_HPP

/// @file     amr_tree3k.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-29
/// @todo     Decide on either k or r as refinement factor variable name
/// @brief    Interface for the Tree3K class

#include "error.hpp"

#include "amr_treek.hpp"

class Tree3K : public TreeK {

  /// @class    Tree3K
  /// @ingroup  Amr
  /// @brief    Generalized k^3-tree

public: // interface

  /// Initialize a Tree3K with given refinement factor k
  Tree3K(int k);

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~Tree3K() 
  { delete root_; };

  /// Copy constructor
  Tree3K(const Tree3K & tree3k) throw()
  { INCOMPLETE_MESSAGE("Tree3K::Tree3K",""); }

  /// Assignment operator
  Tree3K & operator= (const Tree3K & tree3k) throw()
  { INCOMPLETE_MESSAGE("Tree3K::operator =",""); 
    return *this; }

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

  /// Create an image of levels
  void geomview (std::string filename);

private: // attributes

  /// Root of the tree
  Node3K * root_;

};

#endif /* TREE_K_HPP */
