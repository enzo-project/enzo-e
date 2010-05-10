// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_TREEK_HPP
#define MESH_TREEK_HPP

/// @file     mesh_treek.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-30
/// @brief    Include file for mesh_tree[23]k.hpp

#include <string>

class TreeK {

  /// @class    TreeK
  /// @ingroup  Mesh
  /// @brief    Base class for Tree[23]K classes

public: // interface

  /// Initialize a TreeK with refinement factor r
  TreeK() 
  {};

  /// Initialize a TreeK with refinement factor r
  TreeK(int r) : r_(r), levels_(0) 
  {};

  /// Delete a TreeK object
  virtual ~TreeK() 
  {};

  /// Copy constructor
  TreeK(const TreeK & treek) throw()
  {};

  /// Assignment operator
  virtual TreeK & operator= (const TreeK & treek) throw()
  { return *this; };

  /// Return the number of nodes in the tree (@@@ should not be
  /// virtual, but is due to Node2K and Node3K not having a common
  /// base class)
  virtual int num_nodes() = 0;

  /// Refine down to array
  virtual  void refine
    (const int * level_array, 
     int ndx, int ndy, int ndz,
     int max_level, 
     bool full_nodes = true
     ) = 0;

  /// Refine nodes to remove level jumps
  virtual void balance(bool full_nodes = true)= 0;

  /// Replace uniformly-refined patch with single node
  virtual void optimize()= 0;
  
  /// Create an image of levels
  virtual float * create_image (int n, int line_width, int axis=0)= 0;

  /// Create a geomview file
  virtual void geomview (std::string filename) = 0;

  /// Return the number of levels
  int levels() { return levels_; }

protected: // attributes

  /// Refinement factor
  int r_;

  /// Number of levels in the tree
  int levels_;

};

#endif /* TREE_K_HPP */
