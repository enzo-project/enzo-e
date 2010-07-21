// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_TREE2K_HPP
#define MESH_TREE2K_HPP

/// @file     mesh_Tree2K.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-27 
/// @todo     Decide on either k or r as refinement factor variable name
/// @brief    Interface for the Tree2K class

class Tree2K : public TreeK {

  /// @class    Tree2K
  /// @ingroup  Mesh
  /// @brief    Generalized k^2-tree

public: // interface

  /// Initialize a Tree2K with given refinement factor k
  Tree2K(int k);

  /// Destructor
  ~Tree2K() 
  { delete root_; };

  /// Copy constructor
  Tree2K(const Tree2K & tree2k) throw();

  /// Assignment operator
  Tree2K & operator= (const Tree2K & tree2k) throw();

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
  float * create_image (int n, int line_width,int axis=0);

  /// Create a geomview file
  void geomview (std::string filename);

private: // attributes

  /// Root of the tree
  Node2K * root_;
};

#endif /* TREE_K_HPP */
