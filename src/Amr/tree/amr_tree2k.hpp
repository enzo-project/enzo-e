#ifndef TREE2K_HPP
#define TREE2K_HPP

//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      amr_tree2k.hpp
 * @brief     
 * @author    James Bordner (jobordner@ucsd.edu)
 * @date      
 * @bug       
 * @note      
 *
 * 
 *
 * $Id$
 *
 *********************************************************************
 */

class Tree2K {

/** 
 *********************************************************************
 *
 * @class     Tree2K
 * @brief     
 * @ingroup   GROUP
 *
 * 
 *
 *********************************************************************
 */

public:

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

  Tree2K(int k);

  ~Tree2K() { delete root_; };

  /// Refine down to array
  void refine
    (const int * level_array, 
     int ndx, int ndy,
     int max_level, 
     bool full_nodes = true
     );

  /// print levels
  void print_levels();

  /// Refine nodes to remove level jumps
  void balance(bool full_nodes = true);

  /// Refine nodes to remove level jumps
  void fill(bool full_nodes = true);

  /// Replace uniformly-refined patch with single node
  void optimize();

  /// Create an image of levels
  float * create_image (int n, int line_width);

  /// Return the number of levels
  int levels() { return levels_; }

private:

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// Refinement factor
  int r_;

  /// Number of levels in the tree
  int levels_;

  /// Root of the tree
  Node2K * root_;


};

#endif /* TREE_K_HPP */
