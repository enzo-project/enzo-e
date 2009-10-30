#ifndef TREE3K_HPP
#define TREE3K_HPP

//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      amr_tree3k.hpp
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

#include "amr_treek.hpp"

class Tree3K : public TreeK {

/** 
 *********************************************************************
 *
 * @class     Tree3K
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

  Tree3K(int k);

  ~Tree3K() { delete root_; };

  /// Refine down to array
  void refine
    (const int * level_array, 
     int ndx, int ndy, int ndz,
     int max_level, 
     bool full_nodes = true
     );

  /// print levels
  void print_levels();

  /// Refine nodes to remove level jumps
  void balance(bool full_nodes = true);

  /// Replace uniformly-refined patch with single node
  void optimize();

  /// Create an image of levels
  float * create_image (int n, int line_width);

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

private:

  /// Root of the tree
  Node3K * root_;

};

#endif /* TREE_K_HPP */
