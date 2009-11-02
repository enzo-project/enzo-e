#ifndef TREEK_HPP
#define TREEK_HPP

//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      amr_treek.hpp
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

class TreeK {

/** 
 *********************************************************************
 *
 * @class     TreeK
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

  TreeK(int r) : r_(r), levels_(0) {};

  virtual ~TreeK() {};

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
  virtual float * create_image (int n, int line_width)= 0;

  /// Create a geomview file
  virtual void geomview (std::string filename) = 0;

  /// Return the number of levels
  int levels() { return levels_; }

protected:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// Refinement factor
  int r_;

  /// Number of levels in the tree
  int levels_;

};

#endif /* TREE_K_HPP */
