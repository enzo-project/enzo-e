// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef AMR_TREEK_HPP
#define AMR_TREEK_HPP

/// @file
/// @brief     
/// @author    
/// @date      
///
/// Detailed description of file amr_treek.hpp


/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */


/** 
 *********************************************************************
 *
 * @file      amr_treek.hpp
 * @brief     
 * @author    James Bordner (jobordner@ucsd.edu)
 * @date      
 * @note      
 *
 * 
 *
 * $Id$
 *
 *********************************************************************
 */

class TreeK {

  /// @class    Foo
  /// @brief    Brief description of class Foo.
  /// @ingroup  Template
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
  virtual float * create_image (int n, int line_width, int axis=0)= 0;

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
