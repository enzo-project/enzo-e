// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef AMR_TREE2K_HPP
#define AMR_TREE2K_HPP

/// @file     tree2k.hpp
/// @brief    Interface for the Tree2K class
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date      
///
/// Detailed description of file amr_tree2k.hpp


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
 * @file      amr_tree2k.hpp
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

#include "amr_treek.hpp"

class Tree2K : public TreeK {

  /// @class    Foo
  /// @brief    Brief description of class Foo.
  /// @ingroup  Template
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

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

private:

  /// Root of the tree
  Node2K * root_;


};

#endif /* TREE_K_HPP */
