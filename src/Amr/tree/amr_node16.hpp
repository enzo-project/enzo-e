// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef AMR_NODE16_HPP
#define AMR_NODE16_HPP

/// @file
/// @brief     
/// @author    
/// @todo     Remove static for thread safety
/// @date      
///
/// Detailed description of file amr_node16.hpp


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
 * @file      node16.hpp
 * @brief     Node class for 4^2-trees
 * @author    James Bordner (jobordner@ucsd.edu)
 * @date      Tue Oct 27 12:32:07 PDT 2009
 * @note      
 *
 * 
 *
 * $Id$
 *
 *********************************************************************
 */

#include <stdlib.h>

#include "cello.h"

#include "amr_node.hpp"

class Tree16;

class Node16 {

  /// @class    Foo
  /// @brief    Brief description of class Foo.
  /// @ingroup  Template

/** 
 *********************************************************************
 *
 * @class     Node16
 * @brief     Node class for 4^2-trees
 * @ingroup   Amr
 *
 * Node class for 2^2-trees 
 *
 *********************************************************************
 */

public:

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

  /// Create a new leaf node
  Node16( int level_adjust = 0 );

  /// Delete a node and all descedents
  ~Node16();

  /// return the specified child
  Node16 * child (int ix, int iy);

  /// return the specified neighbor
  Node16 * neighbor (face_type face);

  /// make the two nodes neighbors.  friend function since either can be NULL
  friend void make_neighbors 
  (Node16 * node_1, Node16 * node_2, face_type face_1);

  /// get the child's cousin
  Node16 * cousin (face_type face, int ix, int iy);

  /// return the parent
  Node16 * parent ();

  /// Refine if any elements in the array are true and recurse
  /// return the level
  int refine 
    (
     const int * level_array, 
     int nd0,  int nd1,
     int low0, int up0,  
     int low1, int up1,
     int level, 
     int max_level,
     bool is_full = true
     );

  /// Perform a pass of trying to remove level-jumps 
  void balance_pass(bool & refined_tree, bool is_full = true);

  /// Perform a pass of trying to optimize uniformly-refined nodes
  void optimize_pass(bool & refined_tree);

  /// Fill the image region with values
  void fill_image
    (
     float * image,
     int nd0,  int nd1,
     int low0, int up0,  
     int low1, int up1,
     int level,
     int num_levels,
     int line_width
     );

  /// Return whether node has all children
  bool all_children () {
    return 
      ((child_[0][0]) && (child_[1][0]) && (child_[2][0]) && (child_[3][0]) &&
       (child_[0][1]) && (child_[1][1]) && (child_[2][1]) && (child_[3][1]) &&
       (child_[0][2]) && (child_[1][2]) && (child_[2][2]) && (child_[3][2]) &&
       (child_[0][3]) && (child_[1][3]) && (child_[2][3]) && (child_[3][3]));
  };

  /// Return whether node has any children
  bool any_children () { 
    return 
      ((child_[0][0]) || (child_[1][0]) || (child_[2][0]) || (child_[3][0]) ||
       (child_[0][1]) || (child_[1][1]) || (child_[2][1]) || (child_[3][1]) ||
       (child_[0][2]) || (child_[1][2]) || (child_[2][2]) || (child_[3][2]) ||
       (child_[0][3]) || (child_[1][3]) || (child_[2][3]) || (child_[3][3]));
  };

  //-------------------------------------------------------------------
  // STATIC OPERATIONS
  //-------------------------------------------------------------------

  /// Return the number of nodes
  static int num_nodes() { return num_nodes_; };

private:

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------

  /// Create child nodes
  void create_children_();

  /// Connect child nodes
  void update_children_();

  /// Delete children and their descendents
  void delete_children_();

  /// Update neighbors for a child
  void update_child_ (int ix, int iy);

  /// Create a child
  void create_child_ (int ix, int iy);

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// Child nodes in edge_type x edge_type ordering
  Node16 * child_[4][4];

  /// Neighbor nodes in edge_type ordering
  Node16 * neighbor_[4];

  /// Parent node
  Node16 * parent_;

  /// Relative level for coalesced nodes
  int level_adjust_;

  //-------------------------------------------------------------------
  // STATIC ATTRIBUTES
  //-------------------------------------------------------------------

  /// Number of nodes allocated
  static int num_nodes_;

};

#endif /* NODE16_HPP */
