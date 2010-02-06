/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @bug       
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

#ifndef NODE4_HPP
#define NODE4_HPP

/** 
*********************************************************************
*
* @file      node4.hpp
* @brief     Node class for 2^2-trees
* @author    James Bordner (jobordner@ucsd.edu)
* @date      Tue Oct 27 12:32:07 PDT 2009
* @bug       
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

class Tree4;

class Node4 {

/** 
*********************************************************************
*
* @class     Node4
* @brief     Node class for 2^2-trees
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

  // Create a new leaf node
  Node4( int level_adjust = 0 );

  /// Delete a node and all descedents
  ~Node4();

  /// return the specified child
  Node4 * child (corner_type corner);

  /// return the specified neighbor
  Node4 * neighbor (face_type face);

  /// make the two nodes neighbors.  friend function since either can be NULL
  friend 
  void make_neighbors
  (Node4 * node_1, Node4 * node_2, face_type face_1);

  /// get the child's cousin
  Node4 * cousin (face_type face, corner_type corner);

  /// return the parent
  Node4 * parent ();

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
      ((child_[0]) &&
       (child_[1]) &&
       (child_[2]) &&
       (child_[3]));
  };

  /// Return whether node has any children
  bool any_children () { 
    return 
      ((child_[0]) ||
       (child_[1]) ||
       (child_[2]) ||
       (child_[3]));
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
  void update_child_ (corner_type corner);

  /// Create a child
  void create_child_ (corner_type corner);

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// Child nodes in corner_type ordering
  Node4 * child_[4];

  /// Neighbor nodes in face_type ordering
  Node4 * neighbor_[4];

  /// Parent node
  Node4 * parent_;

  /// Relative level for coalesced nodes
  int level_adjust_;

  //-------------------------------------------------------------------
  // STATIC ATTRIBUTES
  //-------------------------------------------------------------------

  /// Number of nodes allocated
  static int num_nodes_;

};

#endif /* NODE4_HPP */
