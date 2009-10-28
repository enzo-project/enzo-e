#ifndef NODE_K_HPP
#define NODE_K_HPP

//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      amr_node-k.hpp
 * @brief     Node class for k^2-trees
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

#include "amr_node_k.hpp"


class Tree_k;

class Node_k {

/** 
 *********************************************************************
 *
 * @class     Node
 * @brief     Node class for K^2-trees
 * @ingroup   Amr
 *
 * Node class for K^2-trees 
 *
 *********************************************************************
 */

public:

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

  /// Create a new leaf node
  Node_K( int k, int level_adjust = 0 );

  /// Delete a node and all descedents
  ~Node_K();

  /// return the specified child
  Node_K * child (int ix, int iy);

  /// return the specified neighbor
  Node_K * neighbor (face_type face);

  /// make the two nodes neighbors.  friend function since either can be NULL
  friend void make_neighbors 
  (Node_K * node_1, Node_K * node_2, face_type face_1);

  /// get the child's cousin
  Node_K * cousin (face_type face, int ix, int iy);

  /// return the parent
  Node_K * parent ();

  /// Refine if any elements in the array are true and recurse
  /// return the level
  int refine 
    (
     const int * level_array, 
     int ndx,  int ndy,
     int lowx, int upx,  
     int lowy, int upy,
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
     int ndx,  int ndy,
     int lowx, int upx,  
     int lowy, int upy,
     int level,
     int num_levels,
     int line_width
     );

  /// Return whether node has all children
  bool all_children () {
    for (int ix=0; ix<k_; ix++) {
      for (int iy=0; iy<k_; iy++) {
	if (! child_[ix][iy]) return false;
      }
    }
    return true;
  };

  /// Return whether node has any children
  bool any_children () { 

    for (int ix=0; ix<k_; ix++) {
      for (int iy=0; iy<k_; iy++) {
	if (child_[ix][iy]) return true;
      }
    }
    return false;

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

  /// Number of cells per node axis, e.g. 2 or 4
  char k_;

  /// Child nodes in edge_type x edge_type ordering
  Node_K ** child_;

  /// Neighbor nodes in edge_type ordering
  Node_K ** neighbor_;

  /// Parent node
  Node_K * parent_;

  /// Relative level for coalesced nodes
  int level_adjust_;

  //-------------------------------------------------------------------
  // STATIC ATTRIBUTES
  //-------------------------------------------------------------------

  /// Number of nodes allocated
  static int num_nodes_;

};

#endif /* NODE_K_HPP */
