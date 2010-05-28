// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_NODE16_HPP
#define MESH_NODE16_HPP

/// @file     mesh_node16.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Oct 27 12:32:07 PDT 2009
/// @brief    Declaration of the Node16 class 

#include <stdlib.h>
#include "cello.hpp"
#include "mesh_node.hpp"

class Tree16;

class Node16 {

  /// @class    Node16
  /// @ingroup  Mesh
  /// @brief    Node class for 4^2-trees

public: // interface

  /// Create a new leaf node
  Node16( int level_adjust = 0 );

  /// Delete a node and all descedents
  ~Node16();

  /// Copy constructor
  Node16(const Node16 & node16) throw();

  /// Assignment operator
  Node16 & operator= (const Node16 & node16) throw();

  /// Number of nodes in subtree rooted at this node
  int num_nodes();

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

  /// Refine if any elements in the array are true and recurse, returning level
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

private: // functions

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

private: // attributes

  /// Child nodes in edge_type x edge_type ordering
  Node16 * child_[4][4];

  /// Neighbor nodes in edge_type ordering
  Node16 * neighbor_[4];

  /// Parent node
  Node16 * parent_;

  /// Relative level for coalesced nodes
  int level_adjust_;

};

#endif /* NODE16_HPP */
