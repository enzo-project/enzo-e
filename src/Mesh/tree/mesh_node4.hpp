// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_NODE4_HPP
#define MESH_NODE4_HPP

/// @file     mesh_node4.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Oct 27 12:32:07 PDT 2009  
/// @brief    [\ref Mesh] Interface for the Node4 class

#include <stdlib.h>

#include "cello.hpp"

#include "mesh_node.hpp"

class Tree4;

class Node4 {

  /// @class    Node4
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Node class for 2^2-trees

public: // interface

  /// Create a new leaf node
  Node4( int level_adjust = 0 );

  /// Delete a node and all descedents
  ~Node4();

  /// Copy constructor
  Node4(const Node4 & node4) throw();

  /// Assignment operator
  Node4 & operator= (const Node4 & node4) throw();

  /// Number of nodes in subtree rooted at this node
  int num_nodes();

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

private: // functions

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

private: // attributes

  /// Child nodes in corner_type ordering
  Node4 * child_[4];

  /// Neighbor nodes in face_type ordering
  Node4 * neighbor_[4];

  /// Parent node
  Node4 * parent_;

  /// Relative level for coalesced nodes
  int level_adjust_;

};

#endif /* NODE4_HPP */
