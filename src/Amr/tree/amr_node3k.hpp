//
// $Id$
//
// See LICENSE_CELLO file for license and copyright information
//

#ifndef AMR_NODE3K_HPP
#define AMR_NODE3K_HPP

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


/** 
 *********************************************************************
 *
 * @file      amr_node3k.hpp
 * @brief     Node class for 3D k^3-trees
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

#include "amr_nodek.hpp"

class Tree3K;

class Node3K {

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
  Node3K( int k, int level_adjust = 0 );

  /// Delete a node and all descedents
  ~Node3K();

  /// return the specified child
  Node3K * child (int ix, int iy, int iz);

  /// return the specified neighbor
  Node3K * neighbor (face_type face);

  /// make the two nodes neighbors.  friend function since either can be NULL
  friend void make_neighbors 
  (Node3K * node_1, Node3K * node_2, face_type face_1);

  /// get the child's cousin
  Node3K * cousin (face_type face, int ix, int iy, int iz);

  /// return the parent
  Node3K * parent ();

  /// Refine if any elements in the array are true and recurse
  /// return the level
  int refine 
    (
     const int * level_array, 
     int ndx,  int ndy, int ndz,
     int lowx, int upx,
     int lowy, int upy,
     int lowz, int upz,
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
     int ndx,  int ndy, int ndz,
     int lowx, int upx,  
     int lowy, int upy,
     int lowz, int upz,
     int level,
     int num_levels,
     int line_width,
     int axis
     );

  /// Write a geomview file
  void geomview
  (
   FILE * fpr,
   double lowx, double upx,  
   double lowy, double upy,
   double lowz, double upz,
   bool full = true);

  /// Return whether node has all children
  bool all_children () {
    for (int iz=0; iz<k_; iz++) {
      for (int iy=0; iy<k_; iy++) {
	for (int ix=0; ix<k_; ix++) {
	  if (! child(ix,iy,iz)) return false;
	}
      }
    }
    return true;
  };

  /// Return whether node has any children
  bool any_children () { 

    if (!child_) return false;

    for (int iz=0; iz<k_; iz++) {
      for (int iy=0; iy<k_; iy++) {
	for (int ix=0; ix<k_; ix++) {
	  if (child(ix,iy,iz)) return true;
	}
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
  void update_child_ (int ix, int iy, int iz);

  /// Create a child
  void create_child_ (int ix, int iy, int iz);

  /// Index into child[] for ix,iy, iz
  int index_(int ix, int iy, int iz) { return ix + k_*(iy + k_*iz); };

  /// Return index of opposite face
  int opposite_face_ (face_type face) { return int(face) ^ 1; };

  /// Return number of faces
  int num_faces_() { return 6; };

  /// Return maximum number of children
  int num_children_() { return k_*k_*k_; };

  /// Return the level increment log2(k_)
  int level_increment_() { 
    switch (k_) {
    case 2:  return ( 1 ); break;
    case 4:  return ( 2 ); break;
    case 8:  return ( 3 ); break;
    case 16: return ( 4 ); break;
    default:
      fprintf (stderr,"Invalid k=%d for Node3K\n",k_);
      exit(1);
      break;
    }
  }

  void allocate_neighbors_ ();
  void deallocate_neighbors_ ();
  void allocate_children_ ();
  void deallocate_children_ ();

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// Number of cells per node axis, e.g. 2, 4, etc.
  char k_;

  /// Child nodes in edge_type x edge_type ordering
  Node3K ** child_;

  /// Neighbor nodes in edge_type ordering
  Node3K ** neighbor_;

  /// Parent node
  Node3K * parent_;

  /// Relative level for coalesced nodes
  int level_adjust_;

  //-------------------------------------------------------------------
  // STATIC ATTRIBUTES
  //-------------------------------------------------------------------

  /// Number of nodes allocated
  static int num_nodes_;

};

#endif /* NODE3K_HPP */
