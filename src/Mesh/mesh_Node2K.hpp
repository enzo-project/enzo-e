// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_NODE2K_HPP
#define MESH_NODE2K_HPP

/// @file     mesh_Node2K.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Oct 27 12:32:07 PDT 2009
/// @brief    Interface for Node2K class

class Tree2K;

class Node2K {

  /// @class    Node2K
  /// @brief    Node class for 2D k^2-trees
  /// @ingroup  Mesh

public: // interface

  /// Create a new leaf node
  Node2K( int k,  int level_adjust = 0  );

  /// Destructor
  ~Node2K();

  /// Copy constructor
  Node2K(const Node2K & node2k) throw();

  /// Assignment operator
  Node2K & operator= (const Node2K & node2k) throw();

  /// Number of nodes in subtree rooted at this node
  int num_nodes();

  /// return the specified child
  Node2K * child (int ix, int iy);

  /// return the specified neighbor
  Node2K * neighbor (face_enum face);

  /// make the two nodes neighbors.  friend function since either can be NULL
  friend void make_neighbors 
  (Node2K * node_1, Node2K * node_2, face_enum face_1);

  /// get the child's cousin
  Node2K * cousin (face_enum face, int ix, int iy);

  /// return the parent
  Node2K * parent ();

  /// set data pointer
  void set_data (void * data)
  { data_ = data; };

  /// get data pointer
  void * data ()
  { return data_; };

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
    for (int iy=0; iy<k_; iy++) {
      for (int ix=0; ix<k_; ix++) {
	if (! child(ix,iy)) return false;
      }
    }
    return true;
  };

  /// Return whether node has any children
  bool any_children () { 
    for (int iy=0; iy<k_; iy++) {
      for (int ix=0; ix<k_; ix++) {
	if (child(ix,iy)) return true;
      }
    }
    return false;

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

  /// Index into child[] for ix,iy, iz
  int index_(int ix, int iy) { return (ix + k_*iy); };

  /// Return index of opposite face
  int opposite_face_ (face_enum face) { return int(face) ^ 1; };

  /// Return number of faces
  int num_faces_() { return 4; };

  /// Return maximum number of children
  int num_children_() { return k_*k_; };

  /// Return the level increment log2(k_)
  int level_increment_() { 
    switch (k_) {
    case 2:  return ( 1 ); break;
    case 4:  return ( 2 ); break;
    case 8:  return ( 3 ); break;
    case 16: return ( 4 ); break;
    default:
      fprintf (stderr,"Invalid k=%d for Node2K\n",k_);
      exit(1);
      break;
    }
  }

  void allocate_neighbors_ ();
  void deallocate_neighbors_ ();
  void allocate_children_ ();
  void deallocate_children_ ();

private: // attributes

  /// Number of cells per node axis, e.g. 2, 4, etc.
  char k_;

  /// Child nodes in edge_enum x edge_enum ordering
  Node2K ** child_;

  /// Neighbor nodes in edge_enum ordering
  Node2K ** neighbor_;

  /// Parent node
  Node2K * parent_;

  /// Relative level for coalesced nodes
  int level_adjust_;

  /// Pointer to data contained in the node
  void * data_;

};

#endif /* MESH_NODE2K_HPP */
