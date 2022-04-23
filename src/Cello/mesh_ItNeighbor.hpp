// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_ItNeighbor.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-06-11
/// @brief    [\ref Mesh] Declaration of the ItNeighbor class
///

#ifndef MESH_IT_NEIGHBOR_HPP
#define MESH_IT_NEIGHBOR_HPP

class ItNeighbor {

  /// @class    ItNeighbor
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  ItNeighbor
  (
   Block * block,
   int min_face_rank,
   int periodic[3],
   int n3[3],
   Index index,
   int neighbor_type,
   int min_level,
   int root_level);

  /// Destructor
  ~ItNeighbor();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    TRACEPUP;
    const bool up = p.isUnpacking();
    if (up) block_ = new Block;
    p | *block_;
    PUParray(p,of3_,3);
    PUParray(p,ic3_,3);
    PUParray(p,ipf3_,3);
    p | rank_;
    p | min_face_rank_;
    PUParray (p,periodic_,3);
    PUParray (p,n3_,3);
    p | index_;
    p | level_;
    p | neighbor_type_;
    p | min_level_;
    p | root_level_;
  }

  /// Reduce another value
  bool next (int of3[3])
  {
    const bool retval = next_();
    if (retval) face_(of3);
    return retval;
  }

  /// Return the level of the current face
  int face_level () const  throw () 
  { return block_->face_level(of3_); }

  void child(int ic3[3]) const ;

  Index index() const ;

  /// Reset the Iterator to the beginning
  void reset();

  bool is_reset() const;

private: // functions

  void face_(int of3[3]) const ;

  /// go to the next face / child
  bool next_();

  /// go to the next face / child
  void increment_();

  /// go to the next child
  void next_child_();

  /// go to the first face
  void set_first_();

  /// go to the first child
  void set_first_child_();

  /// reset child loop
  void reset_child_();

  /// whether child loop is reset
  bool is_reset_child_() const ;

  /// Whether the current face rank is valid
  bool valid_();

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Block that this iterator is associated with
  Block * block_;

  /// Current face
  int of3_[3];

  /// Adjacency child
  int ic3_[3];

  /// Parent face
  int ipf3_[3];

  /// simulation rank
  int rank_;

  /// face rank limit
  int min_face_rank_;

  /// Whether domain is periodic along each face & axis
  int periodic_[3];

  /// Size of the block array
  int n3_[3];

  /// Index
  Index index_;

  /// Level of this block
  int level_;

  /// Neighbor type (neighbor_leaf or neighbor_tree)
  int neighbor_type_;

  /// Minimum level of the Mesh (may be negative)
  int min_level_;
  
  /// Level of coarse grid when neighbor_type_ == neighbor_leaf
  int root_level_;

};

#endif /* MESH_ITNEIGHBOR_HPP */

