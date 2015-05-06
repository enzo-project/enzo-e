// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_ItNeighbor.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-06-11
/// @brief    [\ref Mesh] Declaration of the ItNeighbor class
///

#ifndef MESH_ITNEIGHBOR_HPP
#define MESH_ITNEIGHBOR_HPP

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
   bool periodic[3][2],
   int n3[3],
   Index index);

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
    p | index_;
    p | rank_;
    p | min_face_rank_;
    for (int i=0; i<3; i++) {
      PUParray (p,periodic_[i],2);
    }
    PUParray (p,n3_,3);
    p | index_;
    p | level_;
  }

  /// Reduce another value
  bool next ();

  /// Return the level of the current face
  int face_level () const  throw () 
  { return block_->face_level(of3_); }

  void face(int of3[3]) const ;

  void child(int ic3[3]) const ;

  Index index() const ;

  /// Reset the Iterator to the beginning
  void reset();

  bool is_reset() const;

private: // functions

  /// go to the next face / child
  void next_();

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
  bool periodic_[3][2];

  /// Size of the hierarchy forest
  int n3_[3];

  /// Index
  Index index_;

  /// Level of this block
  int level_;

};

#endif /* MESH_ITNEIGHBOR_HPP */

