// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_ItNeighbor.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-06-11
/// @brief    Implementation of the ItNeighbor class

#include "mesh.hpp"

//----------------------------------------------------------------------

ItNeighbor::ItNeighbor
(
 Block * block,
 int min_face_rank,
 int periodic[3],
 int n3[3],
 Index index,
 int neighbor_type,
 int min_level,
 int root_level)
  : block_(block),
    rank_(cello::rank()),
    min_face_rank_(min_face_rank),
    index_(index),
    level_(index.level()),
    neighbor_type_(neighbor_type),
    min_level_(min_level),
    root_level_(root_level)
{
  if (!block->is_leaf()) {
    WARNING1("ItNeighbor::ItNeighbor",
	   "ItNeighbor assumes Blocks are leaves, but Block %s is not a leaf",
	   block->name().c_str());
  }
  reset();
  for (int axis=0; axis<3; axis++) {
    n3_[axis] = n3[axis];
    periodic_[axis] = periodic[axis];
  }
}

//----------------------------------------------------------------------

ItNeighbor::~ItNeighbor() 
{
}

//----------------------------------------------------------------------

bool ItNeighbor::next_ ()
{
  do {
    increment_();
  } while ( ! valid_() );
  return (! is_reset()) ;
}

//----------------------------------------------------------------------

Index ItNeighbor::index() const
{
  Index index_neighbor = index_.index_neighbor(of3_,n3_);
  int face_level = block_->face_level(of3_);
  if (face_level == level_) {
    return index_neighbor;
  } else if (face_level == level_ + 1) {
    return index_neighbor.index_child(ic3_);
  } else if (face_level == level_ - 1) {
    return index_neighbor.index_parent();
  } else {
    return index_neighbor;
  }
}

//----------------------------------------------------------------------

void ItNeighbor::face_(int of3[3]) const
{
  of3[0] = (rank_ >= 1) ? of3_[0] : 0;
  of3[1] = (rank_ >= 2) ? of3_[1] : 0;
  of3[2] = (rank_ >= 3) ? of3_[2] : 0;
} 

//----------------------------------------------------------------------

void ItNeighbor::child(int ic3[3]) const
{
  ic3[0] = (rank_ >= 1) ? ic3_[0] : 0;
  ic3[1] = (rank_ >= 2) ? ic3_[1] : 0;
  ic3[2] = (rank_ >= 3) ? ic3_[2] : 0;
  if ( ic3_[0]== -2) {
    ic3[0] = 0;
    ic3[1] = 0;
    ic3[2] = 0;
  }
  if (face_level() < level_) {
    index_.child (level_,&ic3[0],&ic3[1],&ic3[2]);
  }
} 

//----------------------------------------------------------------------

void ItNeighbor::reset()
{
  of3_[0] = -2;
  of3_[1] = 0;
  of3_[2] = 0;
  reset_child_();
}

//----------------------------------------------------------------------

void ItNeighbor::reset_child_()
{
  ic3_[0] = -2;
  ic3_[1] = 0;
  ic3_[2] = 0;
}

//----------------------------------------------------------------------

bool ItNeighbor::is_reset() const
{
  return (of3_[0] == -2);
}

//----------------------------------------------------------------------

bool ItNeighbor::is_reset_child_() const 
{
  return (ic3_[0] == -2);
}

//----------------------------------------------------------------------

void ItNeighbor::increment_()
{
  if (is_reset()) {
    set_first_();
  } else {
    if ( face_level() > level_ ) {
      if (is_reset_child_())
	set_first_child_();
      else 
	next_child_();
    } 
    if (is_reset_child_()) {
      if (rank_ >= 1 && of3_[0] < 1) { ++of3_[0]; }
      else {
	of3_[0] = -1;
	if (rank_ >= 2 && of3_[1] < 1) { ++of3_[1]; }
	else {
	  of3_[1] = -1;
	  if (rank_ >= 3 && of3_[2] < 1) { ++of3_[2]; }
	  else reset();
	}
      }
    }
  }
}

//----------------------------------------------------------------------

void ItNeighbor::next_child_()
{
  if (is_reset_child_()) {
    set_first_child_();
  } else {
    if (rank_ >= 1 && ic3_[0] < 1) { ++ic3_[0]; }
    else {
      ic3_[0] = 0;
      if (rank_ >= 2 && ic3_[1] < 1) { ++ic3_[1]; }
      else {
	ic3_[1] = 0;
	if (rank_ >= 3 && ic3_[2] < 1) { ++ic3_[2]; }
	else reset_child_();
      }
    }
  }
}

//----------------------------------------------------------------------

void ItNeighbor::set_first_() 
{
  of3_[0] = (rank_ >= 1) ? -1 : 0;
  of3_[1] = (rank_ >= 2) ? -1 : 0;
  of3_[2] = (rank_ >= 3) ? -1 : 0;
}

//----------------------------------------------------------------------

void ItNeighbor::set_first_child_() 
{
  ic3_[0] = 0;
  ic3_[1] = 0;
  ic3_[2] = 0;
}

//----------------------------------------------------------------------

bool ItNeighbor::valid_()
{
  if (is_reset()) return true;

  // Check that face rank is in range

  int face_rank = rank_;
  for (int axis=0; axis<rank_; axis++) face_rank -= std::abs(of3_[axis]);

  const bool in_range = 
    (min_face_rank_ <= face_rank && face_rank < rank_);

  if (! in_range) return false;

  // Return false if neighbor_tree type and in different root-level tree

  if ((neighbor_type_ == neighbor_tree) &&
      (! index().is_in_same_subtree(index_,min_level_,root_level_))) {
      return false;
  }

  // Return false if on boundary and not periodic

  for (int axis=0; axis<rank_; axis++) {
    for (int face=0; face < 2; face++) {
      const bool is_on_boundary =
	index_.is_on_boundary(axis,of3_[axis],n3_[axis]);
      const bool is_face =
	(face==0 && of3_[axis] == -1) || (face==1 && of3_[axis] == 1);
      const bool is_periodic = periodic_[axis];
      if ( (! is_periodic) && is_on_boundary && is_face ) return false;
    }    
  }

  // Return false if fine child and not adjacent.  Note fine oblique
  // neighbors are not counted as such.
  // 
  //      of3=(1,0)         ic3[]
  //  ----------+      +-----+-----+
  //            |      |     |     |
  //            | ---> | true|false|   valid
  //   self     |      +-----+-----+
  //            | ---> |     |     |         
  //            |      | true|false|   
  //  ----------+      +-----+-----+

  if (face_level() > level_) {

    if (is_reset_child_()) return false;

    const int if3[3] = {-of3_[0],-of3_[1],-of3_[2]};

    bool valid = true;

    for (int axis=0; axis<rank_; axis++) {
      if (if3[axis] == -1 && ic3_[axis] != 0) valid = false;
      if (if3[axis] ==  1 && ic3_[axis] != 1) valid = false;
      if (! valid) return false;
    }

  } else if (face_level() < level_) {

    // Skip coarse oblique neighbors

    int ic3[3] = {0,0,0};
    index_.child (level_,&ic3[0],&ic3[1],&ic3[2]);

    bool valid = true;

    for (int axis=0; axis<rank_; axis++) {
      if  (of3_[axis] == +1 && ic3[axis] == 0) valid = false;
      if  (of3_[axis] == -1 && ic3[axis] == 1) valid = false;
      if (! valid) return false;
    }
  }

  return true;
}
//======================================================================

