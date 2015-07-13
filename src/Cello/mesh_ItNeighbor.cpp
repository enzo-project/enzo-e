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
 bool periodic[3][2],
 int n3[3],
 Index index)
  : block_(block),
    rank_(block->rank()),
    min_face_rank_(min_face_rank),
    index_(index),
    level_(index.level())
{
  if (!block->is_leaf()) {
    WARNING1("ItNeighbor::ItNeighbor",
	   "ItNeighbor assumes Blocks are leaves, but Block %s is not a leaf",
	   block->name().c_str());
  }
  reset();
  for (int axis=0; axis<3; axis++) {
    n3_[axis] = n3[axis];
    for (int face=0; face<2; face++) {
      periodic_[axis][face] = periodic[axis][face];
    }
  }
}

//----------------------------------------------------------------------

ItNeighbor::~ItNeighbor() 
{
}

//----------------------------------------------------------------------

bool ItNeighbor::next ()
{
  do {
    next_();
  } while ( ! valid_() );

  return (! is_reset()) ;
}

//----------------------------------------------------------------------

Index ItNeighbor::index() const
{
  int face_level = block_->face_level(of3_);
  Index index_neighbor = index_.index_neighbor(of3_,n3_);
  if (face_level == level_) {
    return index_neighbor;
  } else if (face_level == level_ + 1) {
    return index_neighbor.index_child(ic3_);
  } else if (face_level == level_ - 1) {
    return index_neighbor.index_parent();
  } else {
    ERROR4("ItNeighbor::index()",
	   "index %s level %d: face_level %d and block level %d differ by more than one",
	   block_->name().c_str(),index_.level(),face_level,level_);
    return index_neighbor;
  }
}

//----------------------------------------------------------------------

void ItNeighbor::face(int of3[3]) const
{
  of3[0] = (rank_ > 0) ? of3_[0] : 0;
  of3[1] = (rank_ > 1) ? of3_[1] : 0;
  of3[2] = (rank_ > 2) ? of3_[2] : 0;
} 

//----------------------------------------------------------------------

void ItNeighbor::child(int ic3[3]) const
{
  ic3[0] = (rank_ > 0) ? ic3_[0] : 0;
  ic3[1] = (rank_ > 1) ? ic3_[1] : 0;
  ic3[2] = (rank_ > 2) ? ic3_[2] : 0;
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

void ItNeighbor::next_()
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
      if (rank_ > 0 && of3_[0] < 1) { ++of3_[0]; }
      else {
	of3_[0] = -1;
	if (rank_ > 1 && of3_[1] < 1) { ++of3_[1]; }
	else {
	  of3_[1] = -1;
	  if (rank_ > 2 && of3_[2] < 1) { ++of3_[2]; }
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
    if (rank_ > 0 && ic3_[0] < 1) { ++ic3_[0]; }
    else {
      ic3_[0] = 0;
      if (rank_ > 1 && ic3_[1] < 1) { ++ic3_[1]; }
      else {
	ic3_[1] = 0;
	if (rank_ > 2 && ic3_[2] < 1) { ++ic3_[2]; }
	else reset_child_();
      }
    }
  }
}

//----------------------------------------------------------------------

void ItNeighbor::set_first_() 
{
  of3_[0] = rank_ > 0 ? -1 : 0;
  of3_[1] = rank_ > 1 ? -1 : 0;
  of3_[2] = rank_ > 2 ? -1 : 0;
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
  for (int i=0; i<rank_; i++) face_rank -= std::abs(of3_[i]);

  const bool in_range = 
    (min_face_rank_ <= face_rank && face_rank < rank_);

  if (! in_range) return false;

  // Return false if on boundary and not periodic

  if (index_.is_on_boundary(of3_,n3_)) {

    for (int axis=0; axis<rank_; axis++) {

      const bool is_lower_face     = (of3_[axis] == -1);
      const bool is_upper_face     = (of3_[axis] == +1);
      const bool is_lower_periodic = periodic_[axis][0];
      const bool is_upper_periodic = periodic_[axis][1];

      if (is_lower_face && (! is_lower_periodic) ) return false;
      if (is_upper_face && (! is_upper_periodic) ) return false;

    }    
  }

  // Return false if fine child and not adjacent

  if (face_level() > level_) {
    if (is_reset_child_()) return false;

    const int if3[3] = {-of3_[0],-of3_[1],-of3_[2]};
    bool valid = true;
    if (rank_ > 0 && if3[0] == -1 && ic3_[0] != 0) valid = false;
    if (rank_ > 0 && if3[0] ==  1 && ic3_[0] != 1) valid = false;
    if (rank_ > 1 && if3[1] == -1 && ic3_[1] != 0) valid = false;
    if (rank_ > 1 && if3[1] ==  1 && ic3_[1] != 1) valid = false;
    if (rank_ > 2 && if3[2] == -1 && ic3_[2] != 0) valid = false;
    if (rank_ > 2 && if3[2] ==  1 && ic3_[2] != 1) valid = false;
    if (valid == false) return false;
  }

  // Skip face if not same as parent
  if (face_level() < level_) {
    int ic3[3] = {0,0,0};
    index_.child (level_,&ic3[0],&ic3[1],&ic3[2]);
    if  (rank_ > 0 && 
  	 ((of3_[0] == +1 && ic3[0] == 0) ||
  	  (of3_[0] == -1 && ic3[0] == 1))) return false;
    if  (rank_ > 1 && 
  	 ((of3_[1] == +1 && ic3[1] == 0) ||
  	  (of3_[1] == -1 && ic3[1] == 1))) return false;
    if  (rank_ > 2 && 
  	 ((of3_[2] == +1 && ic3[2] == 0) ||
  	  (of3_[2] == -1 && ic3[2] == 1))) return false;
  }

  return true;
}
//======================================================================

