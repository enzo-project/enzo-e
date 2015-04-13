// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_ItNeighbor.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-06-11
/// @brief    Implementation of the ItNeighbor class

#include "mesh.hpp"

//----------------------------------------------------------------------

ItNeighbor::ItNeighbor(int rank, 
	       int rank_limit,
	       bool periodic[3][2],
	       int n3[3],
	       Index index,
	       const int * ic3,
	       const int * ipf3) throw()
  : of3_(),
    ic3_(),
    ipf3_(),
    rank_(rank),
    rank_limit_(rank_limit),
    index_(index)
{
  reset();
  if (ic3) {
    ic3_.resize(3);
    for (int i=0; i<rank; i++) ic3_[i] = ic3[i];
    for (int i=rank; i<3; i++) ic3_[i] = 0;
  }
  if (ipf3) {
    ipf3_.resize(3);
    for (int i=0; i<rank; i++)  ipf3_[i] = ipf3[i];
    for (int i=rank; i<3; i++)  ipf3_[i] = ipf3[i];
  }
  for (int axis=0; axis<3; axis++) {
    n3_[axis] = n3[axis];
    for (int face=0; face<2; face++) {
      periodicity_[axis][face] = periodic[axis][face];
    }
  }
}

//----------------------------------------------------------------------

ItNeighbor::~ItNeighbor() throw() 
{
}

//----------------------------------------------------------------------

bool ItNeighbor::next (int of3[3]) throw()
{
  do {
    inc_face_() ;
  } while (!valid_());

  of3[0] = (rank_ >= 1) ? of3_[0] : 0;
  of3[1] = (rank_ >= 2) ? of3_[1] : 0;
  of3[2] = (rank_ >= 3) ? of3_[2] : 0;
  
  return (!is_reset());
}

//----------------------------------------------------------------------

void ItNeighbor::reset() throw()
{
  of3_[0] = -2;
  of3_[1] = 0;
  of3_[2] = 0;
}

//----------------------------------------------------------------------

bool ItNeighbor::is_reset() const 
{
  return of3_[0] == -2; 
}

//----------------------------------------------------------------------

void ItNeighbor::inc_face_()
{
  if (is_reset()) {
    set_first_();
  } else {
    if (rank_ >= 1 && of3_[0] < 1) {
      ++of3_[0];
    }	else {
      of3_[0] = -1;
      if (rank_ >= 2 && of3_[1] < 1) {
	++of3_[1];
      } else {
	of3_[1] = -1;
	if (rank_ >= 3 && of3_[2] < 1) {
	  ++of3_[2];
	} else {
	  reset();
	}
      }
    }
  }
}

//----------------------------------------------------------------------

void ItNeighbor::set_first_() 
{
  of3_[0] = rank_ >= 1 ? -1 : 0;
  of3_[1] = rank_ >= 2 ? -1 : 0;
  of3_[2] = rank_ >= 3 ? -1 : 0;
}

//----------------------------------------------------------------------

bool ItNeighbor::valid_() const
{
  if (is_reset()) return true;
  int rank_face = rank_ - 
    (abs(of3_[0]) + abs(of3_[1]) + abs(of3_[2]));
  bool l_range = (rank_limit_ <= rank_face && rank_face < rank_);
  bool l_face = true;
  bool l_parent = true;
  if (ic3_.size() > 0) {
    if (ipf3_.size() == 0) {
      // Face must be adjacent to child
      if (ic3_.size() >= 1 && rank_ >= 1) {
	if  (of3_[0] == -1 && ic3_[0] != 0) l_face = false;
	if  (of3_[0] ==  1 && ic3_[0] != 1) l_face = false;
      }
      if (ic3_.size() >= 2 && rank_ >= 2) {
	if  (of3_[1] == -1 && ic3_[1] != 0) l_face = false;
	if  (of3_[1] ==  1 && ic3_[1] != 1) l_face = false;
      }
      if (ic3_.size() >= 3 && rank_ >= 3) {
	if  (of3_[2] == -1 && ic3_[2] != 0) l_face = false;
	if  (of3_[2] ==  1 && ic3_[2] != 1) l_face = false;
      }
    } else {

      // Face must be adjacent to same parent's face
      if (ipf3_.size() >= 1) {
	if (ipf3_[0] != 0 && (of3_[0] != ipf3_[0])) l_parent = false; 
	if (ipf3_[0] == 0 && 
	    ((ic3_[0] == 0 && of3_[0] == -1) ||
	     (ic3_[0] == 1 && of3_[0] ==  1))) l_parent = false;
      }
      if (ipf3_.size() >= 2) {

	if (ipf3_[1] != 0 && (of3_[1] != ipf3_[1])) l_parent = false; 
	if (ipf3_[1] == 0 && 
	    ((ic3_[1] == 0 && of3_[1] == -1) ||
	     (ic3_[1] == 1 && of3_[1] ==  1))) l_parent = false;
      }
      if (ipf3_.size() >= 3) {
	if (ipf3_[2] != 0 && (of3_[2] != ipf3_[2]))  l_parent = false; 
	if (ipf3_[2] == 0 && 
	    ((ic3_[2] == 0 && of3_[2] == -1) ||
	     (ic3_[2] == 1 && of3_[2] ==  1))) l_parent = false;
      }
    }
  }

  bool l_periodic = true;
  // Return false if on boundary and not periodic
  if (index_.is_on_boundary(of3_,n3_)) {
    if (of3_[0] == -1 && ! periodicity_[0][0]) l_periodic = false;
    if (of3_[0] == +1 && ! periodicity_[0][1]) l_periodic = false;
    if (of3_[1] == -1 && ! periodicity_[1][0]) l_periodic = false;
    if (of3_[1] == +1 && ! periodicity_[1][1]) l_periodic = false;
    if (of3_[2] == -1 && ! periodicity_[2][0]) l_periodic = false;
    if (of3_[2] == +1 && ! periodicity_[2][1]) l_periodic = false;
  }

  // TEMPORARY
  // static bool warning_displayed = false;
  // if (! warning_displayed) {  
  //   //    WARNING("ItNeighbor::is_valid()", "Setting l_periodic to true for testing");
  //   warning_displayed = true;
  // }

  return (l_face && l_range && l_parent && l_periodic);
}
//======================================================================

