// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_ItFace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-06-11
/// @brief    Implementation of the ItFace class

#include "mesh.hpp"

//----------------------------------------------------------------------

ItFace::ItFace(int rank, 
	       int rank_limit,
	       bool periodic[3][2],
	       int n3[3],
	       Index index,
	       const int * ic3,
	       const int * ipf3) throw()
  : if3_(),
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

ItFace::~ItFace() throw() 
{
}

//----------------------------------------------------------------------

bool ItFace::next () throw()
{
  do {
    increment_() ;
  } while (!valid_());

  return (!is_reset());
}

//----------------------------------------------------------------------

void ItFace::face (int if3[3]) const
{
  if3[0] = rank_ >= 1 ? if3_[0] : 0;
  if3[1] = rank_ >= 2 ? if3_[1] : 0;
  if3[2] = rank_ >= 3 ? if3_[2] : 0;
}
  
//----------------------------------------------------------------------

Index ItFace::index() const
{
  return index_.index_neighbor(if3_,n3_);
}

//----------------------------------------------------------------------

void ItFace::reset() throw()
{
  if3_[0] = -2;
  if3_[1] = 0;
  if3_[2] = 0;
}

//----------------------------------------------------------------------

bool ItFace::is_reset() const 
{
  return if3_[0] == -2; 
}

//----------------------------------------------------------------------

void ItFace::increment_()
{
  if (is_reset()) {
    set_first_();
  } else {
    if (rank_ >= 1 && if3_[0] < 1) {
      ++if3_[0];
    }	else {
      if3_[0] = -1;
      if (rank_ >= 2 && if3_[1] < 1) {
	++if3_[1];
      } else {
	if3_[1] = -1;
	if (rank_ >= 3 && if3_[2] < 1) {
	  ++if3_[2];
	} else {
	  reset();
	}
      }
    }
  }
}

//----------------------------------------------------------------------

void ItFace::set_first_() 
{
  if3_[0] = rank_ >= 1 ? -1 : 0;
  if3_[1] = rank_ >= 2 ? -1 : 0;
  if3_[2] = rank_ >= 3 ? -1 : 0;
}

//----------------------------------------------------------------------

bool ItFace::valid_() const
{
  if (is_reset()) return true;
  int rank_face = rank_ - 
    (abs(if3_[0]) + abs(if3_[1]) + abs(if3_[2]));
  bool l_range = (rank_limit_ <= rank_face && rank_face < rank_);
  bool l_face = true;
  bool l_parent = true;
  if (ic3_.size() > 0) {
    if (ipf3_.size() == 0) {
      // Face must be adjacent to child
      if (ic3_.size() >= 1 && rank_ >= 1) {
	if  (if3_[0] == -1 && ic3_[0] != 0) l_face = false;
	if  (if3_[0] ==  1 && ic3_[0] != 1) l_face = false;
      }
      if (ic3_.size() >= 2 && rank_ >= 2) {
	if  (if3_[1] == -1 && ic3_[1] != 0) l_face = false;
	if  (if3_[1] ==  1 && ic3_[1] != 1) l_face = false;
      }
      if (ic3_.size() >= 3 && rank_ >= 3) {
	if  (if3_[2] == -1 && ic3_[2] != 0) l_face = false;
	if  (if3_[2] ==  1 && ic3_[2] != 1) l_face = false;
      }
    } else {

      // Face must be adjacent to same parent's face
      if (ipf3_.size() >= 1) {
	if (ipf3_[0] != 0 && (if3_[0] != ipf3_[0])) l_parent = false; 
	if (ipf3_[0] == 0 && 
	    ((ic3_[0] == 0 && if3_[0] == -1) ||
	     (ic3_[0] == 1 && if3_[0] ==  1))) l_parent = false;
      }
      if (ipf3_.size() >= 2) {

	if (ipf3_[1] != 0 && (if3_[1] != ipf3_[1])) l_parent = false; 
	if (ipf3_[1] == 0 && 
	    ((ic3_[1] == 0 && if3_[1] == -1) ||
	     (ic3_[1] == 1 && if3_[1] ==  1))) l_parent = false;
      }
      if (ipf3_.size() >= 3) {
	if (ipf3_[2] != 0 && (if3_[2] != ipf3_[2]))  l_parent = false; 
	if (ipf3_[2] == 0 && 
	    ((ic3_[2] == 0 && if3_[2] == -1) ||
	     (ic3_[2] == 1 && if3_[2] ==  1))) l_parent = false;
      }
    }
  }

  bool l_periodic = true;
  // Return false if on boundary and not periodic
  if (index_.is_on_boundary(if3_,n3_)) {
    if (if3_[0] == -1 && ! periodicity_[0][0]) l_periodic = false;
    if (if3_[0] == +1 && ! periodicity_[0][1]) l_periodic = false;
    if (if3_[1] == -1 && ! periodicity_[1][0]) l_periodic = false;
    if (if3_[1] == +1 && ! periodicity_[1][1]) l_periodic = false;
    if (if3_[2] == -1 && ! periodicity_[2][0]) l_periodic = false;
    if (if3_[2] == +1 && ! periodicity_[2][1]) l_periodic = false;
  }

  return (l_face && l_range && l_parent && l_periodic);
}
//======================================================================

