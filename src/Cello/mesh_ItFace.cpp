// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_ItFace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-06-11
/// @brief    Implementation of the ItFace class

#include "mesh.hpp"

//----------------------------------------------------------------------

ItFace::ItFace(int rank, int rank_limit, const int * ic3) throw()
  : rank_(rank),
    rank_limit_(rank_limit),
    ic3_()
{
  reset();
  if (ic3) {
    ic3_.resize(3);
    for (int i=0; i<3; i++) {
      ic3_[i] = i < rank ? ic3[i] : 0;
    }
  }
}

//----------------------------------------------------------------------

ItFace::~ItFace() throw() 
{
}

//----------------------------------------------------------------------

bool ItFace::next (int if3[3]) throw()
{
  do {
    increment_() ;
  } while (!valid_());

  if3[0] = rank_ >= 1 ? if3_[0] : 0;
  if3[1] = rank_ >= 2 ? if3_[1] : 0;
  if3[2] = rank_ >= 3 ? if3_[2] : 0;

  return (!is_reset());
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

void ItFace::value(int if3[3]) const throw()
{
  if3[0] = if3_[0];
  if3[0] = if3_[0];
  if3[0] = if3_[0];
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
  bool in_range = (rank_limit_ <= rank_face && rank_face < rank_);
  bool on_face = true;
  if (ic3_.size() >= 1 && rank_ >= 1) {
    if  (if3_[0] == -1 && ic3_[0] != 0) on_face = false;
    if  (if3_[0] ==  1 && ic3_[0] != 1) on_face = false;
  }
  if (ic3_.size() >= 2 && rank_ >= 2) {
    if  (if3_[1] == -1 && ic3_[1] != 0) on_face = false;
    if  (if3_[1] ==  1 && ic3_[1] != 1) on_face = false;
  }
  if (ic3_.size() >= 3 && rank_ >= 3) {
    if  (if3_[2] == -1 && ic3_[2] != 0) on_face = false;
    if  (if3_[2] ==  1 && ic3_[2] != 1) on_face = false;
  }
  
  return (on_face && in_range);
}
//======================================================================

