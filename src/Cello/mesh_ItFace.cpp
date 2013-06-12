// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_ItFace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-06-11
/// @brief    Implementation of the ItFace class

#include "mesh.hpp"

//----------------------------------------------------------------------

ItFace::ItFace(int rank_simulation, int rank_limit) throw()
  : rank_simulation_(rank_simulation),
    rank_limit_(rank_limit)
{
  reset();
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
  if3[0] = rank_simulation_ >= 1 ? if3_[0] : 0;
  if3[1] = rank_simulation_ >= 2 ? if3_[1] : 0;
  if3[2] = rank_simulation_ >= 3 ? if3_[2] : 0;
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
    if (rank_simulation_ >= 1 && if3_[0] < 1) {
      ++if3_[0];
    }	else {
      if3_[0] = -1;
      if (rank_simulation_ >= 2 && if3_[1] < 1) {
	++if3_[1];
      } else {
	if3_[1] = -1;
	if (rank_simulation_ >= 3 && if3_[2] < 1) {
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
  if3_[0] = rank_simulation_ >= 1 ? -1 : 0;
  if3_[1] = rank_simulation_ >= 2 ? -1 : 0;
  if3_[2] = rank_simulation_ >= 3 ? -1 : 0;
}

//----------------------------------------------------------------------

bool ItFace::valid_() const
{
  if (is_reset()) return true;
  int rank_face = rank_simulation_ - 
    (abs(if3_[0]) + abs(if3_[1]) + abs(if3_[2]));
  return (rank_limit_ <= rank_face && rank_face < rank_simulation_);
}
//======================================================================

