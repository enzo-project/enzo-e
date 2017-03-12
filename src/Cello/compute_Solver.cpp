// See LICENSE_CELLO file for license and copyright information

/// @file     compute_Solver.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-12-06
/// @brief    Implementation of the Solver abstract base class

#include "compute.hpp"

//======================================================================

Solver::~Solver() throw()
{
  for (size_t i=0; i<refresh_list_.size(); i++) {
    delete refresh_list_[i];
    refresh_list_[i] = 0;
  }
}

//----------------------------------------------------------------------

int Solver::add_refresh (int ghost_depth, 
			 int min_face_rank, 
			 int neighbor_type, 
			 int sync_type,
			 int sync_id)
{
  int index=refresh_list_.size();
  refresh_list_.resize(index+1);
  refresh_list_[index] = new Refresh 
    (ghost_depth,min_face_rank,neighbor_type,sync_type,sync_id,true);
  id_sync_ = sync_id;
  return index;
}

//----------------------------------------------------------------------

Refresh * Solver::refresh(size_t index) 
{
  return (index < refresh_list_.size()) ? refresh_list_[index] : NULL;
}

//======================================================================

void Solver::monitor_output_
(Block * block,
 int iter, double rr0,
 double rr_min, double rr, double rr_max,
 bool final) throw()
{
  Monitor * monitor = block->simulation()->monitor();

  monitor->print("Solver", "%s %s iter %04d  err %.16g [%g %g]",
		 this->name().c_str(),
		 final ? "final" : "",
		 iter,
		 (double)(rr    / rr0),
		 (double)(rr_min/ rr0),
		 (double)(rr_max/ rr0));
}

//----------------------------------------------------------------------

void Solver::begin_(Block * block)
{
  block->push_solver(index_);
}

//----------------------------------------------------------------------

void Solver::end_(Block * block)
{
  int index = block->pop_solver();

  ASSERT2("Solver::end_()",
	  "Solver mismatch was %d expected %d",
	  index,index_,(index == index_));
}

//----------------------------------------------------------------------

bool Solver::is_active_(Block * block)
{
  const int level = block->level();
  const bool is_leaf = block->is_leaf();
  const bool is_unigrid = (min_level_ == max_level_);
  const bool in_range = (min_level_ <= level && level <= max_level_);

  return (is_unigrid) ? (in_range) : (is_leaf && in_range);
}

// ----------------------------------------------------------------------

int Solver::neighbor_type_() const throw() {
  return (min_level_ == max_level_) ? neighbor_level : neighbor_leaf;
}

// ----------------------------------------------------------------------

int Solver::sync_type_() const throw() {
  return (min_level_ == max_level_) ? sync_face : sync_neighbor;
}
