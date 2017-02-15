// See LICENSE_CELLO file for license and copyright information

/// @file     compute_Solver.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-12-06
/// @brief    Implementation of the Solver abstract base class

#include "compute.hpp"

//----------------------------------------------------------------------

int Solver::add_refresh (int ghost_depth, 
			 int min_face_rank, 
			 int neighbor_type, 
			 int sync_type,
			 int id)
{
  int index=refresh_list_.size();
  refresh_list_.resize(index+1);
  refresh_list_[index] = new Refresh 
    (ghost_depth,min_face_rank,neighbor_type,sync_type,id,true);
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
