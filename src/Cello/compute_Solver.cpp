// See LICENSE_CELLO file for license and copyright information

/// @file     compute_Solver.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-12-06
/// @brief    Implementation of the Solver abstract base class

#include "compute.hpp"

// #define TRACE_SOLVER

#define CYCLE 0

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
(Block * block, int iter,
 double rr0,
 double rr_min,
 double rr,
 double rr_max,
 bool final) throw()
{
  Monitor * monitor = cello::monitor();

  monitor->print("Solver", "%s %s iter %04d  err %.16g [%g %g]",
		 this->name().c_str(),
		 final ? "final" : "",
		 iter,
		 (rr0 != 0) ? (rr    / rr0) : 0.0,
		 (rr0 != 0) ? (rr_min/ rr0) : 0.0,
		 (rr0 != 0) ? (rr_max/ rr0) : 0.0);
}

//----------------------------------------------------------------------

bool Solver::reuse_solution_ (int cycle) const throw()
{
  // 0   0 1 1 1 1 1    always restart from previous solution (except cycle 0)
  // 1   0 0 0 0 0 0    never restart from previous (default)
  // 2   0 1 0 1 0 1    restart every 2nd
  // 3   0 1 1 0 1 1    restart every 3rd
  //      ...           restart every nth
  return ( ( cycle > 0 ) &&
	   ( ( restart_cycle_ == 0 ) ||
	     ( cycle % restart_cycle_) != 0 ) );
}

//----------------------------------------------------------------------

void Solver::begin_(Block * block)
{
#ifdef TRACE_SOLVER  
  if (block->cycle() >= CYCLE)
    CkPrintf ("%s TRACE_SOLVER %d Solver::begin_(%s)\n",
	    block->name().c_str(),index_,name_.c_str());
#endif  
	    
  block->push_solver(index_);
}

//----------------------------------------------------------------------

void Solver::end_(Block * block)
{
  int index = block->pop_solver();
#ifdef TRACE_SOLVER  
  if (block->cycle() >= CYCLE)
    CkPrintf ("%s TRACE_SOLVER %d Solver::end_(%s)\n",
	      block->name().c_str(),index,name_.c_str());
#endif  

  ASSERT2("Solver::end_()",
	  "Solver mismatch was %d expected %d",
	  index,index_,(index == index_));

}

//----------------------------------------------------------------------

bool Solver::is_active_(Block * block) const
{
  const int level = block->level();
  const bool in_range = (min_level_ <= level && level <= max_level_);

  return (in_range);
}

//----------------------------------------------------------------------

bool Solver::is_finest_ (Block * block) const
{
  const int level = block->level();
  const bool is_leaf = block->is_leaf();
  return is_unigrid_ ? (level == max_level_) : is_leaf;
}
