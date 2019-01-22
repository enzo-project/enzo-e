// See LICENSE_CELLO file for license and copyright information

/// @file     compute_Solver.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-12-06
/// @brief    Implementation of the Solver abstract base class

#include "compute.hpp"

// #define TRACE_SOLVER
// #define  DEBUG_SOLVER_CG

#define CYCLE 0

// NOTE: Update _compute.hpp solve_enum when updating solve_string
const char * solve_string[] = {
  "solve_unknown",
  "solve_leaf",  // Solve on leaf Blocks (default)
  "solve_level", // Solve within a level (e.g. for multigrid smoothers)
  "solve_tree",  // Solve in a root-level octree (e.g. for domain decomposition)
  "solve_block" // Solve in a block (e.g. MG coarse solve on a single block)
};


//======================================================================

Solver::Solver (std::string name,
		std::string field_x,
		std::string field_b,
		int monitor_iter,
		int restart_cycle,
		int solve_type,
		int min_level,
		int max_level) throw()
  : PUP::able(),
  name_(name),
  ix_(-1),ib_(-1),
  refresh_list_(),
  monitor_iter_(monitor_iter),
  restart_cycle_(restart_cycle),
  callback_(0),
  index_(0),
  min_level_(min_level),
  max_level_(max_level),
  id_sync_(0),
  solve_type_(solve_type)
{
  FieldDescr * field_descr = cello::field_descr();
  ix_ = field_descr->field_id(field_x);
  ib_ = field_descr->field_id(field_b);
}

//----------------------------------------------------------------------

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

  monitor->print("Solver", "%s %s %s iter %04d  err %.16g [%g %g]",
		 block->name().c_str(),
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

  CkCallback(callback_,
	     CkArrayIndexIndex(block->index()),
	     block->proxy_array()).send();

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
  switch (solve_enum(solve_type_)) {
  case solve_leaf:
    return block->is_leaf();
    break;
  case solve_level:
#ifdef DEBUG_SOLVER_CG
    CkPrintf ("DEBUG_SOLVER_CG level max_level %d %d\n",
	      block->level() , max_level_);
#endif    
    return (block->level() == max_level_);
    break;
  case solve_tree:
    return block->is_leaf();
    break;
  case solve_block:
    return block->level() == min_level_;
    break;
  default:
    ERROR2("Solver::is_finest_()",
	   "Unexpected solve_type %d in Solver %s",
	   solve_type_,name_.c_str());
    return false;
  }
}
