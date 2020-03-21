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
  monitor_iter_(monitor_iter),
  restart_cycle_(restart_cycle),
  callback_(0),
  index_(0),
  min_level_(min_level),
  max_level_(max_level),
  id_sync_(0),
  solve_type_(solve_type),
  ir_post_(-1)
{
  FieldDescr * field_descr = cello::field_descr();
  ix_ = field_descr->field_id(field_x);
  ib_ = field_descr->field_id(field_b);
  ir_post_ = add_new_refresh_();
  cello::refresh(ir_post_)->set_callback(CkIndex_Block::p_refresh_exit());
}

//----------------------------------------------------------------------

Solver::Solver () throw()
  : PUP::able(),
  name_(""),
  ix_(-1),ib_(-1),
  monitor_iter_(0),
  restart_cycle_(1),
  callback_(0),
  index_(0),
  min_level_(0),
  max_level_(std::numeric_limits<int>::max()),
  id_sync_(0),
  solve_type_(solve_leaf),
  ir_post_(-1)
{
  ir_post_ = add_new_refresh_();
  cello::refresh(ir_post_)->set_callback(CkIndex_Block::p_refresh_exit());
}

//----------------------------------------------------------------------

int Solver::add_new_refresh_ ()
{
  // set Solver::ir_post_

  const int * g3 = cello::config()->field_ghost_depth;
  const int ghost_depth = std::max(g3[0],std::max(g3[1],g3[2]));
  const int min_face_rank = cello::config()->adapt_min_face_rank;

  // Set default refresh object
  Refresh refresh_default
    (ghost_depth,min_face_rank, neighbor_type_(), sync_type_(), 0);

  return cello::simulation()->new_register_refresh(refresh_default);
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
  //
  // if solve_type_ == solve_tree always reuse solution

  return ( solve_type_ == solve_tree ||
	   ( ( cycle > 0 ) &&
	     ( ( restart_cycle_ == 0 ) ||
	       ( cycle % restart_cycle_) != 0 ) ) );
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

  return (in_range || solve_type_ == solve_level);
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
