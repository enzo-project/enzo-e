// See LICENSE_CELLO file for license and copyright information

/// @file     compute_SolverNull.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    

#include "compute.hpp"

SolverNull::SolverNull (std::string name,
			std::string field_x,
			std::string field_b,
			int monitor_iter,
			int restart_cycle,
			int solve_type,
			int min_level,
			int max_level) throw()
  : Solver(name,
	   field_x,
	   field_b,
	   monitor_iter,
	   restart_cycle,
	   solve_type,
	   min_level,
	   max_level)
{
#ifdef NEW_REFRESH  
  cello::simulation()->new_refresh_set_name(ir_post_,name);
#endif  
}
//----------------------------------------------------------------------

void SolverNull::pup (PUP::er &p)
{
  TRACEPUP;
  
  Solver::pup(p);
}

//----------------------------------------------------------------------

void SolverNull::apply ( std::shared_ptr<Matrix> A, Block * block) throw()
{
  Solver::begin_(block);
  Solver::end_(block);
}

//======================================================================

