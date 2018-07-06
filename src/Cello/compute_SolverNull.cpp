// See LICENSE_CELLO file for license and copyright information

/// @file     compute_SolverNull.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    

#include "compute.hpp"

//----------------------------------------------------------------------

void SolverNull::pup (PUP::er &p)
{
  TRACEPUP;
  
  Solver::pup(p);
}

void SolverNull::apply ( std::shared_ptr<Matrix> A,
			 int ix, int ib, Block * block) throw()
{
  Solver::begin_(block);
  Solver::end_(block);
  
  CkCallback(callback_,
	     CkArrayIndexIndex(block->index()),
	     block->proxy_array()).send();
}
//======================================================================

