// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverCg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-1-08
/// @brief    Implements the CG Krylov iterative linear solver

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoSolverCg::EnzoSolverCg() throw ()
{
  INCOMPLETE("EnzoSolverCg::EnzoSolverCg");
}

//----------------------------------------------------------------------

EnzoSolverCg::EnzoSolverCg(const EnzoSolverCg & EnzoSolverCg) throw ()
/// @param     EnzoSolverCg  Object being copied
{
  INCOMPLETE("EnzoSolverCg::EnzoSolverCg(EnzoSolverCg)");
}

//----------------------------------------------------------------------

EnzoSolverCg & EnzoSolverCg::operator= (const EnzoSolverCg & EnzoSolverCg) throw ()
/// @param     EnzoSolverCg  Source object of the assignment
/// @return    The target assigned object
{
  INCOMPLETE("EnzoSolverCg::operator=");
  return *this;
}

//----------------------------------------------------------------------

EnzoSolverCg::~EnzoSolverCg() throw ()
{
  INCOMPLETE("EnzoSolverCg::~EnzoSolverCg");
}

//----------------------------------------------------------------------

void EnzoSolverCg::pup (PUP::er &p)
{
  TRACEPUP;
  // NOTE: change this function whenever attributes change
}

//======================================================================


void EnzoSolverCg::apply ( Matrix * A, int ix, int ib, Block * block) throw()
{
  INCOMPLETE("EnzoSolverCg::apply");
}
