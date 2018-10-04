// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverDd.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2018-10-01
/// @brief    Implements the EnzoSolverDd class
///
/// @brief [\ref Enzo] Multigrid on the root-level grid using Dd, then
/// BiCgStab in overlapping subdomains defined by root-level Blocks.
/// An optional final Jacobi step can be applied to smooth the solution
/// along subdomain boundaries.

#include "cello.hpp"
#include "enzo.hpp"

//======================================================================

EnzoSolverDd::EnzoSolverDd
  (std::string name,
   std::string field_x,
   std::string field_b,
   int index_solve_coarse,
   int index_solve_domain,
   int index_solve_smooth,
   int monitor_iter,
   int restart_cycle)
    : Solver(name,field_x,field_b,monitor_iter,restart_cycle),
    index_solve_coarse_(index_solve_coarse),
    index_solve_domain_(index_solve_domain),
    index_solve_smooth_(index_solve_smooth),
    mx_(0),my_(0),mz_(0),
    gx_(0),gy_(0),gz_(0)
{
  // Initialize temporary fields

  FieldDescr * field_descr = cello::field_descr();

}

//----------------------------------------------------------------------

void EnzoSolverDd::apply ( std::shared_ptr<Matrix> A, Block * block) throw()
{

  Solver::begin_(block);

  A_ = A;

  allocate_temporary_(block);

  end (block);
}

//----------------------------------------------------------------------

void EnzoSolverDd::end (Block* block) throw ()
{
  deallocate_temporary_(block);
  
  Solver::end_(block);
}

//======================================================================
