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
                        int index_prolong,
                        int index_restrict,
                        int min_level,
                        int max_level) throw()
  : Solver(name,
           field_x,
           field_b,
           monitor_iter,
           restart_cycle,
           solve_type,
           index_prolong,
           index_restrict,
           min_level,
           max_level)
{
  cello::simulation()->refresh_set_name(ir_post_,name);
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

