// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodTrace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-11-06
/// @brief    Implementation of the Tracer Particle method

#include "problem.hpp"

//----------------------------------------------------------------------

MethodTrace::MethodTrace 
(
 const FieldDescr * field_descr,
 double courant
 ) throw() 
  : Method (courant)
{
  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier);

  refresh(ir)->add_all_fields(field_descr->field_count());
}

void MethodTrace::compute ( Block * block) throw()
{
  printf ("%s:%d\n",__FILE__,__LINE__);
  block->compute_done(); 
}
