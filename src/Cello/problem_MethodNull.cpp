// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodNull.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-06-13

#include "problem.hpp"

void MethodNull::compute( Block * block) throw()
{  block->compute_done(); }

//======================================================================

void MethodNull::init_refresh_()
{
#ifdef NEW_REFRESH
  Refresh & refresh = new_refresh(ir_post_);
  cello::simulation()->new_refresh_set_name(ir_post_,name());

#else    
  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier,
			     sync_id_method_null);
  Refresh & refresh = *(this->refresh(ir));

#endif    
  //    refresh(ir)->add_all_fields();
  refresh.add_field("density");
  refresh.add_field("velocity_x");
  refresh.add_field("velocity_y");
  refresh.add_field("velocity_z");
  refresh.add_field("total_energy");
  refresh.add_field("internal_energy");
  refresh.add_field("acceleration_x");
  refresh.add_field("acceleration_y");
  refresh.add_field("acceleration_z");
  
}
