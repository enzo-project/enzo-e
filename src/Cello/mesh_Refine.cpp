// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Refine.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-08-18
/// @brief Implementation of the Refine base class for mesh refinement
///        criteria

#include "mesh.hpp"
#include "charm_simulation.hpp"


void Refine::pup (PUP::er &p)
{
  TRACEPUP;
  PUP::able::pup(p);
  // NOTE: change this function whenever attributes change
  p | min_refine_;
  p | max_coarsen_;
  p | max_level_;
  p | include_ghosts_;
  p | schedule_;
  p | output_;
}

//----------------------------------------------------------------------

void Refine::set_schedule (Schedule * schedule) throw()
{ 
  if (schedule_) delete schedule_;
  schedule_ = schedule;
}

//----------------------------------------------------------------------

void * Refine::initialize_output_(FieldData * field_data)
{
  void * output = 0;
  const bool do_output = output_ != "";

  if (do_output) {
    
    Field field (cello::field_descr(),field_data);
    
    const int id_output = field.field_id(output_);
    output = field.values(id_output);
    int mx,my,mz;
    field.dimensions(id_output,&mx,&my,&mz);
    const int m = mx*my*mz;
    precision_type precision = field.precision(id_output);
    if (precision == precision_single) {
      for (int i=0; i<m; i++) ((float*)output)[i] = -1;
    }  else if (precision == precision_double) {
      for (int i=0; i<m; i++) ((double*)output)[i] = -1;
    }  else if (precision == precision_quadruple) {
      for (int i=0; i<m; i++) ((long double*)output)[i] = -1;
    }
  }
  return output;
}
