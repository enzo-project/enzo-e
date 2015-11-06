// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Refine.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-08-18
/// @brief Implementation of the Refine base class for mesh refinement
///        criteria

#include "mesh.hpp"
#include "charm_simulation.hpp"


//----------------------------------------------------------------------

void * Refine::initialize_output_(FieldData * field_data)
{
  void * output = 0;
  const bool do_output = output_ != "";

  if (do_output) {
    
    FieldDescr * field_descr = 
      proxy_simulation.ckLocalBranch()->field_descr();
    Field field (field_descr,field_data);
    const int id_output = field.field_id(output_);
    output = field.values(id_output);
    int nx,ny,nz;
    field.size(&nx,&ny,&nz);
    int gx,gy,gz;
    field.ghost_depth(id_output, &gx,&gy,&gz);
    const int nxd = nx + 2*gx;
    const int nyd = ny + 2*gy;
    const int nzd = nz + 2*gz;
    precision_type precision = field.precision(id_output);
    if (precision == precision_single) {
      for (int i=0; i<nxd*nyd*nzd; i++) ((float*)output)[i] = -1;
    }  else if (precision == precision_double) {
      for (int i=0; i<nxd*nyd*nzd; i++) ((double*)output)[i] = -1;
    }
  }
  return output;
}
