// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineMask.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-07-14
/// @brief    Implementation of RefineMask class

#include "mesh.hpp"

//----------------------------------------------------------------------

RefineMask::RefineMask(Parameters * parameters,
		       const std::string parameter_name,
		       std::string output) throw ()
  : Refine (0.0, 0.0, output),
    value_(new Value(parameters,parameter_name))
{
}

//----------------------------------------------------------------------

void RefineMask::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change
  TRACEPUP;
  Refine::pup(p);
  p | *value_;
}

//----------------------------------------------------------------------

int RefineMask::apply 
(
 Block * block,
 const FieldDescr * field_descr
 ) throw ()
{
  Data * data = block->data();
  FieldData * field_data = data->field_data();

  int nx,ny,nz;
  field_data->size(&nx,&ny,&nz);

  double * x = new double [nx];
  double * y = new double [ny];
  double * z = new double [nz];
  double t = block->time();

  data->field_cells(x,y,z);

  double * v = new double [nx*ny*nz];

  value_->evaluate(v, t,
		   nx,nx,x,
		   ny,ny,y,
		   nz,nz,z);

  double level_want = 0.0;

  // determine desired level
  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
	int i=ix + nx*(iy + ny*iz);
	level_want = std::max(level_want,v[i]);
      }
    }
  }

  const int level_block = block->level();

  if (output_ != "") {
    void * output = initialize_output_(field_data);
    float  * output_float  = (float*) output;
    double * output_double = (double*)output;

    precision_type precision = field_data->precision
      (field_data->field_id(output_));
    for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
	for (int iz=0; iz<nz; iz++) {
	  int i=ix + nx*(iy + ny*iz);
	  if (precision == precision_single) {
	    if (v[i] <  level_block) output_float[i] = -1;
	    if (v[i] == level_block) output_float[i] =  0;
	    if (v[i] >  level_block) output_float[i] = +1;
	  } else if (precision == precision_double) {
	    if (v[i] <  level_block) output_double[i] = -1;
	    if (v[i] == level_block) output_double[i] =  0;
	    if (v[i] >  level_block) output_double[i] = +1;
	  }
	}
      }
    }
  }

  delete [] v;

  delete [] z;
  delete [] y;
  delete [] x;

  return (level_want > level_block) ? adapt_refine :
         (level_want < level_block) ? adapt_coarsen :
                                      adapt_same;
}

//======================================================================

