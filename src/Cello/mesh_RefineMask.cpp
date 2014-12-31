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
 CommBlock * comm_block,
 const FieldDescr * field_descr
 ) throw ()
{
  Block * block = comm_block->block();
  FieldBlock * field_block = block->field_block();

  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);

  double * x = new double [nx];
  double * y = new double [ny];
  double * z = new double [nz];
  double t = comm_block->time();

  block->field_cells(x,y,z);

  double * v = new double [nx*ny*nz];

  value_->evaluate(v, t,
		   nx,nx,x,
		   ny,ny,y,
		   nz,nz,z);

  double level = 0.0;

  // determine level
  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
	int i=ix + nx*(iy + ny*iz);
	level = std::max(level,v[i]);
      }
    }
  }

  if (output_ != "") {
    void * output = initialize_output_(field_block);
    float  * output_float  = (float*) output;
    double * output_double = (double*)output;

    precision_type precision = field_block->precision
      (field_block->field_id(output_));
    for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
	for (int iz=0; iz<nz; iz++) {
	  int i=ix + nx*(iy + ny*iz);
	  if (precision == precision_single) {
	    if (v[i] <  comm_block->level()) output_float[i] = -1;
	    if (v[i] == comm_block->level()) output_float[i] =  0;
	    if (v[i] >  comm_block->level()) output_float[i] = +1;
	  } else if (precision == precision_double) {
	    if (v[i] <  comm_block->level()) output_double[i] = -1;
	    if (v[i] == comm_block->level()) output_double[i] =  0;
	    if (v[i] >  comm_block->level()) output_double[i] = +1;
	  }
	}
      }
    }
  }

  delete [] v;

  delete [] z;
  delete [] y;
  delete [] x;

  return (level > comm_block->level()) ? adapt_refine :
    (level < comm_block->level()) ? adapt_coarsen : adapt_same;
}

//======================================================================

