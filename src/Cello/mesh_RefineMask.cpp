// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineMask.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-07-14
/// @brief    Implementation of RefineMask class

#include "mesh.hpp"

//----------------------------------------------------------------------

RefineMask::RefineMask(Parameters * parameters,
		       const std::string parameter_name) throw ()
  : value_(new Value(parameters,parameter_name))
{
}

//----------------------------------------------------------------------

void RefineMask::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change
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

  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
	int i=ix + nx*(iy + ny*iz);
	level = std::max(level,v[i]);
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

