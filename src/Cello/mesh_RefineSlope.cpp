// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineSlope.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-07-14
/// @brief    Implementation of RefineSlope class

#include "mesh.hpp"

//----------------------------------------------------------------------

RefineSlope::RefineSlope(double slope_min) throw ()
  : slope_min_(slope_min)
{
  TRACE("RefineSlope::RefineSlope");
}

//----------------------------------------------------------------------

int RefineSlope::apply 
(
 FieldBlock * field_block,
 const FieldDescr * field_descr
 ) throw ()
{
  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);

  int gx,gy,gz;
  field_descr->ghosts(0, &gx,&gy,&gz);

  ASSERT4("RefineSlope::apply",
	  "Ghost zone depths for field %d (%d,%d,%d) must be at least 1",
	  0,gx,gy,gz,
	  gx>0 && gy>0 && gz>0);

  precision_type precision = field_descr->precision(0);

  void * void_array = field_block->field_values(0);

  //  int num_fields = field_descr->field_count();

  int count_flagged = 0;

  const int d3[3] = {1,nx,nx*ny};

  switch (precision) {
  case precision_single:
    {
      float * array_float = (float*)void_array;
      float slope_float;
      for (int axis=0; axis<3; axis++) {
	int d = d3[axis];
	for (int ix=0; ix<nx; ix++) {
	  for (int iy=0; iy<ny; iy++) {
	    for (int iz=0; iz<nz; iz++) {
	      int i = (gx+ix) + nx*((gy+iy) + ny*(gz+iz));
	      slope_float = (array_float[i]) ? 
		(array_float[i+d] - array_float[i-d]) / array_float[i] : 0.0;
	      if (slope_float > slope_min_) ++count_flagged;
	    }
	  }
	}
      }
    }
    break;
  case precision_double:
    {
      double * array_double = (double*)void_array;
      double slope_double;

      for (int axis=0; axis<3; axis++) {
	int d = d3[axis];
	for (int ix=0; ix<nx; ix++) {
	  for (int iy=0; iy<ny; iy++) {
	    for (int iz=0; iz<nz; iz++) {
	      int i = (gx+ix) + nx*((gy+iy) + ny*(gz+iz));
	      slope_double = (array_double[i]) ? 
		(array_double[i+d] - array_double[i-d]) / array_double[i] : 0.0;
	      if (slope_double > slope_min_) ++count_flagged;
	    }
	  }
	}
      }
    }
    break;
  default:
    ERROR2("RefineSlope::apply",
	   "Unknown precision %d for field %d",
	   precision,0);
    break;
  }

  return count_flagged;

}


//======================================================================

