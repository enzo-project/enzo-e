// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineSlope.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    
///
/// 

#include "mesh.hpp"

//----------------------------------------------------------------------

RefineSlope::RefineSlope() throw ()
{
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

  precision_enum precision = field_descr->precision(0);

  void * void_array = field_block->field_values(0);

  int num_fields = field_descr->field_count();


  switch (precision) {
  case precision_single:
    {
      float * array = (float*)void_array;
      float * slope = new float [nx*ny*nz];
      const int d3[3] = {1,nx,nx*ny};
      for (int axis=0; axis<3; axis++) {
	int d = d3[axis];
	for (int ix=0; ix<nx; ix++) {
	  for (int iy=0; iy<ny; iy++) {
	    for (int iz=0; iz<nz; iz++) {
	      int i = (gx+ix) + nx*((gy+iy) + ny*(gz+iz));
	      if (array[i]) slope[i] = (array[i+d] - array[i-d]) / array[i];
	    }
	  }
	}
      }
    }
    break;
  case precision_double:
    break;
  default:
    ERROR2("RefineSlope::apply",
	   "Unknown precision %d for field %d",
	   precision,0);
    break;
  }
}


//======================================================================

