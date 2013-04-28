// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineSlope.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-07-14
/// @brief    Implementation of RefineSlope class

#include "mesh.hpp"

//----------------------------------------------------------------------

RefineSlope::RefineSlope(double slope_min_refine,
			 double slope_max_coarsen) throw ()
  : slope_min_refine_ (slope_min_refine),
    slope_max_coarsen_(slope_max_coarsen)
{
  TRACE2("RefineSlope::RefineSlope (%f %f)",
	 slope_min_refine_,slope_max_coarsen_);
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

  // ASSERT4("RefineSlope::apply",
  // 	  "Ghost zone depths for field %d (%d,%d,%d) must be at least 1",
  // 	  0,gx,gy,gz,
  // 	  gx>0 && gy>0 && gz>0);

  precision_type precision = field_descr->precision(0);

  void * void_array = field_block->field_values(0);

  //  int num_fields = field_descr->field_count();

  bool all_coarsen = true;
  bool any_refine = false;

  const int d3[3] = {1,nx,nx*ny};

  // count number of times slope refine and coarsen conditions are satisified
  switch (precision) {
  case precision_single:
    {
      // TODO: use template
      float * array = (float*)void_array;
      float slope;
      for (int axis=0; axis<3; axis++) {
	int d = d3[axis];
	for (int ix=0; ix<nx; ix++) {
	  for (int iy=0; iy<ny; iy++) {
	    for (int iz=0; iz<nz; iz++) {
	      int i = (gx+ix) + nx*((gy+iy) + ny*(gz+iz));
	      slope = (array[i]) ? 
		(array[i+d] - array[i-d]) / array[i] : 0.0;
	      if (slope > slope_min_refine_)  any_refine  = true;
	      if (slope > slope_max_coarsen_) all_coarsen = false;
	    }
	  }
	}
      }
    }
    break;
  case precision_double:
    {
      // TODO: use template
      double * array = (double*)void_array;
      double slope;

      for (int axis=0; axis<3; axis++) {
	int d = d3[axis];
	for (int ix=0; ix<nx; ix++) {
	  for (int iy=0; iy<ny; iy++) {
	    for (int iz=0; iz<nz; iz++) {
	      int i = (gx+ix) + nx*((gy+iy) + ny*(gz+iz));
	      slope = (array[i]) ? 
		(array[i+d] - array[i-d]) / array[i] : 0.0;
	      if (slope > slope_min_refine_)  any_refine  = true;
	      if (slope > slope_max_coarsen_) all_coarsen = false;
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

  return 
    any_refine ?  adapt_refine : (all_coarsen ? adapt_coarsen : adapt_same) ;

}


//======================================================================

