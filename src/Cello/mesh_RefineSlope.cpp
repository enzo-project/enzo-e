// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineSlope.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-07-14
/// @brief    Implementation of RefineSlope class

#include "mesh.hpp"

//----------------------------------------------------------------------

RefineSlope::RefineSlope(FieldDescr * field_descr,
			 double slope_min_refine,
			 double slope_max_coarsen,
			 std::vector<std::string> field_name_list) throw ()
  : slope_min_refine_ (slope_min_refine),
    slope_max_coarsen_(slope_max_coarsen)
{
  if (field_name_list.size() != 0) {
    field_id_list_.resize(field_name_list.size());
    for (size_t i=0; i<field_id_list_.size(); i++) {
      field_id_list_[i] = field_descr->field_id(field_name_list[i]);
      TRACE2("Added field %d %s",field_id_list_[i],field_name_list[i].c_str());
    }
  } else {
    field_id_list_.resize(field_descr->field_count());
    for (size_t i=0; i<field_id_list_.size(); i++) {
      field_id_list_[i] = i;
      TRACE2("Added field %d %s",i,field_descr->field_name(i).c_str());
    }
  }
  TRACE2("RefineSlope::RefineSlope (%f %f)",
	 slope_min_refine_,slope_max_coarsen_);
}

//----------------------------------------------------------------------

int RefineSlope::apply 
(
 CommBlock * comm_block,
 const FieldDescr * field_descr
 ) throw ()
{

  FieldBlock * field_block = comm_block->block()->field_block();

  bool all_coarsen = true;
  bool any_refine = false;

  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);
  TRACE3("n = %d %d %d",nx,ny,nz);

  int rank = nz > 1 ? 3 : (ny > 1 ? 2 : 1);

  double h[3];
  Block * block = comm_block->block();
  double xm[3],xp[3];
  block->lower(&xm[0],&xm[1],&xm[2]);
  block->upper(&xp[0],&xp[1],&xp[2]);
  field_block->cell_width(xm[0],xp[0],&h[0]);
  field_block->cell_width(xm[1],xp[1],&h[1]);
  field_block->cell_width(xm[2],xp[2],&h[2]);

  for (size_t k=0; k<field_id_list_.size(); k++) {

    int id_field = field_id_list_[k];

    int gx,gy,gz;
    field_descr->ghosts(id_field, &gx,&gy,&gz);

    precision_type precision = field_descr->precision(id_field);

    void * void_array = field_block->field_values(id_field);

    const int d3[3] = {1,nx,nx*ny};
    
    // count number of times slope refine and coarsen conditions are satisified
    switch (precision) {
    case precision_single:
      {
	// TODO: use template
	float * array = (float*)void_array;
	float slope;
	for (int axis=0; axis<rank; axis++) {
	  int d = d3[axis];
	  for (int ix=0; ix<nx; ix++) {
	    for (int iy=0; iy<ny; iy++) {
	      for (int iz=0; iz<nz; iz++) {
		int i = (gx+ix) + nx*((gy+iy) + ny*(gz+iz));
		slope = fabs((array[i+d] - array[i-d]) / h[axis]);
		if (slope > slope_min_refine_)  any_refine  = true;
		if (slope > slope_min_refine_)  
		  TRACE6 ("%d: %d %d %d  %g %g",
			  id_field,ix,iy,iz,array[i+d],array[i-d]);
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

	for (int axis=0; axis<rank; axis++) {
	  int d = d3[axis];
	  for (int ix=0; ix<nx; ix++) {
	    for (int iy=0; iy<ny; iy++) {
	      for (int iz=0; iz<nz; iz++) {
		int i = (gx+ix) + nx*((gy+iy) + ny*(gz+iz));
		slope = fabs((array[i+d] - array[i-d]) / h[axis]);
		if (slope > slope_min_refine_)  any_refine  = true;
		if (slope > slope_max_coarsen_) all_coarsen = false;
		if (slope > slope_min_refine_)  
		  TRACE6 ("%d: %d %d %d  %g %g",
			  id_field,ix,iy,iz,array[i+d],array[i-d]);
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
    TRACE3 ("RefineSlope any_refine[%d] = %d all_coarsen = %d",id_field,
	    any_refine,all_coarsen);
  }
  TRACE2 ("RefineSlope any_refine = %d all_coarsen = %d",
	  any_refine,all_coarsen);

  return 
  any_refine ?  adapt_refine : (all_coarsen ? adapt_coarsen : adapt_same) ;

}


//======================================================================

