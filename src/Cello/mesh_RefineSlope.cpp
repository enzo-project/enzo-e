// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineSlope.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-07-14
/// @brief    Implementation of RefineSlope class

#include "mesh.hpp"

//----------------------------------------------------------------------

RefineSlope::RefineSlope(const FieldDescr * field_descr,
			 double slope_min_refine,
			 double slope_max_coarsen,
			 std::vector<std::string> field_name_list) throw ()
  : slope_min_refine_ (slope_min_refine),
    slope_max_coarsen_(slope_max_coarsen),
    debug_(false)
{
  PARALLEL_PRINTF("slope min max %f %f\n",slope_min_refine,slope_max_coarsen);
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
  // TRACE3("n = %d %d %d",nx,ny,nz);

  int rank = nz > 1 ? 3 : (ny > 1 ? 2 : 1);

  double h3[3];
  Block * block = comm_block->block();
  double xm[3],xp[3];
  block->lower(&xm[0],&xm[1],&xm[2]);
  block->upper(&xp[0],&xp[1],&xp[2]);
  field_block->cell_width(xm[0],xp[0],&h3[0]);
  field_block->cell_width(xm[1],xp[1],&h3[1]);
  field_block->cell_width(xm[2],xp[2],&h3[2]);

  debug_ = ((xm[0]==0.375 && xp[0]==0.5   && xm[1]==0.875 && xp[1]==1.0) ||
	    (xm[0]==0.500 && xp[0]==0.625 && xm[1]==0.875 && xp[1]==1.0));


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
      evaluate_block_((float*)void_array, nx,ny,nz,gx,gy,gz,
		      &any_refine,&all_coarsen, rank,h3,d3);
      break;
    case precision_double:
      evaluate_block_((double*)void_array, nx,ny,nz,gx,gy,gz,
		      &any_refine,&all_coarsen, rank,h3,d3);
      break;
    default:
      ERROR2("RefineSlope::apply",
	     "Unknown precision %d for field %d",
	     precision,0);
      break;
    }
    // TRACE3 ("RefineSlope any_refine[%d] = %d all_coarsen = %d",id_field,
    // 	    any_refine,all_coarsen);
  }
  TRACE2 ("REFINE RefineSlope any_refine = %d all_coarsen = %d",
	  any_refine,all_coarsen);

  if (debug_) TRACE7 ("refine-debug %f %f %f  %f %f %f  %d",
		      xm[0],xm[1],xm[2],
		      xp[0],xp[1],xp[2],
		      any_refine);


  return 
  any_refine ?  adapt_refine : (all_coarsen ? adapt_coarsen : adapt_same) ;

}

//----------------------------------------------------------------------

template <class T>
void RefineSlope::evaluate_block_(T * array, 
				  int nx, int ny, int nz,
				  int gx, int gy, int gz,
				  bool *any_refine,
				  bool * all_coarsen, 
				  int rank, 
				  double * h3, const int * d3)
{
  double slope;
  TRACE3("refine-debug nx,ny,nz = %d %d %d",nx,ny,nz);
  for (int axis=0; axis<rank; axis++) {
    int d = d3[axis];
    for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
	for (int iz=0; iz<nz; iz++) {
	  int i = (gx+ix) + (nx+2*gx)*((gy+iy) + (ny+2*gy)*(gz+iz));
	  slope = fabs((array[i+d] - array[i-d]) / (2.0*h3[axis]*array[i]));
	  
	  if (slope > slope_min_refine_)  *any_refine  = true;
	  if (slope > slope_max_coarsen_) *all_coarsen = false;
	  if (debug_) TRACE5 ("refine-debug a%d  [%d %d %d]  %d",
			      axis,ix,iy,iz,*any_refine);
	}
      }
    }
  }
}
//======================================================================

