// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineSlope.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-07-14
/// @brief    Implementation of RefineSlope class

#include "mesh.hpp"

//----------------------------------------------------------------------

RefineSlope::RefineSlope(const FieldDescr * field_descr,
			 double min_refine,
			 double max_coarsen,
			 std::vector<std::string> field_name_list,
			 std::string output) throw ()
  : Refine (min_refine, max_coarsen, output)
{
  if (field_name_list.size() != 0) {
    field_id_list_.resize(field_name_list.size());
    for (size_t i=0; i<field_id_list_.size(); i++) {
      field_id_list_[i] = field_descr->field_id(field_name_list[i]);
    }
  } else {
    field_id_list_.resize(field_descr->field_count());
    for (size_t i=0; i<field_id_list_.size(); i++) {
      field_id_list_[i] = i;
    }
  }
}

//----------------------------------------------------------------------

int RefineSlope::apply 
(
 Block * block,
 const FieldDescr * field_descr
 ) throw ()
{

  FieldData * field_data = block->data()->field_data();

  bool all_coarsen = true;
  bool any_refine = false;

  int nx,ny,nz;
  field_data->size(&nx,&ny,&nz);

  int rank = nz > 1 ? 3 : (ny > 1 ? 2 : 1);

  double h3[3];
  Data * data = block->data();
  double xm[3],xp[3];
  data->lower(&xm[0],&xm[1],&xm[2]);
  data->upper(&xp[0],&xp[1],&xp[2]);
  field_data->cell_width(xm[0],xp[0],&h3[0]);
  field_data->cell_width(xm[1],xp[1],&h3[1]);
  field_data->cell_width(xm[2],xp[2],&h3[2]);

  void * output = initialize_output_(field_data);

  for (size_t k=0; k<field_id_list_.size(); k++) {

    int id_field = field_id_list_[k];

    int gx,gy,gz;
    field_descr->ghosts(id_field, &gx,&gy,&gz);

    int nxd = nx + 2*gx;
    int nyd = ny + 2*gy;
    int nzd = nz + 2*gz;

    precision_type precision = field_descr->precision(id_field);

    void * array = field_data->values(id_field);
   
    // count number of times slope refine and coarsen conditions are satisified
    switch (precision) {
    case precision_single:
      evaluate_block_((float*) array,
		      (float*) output, 
		      nxd,nyd,nzd,nx,ny,nz,gx,gy,gz,
		      &any_refine,&all_coarsen, rank,h3);
      break;
    case precision_double:
      evaluate_block_((double*) array,
		      (double*) output,
		      nxd,nyd,nzd,nx,ny,nz,gx,gy,gz,
		      &any_refine,&all_coarsen, rank,h3);
      break;
    case precision_quadruple:
      evaluate_block_((long double*) array,
		      (long double*) output,
		      nxd,nyd,nzd,nx,ny,nz,gx,gy,gz,
		      &any_refine,&all_coarsen, rank,h3);
      break;
    default:
      ERROR2("RefineSlope::apply",
	     "Unknown precision %d for field %d",
	     precision,0);
      break;
    }
  }

  int refine_result =  any_refine ?  
    adapt_refine : (all_coarsen ? adapt_coarsen : adapt_same) ;

  return refine_result;

}

//----------------------------------------------------------------------

template <class T>
void RefineSlope::evaluate_block_(T * array, T * output ,
				  int ndx, int ndy, int ndz,
				  int nx, int ny, int nz,
				  int gx, int gy, int gz,
				  bool *any_refine,
				  bool * all_coarsen, 
				  int rank, 
				  double * h3 )
{
  // TEMPORARY: evaluate effect of including (some) ghost zones
  T slope;
  const int d3[3] = {1,ndx,ndx*ndy};
  for (int axis=0; axis<rank; axis++) {
    int d = d3[axis];
    for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
	for (int iz=0; iz<nz; iz++) {
	  int i = (gx+ix) + (ndx)*((gy+iy) + (ndy)*(gz+iz));
	  slope = fabs((array[i+d] - array[i-d]) / (2.0*h3[axis]*array[i]));
	  
	  if (slope > min_refine_)  *any_refine  = true;
	  if (slope > max_coarsen_) *all_coarsen = false;
	  if (output) {
	    if (slope > max_coarsen_) output[i] =  0;
	    if (slope > min_refine_)  output[i] = +1;
	  }
	}
      }
    }
  }
}
//======================================================================

