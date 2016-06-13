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
			 int max_level,
			 bool include_ghosts,
			 std::string output) throw ()
  : Refine (min_refine, max_coarsen, max_level, include_ghosts, output)
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

int RefineSlope::apply ( Block * block ) throw ()
{

  Field field = block->data()->field();

  bool all_coarsen = true;
  bool any_refine = false;

  int rank = block->rank();

  double h3[3];
  Data * data = block->data();
  double xm[3],xp[3];
  data->lower(&xm[0],&xm[1],&xm[2]);
  data->upper(&xp[0],&xp[1],&xp[2]);
  field.cell_width(xm[0],xp[0],&h3[0]);
  field.cell_width(xm[1],xp[1],&h3[1]);
  field.cell_width(xm[2],xp[2],&h3[2]);

  void * output = initialize_output_(field.field_data());

  for (size_t k=0; k<field_id_list_.size(); k++) {

    int id_field = field_id_list_[k];

    int gx,gy,gz;
    if (include_ghosts_) {
      gx = (rank >= 1) ? 1 : 0;
      gy = (rank >= 2) ? 1 : 0;
      gz = (rank >= 3) ? 1 : 0;
    } else {
      field.ghost_depth(id_field, &gx,&gy,&gz);
    }

    int mx,my,mz;
    field.dimensions(id_field,&mx,&my,&mz);

    precision_type precision = field.precision(id_field);

    void * array = field.values(id_field);
   
    // count number of times slope refine and coarsen conditions are satisified
    switch (precision) {
    case precision_single:
      evaluate_block_((float*) array,
		      (float*) output, 
		      mx,my,mz,gx,gy,gz,
		      &any_refine,&all_coarsen, rank,h3);
      break;
    case precision_double:
      evaluate_block_((double*) array,
		      (double*) output,
		      mx,my,mz,gx,gy,gz,
		      &any_refine,&all_coarsen, rank,h3);
      break;
    case precision_quadruple:
      evaluate_block_((long double*) array,
		      (long double*) output,
		      mx,my,mz,gx,gy,gz,
		      &any_refine,&all_coarsen, rank,h3);
      break;
    default:
      ERROR2("RefineSlope::apply",
	     "Unknown precision %d for field %d",
	     precision,0);
      break;
    }
  }

  int adapt_result =
    any_refine ? adapt_refine : (all_coarsen ? adapt_coarsen : adapt_same);

  // Don't refine if already at maximum level
  adjust_for_level_ (&adapt_result,block->level());

  return adapt_result;

}

//----------------------------------------------------------------------

template <class T>
void RefineSlope::evaluate_block_(T * array, T * output ,
				  int mx, int my, int mz,
				  int gx, int gy, int gz,
				  bool *any_refine,
				  bool * all_coarsen, 
				  int rank, 
				  double * h3 )
{
  T slope;
  const int d3[3] = {1,mx,mx*my};
  for (int axis=0; axis<rank; axis++) {
    int id = d3[axis];
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
	for (int ix=gx; ix<mx-gx; ix++) {
	  int i = ix + mx*(iy + my*iz);
	  slope = fabs( (array[i+id] - array[i-id]) 
		       / (2.0*h3[axis]*array[i]));
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

