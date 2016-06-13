// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineDensity.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-07-31
/// @brief    Implementation of RefineDensity class

#include "mesh.hpp"

//----------------------------------------------------------------------

RefineDensity::RefineDensity
(
 double min_refine,
 double max_coarsen,
 int max_level,
 bool include_ghosts,
 std::string output) throw ()
  : Refine(min_refine,max_coarsen,max_level,include_ghosts,output)
{
  TRACE("RefineDensity::RefineDensity");
  WARNING ("RefineDensity::RefineDensity()",
	   "Assuming non-Cosmology problem for RefineDensity");

}

//----------------------------------------------------------------------

int RefineDensity::apply ( Block * block ) throw ()
{

  Field field = block->data()->field();

  int id = field.field_id ("density");
  int precision = field.precision(id);

  int mx,my,mz;
  field.dimensions(id,&mx,&my,&mz);
  int gx,gy,gz;
  if (include_ghosts_) {
    gx = gy = gz = 0;
  } else {
    field.ghost_depth(id, &gx,&gy,&gz);
  }
  char * array = field.values(id);

  int adapt_result;

  if (precision == precision_single) {

    adapt_result = apply_ ((const float*)      array,mx,my,mz,gx,gy,gz);

  } else if (precision == precision_double) {

    adapt_result = apply_ ((const double*)     array,mx,my,mz,gx,gy,gz);

  } else if (precision == precision_quadruple) {

    adapt_result = apply_ ((const long double*)array,mx,my,mz,gx,gy,gz);

  } else {
    ERROR1 ("RefineDensity::apply()",
	   "Unrecognized precision %d\n",
	    precision);
    adapt_result = adapt_unknown;
  }

  // Don't refine if already at maximum level
  adjust_for_level_( &adapt_result, block->level() );
    
  return adapt_result;
}

//----------------------------------------------------------------------s    
template <class T>
int RefineDensity::apply_
( const T * array,
  int mx, int my, int mz,
  int gx, int gy, int gz ) const throw ()
{

  bool any_refine  = false;
  bool all_coarsen = true;
  for (int iz=gz; iz<mz-gz; iz++) {
    for (int iy=gy; iy<my-gy; iy++) {
      for (int ix=gx; ix<mx-gx; ix++) {
	int i = ix + mx*(iy + my*iz);
	if (array[i] > min_refine_)  any_refine  = true;
	if (array[i] < max_coarsen_) all_coarsen = false;
      }
    }
  }
  return 
    any_refine ?  adapt_refine :
    (all_coarsen ? adapt_coarsen : adapt_same) ;

}


//======================================================================

