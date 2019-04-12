// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputePressure.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-27 22:37:41
/// @brief    Implements the EnzoComputePressure class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoComputePressure::EnzoComputePressure (double gamma,
					  bool comoving_coordinates)
  : Compute(),
    gamma_(gamma),
    comoving_coordinates_(comoving_coordinates)
{
}

//----------------------------------------------------------------------

void EnzoComputePressure::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Compute::pup(p);

  p | gamma_;
  p | comoving_coordinates_;

}

//----------------------------------------------------------------------

void EnzoComputePressure::compute ( Block * block) throw()
{

  if (!block->is_leaf()) return;

  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  enzo_float * p = field.is_field("pressure") ?
                    (enzo_float*) field.values("pressure", i_hist_) : NULL;

  if (!p){
      ERROR("EnzoComputePressure::compute()",
        " 'pressure' is not defined as a permanent field");
  }

  compute(block, p);
}

//----------------------------------------------------------------------

void EnzoComputePressure::compute (Block * block, enzo_float * p) throw()
{

  if (!block->is_leaf()) return;

  compute_(block, p);
}

//----------------------------------------------------------------------

void EnzoComputePressure::compute_(Block * block,
                                   enzo_float * p
#ifdef CONFIG_USE_GRACKLE
                                 , code_units * grackle_units, /*NULL*/
                                   grackle_field_data * grackle_fields /*NULL*/
#endif
                                   )
{

  if (!block->is_leaf()) return;

  EnzoBlock * enzo_block = enzo::block(block);

  Field field = enzo_block->data()->field();

#ifdef CONFIG_USE_GRACKLE

  code_units grackle_units_;
  grackle_field_data grackle_fields_;

  // setup grackle units if they are not already provided
  if (!grackle_units){
    grackle_units = &grackle_units_;
    EnzoMethodGrackle::setup_grackle_units(enzo_block, grackle_units, i_hist_);
  }

  // if grackle fields are not provided, define them
  bool delete_grackle_fields = false;
  if (!grackle_fields){
    grackle_fields  = &grackle_fields_;
		// NOTE: Add option here to pass history index to setup to
		//       allow for computation of old baryon fields ....
		//       this way we can facilitate computation of
		//       interpolated (in time) fields
    EnzoMethodGrackle::setup_grackle_fields(enzo_block, grackle_fields, i_hist_);
    delete_grackle_fields = true;
  }

  // Compute pressure in Grackle
  if (calculate_pressure(grackle_units, grackle_fields, p) == ENZO_FAIL){
		ERROR("EnzoComputePressure::compute_()",
	        "Error in call to Grackle's calculate_pressure routine.\n");
  }

  if (delete_grackle_fields){
    EnzoMethodGrackle::delete_grackle_fields(grackle_fields);
  }

#else
  const int rank = cello::rank();

  enzo_float * d = (enzo_float*) field.values("density", i_hist_);

  enzo_float * v3[3] =
    { (enzo_float*) (              field.values("velocity_x", i_hist_)),
      (enzo_float*) ((rank >= 2) ? field.values("velocity_y", i_hist_) : NULL),
      (enzo_float*) ((rank >= 3) ? field.values("velocity_z", i_hist_) : NULL) };

  enzo_float * te = (enzo_float*) field.values("total_energy", i_hist_);

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  int gx,gy,gz;
  field.ghost_depth (0,&gx,&gy,&gz);
  if (rank < 2) gy = 0;
  if (rank < 3) gz = 0;

  int m = (nx+2*gx) * (ny+2*gy) * (nz+2*gz);
  enzo_float gm1 = gamma_ - 1.0;
  for (int i=0; i<m; i++) {
    enzo_float e= te[i];
    e -= 0.5*v3[0][i]*v3[0][i];
    if (rank >= 2) e -= 0.5*v3[1][i]*v3[1][i];
    if (rank >= 3) e -= 0.5*v3[2][i]*v3[2][i];
    p[i] = gm1 * d[i] * e;
  }
#endif

  // Place any additional pressure computation here:
	//    Note: if adding more here, may need to also change
	//          location of field pointer declarations above
	//          (inside / outside of Grackle ifdef)

 return;
}
