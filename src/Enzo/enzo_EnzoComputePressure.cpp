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
  compute_(block);
}

//----------------------------------------------------------------------

void EnzoComputePressure::compute_(Block * block
#ifdef CONFIG_USE_GRACKLE
                                    , code_units * grackle_units /*NULL*/,
																	   	grackle_field_data * grackle_fields /*NULL*/
#endif
                                   )
{

  if (!block->is_leaf()) return;

  EnzoBlock * enzo_block = enzo::block(block);

  Field field = enzo_block->data()->field();

  enzo_float * p = (enzo_float*) field.values("pressure");

#ifdef CONFIG_USE_GRACKLE

  code_units grackle_units_;
  grackle_field_data grackle_fields_;

  // setup grackle units if they are not already provided
  if (!grackle_units){
    grackle_units = &grackle_units_;
    EnzoMethodGrackle::setup_grackle_units(enzo_block, grackle_units);
  }

  // if grackle fields are not provided, define them
  if (!grackle_fields){
    grackle_fields  = &grackle_fields_;
    EnzoMethodGrackle::setup_grackle_fields(enzo_block, grackle_fields);
  }

  // Compute pressure in Grackle
  if (calculate_pressure(grackle_units, grackle_fields, p) == ENZO_FAIL){
		ERROR("EnzoComputePressure::compute_()",
	        "Error in call to Grackle's calculate_pressure routine.\n");
	}


#else
  const int rank = cello::rank();

	enzo_float * d = (enzo_float*) field.values("density");

  enzo_float * v3[3] =
    { (enzo_float*) (              field.values("velocity_x")),
      (enzo_float*) ((rank >= 2) ? field.values("velocity_y") : NULL),
      (enzo_float*) ((rank >= 3) ? field.values("velocity_z") : NULL) };

  enzo_float * te = (enzo_float*) field.values("total_energy");

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
