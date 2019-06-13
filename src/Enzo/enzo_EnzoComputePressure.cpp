// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputePressure.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
///           Andrew Emerick (aemerick11@gmail.com)
/// @date     2019-05-07
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

  bool mhd = field.is_field("bfield_x");

  if (enzo::config()->method_grackle_use_grackle){
#ifdef CONFIG_USE_GRACKLE
    code_units grackle_units_;
    grackle_field_data grackle_fields_;

    ASSERT("EnzoComputePressure::compute_",
	   "Not currently equipped to compute MHD with grackle", mhd);

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
    ERROR("EnzoComputePressure::compute_()",
          "Attempting to compute pressure with method Grackle " 
          "but Enzo-E has not been compiled with Grackle (set use_grackle = 1) \n");
#endif

  } else {

    const int rank = cello::rank();

    enzo_float * d = (enzo_float*) field.values("density", i_hist_);

    enzo_float * v3[3] =
      { (enzo_float*) (              field.values("velocity_x", i_hist_)),
        (enzo_float*) ((rank >= 2) ? field.values("velocity_y", i_hist_) : NULL),
        (enzo_float*) ((rank >= 3) ? field.values("velocity_z", i_hist_) : NULL) };

    enzo_float * b3[3] = {NULL, NULL, NULL};
    if (mhd) {
      b3[0]                = (enzo_float*) field.values("bfield_x", i_hist_);
      if (rank >= 2) b3[1] = (enzo_float*) field.values("bfield_y", i_hist_);
      if (rank >= 3) b3[2] = (enzo_float*) field.values("bfield_z", i_hist_);
    }
 

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

      if (mhd)              e -= 0.5*b3[0][i]*b3[0][i]/d[i];
      if (mhd && rank >= 2) e -= 0.5*b3[1][i]*b3[1][i]/d[i];
      if (mhd && rank >= 2) e -= 0.5*b3[2][i]*b3[2][i]/d[i];
      
      p[i] = gm1 * d[i] * e;
    }
  }

  // Place any additional pressure computation here:
	//    Note: if adding more here, may need to also change
	//          location of field pointer declarations above
	//          (inside / outside of Grackle ifdef)

 return;
}
