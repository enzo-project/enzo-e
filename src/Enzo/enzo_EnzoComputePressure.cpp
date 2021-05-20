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
    const EnzoMethodGrackle* grackle_method = enzo::grackle_method();
    grackle_method->calculate_pressure(block, p, grackle_units,
				       grackle_fields, i_hist_);
#else
    ERROR("EnzoComputePressure::compute_()",
          "Attempting to compute pressure with method Grackle " 
          "but Enzo-E has not been compiled with Grackle (set use_grackle = 1) \n");
#endif

  } else {

    const int rank = cello::rank();

    const enzo_float * d = (enzo_float*) field.values("density", i_hist_);

    const enzo_float * vx = (enzo_float*) field.values("velocity_x", i_hist_);
    const enzo_float * vy = (enzo_float*) field.values("velocity_y", i_hist_);
    const enzo_float * vz = (enzo_float*) field.values("velocity_z", i_hist_);

    const enzo_float * bx = mhd ?
      (enzo_float*) field.values("bfield_x", i_hist_) : nullptr;
    const enzo_float * by = mhd ?
      (enzo_float*) field.values("bfield_y", i_hist_) : nullptr;
    const enzo_float * bz = mhd ?
      (enzo_float*) field.values("bfield_z", i_hist_) : nullptr;

    const enzo_float * te = (enzo_float*) field.values("total_energy", i_hist_);
    const enzo_float * ie = (enzo_float*) field.values("internal_energy", i_hist_);

    int nx,ny,nz;
    field.size(&nx,&ny,&nz);

    int gx,gy,gz;
    field.ghost_depth (0,&gx,&gy,&gz);
    if (rank < 2) gy = 0;
    if (rank < 3) gz = 0;

    int m = (nx+2*gx) * (ny+2*gy) * (nz+2*gz);
    enzo_float gm1 = gamma_ - 1.0;

    if (enzo::config()->ppm_dual_energy) {

      for (int i=0; i<m; i++) {
        p[i] = gm1 * d[i] * ie[i];
      }

    } else {
      if (rank == 1) {
        for (int i=0; i<m; i++) {
          enzo_float ke = 0.5*vx[i]*vx[i];
          enzo_float me_den = mhd ? 0.5*bx[i]*bx[i] : 0.;
          p[i] = gm1 * (d[i] * (te[i] - ke) - me_den);
        }
      } else if (rank == 2) {
        for (int i=0; i<m; i++) {
          enzo_float ke = 0.5*(vx[i]*vx[i] + vy[i]*vy[i]);
          enzo_float me_den = mhd ? 0.5*(bx[i]*bx[i] + by[i]*by[i]) : 0.;
          p[i] = gm1 * (d[i] * (te[i] - ke) - me_den);
        }
      } else if (rank == 3) {
        for (int i=0; i<m; i++) {
          enzo_float ke = 0.5*(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
          enzo_float me_den = mhd ?
            0.5*(bx[i]*bx[i] + by[i]*by[i] + bz[i]*bz[i]) : 0.;
          p[i] = gm1 * (d[i] * (te[i] - ke) - me_den);
        }
      }
    }
  }

  // Place any additional pressure computation here:
	//    Note: if adding more here, may need to also change
	//          location of field pointer declarations above
	//          (inside / outside of Grackle ifdef)

 return;
}
