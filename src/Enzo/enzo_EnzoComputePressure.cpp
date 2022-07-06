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

  Field field = block->data()->field();

  ASSERT("EnzoComputePressure::compute()",
         "'pressure' must be defined as a permanent field",
         field.field_id("pressure") >= 0);
  // TODO: possibly check that pressure is cell-centered

  compute(block, (enzo_float*)field.values("pressure", i_hist_));
}

//----------------------------------------------------------------------

void EnzoComputePressure::compute (Block * block, enzo_float *p,
                                   int stale_depth /* 0 */) throw()
{
  if (!block->is_leaf()) return;

  compute_(block, p, stale_depth);
}

//----------------------------------------------------------------------

void EnzoComputePressure::compute_(Block * block,
                                   enzo_float * p,
                                   int stale_depth /* 0 */
#ifdef CONFIG_USE_GRACKLE
                                 , code_units * grackle_units, /*NULL*/
                                   grackle_field_data * grackle_fields /*NULL*/
#endif
                                   )
{
  // make a CelloArray that wraps p
  Field field = block->data()->field();
  int gx,gy,gz;
  field.ghost_depth (field.field_id("density"),&gx,&gy,&gz);
  int nx,ny,nz;
  field.size (&nx,&ny,&nz);
  CelloArray<enzo_float,3> p_arr(p, nz + 2*gz, ny + 2*gy, nx + 2*gx);

  EnzoFieldAdaptor f_adaptor(block, i_hist_);

  bool dual_energy = !enzo::fluid_props()->dual_energy_config().is_disabled();
  bool mhd = field.is_field("bfield_x");

  EnzoComputePressure::compute_pressure(f_adaptor, p_arr, mhd, dual_energy,
                                        gamma_, stale_depth);
}

//----------------------------------------------------------------------

namespace { // local function

  template<class K>
  void exec_loop_(int mz, int my, int mx, int stale_depth, K& kernel){

    const int rank = cello::rank();

    const int ix_start = stale_depth;
    const int ix_stop = mx - stale_depth;

    const int iy_start = (rank > 1) ? stale_depth : 0;
    const int iy_stop  = (rank > 1) ? my - stale_depth : my;

    const int iz_start = (rank > 2) ? stale_depth : 0;
    const int iz_stop  = (rank > 2) ? mz - stale_depth : mz;

    if ((ix_start >= ix_stop) | (iy_start >= iy_stop) | (iz_start >= iz_stop)){
      ERROR("exec_loop_", "stale_depth is too large");
    } else if (stale_depth < 0){
      ERROR("exec_loop_", "stale_depth is negative");
    }

    for (int iz = iz_start; iz < iz_stop; iz++){
      for (int iy = iy_start; iy < iy_stop; iy++){
        for (int ix = ix_start; ix < ix_stop; ix++){
          kernel(iz, iy, ix);
        }
      }
    }
  }

}

//----------------------------------------------------------------------

void EnzoComputePressure::compute_pressure
(const EnzoFieldAdaptor& fadaptor,
 const CelloArray<enzo_float, 3>& p,
 bool mhd,
 bool dual_energy,
 double gamma,
 int stale_depth, /* 0 */
 bool ignore_grackle /*false*/
#ifdef CONFIG_USE_GRACKLE
 , code_units * grackle_units, /*nullptr*/
 grackle_field_data * grackle_fields /*nullptr*/
#endif
 ) throw()
{

  if (enzo::config()->method_grackle_use_grackle & !ignore_grackle){
#ifdef CONFIG_USE_GRACKLE
    // the following assertion is not strictly necessary (the problem will be
    // caught later in EnzoMethodGrackle), but this is more informative...
    ASSERT("EnzoMethodGrackle::calculate_pressure",
           "until PR #106 is merged into grackle, stale_depth must be 0 since "
           "grackle's local_calculate_pressure ignores the grid_start and "
           "grid_stop members of grackle_field_data",
           stale_depth == 0);

    if (!fadaptor.consistent_with_field_strides(p)){
      ERROR("EnzoMethodGrackle::calculate_pressure",
            "When using grackle to compute pressure, the output array must "
            "have identical strides to the fields.");
    }
    const EnzoMethodGrackle* grackle_method = enzo::grackle_method();
    grackle_method->calculate_pressure(fadaptor, p.data(), stale_depth,
                                       grackle_units, grackle_fields);
#else
    ERROR("EnzoComputePressure::compute_()",
          "Attempting to compute pressure with method Grackle " 
          "but Enzo-E has not been compiled with Grackle (set use_grackle = 1) \n");
#endif
  } else {

    using RdOnlyEFltArr = CelloArray<const enzo_float, 3>;

    const RdOnlyEFltArr d = fadaptor.view("density");

    enzo_float gm1 = gamma - 1.0;

    const int mz = d.shape(0);
    const int my = d.shape(1);
    const int mx = d.shape(2);

    ASSERT("EnzoComputePressure::compute_pressure",
           "Shapes must be the same",
           ((mz == p.shape(0)) & (my == p.shape(1)) & (mx == p.shape(2))));

    if (dual_energy) {

      const RdOnlyEFltArr ie = fadaptor.view("internal_energy");

      auto loop_body = [=](int iz, int iy, int ix)
        { p(iz,iy,ix) = gm1 * d(iz,iy,ix) * ie(iz,iy,ix); };
      exec_loop_(mz, my, mx, stale_depth, loop_body);

    } else { // not using dual energy formalism

      const int rank = cello::rank();
      const RdOnlyEFltArr te = fadaptor.view("total_energy");

      // fetch velocity arrays
      const RdOnlyEFltArr vx = fadaptor.view("velocity_x");
      const RdOnlyEFltArr vy = (rank >= 2)
        ? fadaptor.view("velocity_y") : RdOnlyEFltArr();
      const RdOnlyEFltArr vz = (rank >= 3)
        ? fadaptor.view("velocity_z") : RdOnlyEFltArr();

      // fetch bfield arrays
      const RdOnlyEFltArr bx = (mhd)
        ? fadaptor.view("bfield_x") : RdOnlyEFltArr();
      const RdOnlyEFltArr by = (mhd & (rank >= 2))
        ? fadaptor.view("bfield_y") : RdOnlyEFltArr();
      const RdOnlyEFltArr bz = (mhd & (rank >= 3))
        ? fadaptor.view("bfield_z") : RdOnlyEFltArr();

      if (rank == 1) {

        auto loop_body = [=](int iz, int iy, int ix)
          {
            enzo_float ke = 0.5 * (vx(iz,iy,ix) * vx(iz,iy,ix));
            enzo_float me_den = mhd ? 0.5*(bx(iz,iy,ix) * bx(iz,iy,ix)) : 0.;
            p(iz,iy,ix) = gm1 * (d(iz,iy,ix) * (te(iz,iy,ix) - ke) - me_den);
          };
        exec_loop_(mz, my, mx, stale_depth, loop_body);

      } else if (rank == 2) {

        auto loop_body = [=](int iz, int iy, int ix)
          {
            enzo_float ke = 0.5*(vx(iz,iy,ix) * vx(iz,iy,ix) +
                                 vy(iz,iy,ix) * vy(iz,iy,ix));
            enzo_float me_den = mhd ? 0.5*(bx(iz,iy,ix) * bx(iz,iy,ix) +
                                           by(iz,iy,ix) * by(iz,iy,ix)) : 0.;
            p(iz,iy,ix) = gm1 * (d(iz,iy,ix) * (te(iz,iy,ix) - ke) - me_den);
          };
        exec_loop_(mz, my, mx, stale_depth, loop_body);

      } else if (rank == 3) {

        auto loop_body = [=](int iz, int iy, int ix)
          {
            enzo_float ke = 0.5*(vx(iz,iy,ix) * vx(iz,iy,ix) +
                                 vy(iz,iy,ix) * vy(iz,iy,ix) +
                                 vz(iz,iy,ix) * vz(iz,iy,ix));
            enzo_float me_den = mhd ? 0.5*(bx(iz,iy,ix) * bx(iz,iy,ix) +
                                           by(iz,iy,ix) * by(iz,iy,ix) +
                                           bz(iz,iy,ix) * bz(iz,iy,ix)) : 0.;
            p(iz,iy,ix) = gm1 * (d(iz,iy,ix) * (te(iz,iy,ix) - ke) - me_den);
          };
        exec_loop_(mz, my, mx, stale_depth, loop_body);

      }
    }
  }

  // Place any additional pressure computation here:
	//    Note: if adding more here, may need to also change
	//          location of field pointer declarations above
	//          (inside / outside of Grackle ifdef)

 return;
}
