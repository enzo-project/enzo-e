// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialZeldovichPancake.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Sun May 22 2022
/// @brief    Implementation of Zeldovich Pancake initializer

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

namespace{
struct coordinate_pair{ double x_eulerian; double x_lagrangian; };
}

//----------------------------------------------------------------------

void EnzoInitialZeldovichPancake::enforce_block
( Block * block, const Hierarchy  * hierarchy ) throw()
{
  // this function tries to port the initialization routine from
  //   Grid_ZeldovichPancakeInitializeGrid.C
  // The initial conditions seem to come from Ryu+ (1993) - they are further
  // discussed in Collins+ (2010) and Li+ (2008). Unfortunately, it's not
  // obvious to me what the analytic solution for linear collapse is
  const double _TOL = 1e-6;

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  ASSERT("EnzoInitialZeldovichPancake::enforce_block",
         "Cosmology must be in use", enzo::cosmology() != nullptr);
  // todo: add more flexibility or at least be more explicit about the problem
  ASSERT("EnzoInitialZeldovichPancake::enforce_block",
         "invalid cosmological parameters",
         (cosmology->omega_matter_now() == 1.0) &
         (cosmology->omega_baryon_now() == 1.0) &
         (cosmology->omega_cdm_now() == 0.0) &
         (cosmology->omega_lambda_now() == 0.0) );

  const double init_z = cosmology->initial_redshift();
  const double omega_baryon_now = cosmology->omega_baryon_now();

  // todo: make some/most of the following configurable in the future
  const double init_temperature_K = 100.0; // in units of kelvin
  const double pancake_central_offset = 0.0;
  const double collapse_z = 1.0;

  double x_domain_min, x_domain_max;
  cello::hierarchy()->lower(&x_domain_min, nullptr, nullptr);
  cello::hierarchy()->upper(&x_domain_max, nullptr, nullptr);
  const double lambda = x_domain_max - x_domain_min;

  // 1. precompute some relevant quantities:
  const double dbl_pi = 2.0 * cello::pi;
  const double kx = dbl_pi / lambda;
  const double amplitude = (1.0 + collapse_z) / (1.0 + init_z);
  const double amplitude_vel =
    -std::sqrt(2.0/3.0) * (1.0 + collapse_z) / ((1.0 + init_z) * dbl_pi);

  // 2. Construct a function that returns lagrangian and eularian positions

  Data* data = block->data();
  Field field = data->field();
  int mx,my,mz;
  int gx,gy,gz;
  field.dimensions (field.field_id("density"),&mx,&my,&mz);
  field.ghost_depth(field.field_id("density"),&gx,&gy,&gz);

  std::vector<double> xc(mx);
  std::vector<double> dummy(std::max(my, mz));
  data->field_cells (xc.data(), dummy.data(), dummy.data(), gx, gy, gz);
  const double pancake_center = (0.5*(x_domain_min + x_domain_max) +
                                 pancake_central_offset);

  auto get_pos_pair = [&](int ix)
    {
      // I don't know why we aren't starting with eulerian coordinates and then
      // computing lagrangian coordinates... (we're following what Enzo does)
      double x_lagrange = xc[ix] - pancake_center;

      double x_euler = x_lagrange;
      double x_euler_old = std::numeric_limits<double>::max();
      while (fabs((x_euler-x_euler_old)/x_euler) > _TOL) {
        x_euler_old = x_euler;
        x_euler += (x_lagrange - x_euler + amplitude*std::sin(kx*x_euler)/kx)
                  /(1                    - amplitude*std::cos(kx*x_euler)   );
      }
      coordinate_pair out = {x_euler, x_lagrange};
      return out;
    };

  // todo: use CelloArrays (we aren't using them now to avoid breaking stuff
  // when we merge in another PR)
  enzo_float * density = (enzo_float *) field.values("density");
  enzo_float * e_int = (enzo_float *) field.values("internal_energy");
  const bool idual = (e_int != nullptr);
  enzo_float * e_tot = (enzo_float *) field.values("total_energy");
  enzo_float * vx = (enzo_float *) field.values("velocity_x");
  enzo_float * vy = (enzo_float *) field.values("velocity_y");
  enzo_float * vz = (enzo_float *) field.values("velocity_z");

  EnzoUnits * units = enzo::units();

  const double gm1 = gamma_ - 1.0;
  // convert to problem units
  const double bulkv = 0.0;
  const double kelvin_per_energy_units = units->kelvin_per_energy_units();

  for (int iz=gz; iz<mz-gz; iz++) {
    for (int iy=gy; iy<my-gy; iy++) {
      for (int ix=gx; ix<mx-gx; ix++) {
        coordinate_pair coord_pair = get_pos_pair(ix);
        double x_euler = coord_pair.x_eulerian;

	int i = ix + mx*(iy + my*iz);

        double cur_density =
          omega_baryon_now / (1 - amplitude*std::cos(kx*x_euler));
        double cur_vx = amplitude_vel * sin(kx*x_euler) + bulkv;
        double cur_eint =
          ((init_temperature_K / kelvin_per_energy_units) *
           std::pow(cur_density / omega_baryon_now, gm1) / gm1);
        double cur_etot = cur_eint + 0.5 * cur_vx * cur_vx;

        density[i] = (enzo_float)cur_density;
        vx[i] = cur_vx;
        vy[i] = 0.0;
        vz[i] = 0.0;
        e_tot[i] = cur_etot;
        if (idual) { e_int[i] = cur_eint; }

      }
    }
  }

  // set initial time!
  EnzoInitialCosmology::init_cosmology(block, -1, -1);

  block->initial_done();
}

