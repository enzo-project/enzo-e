// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialZeldovichPancake.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Sun May 22 2022
/// @brief    Implementation of Zeldovich Pancake initializer

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoInitialZeldovichPancake::EnzoInitialZeldovichPancake
(int cycle, double time, std::string aligned_ax_name)
  : Initial(cycle, time), aligned_ax_(0)
{
  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  ASSERT("EnzoInitialZeldovichPancake::enforce_block",
         "Cosmology must be in use", enzo::cosmology() != nullptr);

  ASSERT("EnzoInitialZeldovichPancake::enforce_block",
         "invalid cosmological parameters. Currently, the problem only "
         "supports a flat universe with omega_baryon_now = 1.0.",
         (cosmology->omega_matter_now() == 1.0) &
         (cosmology->omega_baryon_now() == 1.0) &
         (cosmology->omega_cdm_now() == 0.0) &
         (cosmology->omega_lambda_now() == 0.0) );

  if (aligned_ax_name == "x") {
    aligned_ax_ = 0;
  } else if (aligned_ax_name == "y") {
    aligned_ax_ = 1;
  } else if (aligned_ax_name == "y") {
    aligned_ax_ = 2;
  } else {
    ERROR1("EnzoInitialZeldovichPancake::EnzoInitialZeldovichPancake",
           "aligned_ax_name must be x, y, or z. It can't be \"%s\"",
           aligned_ax_name.c_str());
  }

  // I suspect that things will work correctly, but I have not tried yet
  ASSERT("EnzoInitialZeldovichPancake::enforce_block",
         "initializer hasn't been tested with a rank of 1 or 2",
         cello::rank() == 3);
}

//----------------------------------------------------------------------

namespace{ // things within this anonymous namespace are local to this file

struct PressureFreeFluidProps{ double density; double velocity; };

class ZeldovichPancakeSoln {

  /// @class    ZeldovichPancakeSoln
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Functor that computes the solution for the
  ///     Zeldovich Pancake problem (for a pressure-free fluid) and a given
  ///     redshift

public:

  ZeldovichPancakeSoln(double collapse_z, double target_z,
                       double wavelength, double rho_bkg)
  {
    k_ = (2.0 * cello::pi) / wavelength;

    // Section 7.2 of Anninos+ (1994), defines an amplitude, A', such that the
    // caustic forms at collapse_z. This is given by:
    //       A' = -5*(1+collapse_z) / (2 * k).
    //  in this equation, k is the wavenumber of the perturbation
    //
    // Like in the original Enzo implementation, we use an alternative, more
    // convenient, definition for amplitude, of
    //       A' * (2/5) * a * k
    // where a is the scale factor of interest (the value at target_z).
    // (Note: Unlike the version in Enzo, we preserve the minus sign in A'
    // as part of our new amplitude definiton). This is given by:
    my_amp_ = -(1.0 + collapse_z) / (1.0 + target_z);

    // -> We also define the amplitude for the velocity eqn as (2/5) * adot * A'
    my_amp_vel_ =
      -std::sqrt(2.0/3.0) * (1.0 + collapse_z) / ((1.0 + target_z) * k_);

    rho_bkg_ = rho_bkg;
  }

  ZeldovichPancakeSoln(const ZeldovichPancakeSoln&) = default;

  PressureFreeFluidProps operator()(double x_euler) const {
    // the argument is the position relative to center of pancake

    // The analytic functions that describe the evolution of the
    // pressure-free fluid are written in terms of the initial position of
    // the fluid-elements (i.e. when the perturbation is introduced). These
    // "initial positions" are the Lagrangian positions.
    //
    // We will use the a Newton-Raphson method to convert the given position
    // into the Lagrangian position
    //
    // NOTE: the implementation in Enzo reverses Eulerian & Lagrangian names:
    //  * it considers the position on the grid at init_z to be Lagrangian
    //  * it considers the pressure-free fluid element's position at the time
    //    when the perturbation got introduced to be Eulerian
    // That convention is confusing

    double x_lagrange = x_euler;
    double x_lagrange_old = std::numeric_limits<double>::max();

    const double _TOL = 1e-6;
    while (fabs((x_lagrange-x_lagrange_old)/x_lagrange) > _TOL) {
      x_lagrange_old = x_lagrange;
      // unlike the version in Enzo, my_amp includes the negative sign
      x_lagrange +=(x_euler - x_lagrange - my_amp_*std::sin(k_*x_lagrange)/k_)
                  /(1                    + my_amp_*std::cos(k_*x_lagrange));
    }

    // now that we have x_lagrange, compute the density and velocity

    PressureFreeFluidProps out;
    // unlike the version in Enzo, my_amp includes the negative sign
    out.density = rho_bkg_ / (1 + my_amp_ * std::cos(k_ * x_lagrange));
    out.velocity = my_amp_vel_ * std::sin(k_ * x_lagrange);
    return out;
  }

private:

  /// wavenumber
  double k_;
  /// revised amplitude
  double my_amp_;
  /// revised ampliude
  double my_amp_vel_;
  /// background density (usually taken to be the critical density)
  double rho_bkg_;
};

}


//----------------------------------------------------------------------

void EnzoInitialZeldovichPancake::enforce_block
( Block * block, const Hierarchy  * hierarchy ) throw()
{
  // port of initialization routine from Grid_ZeldovichPancakeInitializeGrid.C
  // The initial conditions come from section 7.2 of Anninos+ (1994).
  //
  // The following features were supported by the original Enzo initializer,
  // but are not currently implemented
  //   - initializing Bfields (I imagine this is something like what the way in
  //     which Li+ 2008 & Collins+ 2010 extended the variant of the problem
  //     introduced in Ryu+ 1993)
  //   - initializing a problem where omega_cdm is non-zero

  // todo: make some/most of the following configurable in the future
  const double init_temperature_K = 100.0; // in units of kelvin
  const double pancake_central_offset = 0.0;
  const double collapse_z = 1.0;

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();
  const double omega_baryon_now = cosmology->omega_baryon_now();
  const double init_z = cosmology->initial_redshift();

  ASSERT("EnzoInitialZeldovichPancake::enforce_block",
         "The initial conditions are only meaningful when the initial "
         "is (considerably) larger than the collapse redshift",
         init_z > collapse_z);

  double domain_min_xyz[3], domain_max_xyz[3];
  cello::hierarchy()->lower(domain_min_xyz, domain_min_xyz+1, domain_min_xyz+2);
  cello::hierarchy()->upper(domain_max_xyz, domain_max_xyz+1, domain_max_xyz+2);
  const double lambda =
    domain_max_xyz[aligned_ax_] - domain_min_xyz[aligned_ax_];

  // In the original version of the code (from enzo), some multiplications &
  // divisions by factors of lambda were dropped since it was defined to be 1.
  // I've tried to make things more explicit, but I may have missed some spots
  ASSERT("EnzoInitialZeldovichPancake::enforce_block",
         "This has only been tested in the case where lambda is equal to the "
         "entire domain length, which is equal to 1.0 in code units.",
         lambda == 1.0);

  // Imagine a universe filled with a pressure-free fluid (like dark matter),
  // which has a sinusoidal perturbation. That universe evolves forward in
  // time to some redshift `init_z`. We are starting our simulation at that
  // redshift and we are initializing our (ordinary) gas to follow the
  // density and velocity profile (we're replacing the pressure-free fluid).

  // 1. Define a function to compute the properties of that pressure-free fluid
  //    at a given position on the grid at `init_z`.

  // the background-density used in the the ZeldovichPancake problem is the
  // critical mass density of the universe. Density code-units are defined such
  // that when cosmology->omega_matter_now() == 1, this density is equal to 1.0
  const double rho_bkg = 1.0;
  const ZeldovichPancakeSoln get_fluid_props(collapse_z, init_z,
                                             lambda, rho_bkg);

  // 2. load some other information (about the grid):
  Data* data = block->data();
  Field field = data->field();
  int mx,my,mz;
  int gx,gy,gz;
  field.dimensions (field.field_id("density"),&mx,&my,&mz);
  field.ghost_depth(field.field_id("density"),&gx,&gy,&gz);

  std::vector<double> xc(mx);
  std::vector<double> yc(my);
  std::vector<double> zc(mz);
  data->field_cells (xc.data(), yc.data(), zc.data(), gx, gy, gz);
  const double pancake_center =
    (0.5*(domain_max_xyz[aligned_ax_] + domain_min_xyz[aligned_ax_]) +
     pancake_central_offset);

  // ensure consistency between cello:rank and the defined vel fields:
  enzo_utils::assert_rank_velocity_field_consistency(*field.field_descr());

  const bool idual = enzo::fluid_props()->dual_energy_config().any_enabled();

  CelloView<enzo_float,3> density = field.view<enzo_float>("density");
  CelloView<enzo_float,3> e_tot = field.view<enzo_float>("total_energy");
  CelloView<enzo_float,3> e_int = (idual) ?
    field.view<enzo_float>("internal_energy") : CelloView<enzo_float,3>();

  // load the velocities
  ASSERT("EnzoInitialZeldovichPancake::enforce_block",
         "rank is incompatible with the axis-direction",
         cello::rank() >= aligned_ax_);
  const EnzoPermutedCoordinates coord(aligned_ax_);
  const std::string velocities[3] = {"velocity_x", "velocity_y", "velocity_z"};

  CelloView<enzo_float,3> aligned_vel =
    field.view<enzo_float>(velocities[coord.i_axis()]);
  std::array<CelloView<enzo_float,3>, 2> transverse_vel =
    { field.is_field(velocities[coord.j_axis()])
      ? field.view<enzo_float>(velocities[coord.j_axis()])
      : CelloView<enzo_float,3>(),
      field.is_field(velocities[coord.k_axis()])
      ? field.view<enzo_float>(velocities[coord.k_axis()])
      : CelloView<enzo_float,3>()
    };

  EnzoUnits * units = enzo::units();

  const double gm1 = enzo::fluid_props()->gamma() - 1.0;
  const double bulkv = 0.0; // convert to problem units?
  const double kelvin_per_energy_units = units->kelvin_per_energy_units();

  // 3. Perform initialization:
  for (int iz=gz; iz<mz-gz; iz++) {
    for (int iy=gy; iy<my-gy; iy++) {
      for (int ix=gx; ix<mx-gx; ix++) {

        // the following line is not very efficient, but that's probably fine
        // (since this is only an initializer)
        double pos = (aligned_ax_ == 0) ? xc[ix]
          : (aligned_ax_ == 1) ? yc[iy] : zc[iz];
        PressureFreeFluidProps props = get_fluid_props(pos - pancake_center);

        // We multiply by omega_baryon_now to be explicit. But it's not really
        // necessary since we currently force it to be 1.0
        double cur_density = omega_baryon_now * props.density;
        double cur_vel = props.velocity + bulkv;

        // we are initializing the internal energy such that the entropy is
        // initialy constant
        double cur_eint =
          ((init_temperature_K / kelvin_per_energy_units) *
           std::pow(cur_density / omega_baryon_now, gm1) / gm1);
        double cur_etot = cur_eint + 0.5 * (cur_vel * cur_vel);

        density(iz,iy,ix) = cur_density;
        aligned_vel(iz,iy,ix) = cur_vel;
        e_tot(iz,iy,ix) = cur_etot;
        if (idual) { e_int(iz,iy,ix) = cur_eint; }
      }
    }
  }

  // initialize the other velocity components:
  for (int i = 0; i <2; i++) {
    if (transverse_vel[i].is_null()) { continue; }
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
        for (int ix=gx; ix<mx-gx; ix++) {
          transverse_vel[i](iz,iy,ix) = 0.0;
        }
      }
    }
  }

  // set initial time!
  EnzoInitialCosmology::init_cosmology(block, -1, -1);

  block->initial_done();
}

