// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialZeldovichPancake.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Sun May 22 2022
/// @brief    Implementation of Zeldovich Pancake initializer

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoInitialZeldovichPancake::EnzoInitialZeldovichPancake
(double gamma, int cycle, double time)
  : Initial(cycle, time), gamma_(gamma)
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
}

//----------------------------------------------------------------------

namespace{

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
  // ports of initialization routine from Grid_ZeldovichPancakeInitializeGrid.C
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

  double x_domain_min, x_domain_max;
  cello::hierarchy()->lower(&x_domain_min, nullptr, nullptr);
  cello::hierarchy()->upper(&x_domain_max, nullptr, nullptr);
  const double lambda = x_domain_max - x_domain_min;

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
  std::vector<double> dummy(std::max(my, mz));
  data->field_cells (xc.data(), dummy.data(), dummy.data(), gx, gy, gz);
  const double pancake_center = (0.5*(x_domain_min + x_domain_max) +
                                 pancake_central_offset);

  const bool idual = enzo::fluid_props()->dual_energy_config().any_enabled();

  CelloView<enzo_float,3> density = field.view<enzo_float>("density");
  CelloView<enzo_float,3> e_tot = field.view<enzo_float>("total_energy");
  CelloView<enzo_float,3> e_int = (idual) ?
    field.view<enzo_float>("internal_energy") : CelloView<enzo_float,3>();
  CelloView<enzo_float,3> vx = field.view<enzo_float>("velocity_x");
  CelloView<enzo_float,3> vy = (cello::rank() >= 2) ?
    field.view<enzo_float>("velocity_y") : CelloView<enzo_float,3>();
  CelloView<enzo_float,3> vz = (cello::rank() >= 3) ?
    field.view<enzo_float>("velocity_z") : CelloView<enzo_float,3>();

  ASSERT("EnzoInitialZeldovichPancake::enforce_block",
         "current implementation only supports 3D sims", cello::rank() == 3);

  EnzoUnits * units = enzo::units();

  const double gm1 = gamma_ - 1.0;
  const double bulkv = 0.0; // convert to problem units?
  const double kelvin_per_energy_units = units->kelvin_per_energy_units();

  // 3. Perform initialization:
  for (int iz=gz; iz<mz-gz; iz++) {
    for (int iy=gy; iy<my-gy; iy++) {
      for (int ix=gx; ix<mx-gx; ix++) {

        PressureFreeFluidProps props = get_fluid_props(xc[ix] - pancake_center);

        // We multiply by omega_baryon_now to be explicit. But it's not really
        // necessary since we currently force it to be 1.0
        double cur_density = omega_baryon_now * props.density;
        double cur_vx = props.velocity + bulkv;

        // we are initializing the internal energy such that the entropy is
        // initialy constant
        double cur_eint =
          ((init_temperature_K / kelvin_per_energy_units) *
           std::pow(cur_density / omega_baryon_now, gm1) / gm1);
        double cur_etot = cur_eint + 0.5 * (cur_vx * cur_vx);

        density(iz,iy,ix) = cur_density;
        vx(iz,iy,ix) = cur_vx;
        vy(iz,iy,ix) = 0.0;
        vz(iz,iy,ix) = 0.0;
        e_tot(iz,iy,ix) = cur_etot;
        if (idual) { e_int(iz,iy,ix) = cur_eint; }
      }
    }
  }

  // set initial time!
  EnzoInitialCosmology::init_cosmology(block, -1, -1);

  block->initial_done();
}

