/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodStarMaker.cpp
/// @author     Andrew Emerick (aemerick11@gmail.com)
/// @date
/// @brief  Implements a star maker class
///
///       This is supposed to be a general class for star formation routines
///       where individual SF routines can be created as derived classes.
///       The intention is to try and improve upon the clutter / mess
///       in Enzo's Grid_StarParticleHandler. This will do this, but
///       will still require quite a bit of repeated code across
///       individual SF (the derived classes) routines... so not perfect...

#include "cello.hpp"
#include "enzo.hpp"

//-------------------------------------------------------------------

EnzoMethodStarMaker::EnzoMethodStarMaker
()
  : Method()
{

  const EnzoConfig * enzo_config = enzo::config();

  // AJE: This was the old way this was done
  // Initialize default Refresh object
  // const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier,
  //                           enzo_sync_id_method_star_maker);
  // refresh(ir)->add_all_fields();

  Refresh & refresh = new_refresh(ir_post_);
  cello::simulation()->new_refresh_set_name(ir_post_,name());

  refresh.add_all_fields();

  // Copy over parameters from config to local names here for convenience
  use_density_threshold_     = enzo_config->method_star_maker_use_density_threshold;
  use_velocity_divergence_   = enzo_config->method_star_maker_use_velocity_divergence;
  use_dynamical_time_        = enzo_config->method_star_maker_use_dynamical_time;
  use_self_gravitating_      = enzo_config->method_star_maker_use_self_gravitating;
  use_h2_self_shielding_     = enzo_config->method_star_maker_use_h2_self_shielding;
  use_jeans_mass_            = enzo_config->method_star_maker_use_jeans_mass;
  number_density_threshold_  = enzo_config->method_star_maker_number_density_threshold;
  efficiency_                = enzo_config->method_star_maker_efficiency;
  maximum_star_fraction_     = enzo_config->method_star_maker_maximum_mass_fraction;
  star_particle_min_mass_    = enzo_config->method_star_maker_minimum_star_mass;
  star_particle_max_mass_    = enzo_config->method_star_maker_maximum_star_mass;
}

//-------------------------------------------------------------------

void EnzoMethodStarMaker::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | use_density_threshold_;
  p | use_velocity_divergence_;
  p | use_dynamical_time_;
  p | number_density_threshold_;
  p | efficiency_;
  p | maximum_star_fraction_;
  p | star_particle_min_mass_;
  p | star_particle_max_mass_;
  p | use_self_gravitating_;
  p | use_h2_self_shielding_;
  p | use_jeans_mass_;

  return;
}

//------------------------------------------------------------------
//   This does nothing at the moment - business is done in derived
//   class (Currently EnzoMethodStarMakerStochasticSF)
void EnzoMethodStarMaker::compute ( Block *block) throw()
{

  if (! block->is_leaf()) return;

  block->compute_done();

  return;
}

// Required
double EnzoMethodStarMaker::timestep ( Block *block) const throw()
{
  return std::numeric_limits<double>::max();
}

void EnzoMethodStarMaker::rescale_densities(EnzoBlock * enzo_block,
                                            const int index,
                                            const double density_ratio) throw() {

  // Loop through all passive scalars (colour fields)
  // which are mass fractions stored as densities, and rescale
  // to the new density after star formation.
  //
  // AE: NOTE: Change this routine if / whenever there needs to be
  //           fraction fields that are not labelled as colour
  //           Obviously requires these fields to be declared as colour
  //           in input file to work.
  //           This can / should likely get moved to be a more general
  //           function of the block / data / field class (one of those)
  //
  //    density_ratio = new_density / old_density
  //

  Field field = enzo_block->data()->field();

  Grouping * field_groups = field.groups();
  int nc = field_groups->size("colour");

  for (int ic = 0; ic < nc; ic++){
    enzo_float * cfield = (enzo_float *)
      field.values(field_groups->item("colour",ic));

    cfield[index] *= density_ratio;

  }

  return;
}

/*
void EnzoMethodStarMaker::convert_densities_to_fraction(EnzoBlock * enzo_block,
                                                        int direction) throw() {

  // Actually, I don't really think we need this...
  //   this only needs to be done with cells that either have
  //   star formation or get gas removed for star formation. This is
  //   likely to be a small number of cells on a given grid, so
  //   really there is no need to do this conversion for every cell...


  Field field = enzo_block->data()->field();

  int gx,gy,gz;
  field.ghost_depth (0, &gx, &gy, &gz);

  int mx, my, mz;
  field.dimensions (0, &mx, &my, &mz);

  int nx, ny, nz;
  field.size ( &nx, &ny, &nz);

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;
  int ngz = nz + 2*gz;

  enzo_float * d = (enzo_float *) field.values("density");

  if (direction == 0){ // convert density to fraction

    for (int iz = 0; iz < ngz; iz++){
      for (int iy = 0; iy < ngy; iy++){
        for (int ix = 0; ix < ngx; ix++){
          int i = INDEX(ix,iy,iz,ngx,ngy);

          double inv_dens = 1.0 / d[i];

          if (metal) metal[i] = metal[i] * inv_dens;
        }
      }
    }


  } else { // convert fraction to density

    if (metal) metal[i] *= d[i];
  }

  return;
}
*/

// ---------------------------------------------------------

int EnzoMethodStarMaker::check_number_density_threshold(
                                                       const double &d
                                                        ){

  ///  Apply the criteria that the local number density be greater
  ///  than the provided number density if use_density_threshold_ is
  ///  desired by the user.

  return !(use_density_threshold_) +
          (d >= number_density_threshold_);
}

int EnzoMethodStarMaker::check_self_gravitating(
                const double mean_particle_mass, const double rho_cgs, const enzo_float temperature,
                enzo_float *vx, enzo_float *vy, enzo_float *vz,
                const double lunit, const double vunit,
                const int &index, const int &dix, const int &diy, const int &diz,
                const double dx, const double dy, const double dz)
{

  if (!use_self_gravitating_)
    return 1;

  // Hopkins et al. (2013). Virial parameter: alpha < 1 -> self-gravitating
  
  double div_v_norm2, cs2, alpha;
  double dx2 = dx*dx * lunit*lunit;
  double dy2 = dy*dy * lunit*lunit;
  double dz2 = dz*dz * lunit*lunit;

  // Frobenius norm of the velocity gradient tensor
  div_v_norm2 = (pow(vx[index+dix] - vx[index-dix], 2) +
                 pow(vy[index+dix] - vy[index-dix], 2) +
                 pow(vz[index+dix] - vz[index-dix], 2)) / dx2 +
                (pow(vx[index+diy] - vx[index-diy], 2) +
                 pow(vy[index+diy] - vy[index-diy], 2) +
                 pow(vz[index+diy] - vz[index-diy], 2)) / dy2 +
                (pow(vx[index+diz] - vx[index-diz], 2) +
                 pow(vy[index+diz] - vy[index-diz], 2) +
                 pow(vz[index+diz] - vz[index-diz], 2)) / dz2;
  div_v_norm2 *= (vunit*vunit);

  // constant for testing. TODO: change to variable
  const double gamma = 5.0 / 3.0;
  cs2 = (gamma * cello::kboltz * temperature) / mean_particle_mass;

  alpha = (div_v_norm2 + cs2/dx2) / (8 * cello::pi * cello::grav_constant * rho_cgs);
  return (alpha < 1);

}

double EnzoMethodStarMaker::h2_self_shielding_factor(
                enzo_float *rho, const double rho_cgs,
                const double dunit, const double lunit,
                const int &index, const int &dix, const int &diy, const int &diz,
                const double dx, const double dy, const double dz,
                const double metallicity)
{

  if (!use_h2_self_shielding_)
    return 1;

  // Hopkins et al. (2017) and Krumholz & Gnedin (2011). Constant numbers come from their models and fits.
  // Mass fraction that is self-shielded and able to cool. f_shield > 0

  double tau, phi, psi, f_shield, grad_rho;

  grad_rho = sqrt(pow((rho[index+dix] - rho[index-dix]) / dx, 2) +
                  pow((rho[index+diy] - rho[index-diy]) / dy, 2) +
                  pow((rho[index+diz] - rho[index-diz]) / dz, 2));
  grad_rho *= dunit / lunit;
  tau = 434.8 * rho_cgs * (dx + rho_cgs / grad_rho);  // 434.8 cm^2 / g
  phi = 0.756 * pow(1.0 + 3.1 * metallicity, 0.365);
  psi = (0.6 * tau * (0.01 + metallicity)) / (log(1.0 + 0.6*phi + 0.01*phi*phi));
  f_shield = 1.0 - 3.0 / (1.0 + 4.0*psi);
  return f_shield;

}

int EnzoMethodStarMaker::check_jeans_mass(
  const double temperature, const double mean_particle_mass, 
  const double rho_cgs, const double mass
)
{
  if (!use_jeans_mass_)
    return 1;

  const double gamma = 5.0 / 3.0;
  const double minimum_jeans_mass = 1000 * cello::mass_solar;
  double cs2 = (gamma * cello::kboltz * temperature) / mean_particle_mass;
  double m_jeans = (cello::pi/6) * pow(cs2, 1.5) / (pow(cello::grav_constant, 1.5) * sqrt(rho_cgs));
  double m_jcrit = MAX(minimum_jeans_mass, m_jeans);
  return (mass < m_jcrit);
}

int EnzoMethodStarMaker::check_velocity_divergence(
                enzo_float *vx, enzo_float *vy, enzo_float *vz,
                const int &index, const int &dix, const int &diy,
                const int &diz){

    ///  Apply the criteria that the divergence of the velocity
    ///  be negative, if so desired by user (use_velocity_divergence).

    if (!(this->use_velocity_divergence_)){
      return 1.0;
    }

    int result = 0;

    if(vx){
      result = (vx[index+dix] - vx[index-dix] < 0) ? 1 : 0;
    }

    if(vy && result){
      result = (vy[index+diy] - vy[index-diy] < 0) ? 1 : 0;
    }

    if(vz && result){
      result = (vz[index+diz] - vz[index-diz] < 0) ? 1 : 0;
    }

    return result;
}

int EnzoMethodStarMaker::check_mass(const double &m){
  /// Apply the condition that the mass of gas converted into
  /// stars in a single cell cannot exceed a certain fraction
  /// of that cell's mass
  int minlimit = ((maximum_star_fraction_ * m) > star_particle_min_mass_);
  int maxlimit = ((maximum_star_fraction_ * m) < star_particle_max_mass_);
  return minlimit && maxlimit;

}
