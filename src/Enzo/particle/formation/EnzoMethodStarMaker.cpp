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

#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"
#include "Enzo/particle/particle.hpp"

//#define DEBUG_SF

//-------------------------------------------------------------------

EnzoMethodStarMaker::EnzoMethodStarMaker(ParameterGroup p)
  : Method()
{
  cello::particle_descr()->check_particle_attribute("star","mass");

  cello::simulation()->refresh_set_name(ir_post_,name());

  Refresh * refresh = cello::refresh(ir_post_);

  refresh->add_all_fields();

  // parameter parsing:
  use_density_threshold_     = p.value_logical("use_density_threshold",false);
  use_velocity_divergence_   = p.value_logical("use_velocity_divergence",false);
  use_self_gravitating_      = p.value_logical("use_self_gravitating", false);
  use_altAlpha_              = p.value_logical("use_altAlpha",false);
  use_h2_self_shielding_     = p.value_logical("use_h2_self_shielding", false);
  use_jeans_mass_            = p.value_logical("use_jeans_mass", false);
  number_density_threshold_  = p.value_float("number_density_threshold",0.0);
  efficiency_                = p.value_float("efficiency",0.01);
  maximum_star_fraction_     = p.value_float("maximum_mass_fraction",0.05);
  star_particle_min_mass_    = p.value_float("minimum_star_mass",0.0);
  star_particle_max_mass_    = p.value_float("maximum_star_mass",-1.0);
  use_dynamical_time_        = p.value_logical("use_dynamical_time",false);

  use_overdensity_threshold_ = p.value_logical("use_overdensity_threshold",
                                               false);
  overdensity_threshold_     = p.value_float("overdensity_threshold",0.0);
  use_critical_metallicity_  = p.value_logical("use_critical_metallicity",
                                               false);
  critical_metallicity_      = p.value_float("critical_metallicity",0.0);
  use_cooling_time_          = p.value_logical("use_cooling_time",false);
  use_temperature_threshold_ = p.value_logical("use_temperature_threshold",
                                               false);
  temperature_threshold_     = p.value_float("temperature_threshold",1.0E4);
}

//-------------------------------------------------------------------

void EnzoMethodStarMaker::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | use_density_threshold_;
  p | use_velocity_divergence_;
  p | use_self_gravitating_;
  p | use_altAlpha_;
  p | use_h2_self_shielding_;
  p | use_jeans_mass_;
  p | use_overdensity_threshold_;
  p | use_critical_metallicity_;
  p | use_temperature_threshold_;
  p | use_cooling_time_;
  p | use_dynamical_time_;
  p | overdensity_threshold_;
  p | critical_metallicity_;
  p | number_density_threshold_;
  p | efficiency_;
  p | maximum_star_fraction_;
  p | star_particle_min_mass_;
  p | star_particle_max_mass_;
  p | temperature_threshold_;
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

//----------------------------------------------------------------------

void EnzoMethodStarMaker::rescale_densities(EnzoBlock * enzo_block,
                                            const int index,
                                            const double density_ratio) throw() {

  // Loop through all passive scalars (color fields)
  // which are mass fractions stored as densities, and rescale
  // to the new density after star formation.
  //
  // AE: NOTE: Change this routine if / whenever there needs to be
  //           fraction fields that are not labelled as color
  //           Obviously requires these fields to be declared as color
  //           in input file to work.
  //           This can / should likely get moved to be a more general
  //           function of the block / data / field class (one of those)
  //
  //    density_ratio = new_density / old_density
  //

  Field field = enzo_block->data()->field();

  Grouping * field_groups = field.groups();
  int nc = field_groups->size("color");

  for (int ic = 0; ic < nc; ic++){
    enzo_float * cfield = (enzo_float *)
      field.values(field_groups->item("color",ic));
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

  return !(this->use_density_threshold_) +
          (d >= this->number_density_threshold_);
}

int EnzoMethodStarMaker::check_overdensity_threshold(const double &rho)
{
  //#ifdef DEBUG_SF
  //  CkPrintf("MethodStarMaker -- overdensity = %f\n", rho);
  //#endif
  return !(this->use_overdensity_threshold_) +  
          (rho >= this->overdensity_threshold_);
}


int EnzoMethodStarMaker::check_self_gravitating_alt(const double total_energy, const double potential)
{

  if (!this->use_self_gravitating_)
    return 1;

  #ifdef DEBUG_SF
    CkPrintf("MethodStarMaker -- alpha = %f\n",total_energy/potential); 
  #endif

  return (total_energy/potential < 1);

}

int EnzoMethodStarMaker::check_self_gravitating(
                const double mean_particle_mass, const double density, const enzo_float temperature,
                enzo_float *vx, enzo_float *vy, enzo_float *vz,
                const double lunit, const double vunit, const double rhounit,
                const int &index, const int &dix, const int &diy, const int &diz,
                const double dx, const double dy, const double dz)
{

  if (!this->use_self_gravitating_)
    return 1;

  // Hopkins et al. (2013). Virial parameter: alpha < 1 -> self-gravitating

  double div_v_norm2, cs2, alpha;
  double dx2 = dx*dx * lunit*lunit;
  double dy2 = dy*dy * lunit*lunit;
  double dz2 = dz*dz * lunit*lunit;

  // Frobenius norm of the velocity gradient tensor
  div_v_norm2 = (pow(vx[index+dix] - vx[index-dix], 2) +
                 pow(vy[index+dix] - vy[index-dix], 2) +
                 pow(vz[index+dix] - vz[index-dix], 2)) / (4*dx2) +
                (pow(vx[index+diy] - vx[index-diy], 2) +
                 pow(vy[index+diy] - vy[index-diy], 2) +
                 pow(vz[index+diy] - vz[index-diy], 2)) / (4*dy2) +
                (pow(vx[index+diz] - vx[index-diz], 2) +
                 pow(vy[index+diz] - vy[index-diz], 2) +
                 pow(vz[index+diz] - vz[index-diz], 2)) / (4*dz2);
  div_v_norm2 *= (vunit*vunit); 

  // constant for testing. TODO: change to variable
  const double gamma = 5.0 / 3.0;
  cs2 = (gamma * enzo_constants::kboltz * temperature) / mean_particle_mass;

  alpha = (div_v_norm2 + cs2/dx2) / (8 * cello::pi * enzo::grav_constant_cgs() * density*rhounit);
  #ifdef DEBUG_SF
    CkPrintf("MethodStarMaker -- alpha = %f\n",alpha); 
  #endif

  return (alpha < 1);

}

double EnzoMethodStarMaker::h2_self_shielding_factor(
                enzo_float *rho, const double metallicity,
                const double dunit, const double lunit,
                const int &index, const int &dix, const int &diy, const int &diz,
                const double dx, const double dy, const double dz)
{

  if (!this->use_h2_self_shielding_)
    return 1;

  // Hopkins et al. (2017) and Krumholz & Gnedin (2011). Constant numbers come from their models and fits.
  // Mass fraction that is self-shielded and able to cool. f_shield > 0

  double tau, phi, psi, f_shield, grad_rho;

  const double rho_cgs = rho[index] * dunit;

  grad_rho = sqrt(pow((rho[index+dix] - rho[index-dix]) / (2*dx), 2) +
                  pow((rho[index+diy] - rho[index-diy]) / (2*dy), 2) +
                  pow((rho[index+diz] - rho[index-diz]) / (2*dz), 2));
  grad_rho *= dunit / lunit;
  tau = 434.8 * rho_cgs * (dx*lunit + rho_cgs / grad_rho);  // 434.8 cm^2 / g
  phi = 0.756 * pow(1.0 + 3.1 * metallicity, 0.365);
  psi = (0.6 * tau * (0.01 + metallicity)) / (log(1.0 + 0.6*phi + 0.01*phi*phi));
  f_shield = 1.0 - 3.0 / (1.0 + 4.0*psi);
  #ifdef DEBUG_SF
    CkPrintf("MethodStarMaker -- f_shield = %f\n",f_shield); 
  #endif
  return f_shield;

}

int EnzoMethodStarMaker::check_jeans_mass(
  const double temperature, const double mean_particle_mass,
  const double density, const double mass,
  const double munit, const double rhounit
)
{
  if (!use_jeans_mass_)
    return 1;

  const double gamma = 5.0 / 3.0;
  const double minimum_jeans_mass = 1000 * enzo_constants::mass_solar;
  double cs2 = (gamma * enzo_constants::kboltz * temperature) / mean_particle_mass;

  double m_jeans = (cello::pi/6) * pow(cs2, 1.5) / 
                   (pow(enzo::grav_constant_cgs(), 1.5) * sqrt(density*rhounit));

  double m_jcrit = MAX(minimum_jeans_mass, m_jeans);
  #ifdef DEBUG_SF
    CkPrintf("MethodStarMaker -- jeans_mass = %f\n",m_jeans); 
  #endif
  return (mass*munit > m_jcrit);
}

int EnzoMethodStarMaker::check_velocity_divergence(
                enzo_float *vx, enzo_float *vy, enzo_float *vz,
                const int &index, const int &dix, const int &diy,
                const int &diz, const double dx, const double dy, const double dz) {

    ///  Apply the criteria that the divergence of the velocity
    ///  be negative, if so desired by user (use_velocity_divergence).

    if (!(this->use_velocity_divergence_)){
      return 1.0;
    }

   double div = 0.0;
   if (vx) div += 0.5 * (vx[index+dix] - vx[index-dix]) / dx; // in units of dx
   if (vy) div += 0.5 * (vy[index+diy] - vy[index-diy]) / dy; // in units of dy
   if (vz) div += 0.5 * (vz[index+diz] - vz[index-diz]) / dz; // in units of dz

   #ifdef DEBUG_SF
     CkPrintf("MethodStarMaker -- velocity_divergence = %f\n",div); 
   #endif

   return div < 0;
}

int EnzoMethodStarMaker::check_mass(const double &m){
  /// Apply the condition that the mass of gas converted into
  /// stars in a single cell cannot exceed a certain fraction
  /// of that cell's mass. There does not need to be a check on
  /// the maximum particle mass.

  int minlimit = ((maximum_star_fraction_ * m) > star_particle_min_mass_);
  return minlimit;

}

int EnzoMethodStarMaker::check_cooling_time(const double &cooling_time,const double &total_density,
                         const double tunit, const double rhounit)
{
  // Check whether cooling_time < dynamical_time
  if (!(this->use_cooling_time_)) {
    return 1;
  }

  double dynamical_time = pow(3.0*cello::pi/32.0/enzo::grav_constant_cgs()/(total_density*rhounit),0.5); //s
  #ifdef DEBUG_SF
    CkPrintf("MethodStarMaker -- cooling_time = %f, dynamical_time = %f\n",cooling_time*tunit, dynamical_time); 
  #endif
  return cooling_time*tunit < dynamical_time;
  
}

int EnzoMethodStarMaker::check_metallicity(const double &Z) 
{
  // Enforce a critical metallicity for star formation
  #ifdef DEBUG_SF
    CkPrintf("MethodStarMaker -- Z = %f, Zcrit = %f\n",Z, critical_metallicity_); 
  #endif
  return !(this->use_critical_metallicity_) +
          (Z >= critical_metallicity_);
}

int EnzoMethodStarMaker::check_temperature(const double &T)
{
  #ifdef DEBUG_SF
    CkPrintf("MethodStarMaker -- T = %f, Tcrit = %f\n", T, temperature_threshold_); 
  #endif
  return !(this->use_temperature_threshold_) +
          (T < temperature_threshold_);
}
