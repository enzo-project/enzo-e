// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoInitialIsolatedGalaxy.cpp
/// @author	Andrew Emerick (aemerick11@gmail.com)
/// @date	Tue May 8
/// @brief	Implementation of an isolated galaxy simulation
///
/// More explanation here
///

#include <fstream>

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------
#define DEBUG_PERFORMANCE
//----------------------------------------------------------------------

// convenience flag when initializing using MakeDiskGalaxy gas particles
#define GAS_PARTICLE_FLAG -999

int nlines(std::string fname) {

  int n = 0;
  std::string line;
  std::ifstream inFile(fname);

  while (std::getline(inFile, line))
    ++n;

  inFile.close();

  return n;
}

EnzoInitialIsolatedGalaxy::EnzoInitialIsolatedGalaxy
(const EnzoConfig * config) throw()
: Initial(config->initial_cycle, config->initial_time)
{
  // read in parameter settings from config and
  // set corresponding member variables
  for(int i = 0; i < 3; i ++){
    center_position_[i] = config->initial_IG_center_position[i];
    bfield_[i] = config->initial_IG_bfield[i];
  }

  EnzoUnits * enzo_units = enzo::units();
  // Store variables locally with code units

  this->scale_height_         = config->initial_IG_scale_height * cello::kpc_cm /
                                enzo_units->length();
  this->scale_length_	      = config->initial_IG_scale_length * cello::kpc_cm /
                                enzo_units->length();
  this->disk_mass_            = config->initial_IG_disk_mass * cello::mass_solar /
                                enzo_units->mass();
  this->gas_fraction_         = config->initial_IG_gas_fraction;
  this->disk_temperature_     = config->initial_IG_disk_temperature /
                                enzo_units->temperature();
  this->disk_metal_fraction_  = config->initial_IG_disk_metal_fraction;
  this->gas_halo_mass_        = config->initial_IG_gas_halo_mass * cello::mass_solar /
                                enzo_units->mass();
  this->gas_halo_temperature_ = config->initial_IG_gas_halo_temperature /
                                enzo_units->temperature();
  this->gas_halo_metal_fraction_ = config->initial_IG_gas_halo_metal_fraction;
  this->gas_halo_density_     = config->initial_IG_gas_halo_density /
                                enzo_units->density();
  this->gas_halo_radius_      = config->initial_IG_gas_halo_radius * cello::kpc_cm /
                                enzo_units->length();

  this->use_gas_particles_ = config->initial_IG_use_gas_particles;

  this->live_dm_halo_         = config->initial_IG_live_dm_halo;
  this->stellar_disk_         = config->initial_IG_stellar_disk;
  this->stellar_bulge_        = config->initial_IG_stellar_bulge;

  this->analytic_velocity_    = config->initial_IG_analytic_velocity;

  if ((this->gas_halo_density_ == 0.0) && this->gas_halo_mass_ > 0)
  {
    if (this->gas_halo_radius_ > 0.0)
    {
      this->gas_halo_density_ = this->gas_halo_mass_ /
                    (4.0 * cello::pi / 3.0 * pow( this->gas_halo_radius_, 3));
    } else { // assume full box of LUxLUxLU
      this->gas_halo_density_ = this->gas_halo_mass_; // In code density  mass / 1**3
    }
  } // AE: Need an error check to make sure these are set correctly


  this->uniform_density_    = config->field_uniform_density;
  this->dual_energy_        = config->ppm_dual_energy; // or ppm? which one?
  this->gamma_              = config->field_gamma;
  this->mu_                 = config->ppm_mol_weight;

  // read in data for initialization
  this->ntypes_            = 0;
  this->ndim_              = cello::rank();
  this->nparticles_        = 0;

  this->ReadInVcircData();        // circular velocity curve
  this->ReadParticles();     // IC particles - only if available
}

//----------------------------------------------------------------------

void EnzoInitialIsolatedGalaxy::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP; //

  Initial::pup(p); //

  PUParray(p,center_position_,3);
  PUParray(p, bfield_,3);
  p | scale_length_;
  p | disk_mass_;
  p | gas_fraction_;
  p | disk_temperature_;
  p | disk_metal_fraction_;
  p | gas_halo_mass_;
  p | gas_halo_temperature_;
  p | gas_halo_metal_fraction_;
  p | gas_halo_density_;
  p | gas_halo_radius_;

  p | uniform_density_;
  p | dual_energy_;
  p | gamma_;
  p | mu_;

  p | use_gas_particles_;
  p | live_dm_halo_;
  p | stellar_bulge_;
  p | stellar_disk_;
  p | analytic_velocity_;

  p | ntypes_;
  p | ndim_;
  p | nparticles_;

  if (p.isUnpacking()){
    // we know aboutntype, ndim, and nparticles..
    // allocate particle IC arrays
    allocateParticles();
  }

  // Pup grid values element-by-element
  for (int k = 0; k < ntypes_; k++){
    for (int j = 0; j < ndim_; j ++){
      PUParray(p, particleIcPosition[k][j], nparticles_);
      PUParray(p, particleIcVelocity[k][j], nparticles_);
    }
    PUParray(p, particleIcMass[k], nparticles_);
  }
  PUParray(p, particleIcTypes, ntypes_);

  p | particleIcFileNames;

  return;
}

void EnzoInitialIsolatedGalaxy::enforce_block
(
 Block * block,
 const Hierarchy * hierarchy
 ) throw()
{

  //
  // Make sure we can operate on this block
  //
  if (!block->is_leaf()) return;

  Timer timer;
  timer.start();

  ASSERT("EnzoInitialIsolatedGalaxy",
         "Block does not exist",
         block != NULL);

#ifdef CONFIG_USE_GRACKLE
   grackle_field_data grackle_fields_;
   EnzoMethodGrackle::setup_grackle_fields(block, &grackle_fields_);
#endif

  if (this->use_gas_particles_){
    this->InitializeGasFromParticles(block);
  } else {
    this->InitializeExponentialGasDistribution(block);
  }


#ifdef CONFIG_USE_GRACKLE
  /* Assign grackle chemistry fields to default fractions based on density */
  EnzoMethodGrackle::update_grackle_density_fields(&grackle_data);
#endif

  // Update temperature field if it exists
  Field field = block->data()->field();
  const EnzoConfig * enzo_config = enzo::config();
  enzo_float * temperature = field.is_field("temperature") ?
                   (enzo_float*) field.values("temperature") : NULL;

  if (temperature) {
    EnzoComputeTemperature compute_temperature
      (enzo_config->ppm_density_floor,
       enzo_config->ppm_temperature_floor,
       enzo_config->ppm_mol_weight,
       enzo_config->physics_cosmology);

    compute_temperature.compute(block);
  }

  // now initialize particles
  Particle particle = block->data()->particle();
  InitializeParticles(block, &particle);

  return;
}

void EnzoInitialIsolatedGalaxy::InitializeExponentialGasDistribution(Block * block){

  // Initialize gas distribution using a double exponential

  EnzoUnits * enzo_units = enzo::units();
  const EnzoConfig * enzo_config = enzo::config();
  Field field = block->data()->field();

  //
  // Get Grid and Field parameters
  //

  // Block sizes (excluding ghost zones)
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  // Ghost zone depth (gx + gx total ghost zones in x dimension)
  int gx, gy, gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  // total number of cells in each dimension
  int mx = nx + 2*gx;
  int my = ny + 2*gy;
  int mz = nz + 2*gz;

  // Grab information about grid properties
  //     Cell min coords (xm,ym,zm)
  //     Cell max coords (xp,yp,zp)
  //     Cell widths (hx,hy,hz)
  double xm, ym, zm, xp, yp, zp, hx, hy, hz;
  block->data()->lower(&xm,&ym,&zm);
  block->data()->upper(&xp,&yp,&zp);
  block->cell_width(&hx,&hy,&hz);

  // Get Fields
  enzo_float * d = (enzo_float *) field.values ("density");
  enzo_float * temperature = field.is_field("temperature") ?
                   (enzo_float*) field.values("temperature") : NULL;
  enzo_float * p = field.is_field("pressure") ?
                   (enzo_float *) field.values ("pressure") : NULL;
  enzo_float * a3[3] = { (enzo_float *) field.values("acceleration_x"),
                         (enzo_float *) field.values("acceleration_y"),
                         (enzo_float *) field.values("acceleration_z")};
  enzo_float * v3[3] = { (enzo_float *) field.values("velocity_x"),
                         (enzo_float *) field.values("velocity_y"),
                         (enzo_float *) field.values("velocity_z")};

  enzo_float * te  = (enzo_float *) field.values("total_energy");
  enzo_float * ge  = (enzo_float *) field.values("internal_energy");
  enzo_float * pot = (enzo_float *) field.values("potential");

  enzo_float * metal = field.is_field("metal_density") ?
                      (enzo_float *) field.values("metal_density") : NULL;


  // double velocity_units = this->length_units_ / this->time_units_;
  // double temperature_units = cello::mass_hydrogen *
  //                           pow(velocity_units, 2) / cello::k;

  //
  // Now lets calculate some physical properties of the galaxy and halo
  //
  double rho_zero = this->disk_mass_ * this->gas_fraction_ / (4.0 * cello::pi)/
      (pow((this->scale_length_),2)*(this->scale_height_));

  double halo_gas_energy = this->gas_halo_temperature_ / this->mu_ / (this->gamma_ -1);

  double disk_gas_energy = this->disk_temperature_ / this->mu_ / (this->gamma_ -1) ;

  double tiny_number = 1.0E-10;

  // initialize fields to something
  for (int iz=0; iz<mz; iz++){
    for (int iy=0; iy<my; iy++){
      for (int ix=0; ix<mx; ix++){
        int i = INDEX(ix,iy,iz,mx,my);
        d[i]   = this->uniform_density_;
        te[i]  = halo_gas_energy;
        p[i]   = (this->gamma_ - 1.0) * te[i] * d[i];

        if(pot) pot[i] = 0.0;

        if (this->dual_energy_)
        {
          ge[i] = halo_gas_energy;
        }

        for (int dim = 0; dim < 3; dim++){
          a3[dim][i] = 0.0;
          v3[dim][i] = 0.0;
        }

        if(metal) metal[i] = tiny_number * d[i];

      }
    }
  } // end loop over all cells for background values


  for (int iz=0; iz<mz; iz++){
    // compute z coordinate of cell center
    double z = zm + (iz - gz + 0.5)*hz - this->center_position_[2];
    z *= enzo_units->length(); // convert to cgs units

    for (int iy=0; iy<my; iy++){
      double y = ym + (iy - gy + 0.5)*hy - this->center_position_[1];
      y *= enzo_units->length();

      for (int ix=0; ix<mx; ix++){
        double x = xm + (ix - gx + 0.5)*hx - this->center_position_[0];
        x *= enzo_units->length();

        // 1D index of current cell
        int i = INDEX(ix,iy,iz,mx,my);

        // compute spherical and cylindrical radii (in cgs)
        double radius = sqrt(x*x + y*y + z*z);
        double r_cyl  = sqrt(x*x + y*y);

        // compute the disk density (in code units)
        double disk_density = this->gauss_mass(rho_zero, x/enzo_units->length(), y/enzo_units->length(),
                                              z/enzo_units->length(), hx) / (hx*hx*hx);


        if ((this->gas_halo_density_ * this->gas_halo_temperature_ > disk_density*this->disk_temperature_) &&
            (radius < this->gas_halo_radius_*enzo_units->length())){
          // in the halo, set the halo properties
          d[i]  = this->gas_halo_density_;
          te[i] = halo_gas_energy;
          p[i]  = (this->gamma_ - 1.0) * te[i] * d[i];

          if (this->dual_energy_) {
            ge[i] = halo_gas_energy;
          }

          // set metal fraction here
          if(metal) metal[i] = this->gas_halo_metal_fraction_ * d[i];

        }
        else if ( radius < this->gas_halo_radius_*enzo_units->length() ) // we are in the disk
        {
          // in the disk, set the disk properties
          d[i]   = disk_density;

          double vcirc = 0.0;
          if (this->analytic_velocity_){
//            double rhodm = enzo_config->method_background_acceleration_DM_density;
            double rcore = enzo_config->method_background_acceleration_core_radius;
            double rvir  = enzo_config->method_background_acceleration_DM_mass_radius;

            double Mvir  = enzo_config->method_background_acceleration_DM_mass;

            rcore = rcore * cello::kpc_cm;
            rvir  = rvir  * cello::kpc_cm;
            Mvir  = Mvir  * cello::mass_solar;

            double conc = rvir  / rcore;
            double    x = r_cyl / rvir;
            

            vcirc = (std::log(1.0 + conc*x) - (conc*x)/(1.0+conc*x))/
                      (std::log(1.0 + conc) - (conc / (1.0 + conc))) / x;
            vcirc = std::sqrt(vcirc * cello::grav_constant * Mvir / rvir);

          } else {
            vcirc = this->InterpolateVcircTable(r_cyl);
          }

          v3[0][i] = -(vcirc*(y/r_cyl))/enzo_units->velocity();
          v3[1][i] = -(vcirc*(x/r_cyl))/enzo_units->velocity();
          v3[2][i] = 0.0;

          te[i]  = disk_gas_energy;
          p[i]   = (this->gamma_ - 1.0) * te[i] * d[i];

          for (int dim = 0; dim < 3; dim++) // AE: do I need a check for not ppm?
          {
            te[i] += 0.5*v3[dim][i]*v3[dim][i];
          }

          if (this->dual_energy_)
          {
            ge[i] = disk_gas_energy;
          }

          if(metal) metal[i] = this->disk_metal_fraction_ * d[i];

        } // end disk / halo check

      }
    }
  } // end loop over all cells

  return;
}

void EnzoInitialIsolatedGalaxy::InitializeGasFromParticles(Block * block){

  EnzoUnits * enzo_units = enzo::units();
  const EnzoConfig * enzo_config = enzo::config();
  Field field = block->data()->field();

  // Figure out which index corresponds to the gas particles
  int ipt = 0;
  for (ipt = 0; ipt < ntypes_; ipt++){
      if ( particleIcTypes[ipt] == GAS_PARTICLE_FLAG) break;
  }
  ASSERT("EnzoInitialIsolatedGalaxy",
         "Using MakeDiskGalaxy particles but gas particles not found",
         ipt < ntypes_);

  //
  // Get Grid and Field parameters
  //

  // Block sizes (excluding ghost zones)
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  // Ghost zone depth (gx + gx total ghost zones in x dimension)
  int gx, gy, gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  // total number of cells in each dimension
  int mx = nx + 2*gx;
  int my = ny + 2*gy;
  int mz = nz + 2*gz;

  // Grab information about grid properties
  //     Cell min coords (xm,ym,zm)
  //     Cell max coords (xp,yp,zp)
  //     Cell widths (hx,hy,hz)
  double xm, ym, zm, xp, yp, zp, hx, hy, hz;
  block->data()->lower(&xm,&ym,&zm);
  block->data()->upper(&xp,&yp,&zp);
  block->cell_width(&hx,&hy,&hz);

  // now loop over all gas particles and deposit
  // Get Fields
  enzo_float * d = (enzo_float *) field.values ("density");
  enzo_float * temperature = field.is_field("temperature") ?
                   (enzo_float*) field.values("temperature") : NULL;
  enzo_float * p = field.is_field("pressure") ?
                   (enzo_float *) field.values ("pressure") : NULL;
  enzo_float * a3[3] = { (enzo_float *) field.values("acceleration_x"),
                         (enzo_float *) field.values("acceleration_y"),
                         (enzo_float *) field.values("acceleration_z")};
  enzo_float * v3[3] = { (enzo_float *) field.values("velocity_x"),
                         (enzo_float *) field.values("velocity_y"),
                         (enzo_float *) field.values("velocity_z")};

  enzo_float * te  = (enzo_float *) field.values("total_energy");
  enzo_float * ge  = (enzo_float *) field.values("internal_energy");
  enzo_float * pot = (enzo_float *) field.values("potential");

  enzo_float * metal = field.is_field("metal_density") ?
                      (enzo_float *) field.values("metal_density") : NULL;

  double disk_gas_energy = this->disk_temperature_ / this->mu_ / (this->gamma_ - 1);
  double halo_gas_energy = this->gas_halo_temperature_ / this->mu_ / (this->gamma_ - 1);
  double tiny_number = 1.0E-10;

  // Initialize the full grid to either CGM properties or
  // IGM properties (if beyond virial radius)
  for (int iz=0; iz<mz; iz++){
    double z = zm + (iz - gz + 0.5)*hz - this->center_position_[2];
    z *= enzo_units->length();

    for (int iy=0; iy<my; iy++){
      double y = ym + (iy - gy + 0.5)*hy - this->center_position_[1];
      y *= enzo_units->length();

      for (int ix=0; ix<mx; ix++){
        double x = xm + (ix - gx + 0.5)*hx - this->center_position_[0];
        x *= enzo_units->length();

        // 1D index of current cell
        int i = INDEX(ix,iy,iz,mx,my);

        // compute spherical position (in cgs)
        double  radius = std::sqrt(x*x + y*y + z*z);

        // set everywhere internal to Rvir to CGM properties
        double metal_fraction = tiny_number;
        if (radius < this->gas_halo_radius_*enzo_units->length()){
          // in the halo
          d[i]  = this->gas_halo_density_;
          metal_fraction = this->gas_halo_metal_fraction_;
        } else {
          d[i]  = this->uniform_density_;
        }

        te[i] = halo_gas_energy;

        p[i]  = (this->gamma_ - 1.0) * te[i] * d[i];

        if(pot) pot[i] = 0.0;

        if (this->dual_energy_) ge[i] = halo_gas_energy;

        for (int dim=0; dim < 3; dim++){
          a3[dim][i] = 0.0;
          v3[dim][i] = 0.0; // static
        }

        if (metal) metal[i] = metal_fraction * d[i];

      }
    }
  } // -- end loop over grid for background conditions

  // flagging field to flag cells with particle deposit
  //   this will be used later to do some smoothing
  //   of cells adjacent to cells with particle deposits
  int *iflag = new int[mx*my*mz];
  for (int i = 0; i < mx*my*mz; i++) iflag[i] = 0;

  // Now loop over all particles and deposit
  int np = nlines(particleIcFileNames[ipt]);

  for (int ip = 0; ip < np; ip++){

    if ( !(block->check_position_in_block(particleIcPosition[ipt][0][ip],
                                          particleIcPosition[ipt][1][ip],
                                          particleIcPosition[ipt][2][ip]))){
      continue;
    }

    // get corresponding grid position
    double xp = (particleIcPosition[ipt][0][ip] - xm) / hx;
    double yp = (particleIcPosition[ipt][1][ip] - ym) / hy;
    double zp = (particleIcPosition[ipt][2][ip] - zm) / hz;

    int ix = ((int) std::floor(xp));
    int iy = ((int) std::floor(yp));
    int iz = ((int) std::floor(zp));

    int i  = INDEX(ix,iy,iz,mx,my);

    d[i]    = d[i]*iflag[i] + (particleIcMass[ipt][ip] / (hx*hx*hx));

    // add momentum
    for (int dim = 0; dim < 3; dim++){
      v3[dim][i] = particleIcMass[ipt][ip] * particleIcVelocity[ipt][dim][ip] / (hx*hx*hx);
    }

    iflag[i] = 1;
  }

  // re-adjust velocities and set energies
  for (int iz=0; iz<mz; iz++){
    for (int iy=0; iy<my; iy++){
      for (int ix=0; ix<mx; ix++){
        // 1D index of current cell
        int i = INDEX(ix,iy,iz,mx,my);

        if (iflag[i] == 0) continue;

        te[i] = disk_gas_energy;
        for (int dim = 0; dim < 3; dim++){
          v3[dim][i] /= d[i];
          te[i]     += 0.5*v3[dim][i]*v3[dim][i];
        }

        if (this->dual_energy_) ge[i] = disk_gas_energy;

        if (metal) metal[i] = this->disk_metal_fraction_ * d[i];

      }
    }
  }

  // finally, average velocity of cells adjcent to deposit cells
  // to reduce shocks. 
  int xo = 1;
  int yo = mx;
  int zo = mx*my;
  for (int iz=0; iz<mz; iz++){
    for (int iy=0; iy<my; iy++){
      for (int ix=0; ix<mx; ix++){
        // 1D index of current cell
        int i = INDEX(ix,iy,iz,mx,my);

        if (iflag[i]) continue; // skip deposit cells

        // get indeces of all adjecent cells
        int xl, xh, yl, yh, zl, zh;
        xl = std::max(i - xo, 0);
        xh = std::min(i + xo, mx*my*mz - 1);
        yl = std::max(i - yo, 0);
        yh = std::min(i + yo, mx*my*mz - 1);
        zl = std::max(i - zo, 0);
        zh = std::min(i + zo, mx*my*mz - 1);

        // check if adjacent to a deposit cell and compute
        // mass weighting
        double mass_deposit = 0.0;
        mass_deposit =  iflag[xh]*d[xh] + iflag[xl]*d[xl] +
                        iflag[yh]*d[yh] + iflag[yl]*d[yl] +
                        iflag[zh]*d[zh] + iflag[zl]*d[zl];

        // if adjcent to deposit cells, set velocity to
        // mass-weighted average
        if (mass_deposit > 0){
          for (int dim=0; dim<3; dim++){
            v3[dim][i] = ( iflag[xh]*v3[dim][xh]*d[xh] + iflag[xl]*v3[dim][xl]*d[xl] +
                          iflag[yh]*v3[dim][yh]*d[yh] + iflag[yl]*v3[dim][yl]*d[yl] +
                          iflag[zh]*v3[dim][zh]*d[zh] + iflag[zl]*v3[dim][zl]*d[zl]
                             ) / mass_deposit ;
          }
        }

        // 

      }
    }
  }

  delete [] iflag;
  iflag = NULL;

  return;
}

void EnzoInitialIsolatedGalaxy::InitializeParticles(Block * block,
                                                    Particle * particle){

  if (this->ntypes_ == 0) return;

  int rank = cello::rank();

  // Loop over all particle types and initialize
  for(int ipt = 0; ipt < ntypes_; ipt++){

    int it   = particleIcTypes[ipt];

    if (it == GAS_PARTICLE_FLAG) continue; // do not make these actual particles

    int np   = nlines(particleIcFileNames[ipt]);

    //
    // Now create the particles and assign values
    //    There is likely a better way of doing this that is both cleaner
    //    and more efficient, but this is likely the most readable
    //

    // obtain particle attribute indexes for this type
    int ia_m = particle->attribute_index (it, "mass");
    int ia_x = particle->attribute_index (it, "x");
    int ia_y = particle->attribute_index (it, "y");
    int ia_z = particle->attribute_index (it, "z");
    int ia_vx = particle->attribute_index (it, "vx");
    int ia_vy = particle->attribute_index (it, "vy");
    int ia_vz = particle->attribute_index (it, "vz");

    int ib  = 0; // batch counter
    int ipp = 0; // counter

    // this will point to the particular value in the
    // particle attribute array
    enzo_float * pmass = 0;
    enzo_float * px   = 0;
    enzo_float * py   = 0;
    enzo_float * pz   = 0;
    enzo_float * pvx  = 0;
    enzo_float * pvy  = 0;
    enzo_float * pvz  = 0;

    // now loop over all particles
    for (int i = 0; i < np; i ++){
      ASSERT("EnzoInitialIsolatedGalaxy",
             "Attempting to initialize a particle with negative mass",
              particleIcMass[ipt][i] > 0);

      // make sure particle exists on this grid before depositing
      //    AE:  Note, I'm not sure if this is needed (it is in Enzo),
      //         but there may be a smarter way to do this for Enzo-E
      //
      if (!(block->check_position_in_block(particleIcPosition[ipt][0][i],
                                           particleIcPosition[ipt][1][i],
                                           particleIcPosition[ipt][2][i]))){
        continue;
      }

      int new_particle = particle->insert_particles(it, 1);
      particle->index(new_particle,&ib,&ipp);

      pmass = (enzo_float *) particle->attribute_array(it, ia_m, ib);
      px    = (enzo_float *) particle->attribute_array(it, ia_x, ib);
      py    = (enzo_float *) particle->attribute_array(it, ia_y, ib);
      pz    = (enzo_float *) particle->attribute_array(it, ia_z, ib);
      pvx   = (enzo_float *) particle->attribute_array(it, ia_vx, ib);
      pvy   = (enzo_float *) particle->attribute_array(it, ia_vy, ib);
      pvz   = (enzo_float *) particle->attribute_array(it, ia_vz, ib);

      pmass[ipp] = particleIcMass[ipt][i];
      px[ipp]    = particleIcPosition[ipt][0][i];
      py[ipp]    = particleIcPosition[ipt][1][i];
      pz[ipp]    = particleIcPosition[ipt][2][i];
      pvx[ipp]   = particleIcVelocity[ipt][0][i];
      pvy[ipp]   = particleIcVelocity[ipt][1][i];
      pvz[ipp]   = particleIcVelocity[ipt][2][i];
    } // end loop over particles

  } // end loop over particle types


  return;
}

void EnzoInitialIsolatedGalaxy::ReadParticles(void){

  //
  // Check if particles are being used in initialization
  // and read in particle ICs if so.
  //
  //


  if (!(this->live_dm_halo_ || this->stellar_disk_ ||
       this->stellar_bulge_ || this->use_gas_particles_)){
    // do not load initial particles
    return ;
  } {
    ntypes_ = int(this->live_dm_halo_) + int(this->stellar_disk_) +
             int(this->stellar_bulge_) + int(this->use_gas_particles_);
  }

  ParticleDescr * particle_descr = cello::particle_descr();

  particleIcTypes = new int[ntypes_];

  int ipt = 0;
  if (this->live_dm_halo_){
    particleIcTypes[ipt] = particle_descr->type_index("dark");
    particleIcFileNames.push_back("halo.dat");
    nparticles_ = std::max(nparticles_, nlines("halo.dat"));
    ipt++;
  }
  if (this->stellar_disk_){
    particleIcTypes[ipt] = particle_descr->type_index("star");
    particleIcFileNames.push_back("disk.dat");
    nparticles_ = std::max(nparticles_, nlines("disk.dat"));
    ipt++;
  }
  if (this->stellar_bulge_){
    particleIcTypes[ipt] = particle_descr->type_index("star");
    particleIcFileNames.push_back("bulge.dat");
    nparticles_ = std::max(nparticles_, nlines("bulge.dat"));
    ipt++;
  }

  if (this->use_gas_particles_){
    // not actually using these as particles, but uses same IC reader
    // set particle type to specificied flag to ensure these won't get
    // initialized as actual particles
    particleIcTypes[ipt] = GAS_PARTICLE_FLAG;
    particleIcFileNames.push_back("gas.dat");
    nparticles_ = std::max(nparticles_, nlines("gas.dat"));
    ipt++;
  }

  // allocate particle IC arrays
  allocateParticles();

  // Read in data from files to arrays
  for (ipt = 0; ipt < ntypes_; ipt++){
      int num_lines = nlines(particleIcFileNames[ipt]);

      ReadParticlesFromFile(num_lines, particleIcPosition[ipt],
                                particleIcVelocity[ipt],
                                particleIcMass[ipt],
                                particleIcFileNames[ipt]);
  }

  return;
}


void EnzoInitialIsolatedGalaxy::ReadParticlesFromFile(const int& nl,
                                                      enzo_float ** position,
                                                      enzo_float ** velocity,
                                                      enzo_float * mass,
                                                      const std::string& filename){
  //
  // Generic particle IC reader. Position and velocity
  // are assumed to be 2D (dim, num_particles). Files must have columns
  // of position (x,y,z),  velocity (x,y,z), and mass. Units
  // are assumed to be pc, km/s, and Msun respectively
  //


  EnzoUnits * enzo_units = enzo::units();

  std::fstream inFile;
  inFile.open(filename, std::ios::in);

  ASSERT("EnzoInitialIsolatedGalaxy",
         "ParticleFile not found", inFile.is_open());

  int i = 0;
  while(inFile >>
        position[0][i] >> position[1][i] >> position[2][i] >>
        velocity[0][i] >> velocity[1][i] >> velocity[2][i] >>
        mass[i]){
    ASSERT("EnzoInitialIsolatedGalaxy",
           "Too many lines in particle file",
            i < nl);
    // positions are in pc
    // velocities in km / s
    // mass in Msun

    for (int dim = 0; dim < cello::rank(); dim++){

      position[dim][i] = position[dim][i] * cello::pc_cm /
                               enzo_units->length() + this->center_position_[dim];
      velocity[dim][i] = velocity[dim][i] * 1000.0 /
                               enzo_units->velocity();
    }

    mass[i]          = mass[i] * cello::mass_solar /
                       enzo_units->mass();

    i++;
  }

  inFile.close();

  return;
 }

void EnzoInitialIsolatedGalaxy::ReadInVcircData(void)
{
  //
  // Read in circular velocity data from file. Units
  // are assumed to be kpc (radius) and km/s (velocity)
  //
  std::fstream inFile;
  inFile.open("vcirc.dat", std::ios::in);

  ASSERT("EnzoInitialIsolatedGalaxy",
         "Circular velocity file failed to open",
         inFile.is_open());

  int i = 0;
  while(inFile >> this->vcirc_radius[i] >> this->vcirc_velocity[i])
  {
    ASSERT("EnzoInitialIsolatedGalaxy",
           "Too many lines in circular velocity file",
           i < this->VCIRC_TABLE_LENGTH);

    this->vcirc_radius[i]   *= cello::kpc_cm;   // kpc  -> cm
    this->vcirc_velocity[i] *= 1.0E5; // km/s -> cm/s
    i++;
  }

  inFile.close();

  return;
}

double EnzoInitialIsolatedGalaxy::InterpolateVcircTable(double radius)
{
  //
  // Interpolate the circular velocity from the read-in
  // circular velocity table. Throws error if position is
  // past maximum radius. Takes first bin value if radius is
  // below lower radius limit
  //
  int i;
  double vcirc;

  for (i = 0; i < VCIRC_TABLE_LENGTH; i++){
    if (radius < this->vcirc_radius[i])
      break;
  }

  if (i == 0){
    vcirc = (this->vcirc_velocity[i]) *
            (radius - this->vcirc_radius[0]) / this->vcirc_radius[0];
  } else{
    ASSERT("EnzoInitialIsolatedGalaxy",
           "Radius is outside of circular velocity table",
           i <= VCIRC_TABLE_LENGTH);

    vcirc = this->vcirc_velocity[i-1] +
        (this->vcirc_velocity[i] - this->vcirc_velocity[i-1]) *
        (radius - this->vcirc_radius[i-1])   /
        (this->vcirc_radius[i] - this->vcirc_radius[i-1]);
  }

  return vcirc;
}

double EnzoInitialIsolatedGalaxy::gauss_mass(
             const double rho_zero,
             const double xpos,
             const double ypos,
             const double zpos,
             const double dx)
{

    // Computes the total mass in a given cell by integrating the density
    // profile using 5-point Gaussian quadrature.
    // http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
    double evaluation_points [5] = {-0.90617985,-0.53846931,0.0,0.53846931,0.90617985};
    double weights [5] = {0.23692689,0.47862867,0.56888889,0.47862867,0.23692689};
    double x_result [5];
    double y_result [5];
    double r, z;
    double mass = 0;

    for (int i=0;i<5;i++)
    {
      x_result[i] = 0.0;
      for (int j=0;j<5;j++)
      {
        y_result[j] = 0.0;
        for (int k=0;k<5;k++)
        {
          r = sqrt((pow(xpos+evaluation_points[i]*dx/2.0, 2.0) +
                    pow(ypos+evaluation_points[j]*dx/2.0, 2.0) ) );
          z = fabs(zpos+evaluation_points[k]*dx/2.0);
          y_result[j] +=
            dx/2.0 * weights[k] * rho_zero *
            exp(-r/this->scale_length_) *
            exp(-fabs(z)/this->scale_height_);
        }
        x_result[i] += dx/2.0*weights[j]*y_result[j];
      }
      mass += dx/2.0*weights[i]*x_result[i];
    }
  return mass;
}
