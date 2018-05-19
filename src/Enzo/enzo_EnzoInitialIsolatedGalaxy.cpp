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
  this->scale_height_         = config->initial_IG_scale_height;
  this->scale_length_	      = config->initial_IG_scale_length;
  this->disk_mass_            = config->initial_IG_disk_mass;
  this->gas_fraction_         = config->initial_IG_gas_fraction;
  this->disk_temperature_     = config->initial_IG_disk_temperature;
  this->disk_metallicity_     = config->initial_IG_disk_metallicity;
  this->gas_halo_mass_        = config->initial_IG_gas_halo_mass;
  this->gas_halo_temperature_ = config->initial_IG_gas_halo_temperature;
  this->gas_halo_metallicity_ = config->initial_IG_gas_halo_metallicity;
  this->gas_halo_density_     = config->initial_IG_gas_halo_density;
  this->gas_halo_radius_      = config->initial_IG_gas_halo_radius;


  this->uniform_density_    = config->field_uniform_density;
  this->dual_energy_        = config->ppm_dual_energy; // or ppm? which one?
  this->gamma_              = config->field_gamma;
  this->mu_                 = config->ppm_mol_weight;

  // do I have to do this this way?
  this->mass_units_         = config->units_mass;
  this->length_units_       = config->units_length;
  this->density_units_      = config->units_density;
  this->time_units_         = config->units_time;

  // read in data for initialization
  this->ReadInVcircData();

  this->ReadInParticleData(); // AE: Does nothing at the moment
}

//----------------------------------------------------------------------

/* I do not know what this does: */
void EnzoInitialIsolatedGalaxy::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP; // ?

  Initial::pup(p); // ?

  PUParray(p,center_position_,3);
  PUParray(p, bfield_,3);
  p | scale_length_;
  p | disk_mass_;
  p | gas_fraction_;
  p | disk_temperature_;
  p | disk_metallicity_;
  p | gas_halo_mass_;
  p | gas_halo_temperature_;
  p | gas_halo_metallicity_;
  p | gas_halo_density_;
  p | gas_halo_radius_;

  p | uniform_density_;
  p | dual_energy_;
  p | gamma_;
  p | mu_;
  // not actually sure what this function does and how to do things here

  p | mass_units_;
  p | length_units_;
  p | density_units_;
  p | time_units_;
}

void EnzoInitialIsolatedGalaxy::enforce_block
(
 Block * block,
 const FieldDescr * field_descr,
 const ParticleDescr * particle_descr,
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
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);

  // Get Fields
  enzo_float * d = (enzo_float *) field.values ("density");
  enzo_float * p = (enzo_float *) field.values ("pressure");
  enzo_float * a3[3] = { (enzo_float *) field.values("acceleration_x"),
                         (enzo_float *) field.values("acceleration_y"),
                         (enzo_float *) field.values("acceleration_z")};
  enzo_float * v3[3] = { (enzo_float *) field.values("velocity_x"),
                         (enzo_float *) field.values("velocity_y"),
                         (enzo_float *) field.values("velocity_z")};

  enzo_float * te  = (enzo_float *) field.values("total_energy");
  enzo_float * ge  = (enzo_float *) field.values("internal_energy");
//  enzo_float * pot = (enzo_float *) field.values("potential");


  double velocity_units = this->length_units_ / this->time_units_;
  double temperature_units = cello::mass_hydrogen *
                             pow(velocity_units, 2) / cello::k;

  //
  // Now lets calculate some physical properties of the galaxy and halo
  //
  double rho_zero = this->disk_mass_ * this->gas_fraction_ / (4.0 * cello::pi)/
      (pow((this->scale_length_),2)*(this->scale_height_));

  double halo_gas_energy = this->gas_halo_temperature_ / this->mu_ / (this->gamma_ -1) /
        temperature_units;

  double disk_gas_energy = this->disk_temperature_ / this->mu_ / (this->gamma_ -1)/
        temperature_units;

  // initialize fields to something
  for (int iz=0; iz<mz; iz++){
    for (int iy=0; iy<my; iy++){
      for (int ix=0; ix<mx; ix++){
        int i = INDEX(ix,iy,iz,mx,my);
        d[i] = this->uniform_density_;
        te[i] = halo_gas_energy;
        p[i]  = (this->gamma_ - 1.0) * te[i] * d[i];

        if (this->dual_energy_)
        {
          ge[i] = halo_gas_energy;
        }

        for (int dim = 0; dim < 3; dim++){
          a3[dim][i] = 0.0;
          v3[dim][i] = 0.0;
        }
      }
    }
  } // end loop over all cells for background values


  //
  int outcount = 0;
  for (int iz=0; iz<mz; iz++){
    // compute z coordinate of cell center
    double z = zm + (iz - gz + 0.5)*hz - this->center_position_[2];
    z *= this->length_units_; // convert to cgs units

    for (int iy=0; iy<my; iy++){
      double y = ym + (iy - gy + 0.5)*hy - this->center_position_[1];
      y *= this->length_units_;

      for (int ix=0; ix<mx; ix++){
        double x = xm + (ix - gx + 0.5)*hx - this->center_position_[0];
        x *= this->length_units_;

        // 1D index of current cell
        int i = INDEX(ix,iy,iz,mx,my);

        // compute spherical and cylindrical radii (in cgs)
        double radius = sqrt(x*x + y*y + z*z);
        double r_cyl  = sqrt(x*x + y*y);

        // compute the disk density (in code units)
        double disk_density = this->gauss_mass(rho_zero, x/this->length_units_, y/this->length_units_,
                                              z/this->length_units_, hx) / (hx*hx*hx);

        if (outcount < 100){
            std::cout << disk_density << "   " << radius << "   " << r_cyl << "\n";
            outcount++;
        }
        if ((this->gas_halo_density_ * this->gas_halo_temperature_ > disk_density*this->disk_temperature_) &&
            (radius < this->gas_halo_radius_)){
          // in the halo, set the halo properties
          d[i]  = this->gas_halo_density_;
          te[i] = halo_gas_energy;
          p[i]  = (this->gamma_ - 1.0) * te[i] * d[i];

          if (this->dual_energy_) {
            ge[i] = halo_gas_energy;
          }

          // set metal fraction here

          //

        }
        else
        {
          // in the disk, set the disk properties
          d[i]   = disk_density;

          double vcirc = this->InterpolateVcircTable(r_cyl);

          v3[0][i] = -(vcirc*(y/r_cyl))/velocity_units;
          v3[1][i] = -(vcirc*(x/r_cyl))/velocity_units;
          v3[2][i] = 0.0;

          te[i]  = disk_gas_energy;

          for (int dim = 0; dim < 3; dim++) // AE: do I need a check for not ppm?
          {
            te[i] += 0.5*v3[dim][i]*v3[dim][i];
          }

          if (this->dual_energy_)
          {
            ge[i] = disk_gas_energy;
          }

        } // end disk / halo check

        /* Set multispecies information here */

        /*                                   */
      }
    }
  } // end loop over all cells

  // do particle initialization here (not sure what this does)
  Particle particle = block->data()->particle();

}

void EnzoInitialIsolatedGalaxy::ReadInParticleData(void)
{
  return;
}

void EnzoInitialIsolatedGalaxy::ReadInVcircData(void)
{
  //
  // Read in circular velocity date from file
  //
  std::fstream inFile;
  inFile.open("vcirc.dat", std::ios::in);

  ASSERT("EnzoInitialIsolatedGalaxy",
         "Circular velocity file failed to open",
         inFile.is_open());

  int i = 0;
  while(!inFile.eof())
  {
    ASSERT("EnzoInitialIsolatedGalaxy",
           "Too many lines in circular velocity file",
           i < this->VCIRC_TABLE_LENGTH);

    inFile >> this->vcirc_radius[i] >> this->vcirc_velocity[i];

    this->vcirc_radius[i]   *= cello::kpc;   // kpc  -> cm
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
