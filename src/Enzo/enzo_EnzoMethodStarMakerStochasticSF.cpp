/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodStarMakerStochasticSF.cpp
/// @author     Andrew Emerick (aemerick11@gmail.com)
/// @date
/// @brief  Implements the star maker stochastic star formation clas
///
///     Derived star maker class that actually makes stars. This is
///     adapted after the star_maker_ssn method from Enzo

#include "cello.hpp"
#include "enzo.hpp"
#include <time.h>

//-------------------------------------------------------------------

EnzoMethodStarMakerStochasticSF::EnzoMethodStarMakerStochasticSF
()
  : EnzoMethodStarMaker()
{
  // To Do: Make the seed an input parameter
  srand(time(NULL)); // need randum number generator for later
  return;
}

//-------------------------------------------------------------------

void EnzoMethodStarMakerStochasticSF::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  EnzoMethodStarMaker::pup(p); // call parent class pup

  return;
}

//------------------------------------------------------------------

void EnzoMethodStarMakerStochasticSF::compute ( Block *block) throw()
{

  // Loop through the grid and check star formation criteria
  // stochastically form stars if zone meets these criteria

  int count = 0;

  // Are we at the highest level?
  if (! block->is_leaf()){
    block->compute_done();
    return;
  }

  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();


  Particle particle = enzo_block->data()->particle();
  Field field = enzo_block->data()->field();

  double dx, dy, dz;
  block->cell_width(&dx, &dy, &dz);

  double lx, ly, lz;
  block->lower(&lx,&ly,&lz);

  // declare particle position arrays
  //  default particle type is "star", but this will default
  //  to subclass particle_type
  const int it   = particle.type_index (this->particle_type());

  const int ia_m = particle.attribute_index (it, "mass");
  const int ia_x = particle.attribute_index (it, "x");
  const int ia_y = particle.attribute_index (it, "y");
  const int ia_z = particle.attribute_index (it, "z");
  const int ia_vx = particle.attribute_index (it, "vx");
  const int ia_vy = particle.attribute_index (it, "vy");
  const int ia_vz = particle.attribute_index (it, "vz");

  // additional particle attributes
  const int ia_metal = particle.attribute_index (it, "metal_fraction");
  const int ia_to    = particle.attribute_index (it, "creation_time");
  const int ia_l     = particle.attribute_index (it, "lifetime");

  int ib  = 0; // batch counter
  int ipp = 0; // counter

  /// pointers for particle attribut arrays (later)
  enzo_float * pmass = 0;
  enzo_float * px   = 0;
  enzo_float * py   = 0;
  enzo_float * pz   = 0;
  enzo_float * pvx  = 0;
  enzo_float * pvy  = 0;
  enzo_float * pvz  = 0;
   ///
  enzo_float * pmetal = 0;
  enzo_float * pform  = 0;
  enzo_float * plifetime = 0;

  // obtain the particle stride length
  const int ps = particle.stride(it, ia_m);

  int rank = cello::rank();

  int gx,gy,gz;
  field.ghost_depth (0, &gx, &gy, &gz);

  int nx, ny, nz;
  field.size ( &nx, &ny, &nz);

  int mx = nx + 2*gx;
  int my = ny + 2*gy;
  int mz = nz + 2*gz;


  enzo_float * density     = (enzo_float *) field.values("density");
  enzo_float * temperature = (enzo_float *) field.values("temperature");

  enzo_float * velocity_x = (rank >= 1) ?
    (enzo_float *)field.values("velocity_x") : NULL;
  enzo_float * velocity_y = (rank >= 2) ?
    (enzo_float *)field.values("velocity_y") : NULL;
  enzo_float * velocity_z = (rank >= 3) ?
    (enzo_float *)field.values("velocity_z") : NULL;

  enzo_float * metal = field.is_field("metal_density") ?
    (enzo_float *) field.values("metal_density") : NULL;

  const double Zsolar = 0.02;  // TODO: Update to more accurate value

  // Idea for multi-metal species - group these using 'group'
  // class in IC parameter file and in SF / Feedback routines simply
  // check if this group exists, and if it does, loop over all of these
  // fields to assign particle chemical tags and deposit yields

  // compute the temperature (we need it here)
  EnzoComputeTemperature compute_temperature
    (enzo_config->ppm_density_floor,
     enzo_config->ppm_temperature_floor,
     enzo_config->ppm_mol_weight,
     enzo_config->physics_cosmology);

  compute_temperature.compute(enzo_block);

  // iterate over all cells (not including ghost zones)
  //
  //   To Do: Allow for multi-zone star formation by adding mass in
  //          surrounding cells if needed to accumulte enough mass
  //          to hit target star particle mass ()
  for (int iz=gz; iz<nz+gz; iz++){
    for (int iy=gy; iy<ny+gy; iy++){
      for (int ix=gx; ix<nx+gx; ix++){

        int i = ix + mx*(iy + my*iz);

        // need to compute this better for Grackle fields (on to-do list)
        double rho_cgs = density[i] * enzo_units->density();
        double mean_particle_mass = enzo_config->ppm_mol_weight * cello::mass_hydrogen;
        double ndens = rho_cgs / mean_particle_mass;

        double mass  = density[i] *dx*dy*dz * enzo_units->mass() / cello::mass_solar;
        double metallicity = (metal) ? metal[i]/density[i]/Zsolar : 0.0;

        //
        // Apply the criteria for star formation
        //
        if (! this->check_number_density_threshold(ndens)) continue;
        if (! this->check_self_gravitating(mean_particle_mass, rho_cgs, temperature[i],
                                           enzo_units->length(), enzo_units->density(),
                                           velocity_x, velocity_y, velocity_z, 
                                           i, 1, my, my*mz, dx)) continue;
        if (! this->check_h2_self_shielding(density, metallicity, i, 1, my, my*mz, dx)) continue;
        if (! this->check_velocity_divergence(velocity_x, velocity_y,
                                              velocity_z, i,
                                              1, my, my*mz)) continue;
        // Check whether mass in [min_mass, max_range] range and if specified, Jeans unstable
        if (! this->check_mass(mass)) continue;

        double tdyn = sqrt(3.0 * cello::pi / 32.0 / cello::grav_constant /
                      (density[i] * enzo_units->density()));

        //
        // compute fraction that can / will be converted to stars this step
        // (just set to efficiency if dynamical time is ignored)
        //
        double star_fraction =  this->use_dynamical_time_ ?
                                std::min(this->efficiency_ * enzo_block->dt * enzo_units->time() / tdyn, 1.0) :
                                         this->efficiency_ ;

        // if this is less than the mass of a single particle,
        // use a random number draw to generate the particle
        if ( star_fraction * mass < this->star_particle_mass_){
          // get a random number
          double rnum = (double(rand())) / (double(RAND_MAX));
          double probability = this->efficiency_*mass/this->star_particle_mass_;
          if (rnum > probability){
              continue; // do not form stars
          } else{
            star_fraction = this->star_particle_mass_ / mass;
          }
        } else {
          // (else, leave star fraction alone )
          star_fraction = this->star_particle_mass_ / mass;
        }

        count++; //

        // now create a star particle
        //    insert_particles( particle_type, number_of_particles )
        int my_particle = particle.insert_particles(it, 1);

        // For the inserted particle, obtain the batch number (ib)
        //  and the particle index (ipp)
        particle.index(my_particle, &ib, &ipp);

        int io = ipp; // ipp*ps
        // pointer to mass array in block
        pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);

        pmass[io] = star_fraction * (density[i] * dx * dy * dz);
        px = (enzo_float *) particle.attribute_array(it, ia_x, ib);
        py = (enzo_float *) particle.attribute_array(it, ia_y, ib);
        pz = (enzo_float *) particle.attribute_array(it, ia_z, ib);

        // need to double check that these are correctly handling ghost zones
        //   I believe lx is lower coordinates of active region, but
        //   ix is integer index of whole grid (active + ghost)
        //
        px[io] = lx + (ix - gx + 0.5) * dx;
        py[io] = ly + (iy - gy + 0.5) * dy;
        pz[io] = lz + (iz - gz + 0.5) * dz;

        pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
        pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
        pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ib);

        pvx[io] = velocity_x[i];
        if (velocity_y) pvy[io] = velocity_y[i];
        if (velocity_z) pvz[io] = velocity_z[i];

        // finalize attributes
        plifetime = (enzo_float *) particle.attribute_array(it, ia_l, ib);
        pform     = (enzo_float *) particle.attribute_array(it, ia_to, ib);

        pform[io]     =  enzo_block->time();   // formation time
        plifetime[io] =  10.0 * cello::Myr_s / enzo_units->time() ; // lifetime

        if (metal){
          pmetal     = (enzo_float *) particle.attribute_array(it, ia_metal, ib);
          pmetal[io] = metal[i] / density[i];
        }

        // Remove mass from grid and rescale fraction fields
        density[i] = (1.0 - star_fraction) * density[i];
        double scale = (1.0 - star_fraction) / 1.0;

        // rescale tracer fields to maintain constant mass fraction
        // with the corresponding new density...
        //    scale = new_density / old_density
        rescale_densities(enzo_block, i, scale);
      }
    }
  } // end loop iz

  if (count > 0){
      std::cout << "Number of particles formed:   " << count << "\n";
  }


  block->compute_done();

  return;
}

/*
   Defaults to parent class timestep if nothing declared here
double EnzoMethodStarMakerStochasticSF::timestep ( Block *block) const throw()
{
  return std::numeric_limits<double>::max();
}
*/
