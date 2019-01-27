/// See LICENSE_CELLO file for license and copyright information

///
///
///
///

#include "cello.hpp"
#include "enzo.hpp"
#include <time.h>

//-------------------------------------------------------------------

EnzoMethodStarMakerStochasticSF::EnzoMethodStarMakerStochasticSF
()
  : EnzoMethodStarMaker()
{
  srand(time(NULL));
  return;
}

//-------------------------------------------------------------------

void EnzoMethodStarMakerStochasticSF::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  EnzoMethodStarMaker::pup(p); // call parent class pup

}

//------------------------------------------------------------------

void EnzoMethodStarMakerStochasticSF::compute ( Block *block) throw()
{


  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

  int count = 0;

  // Are we at the highest level?
  if (block->is_leaf()) {

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

    int ib  = 0; // batch counter
    int ipp = 0; // counter

    /// should these be set to double????
    enzo_float * pmass = 0;
    enzo_float * px   = 0;
    enzo_float * py   = 0;
    enzo_float * pz   = 0;
    enzo_float * pvx  = 0;
    enzo_float * pvy  = 0;
    enzo_float * pvz  = 0;

    //std::mt19937 generator(12345);
    //std::uniform_real_distribution<double> generate_rnum(0.0, 1.0);

    // obtain the particle stride length
    const int ps = particle.stride(it, ia_m);

    int rank = cello::rank();

    int gx,gy,gz;
    field.ghost_depth (0, &gx, &gy, &gz);

    int mx, my, mz;
    field.dimensions (0, &mx, &my, &mz);

    int nx, ny, nz;
    field.size ( &nx, &ny, &nz);


    enzo_float * density     = (enzo_float *) field.values("density");
    enzo_float * temperature = (enzo_float *) field.values("temperature");

    enzo_float * velocity_x = (rank >= 1) ?
      (enzo_float *)field.values("velocity_x") : NULL;
    enzo_float * velocity_y = (rank >= 2) ?
      (enzo_float *)field.values("velocity_y") : NULL;
    enzo_float * velocity_z = (rank >= 3) ?
      (enzo_float *)field.values("velocity_z") : NULL;

    // compute the temperature
    EnzoComputeTemperature compute_temperature
      (enzo_config->ppm_density_floor,
       enzo_config->ppm_temperature_floor,
       enzo_config->ppm_mol_weight,
       enzo_config->physics_cosmology);

    compute_temperature.compute(enzo_block);

    // iterate over all cells
    for (int iz=gz; iz<nz+gz; iz++){
      for (int iy=gy; iy<ny+gy; iy++){
        for (int ix=gx; ix<nx+gx; ix++){

          int i = ix + mx*(iy + my*iz);


          double ndens = density[i] * enzo_units->density() /
                         (enzo_config->ppm_mol_weight * cello::mass_hydrogen);

          double mass  = density[i] *dx*dy*dz * enzo_units->mass() / cello::mass_solar;

          // Apply the criteria for star formation
          if (! this->check_number_density_threshold(ndens)) continue;


          if (! this->check_velocity_divergence(velocity_x, velocity_y,
                                               velocity_z, i,
                                               1, my, my*mz)) continue;


          if (! this->check_minimum_mass(mass)) continue;

          double tdyn = sqrt(3.0 * cello::pi / 32.0 / cello::grav_constant / (density[i] * enzo_units->density()));

//          double isothermal_sound_speed2 = 1.3095E8 * temperature[i] * enzo_units->temperature();

          // compute SF
          double star_fraction =  use_dynamical_time_ ?
                                         std::min(this->efficiency_ * enzo_block->dt * enzo_units->time() / tdyn, 1.0) :
                                         this->efficiency_ ;

          if ( star_fraction * mass < this->star_particle_mass_){
            // get a random number
            double rnum = (double(rand())) / (double(RAND_MAX));
            double probability = this->efficiency_*mass/this->star_particle_mass_;
            if (rnum > probability){
                continue; // do not form stars
            } else{
              star_fraction = this->star_particle_mass_ / mass;
            }
          }

          count++;
          // now create a star particle
          //    insert_particles( particle_type, number_of_particles )
          int my_particle = particle.insert_particles(it, 1);

          // For the inserted particle, obtain the batch number (ib)
          //  and the particle index (ipp)
          particle.index(my_particle, &ib, &ipp);

          int io = ipp; // ipp*ps
          // pointer to mass array in block
          pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);

          pmass[io] = star_fraction * density[i];
          px = (enzo_float *) particle.attribute_array(it, ia_x, ib);
          py = (enzo_float *) particle.attribute_array(it, ia_y, ib);
          pz = (enzo_float *) particle.attribute_array(it, ia_z, ib);

          px[io] = lx + (ix + 0.5) * dx;
          py[io] = ly + (iy + 0.5) * dy;
          pz[io] = lz + (iz + 0.5) * dz;

          pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
          pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
          pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ib);

          pvx[io] = velocity_x[i];
          if (velocity_y) pvy[io] = velocity_y[i];
          if (velocity_z) pvz[io] = velocity_z[i];

          // Adjust density and density dependent fields

          density[i] = (1.0 - star_fraction) * density[i];
          double scale = 1.0 / (1.0 - star_fraction);

          // rescale tracer fields to maintain constant mass fraction
          // with the corresponding new density
          rescale_densities(enzo_block, i, scale);
        }
      }
    } // end loop iz

  } // end leaf check

  if (count > 0){
      std::cout << "Number of particles formed:   " << count << "\n";
  }

  enzo_block->compute_done();

  return;
}
/*
double EnzoMethodStarMakerStochasticSF::timestep ( Block *block) const throw()
{
  return std::numeric_limits<double>::max();
}
*/
