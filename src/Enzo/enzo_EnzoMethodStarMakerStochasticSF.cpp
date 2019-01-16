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

  Method::pup(p);

}

//------------------------------------------------------------------

void EnzoMethodStarMakerStochasticSF::compute ( Block *block) throw()
{


  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

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

    double * pmass = 0;
    double * px   = 0;
    double * py   = 0;
    double * pz   = 0;
    double * pvx  = 0;
    double * pvy  = 0;
    double * pvz  = 0;

    //std::mt19937 generator(12345);
    //std::uniform_real_distribution<double> generate_rnum(0.0, 1.0);

    // obtain the particle stride length
    const int ps = particle.stride(it, ia_m);

    int rank = cello::rank();

    enzo_float * density     = (enzo_float *) field.values("density");
    enzo_float * velocity_x = (rank >= 1) ?
      (enzo_float *)field.values("velocity_x") : NULL;
    enzo_float * velocity_y = (rank >= 2) ?
      (enzo_float *)field.values("velocity_y") : NULL;
    enzo_float * velocity_z = (rank >= 3) ?
      (enzo_float *)field.values("velocity_z") : NULL;

    int gx,gy,gz;
    field.ghost_depth (0, &gx, &gy, &gz);

    int mx, my, mz;
    field.dimensions (0, &mx, &my, &mz);

    int nx, ny, nz;
    field.size ( &nx, &ny, &nz);

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

          // compute SF
          double star_fraction = 0.0;
          if ( this->efficiency_*mass < this->star_particle_mass_){
            // get a random number
            double rnum = (double(rand())) / (double(RAND_MAX));
            double probability = this->efficiency_*mass/this->star_particle_mass_;
            if (rnum > probability){
                continue; // do not form stars
            } else{
              star_fraction = this->star_particle_mass_ / mass;
            }
          }

          // now create a star particle
          int my_particle = particle.insert_particles(it, 1);

          //
          particle.index(my_particle, &ib, &ipp);

          // pointer to mass array in block
          pmass = (double *) particle.attribute_array(it, ia_m, ib);

          pmass[ipp*ps] = star_fraction * density[i];
          px = (double *) particle.attribute_array(it, ia_x, ib);
          py = (double *) particle.attribute_array(it, ia_y, ib);
          pz = (double *) particle.attribute_array(it, ia_z, ib);

          px[ipp*ps] = lx + (ix + 0.5) * dx;
          py[ipp*ps] = ly + (iy + 0.5) * dy;
          pz[ipp*ps] = lz + (iz + 0.5) * dz;

          pvx = (double *) particle.attribute_array(it, ia_vx, ib);
          pvy = (double *) particle.attribute_array(it, ia_vy, ib);
          pvz = (double *) particle.attribute_array(it, ia_vz, ib);

          pvx[ipp*ps] = velocity_x[i];
          if (velocity_y) pvy[ipp*ps] = velocity_y[i];
          if (velocity_z) pvz[ipp*ps] = velocity_z[i];

        }
      }
    } // end loop iz

  } // end leaf check

  enzo_block->compute_done();

  return;
}
/*
double EnzoMethodStarMakerStochasticSF::timestep ( Block *block) const throw()
{
  return std::numeric_limits<double>::max();
}
*/
