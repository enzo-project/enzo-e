/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodStarMakerSTARSS.cpp
/// @author     Will Hicks (whicks@ucsd.edu)
/// @date
/// @brief  Implements the star maker class using STARSS algorithm developed
//          by Azton Wells (TODO: link paper) 
///
///     Derived star maker class that actually makes stars. This is
///     adapted from the algorithm described in Hopkins et al. 
///     (2018, MNRAS, 480, 800)
//
//      Checks the following star formation criteria by default:
//           1. overdensity > critical_overdensity
//           2. div(V) < 0
//           3. alpha < 0
//           4. T < 10^4 K or t_cool < t_dynamical
//           5. mass[i] > jeans mass
//           6. f_H2_shielded > 0

#include "cello.hpp"
#include "enzo.hpp"
#include <time.h>

// TODO (JHW, 27 Feb 2020)
// Copied from EnzoMethodStarMakerStochasticSF with no changes whatsoever.
// Plan is to make this a separate class that inherits from the stochastic
// algorithm.

//#define DEBUG_SF_CRITERIA
//-------------------------------------------------------------------

EnzoMethodStarMakerSTARSS::EnzoMethodStarMakerSTARSS
()
  : EnzoMethodStarMaker()
{
  // To Do: Make the seed an input parameter
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_fields();
  refresh->add_all_particles();

  ParticleDescr * particle_descr = cello::particle_descr();

  //
  // Refresh copies of all star particles on neighboring grids
  //
  //const int it = particle_descr->type_index("star");
  //refresh->add_particle(it,true);
  //refresh->all_particles_copy(true);
  return;
}

//-------------------------------------------------------------------

void EnzoMethodStarMakerSTARSS::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  EnzoMethodStarMaker::pup(p); // call parent class pup

  return;
}

//------------------------------------------------------------------

void EnzoMethodStarMakerSTARSS::compute ( Block *block) throw()
{

  // Loop through the grid and check star formation criteria
  // stochastically form stars if zone meets these criteria

  std::mt19937 mt(std::time(nullptr));
  int count = 0;

  const EnzoConfig * enzo_config = enzo::config();

  // Are we at the highest level?
  // Can we form stars at this level?
  if ( (! block->is_leaf() ) || (block->level() < enzo_config->method_star_maker_min_level) ){
    block->compute_done();
    return;
  }


  EnzoBlock * enzo_block = enzo::block(block);
  EnzoUnits * enzo_units = enzo::units();
  double lunit = enzo_units->length();
  double tunit = enzo_units->time();
  double vunit = enzo_units->velocity();
  double rhounit = enzo_units->density();
  double munit = enzo_units->mass();
  double munit_solar = munit / cello::mass_solar;
  double Tunit = enzo_units->temperature();


  Particle particle = enzo_block->data()->particle();
  Field field = enzo_block->data()->field();

  double dx, dy, dz;
  block->cell_width(&dx, &dy, &dz);

  const int rank = cello::rank();

  double dt    = enzo_block->dt; // timestep in code units
  enzo_float cosmo_a = 1.0;
  enzo_float cosmo_dadt = 0.0;
  EnzoPhysicsCosmology * cosmology = enzo::cosmology();
/*
  if (cosmology) {
    double current_time = enzo_block->time();
    cosmology->compute_expansion_factor(&cosmo_a,&cosmo_dadt,current_time+0.5*dt);
    if (rank >= 1) dx *= cosmo_a;
    if (rank >= 2) dy *= cosmo_a;
    if (rank >= 3) dz *= cosmo_a;
  }
*/
  double cell_volume = dx*dy*dz;
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

  // obtain the particle stride length (TODO: remove--don't need?)
  //const int ps = particle.stride(it, ia_m);

  int gx,gy,gz;
  field.ghost_depth (0, &gx, &gy, &gz);

  int nx, ny, nz;
  field.size ( &nx, &ny, &nz);

  int mx = nx + 2*gx;
  int my = ny + 2*gy;
  int mz = nz + 2*gz;

  // array increments
  const int idx = 1;
  const int idy = mx;
  const int idz = mx*my;


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

  enzo_float * cooling_time = field.is_field("cooling_time") ?
    (enzo_float *) field.values("cooling_time") : NULL;
  
  enzo_float * density_particle_accumulate = field.is_field("density_particle_accumulate") ?
    (enzo_float *) field.values("density_particle_accumulate") : NULL;


  //const double Zsolar = 0.02;  // TODO: Update to more accurate value

  // Idea for multi-metal species - group these using 'group'
  // class in IC parameter file and in SF / Feedback routines simply
  // check if this group exists, and if it does, loop over all of these
  // fields to assign particle chemical tags and deposit yields

  // compute the temperature (we need it here)
  // TODO: Calling compute_temperature like this
  // returns temperature in Kelvin--not code_temperature??
  // EnzoMethodGrackle::claculate_temperature called
  // without passing in grackle_units or grackle_fields.
  // How does Grackle calculate temperatures? Tabulated?
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

        int i = INDEX(ix,iy,iz,mx,my);//ix + mx*(iy + my*iz);

        // need to compute this better for Grackle fields (on to-do list)
        double rho_cgs = density[i] * enzo_units->density();
        double mean_particle_mass = enzo_config->ppm_mol_weight * cello::mass_hydrogen;
        double ndens = rho_cgs / mean_particle_mass;

        double cell_mass  = density[i] * cell_volume;
        double metallicity = (metal) ? metal[i]/density[i]/cello::metallicity_solar : 0.0;

        //
        // Apply the criteria for star formation
        //
        
        // calculate baryon overdensity
        // take a local mean, but weight the central cell more
        
        // either use number density threshold or overdensity threshold
        // defaults to number density if both density and overdensity are set to 1
        if (this->use_density_threshold_) { 
          if (! this->check_number_density_threshold(ndens) ) continue;
        }
        else if (this->use_overdensity_threshold_){
          double dmean = (10*density[i] 
                    + density[i-idx] + density[i+idx] 
                    + density[i-idy] + density[i+idy]
                    + density[i-idz] + density[i+idz]) / 17.0;
          if (! this->check_overdensity_threshold(dmean) ) continue;
          //#ifdef DEBUG_SF_CRITERIA 
          //  CkPrintf("MethodFeedbackSTARSS -- overdensity threshold passed! rho=%f\n",dmean);
          //#endif
         }
        // check velocity divergence < 0
        if (! this->check_velocity_divergence(velocity_x, velocity_y,
                                              velocity_z, i,
                                              idx, idy, idz)) continue;

        // check that alpha < 0
        if (! this->check_self_gravitating(mean_particle_mass, density[i], temperature[i],
                                            velocity_x, velocity_y, velocity_z,
                                            lunit, vunit, rhounit, Tunit,
                                            i, idx, idy, idz, dx, dy, dz)) continue;

        // check that (T<10^4 K) or (dynamical_time < cooling_time)
        // In order to check cooling time, must have use_temperature_threshold=true;
        double total_density = density[i] + density_particle_accumulate[i];

        if (! this->check_temperature(temperature[i], Tunit)) { // if T > 10^4 K
           if (enzo_config->method_grackle_chemistry) continue; //no hot gas forming stars!
           if (cooling_time){
             if (! this->check_cooling_time(cooling_time[i], total_density, tunit, rhounit)) continue;
           }
        }
        // check that M > Mjeans
        if (! check_jeans_mass(temperature[i], mean_particle_mass, density[i], cell_mass,
                               Tunit,munit,rhounit )) continue;     
        
        // check that H2 self shielded fraction f_shield > 0
        double f_shield = this->h2_self_shielding_factor(density,metallicity,
                     rhounit,lunit,i,idx,idy,idz,dx,dy,dz); 
        if (f_shield < 0) continue;
        // check that Z > Z_crit
        if (! check_metallicity(metallicity)) continue;

//-----------------------------------CREATION ROUTINE-------------------------------
      
        #ifdef DEBUG_SF_CRITERIA
           CkPrintf("MethodStarMakerSTARSS -- SF criteria passed in cell %d\n", i);
        #endif 
        // If cell passed all of the tests, form a star
        double tdyn = sqrt(3.0 * cello::pi / 32.0 / cello::grav_constant /
                      (density[i] * rhounit)) / tunit; //dynamical time in code units (not used anywhere)
 
        //free fall time in code units
        double tff = sqrt(3*cello::pi/(32*cello::grav_constant*density[i]*rhounit))/tunit;        
       /* Determine Mass of new particle
                WARNING: this removes the mass of the formed particle from the
                         host cell.  If your simulation has very small (>15 Msun) baryon mass
                         per cell, it will break your sims! - AIW
       */
        double divisor = std::max(1.0, tff * tunit/cello::Myr_s);
        double maximum_star_mass = enzo_config->method_star_maker_maximum_star_mass;
        double minimum_star_mass = enzo_config->method_star_maker_minimum_star_mass;
         
        if (maximum_star_mass < 0){
            maximum_star_mass = this->maximum_star_fraction_ * cell_mass * munit_solar; //Msun
        }

        double bulk_SFR = f_shield * cell_mass*munit_solar/divisor;
        double mass_should_form = bulk_SFR * dt*tunit/cello::Myr_s; //proposed stellar mass in Msun
        //double mass_should_form = std::min(f_shield/divisor * dt*tunit/cello::Myr_s,
        //                              this->maximum_star_fraction_) * cell_mass*munit_solar;                            
        
        // Probability has the last word
        // FIRE-2 uses p = 1 - exp (-MassShouldForm*dt / M_gas_particle) to convert a whole particle to star particle
        //  We convert a fixed portion of the baryon mass (or the calculated amount)
        //TODO: Difference between dt and dtFixed???
       
        double p_form = 1.0 - std::exp(-mass_should_form/
                (this->maximum_star_fraction_*cell_mass*munit_solar));

        //double p_form = 1.0 - std::exp(-mass_should_form/maximum_star_mass);
        if (enzo_config->method_star_maker_turn_off_probability) p_form = 1.0;

        #ifdef DEBUG_SF_CRITERIA
          CkPrintf("MethodStarMakerSTARSS -- mass_should_form = %f; p_form = %f\n", mass_should_form, p_form);
          CkPrintf("MethodStarMakerSTARSS -- cell_mass = %f Msun; divisor = %f\n", cell_mass*munit_solar,divisor);
          CkPrintf("MethodStarMakerSTARSS -- (ix, iy, iz) = (%d, %d, %d)\n", ix,iy,iz);
        #endif

        double random = double(mt()) / double(mt.max()); 

        /* New star is mass_should_form up to `conversion_fraction` * baryon mass of the cell, but at least 15 msun */       
        //double new_mass = std::min(mass_should_form/munit_solar, maximum_star_mass/munit_solar);
        double new_mass = std::min(mass_should_form/munit_solar, maximum_star_mass/munit_solar);
        if (
                 (new_mass * munit_solar < minimum_star_mass) // too small
                 || (random > p_form) // too unlikely
                 || (new_mass > cell_mass) // too big compared to cell    
           ) 
           {
           #ifdef DEBUG_SF_CRITERIA
             CkPrintf("MethodStarMakerSTARSS -- star mass is either too big, too small, or failed the dice roll\n");
           #endif
             continue;
           }

        count++; // time to form a star! 

        // now create a star particle
        //    insert_particles( particle_type, number_of_particles )
        int my_particle = particle.insert_particles(it, 1);

        // For the inserted particle, obtain the batch number (ib)
        //  and the particle index (ipp)
        particle.index(my_particle, &ib, &ipp);

        int io = ipp; // ipp*ps
        // pointer to mass array in block
        pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);

        pmass[io] = new_mass;
        px = (enzo_float *) particle.attribute_array(it, ia_x, ib);
        py = (enzo_float *) particle.attribute_array(it, ia_y, ib);
        pz = (enzo_float *) particle.attribute_array(it, ia_z, ib);

        // give it position at center of host cell
        px[io] = lx + (ix - gx + 0.5) * dx;
        py[io] = ly + (iy - gy + 0.5) * dy;
        pz[io] = lz + (iz - gz + 0.5) * dz;

        pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
        pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
        pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ib);

        // average particle velocity over many cells to prevent runaway
        double rhosum = 0.0;
        double vx = 0.0;
        double vy = 0.0;
        double vz = 0.0;
        for (int ix_=std::max(0,ix-3); ix_ <= std::min(ix+3,mx); ix_++) {
            for (int iy_=std::max(0,iy-3); iy_ <= std::min(iy+3,my); iy_++) {
                for (int iz_=std::max(0,iz-3); iz_ <= std::min(iz+3,mz); iz_++) {
                    int i_ = INDEX(ix_,iy_,iz_,mx,my);
                    vx += velocity_x[i_]*density[i_];
                    vy += velocity_y[i_]*density[i_];
                    vz += velocity_z[i_]*density[i_];
                    rhosum += density[i_];
                }
            } 
        }
        vx /= rhosum;
        vy /= rhosum;
        vz /= rhosum;

        // TODO: Make this an input parameter
        double max_velocity = 150e5/vunit; 

        pvx[io] = vx;
        pvy[io] = vy;
        pvz[io] = vz;

        // finalize attributes
        plifetime = (enzo_float *) particle.attribute_array(it, ia_l, ib);
        pform     = (enzo_float *) particle.attribute_array(it, ia_to, ib);

        pform[io]     =  enzo_block->time();   // formation time
        //TODO: Need to have some way of calculating lifetime based on particle mass
        plifetime[io] =  25.0 * cello::Myr_s / enzo_units->time() ; // lifetime

        if (metal){
          pmetal     = (enzo_float *) particle.attribute_array(it, ia_metal, ib);
          pmetal[io] = metal[i] / density[i];
        }

        // Remove mass from grid and rescale fraction fields
        double scale = (1.0 - new_mass / cell_mass);
        density[i] *= scale;
        // rescale color fields too 
        this->rescale_densities(enzo_block, i, scale);
      }
    }
  } // end loop iz

  if (count > 0){
    CkPrintf("MethodStarMakerSTARSS -- Number of particles formed: %d\n", count);
  }

// TODO: Set max_number_of_new_particles parameter????? 

/* TODO: Add this part in once Grid_MechStarsSeedSupernova.C is implemented
 
    if (gridShouldFormStars && MechStarsSeedField && (nCreated == 0)){
        // set off a p3 supernova at at the last cell that could 
        // host star formation in the grid if the
        // grid can host star formation but has no appreciable metals
        fprintf(stdout, "\n\n\n[%d] %d %d %d Creating seed field!\n\n\n\n", 
                ID,seedIndex[0], seedIndex[1], seedIndex[2]) ;
        MechStars_SeedSupernova(&totalMetal[0], Temperature, seedIndex);
        
    }
*/

  block->compute_done();

  return;
}
