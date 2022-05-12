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

#define DEBUG_SF_CRITERIA
// #define DEBUG_SF_CRITERIA_EXTRA
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
  double eunit = vunit*vunit; // specific energy units

  Particle particle = enzo_block->data()->particle();
  Field field = enzo_block->data()->field();

  double dx, dy, dz;
  block->cell_width(&dx, &dy, &dz);

  const int rank = cello::rank();

  double dt    = enzo_block->dt; // timestep in code units
  enzo_float cosmo_a = 1.0;
  enzo_float cosmo_dadt = 0.0;
  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  double cell_volume = dx*dy*dz;
  double lx, ly, lz;
  block->lower(&lx,&ly,&lz);

  // declare particle position arrays
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
  const int ia_lev   = particle.attribute_index (it, "creation_level");

  int ib  = 0; // batch counter
  int ipp = 0; // counter

  /// pointers for particle attribute arrays (later)
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
  enzo_float * plevel = 0;

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

  enzo_float * total_energy = (enzo_float *) field.values("total_energy");
  //enzo_float * potential    = (enzo_float *) field.values("potential");

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

  enzo_float * dHI    = field.is_field("HI_density") ? 
          (enzo_float*) field.values("HI_density") : NULL;
  enzo_float * dHII   = field.is_field("HII_density") ? 
          (enzo_float*) field.values("HII_density") : NULL;
  enzo_float * dHeI   = field.is_field("HeI_density") ? 
          (enzo_float*) field.values("HeI_density") : NULL;
  enzo_float * dHeII  = field.is_field("HeII_density") ? 
          (enzo_float*) field.values("HeII_density") : NULL;
  enzo_float * dHeIII = field.is_field("HeIII_density") ? 
          (enzo_float*) field.values("HeIII_density") : NULL;
  enzo_float * d_el   = field.is_field("e_density") ?
          (enzo_float*) field.values("e_density") : NULL;
 
  enzo_float * dH2I   = field.is_field("H2I_density") ? 
          (enzo_float*) field.values("H2I_density") : NULL;
  enzo_float * dH2II  = field.is_field("H2II_density") ? 
          (enzo_float*) field.values("H2II_density") : NULL;
  enzo_float * dHM    = field.is_field("HM_density") ? 
          (enzo_float*) field.values("HM_density") : NULL;

  enzo_float * dDI    = field.is_field("DI_density") ? 
         (enzo_float *) field.values("DI_density") : NULL;
  enzo_float * dDII   = field.is_field("DII_density") ? 
         (enzo_float *) field.values("DII_density") : NULL;
  enzo_float * dHDI   = field.is_field("HDI_density") ? 
         (enzo_float *) field.values("HDI_density") : NULL;



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

  chemistry_data * grackle_chemistry =
    enzo::config()->method_grackle_chemistry;

  int primordial_chemistry = grackle_chemistry->primordial_chemistry;
  double mu = enzo_config->ppm_mol_weight;
  // iterate over all cells (not including ghost zones)
  for (int iz=gz; iz<nz+gz; iz++){
    for (int iy=gy; iy<ny+gy; iy++){
      for (int ix=gx; ix<nx+gx; ix++){

        int i = INDEX(ix,iy,iz,mx,my);//ix + mx*(iy + my*iz);
        
     /*
        // compute MMW -- TODO: Make EnzoComputeMeanMolecularWeight class and reference
        // mu_field here
        if (primordial_chemistry == 0) mu = enzo_config->ppm_mol_weight;
        else {
          mu = d_el[i] + dHI[i] + dHII[i] + 0.25*(dHeI[i]+dHeII[i]+dHeIII[i]);

          if (primordial_chemistry > 1) {
            mu += dHM[i] + 0.5*(dH2I[i]+dH2II[i]);
          }
          if (primordial_chemistry > 2) {
            mu += 0.5*(dDI[i] + dDII[i]) + dHDI[i]/3.0;
          }
        }
        mu /= density[i]; 
        mu = 1.0/mu;
    */

        // TODO: need to compute this better for Grackle fields (on to-do list)
        double rho_cgs = density[i] * enzo_units->density();
        double mean_particle_mass = mu * cello::mass_hydrogen;
        double ndens = rho_cgs / mean_particle_mass;

        double cell_mass  = density[i] * cell_volume;
        double metallicity = (metal) ? metal[i]/density[i]/cello::metallicity_solar : 0.0;

        //
        // Apply the criteria for star formation
        //
        
        // either use number density threshold or overdensity threshold
        // defaults to number density if both density and overdensity are set to 1
        if (this->use_density_threshold_) { 
          if (! this->check_number_density_threshold(ndens) ) continue;
        }
        else if (this->use_overdensity_threshold_){
          if (! this->check_overdensity_threshold(density[i]) ) continue;
          //#ifdef DEBUG_SF_CRITERIA 
          //  CkPrintf("MethodFeedbackSTARSS -- overdensity threshold passed! rho=%f\n",dmean);
          //#endif
         }
        // check velocity divergence < 0
        if (! this->check_velocity_divergence(velocity_x, velocity_y,
                                              velocity_z, i,
                                              idx, idy, idz)) continue;
      
        #ifdef DEBUG_SF_CRITERIA_EXTRA
           CkPrintf("MethodStarMakerSTARSS -- div(v) < 0 in cell %d\n", i);
        #endif 

        // check that alpha < 1
        if (enzo_config->method_star_maker_use_altAlpha) {
          // approximate grav potential of cell, taking r = dx
          // NOTE: potential field is currently cleared at the end of EnzoMethodGravity.
          //       Can just access it here if we decide to either not clear or 
          //       make a copy (DEBUG_COPY_POTENTIAL).
          double potential_i = cello::grav_constant * cell_mass*munit / (dx*lunit) / eunit;
          if (! this->check_self_gravitating_new(total_energy[i], potential_i)) continue;
        }

        else {
          if (! this->check_self_gravitating(mean_particle_mass, density[i], temperature[i],
                                             velocity_x, velocity_y, velocity_z,
                                             lunit, vunit, rhounit,
                                             i, idx, idy, idz, dx, dy, dz)) continue;
        }
        #ifdef DEBUG_SF_CRITERIA_EXTRA
           CkPrintf("MethodStarMakerSTARSS -- alpha < 1 in cell %d\n", i);
        #endif 

        // check that (T<Tcrit) or (dynamical_time < cooling_time)
        // In order to check cooling time, must have use_temperature_threshold=true;
        double total_density = density[i] + density_particle_accumulate[i];

        if (! this->check_temperature(temperature[i])) { // if T > Tcrit
           if (enzo_config->method_grackle_chemistry) continue; //no hot gas forming stars!
           if (cooling_time){
             if (! this->check_cooling_time(cooling_time[i], total_density, tunit, rhounit)) continue;
           }
        }
        // check that M > Mjeans
        if (! check_jeans_mass(temperature[i], mean_particle_mass, density[i], cell_mass,
                               munit,rhounit )) continue;

        #ifdef DEBUG_SF_CRITERIA_EXTRA
           CkPrintf("MethodStarMakerSTARSS -- M > M_jeans in cell %d\n", i);
        #endif     
        
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

        double bulk_SFR = f_shield * this->maximum_star_fraction_ * cell_mass*munit_solar/divisor;
        
        // Probability has the last word
        // FIRE-2 uses p = 1 - exp (-MassShouldForm*dt / M_gas_particle) to convert a whole particle to star particle
        //  We convert a fixed portion of the baryon mass (or the calculated amount)
        //TODO: Difference between dt and dtFixed???
       
        double p_form = 1.0 - std::exp(-bulk_SFR*dt*(tunit/cello::Myr_s)/
                (this->maximum_star_fraction_*cell_mass*munit_solar));

        if (enzo_config->method_star_maker_turn_off_probability) p_form = 1.0;

        double random = double(mt()) / double(mt.max()); 

        /* New star is mass_should_form up to f_shield*maximum_star_fraction_ * baryon mass of the cell,
           but at least 15 msun */       
        double new_mass = std::min(f_shield*maximum_star_fraction_*cell_mass, maximum_star_mass/munit_solar);

        #ifdef DEBUG_SF_CRITERIA
          CkPrintf("MethodStarMakerSTARSS -- new_mass = %f; p_form = %f\n", new_mass*munit_solar, p_form);
          CkPrintf("MethodStarMakerSTARSS -- cell_mass = %f Msun; divisor = %f\n", cell_mass*munit_solar,divisor);
          CkPrintf("MethodStarMakerSTARSS -- (ix, iy, iz) = (%d, %d, %d)\n", ix,iy,iz);
        #endif

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

        int n_newStars = 1; // track how many stars to form
        double mFirmed = 0.0; // track total mass formed thus far
        double massPerStar = new_mass;
        double max_massPerStar = 5e3; // TODO: either make this a new parameter or replace 'maximum_star_mass'
        double mass_split = 1e3; 
        if (new_mass * munit_solar > max_massPerStar) { //
          // for large particles, split them into several 1e3-ish Msun particles 
          // that have slightly different birth times,
          // spread over three dynamical times
          // TODO: Enforce that if dt < 3 dynamical times, new star particles form in the next timestep
          //       Can do this by creating the particle so we can still access it at the next timestep,
          //       but not assigning it a mass until it's actually supposed to form
         
          n_newStars = std::floor(new_mass * munit_solar / mass_split);
          massPerStar = new_mass / n_newStars;
        #ifdef DEBUG_SF_CRITERIA
          CkPrintf("MethodStarMakerSTARSS -- Predicted cluster mass %1.3e Msun > %1.3e Msun;\n" 
                   "                         splitting into %d particles with mass %1.3e Msun\n",
                                             new_mass*munit_solar, max_massPerStar, n_newStars, massPerStar*munit_solar);
        #endif
        } 
          for (int n=0; n<n_newStars; n++) {
            double ctime = enzo_block->time(); 
            if (n > 0) {
              double mod = n * 3.0 * tff/tunit/n_newStars;
              ctime += mod;
            }
            count++; // time to form a star!
      
          #ifdef DEBUG_SF_CRITERIA
            CkPrintf("MethodStarMakerSTARSS -- Forming star in gas with number density %f cm^-3\n", ndens);
          #endif

          // now create a star particle
          int my_particle = particle.insert_particles(it, 1);

          // For the inserted particle, obtain the batch number (ib)
          // and the particle index (ipp)
          particle.index(my_particle, &ib, &ipp);

          int io = ipp; // ipp*ps
          // pointer to mass array in block
          pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);

          // TODO: saving particle mass as density. May need to update this in the future
          //       when PR #89 passes 
          //pmass[io] = new_mass / cell_volume;

          pmass[io] = massPerStar;
          px = (enzo_float *) particle.attribute_array(it, ia_x, ib);
          py = (enzo_float *) particle.attribute_array(it, ia_y, ib);
          pz = (enzo_float *) particle.attribute_array(it, ia_z, ib);

          // give it position at center of host cell
          // TODO: Calculate CM instead?
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
          plevel    = (enzo_float *) particle.attribute_array(it, ia_lev, ib);

          pform[io]     =  ctime;   // formation time

          //TODO: Need to have some way of calculating lifetime based on particle mass
          plifetime[io] =  25.0 * cello::Myr_s / enzo_units->time() ; // lifetime (not accessed for STARSS FB)

          plevel[io] = enzo_block->level(); // formation level

          if (metal){
            pmetal     = (enzo_float *) particle.attribute_array(it, ia_metal, ib);
            pmetal[io] = metal[i] / density[i]; // in ABSOLUTE units
          }

          // Remove mass from grid and rescale fraction fields
          // TODO: Remove mass using CiC instead? Should do this if particle's initial position
          //       isn't cell-centered
          double scale = (1.0 - pmass[io] / cell_mass);
          density[i] *= scale;
          // rescale color fields too 
          this->rescale_densities(enzo_block, i, scale);

        } // end loop through particles created in this cell


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
