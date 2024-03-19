/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodStarMakerSTARSS.cpp
/// @author     Will Hicks (whicks@ucsd.edu)
/// @date
/// @brief  Implements the star maker class using STARSS algorithm developed
//          by Azton Wells (TODO: link paper when it's published) 
///
///     Derived star maker class that actually makes stars. This is
///     adapted from the algorithm described in Hopkins et al. 
///     (2018, MNRAS, 480, 800)

#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"
#include "Enzo/particle/particle.hpp"

#include <time.h>

// #define DEBUG_SF_CRITERIA
// #define DEBUG_STORE_INITIAL_PROPERTIES
//-------------------------------------------------------------------

EnzoMethodStarMakerSTARSS::EnzoMethodStarMakerSTARSS
()
  : EnzoMethodStarMaker()
{
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_fields();
  refresh->add_all_particles();

  ParticleDescr * particle_descr = cello::particle_descr();

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

  // TODO: Make random seed an input parameter
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
  double munit_solar = munit / enzo_constants::mass_solar;
  double eunit = vunit*vunit; // specific energy units

  Particle particle = enzo_block->data()->particle();
  Field field = enzo_block->data()->field();

  double dx, dy, dz;
  block->cell_width(&dx, &dy, &dz);

  const int rank = cello::rank();

  auto dt = enzo_block->state()->dt(); // timestep in code units
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

  #ifdef DEBUG_STORE_INITIAL_PROPERTIES 
    const int ia_m_0 = particle.attribute_index (it, "mass_0");
    const int ia_x_0 = particle.attribute_index (it, "x_0");
    const int ia_y_0 = particle.attribute_index (it, "y_0");
    const int ia_z_0 = particle.attribute_index (it, "z_0");
    const int ia_vx_0 = particle.attribute_index (it, "vx_0");
    const int ia_vy_0 = particle.attribute_index (it, "vy_0");
    const int ia_vz_0 = particle.attribute_index (it, "vz_0");
  #endif

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
  
  enzo_float * potential    = (enzo_float *) field.values("potential");

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


  if (use_altAlpha_) {
    // "potential" field holds potential in proper coordinates,
    //  need to convert back to comoving. 
    if (cosmology) {
      enzo_float cosmo_a = 1.0;
      enzo_float cosmo_dadt = 0.0;
      double dt   = enzo_block->state()->dt();
      double time = enzo_block->state()->time();
      cosmology-> compute_expansion_factor (&cosmo_a,&cosmo_dadt,time+0.5*dt);
      for (int i=0; i<mx*my*mz; i++) potential[i] *= cosmo_a;
    }
  }

  if (use_overdensity_threshold_) {
    ASSERT("EnzoMethodStarMakerSTARSS::compute()",
           "parameter use_overdensity_threshold is only valid for cosmology!",
            cosmology);
  }

  // compute the temperature
  EnzoComputeTemperature compute_temperature(enzo::fluid_props(),
                                             enzo_config->physics_cosmology);

  compute_temperature.compute(enzo_block);

  const GrackleChemistryData * grackle_chem = enzo::grackle_chemistry();
  const int primordial_chemistry = (grackle_chem == nullptr) ?
    0 : grackle_chem->get<int>("primordial_chemistry");

  const double dflt_mu = static_cast<double>(enzo::fluid_props()->mol_weight());
  // iterate over all cells (not including ghost zones)
  for (int iz=gz; iz<nz+gz; iz++){
    for (int iy=gy; iy<ny+gy; iy++){
      for (int ix=gx; ix<nx+gx; ix++){

        int i = INDEX(ix,iy,iz,mx,my);//ix + mx*(iy + my*iz);

        double mu;
        // compute MMW -- TODO: Make EnzoComputeMeanMolecularWeight class and reference
        // mu_field here?
        if (primordial_chemistry > 0) {
          // use species fields to get number density times mass_Hydrogen
          // (note: "e_density" field tracks ndens_electron * mass_Hydrogen)
          double ndens_times_mH
            =  d_el[i] + dHI[i] + dHII[i] + 0.25*(dHeI[i]+dHeII[i]+dHeIII[i]);

          if (primordial_chemistry > 1) {
            ndens_times_mH += dHM[i] + 0.5*(dH2I[i]+dH2II[i]);
          }
          if (primordial_chemistry > 2) {
            ndens_times_mH += 0.5*(dDI[i] + dDII[i]) + dHDI[i]/3.0;
          }
          // MA: NOTE previous versions of the code did NOT include the effect
          //     of metals on the mmw. To include it, uncomment next line:
          // if (metal) { ndens_times_mH += metal[i]/16.0; }

          mu = density[i] / ndens_times_mH;
        } else {
          mu = dflt_mu;
        }

        double rho_cgs = density[i] * enzo_units->density();
        double mean_particle_mass = mu * enzo_constants::mass_hydrogen;
        double ndens = rho_cgs / mean_particle_mass;

        double cell_mass  = density[i] * cell_volume; // code units
        double metallicity = (metal) ? metal[i]/density[i]/enzo_constants::metallicity_solar : 0.0;

        //
        // Apply the criteria for star formation
        //
        
        // check number density threshold
        if (! this->check_number_density_threshold(ndens) ) continue;

        // check overdensity threshold
        // In cosmology, units are scaled such that mean(density) = 1,
        // so density IS overdensity in these units  

        if (! this->check_overdensity_threshold(density[i]) ) continue;

        #ifdef DEBUG_SF_CRITERIA 
           CkPrintf("MethodFeedbackSTARSS -- overdensity threshold passed! rho=%f\n",dmean);
        #endif
        
        // check velocity divergence < 0
        if (! this->check_velocity_divergence(velocity_x, velocity_y,
                                              velocity_z, i,
                                              idx, idy, idz, dx, dy, dz)) continue;
      
        #ifdef DEBUG_SF_CRITERIA
           CkPrintf("MethodStarMakerSTARSS -- div(v) < 0 in cell %d\n", i);
        #endif 

        // check that alpha < 1
        if (use_altAlpha_) {
          if (! this->check_self_gravitating_alt(total_energy[i], potential[i])) continue;
        }

        else {
          if (! this->check_self_gravitating(mean_particle_mass, density[i], temperature[i],
                                             velocity_x, velocity_y, velocity_z,
                                             lunit, vunit, rhounit,
                                             i, idx, idy, idz, dx, dy, dz)) continue;
        }

        #ifdef DEBUG_SF_CRITERIA
           CkPrintf("MethodStarMakerSTARSS -- alpha < 1 in cell %d\n", i);
        #endif 

        // check that (T<Tcrit) or (dynamical_time < cooling_time)
        // In order to check cooling time, must have use_temperature_threshold=true;
        double total_density = density[i] + density_particle_accumulate[i];

        // if T > Tcrit
        if (! this->check_temperature(temperature[i])) continue;

        if (cooling_time){ // if we are evolving a "cooling_time" field
           if (! this->check_cooling_time(cooling_time[i], total_density, tunit, rhounit)) continue;
        }
        
        // check that M > Mjeans
        if (! check_jeans_mass(temperature[i], mean_particle_mass, density[i], cell_mass,
                               munit,rhounit )) continue;

        #ifdef DEBUG_SF_CRITERIA
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
        double tff = sqrt(3*cello::pi/(32*enzo::grav_constant_cgs()*density[i]*rhounit))/tunit;        
       /* Determine Mass of new particle
                WARNING: this removes the mass of the formed particle from the
                         host cell.  If your simulation has very small (>15 Msun) baryon mass
                         per cell, it will break your sims! - AIW
       */
        double divisor = std::max(1.0, tff * tunit/enzo_constants::Myr_s);
        double maximum_star_mass = enzo_config->method_star_maker_maximum_star_mass;
        double minimum_star_mass = enzo_config->method_star_maker_minimum_star_mass;
         
        if (maximum_star_mass < 0){
            maximum_star_mass = this->maximum_star_fraction_ * cell_mass * munit_solar; //Msun
        }

        double bulk_SFR = f_shield * this->maximum_star_fraction_ * cell_mass*munit_solar/divisor;
        
        // Probability has the last word
        // FIRE-2 uses p = 1 - exp (-MassShouldForm*dt / M_gas_particle) to convert a whole particle to star particle
        //  We convert a fixed portion of the baryon mass (or the calculated amount)
       
        double p_form = 1.0 - std::exp(-bulk_SFR*dt*(tunit/enzo_constants::Myr_s)/
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
          // spread over three dynamical times. This is to prevent huge particles suddenly dumping a HUGE
          // amount of ionizing radiation at once.
          // TODO: Enforce that if dt < 3 dynamical times, new star particles form in the next timestep
          //       Can do this by creating the particle now so we can still access it at the next timestep,
          //       but not assigning it any attributes until it's actually supposed to form
         
          n_newStars = std::floor(new_mass * munit_solar / mass_split);
          massPerStar = new_mass / n_newStars;
        #ifdef DEBUG_SF_CRITERIA
          CkPrintf("MethodStarMakerSTARSS -- Predicted cluster mass %1.3e Msun > %1.3e Msun;\n" 
                   "                         splitting into %d particles with mass %1.3e Msun\n",
                                             new_mass*munit_solar, max_massPerStar, n_newStars, massPerStar*munit_solar);
        #endif
        } 
          for (int n=0; n<n_newStars; n++) {
            double ctime = enzo_block->state()->time(); 
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
          if (std::abs(vx) > max_velocity) vx = vx/std::abs(vx) * max_velocity; 
          if (std::abs(vy) > max_velocity) vy = vy/std::abs(vy) * max_velocity;
          if (std::abs(vz) > max_velocity) vz = vz/std::abs(vz) * max_velocity;

          pvx[io] = vx;
          pvy[io] = vy;
          pvz[io] = vz;

          // finalize attributes
          plifetime = (enzo_float *) particle.attribute_array(it, ia_l, ib);
          pform     = (enzo_float *) particle.attribute_array(it, ia_to, ib);
          plevel    = (enzo_float *) particle.attribute_array(it, ia_lev, ib);

          pform[io]     =  ctime;   // formation time

          //TODO: Need to have some way of calculating lifetime based on particle mass
          plifetime[io] =  25.0 * enzo_constants::Myr_s / enzo_units->time() ; // lifetime (not accessed for STARSS FB)

          plevel[io] = enzo_block->level(); // formation level

          if (metal){
            pmetal     = (enzo_float *) particle.attribute_array(it, ia_metal, ib);
            pmetal[io] = metal[i] / density[i]; // in ABSOLUTE units
          }

          // Remove mass from grid and rescale fraction fields
          // TODO: If particle position is updated to CM instead of being cell-centered, will have to 
          //       remove mass using CiC, which could complicate things because that CiC cloud could
          //       leak into the ghost zones. Would have to use same refresh+accumulate machinery
          //       as MethodFeedbackSTARSS to account for this.
          double scale = (1.0 - pmass[io] / cell_mass);
          density[i] *= scale;
          // rescale color fields too 
          this->rescale_densities(enzo_block, i, scale);

          #ifdef DEBUG_STORE_INITIAL_PROPERTIES 
            enzo_float * pmass0 = (enzo_float *) particle.attribute_array(it, ia_m_0 , ib);
            enzo_float * px0    = (enzo_float *) particle.attribute_array(it, ia_x_0 , ib);
            enzo_float * py0    = (enzo_float *) particle.attribute_array(it, ia_y_0 , ib);
            enzo_float * pz0    = (enzo_float *) particle.attribute_array(it, ia_z_0 , ib);
            enzo_float * pvx0   = (enzo_float *) particle.attribute_array(it, ia_vx_0, ib);
            enzo_float * pvy0   = (enzo_float *) particle.attribute_array(it, ia_vy_0, ib);
            enzo_float * pvz0   = (enzo_float *) particle.attribute_array(it, ia_vz_0, ib);

            pmass0[io] = pmass[io];
            px0 [io] = px[io];
            py0 [io] = py[io];
            pz0 [io] = pz[io];
            pvx0[io] = pvx[io];
            pvy0[io] = pvy[io];
            pvz0[io] = pvz[io];
          #endif

        } // end loop through particles created in this cell


      }
    }
  } // end loop iz

  #ifdef DEBUG_SF_CRITERIA
    if (count > 0){
      CkPrintf("MethodStarMakerSTARSS -- Number of particles formed: %d\n", count);
    }
  #endif

/* TODO: Add this part in once Pop III SF/FB is implemented
 
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
