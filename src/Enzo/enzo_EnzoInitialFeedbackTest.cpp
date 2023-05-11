// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialFeedbackTest.cpp
/// @author   Andrew Emerick (emerick@astro.columbia.edu)
/// @date
/// @brief    [\ref Enzo] Initialization routine for Feedback test problem

#include <fstream>

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

// #define DEBUG_INITIALFB

int nlines(std::string fname);


EnzoInitialFeedbackTest::EnzoInitialFeedbackTest
(const EnzoConfig * config) throw ()
  : Initial(config->initial_cycle, config->initial_time)
{

  cello::particle_descr()->check_particle_attribute("star","mass");

  if (config->initial_feedback_test_from_file){
    this->num_particles = nlines("initial_feedback_stars.in");

    for (int dim = 0; dim < 3; dim++){
      this->position[dim].resize(this->num_particles);
    }
    this->mass.resize(this->num_particles);
    this->luminosity.resize(this->num_particles);

    std::fstream inFile;
    inFile.open("initial_feedback_stars.in", std::ios::in);

    ASSERT("EnzoInitialFeedbackTest", "initial_feedback_stars.in failed to open",
           inFile.is_open());

    int i = 0;
    while(inFile >> this->mass[i] >> this->position[0][i] >> this->position[1][i] >> this->position[2][i] >> this->luminosity[i]){
      i++;
    }

    inFile.close();

  } else {
    this->num_particles = 1;
    this->mass.resize(this->num_particles);
    this->luminosity.resize(this->num_particles);

    for (int dim = 0; dim < 3; dim++){
      this->position[dim].resize(this->num_particles);

      this->position[dim][0] = config->initial_feedback_test_position[dim];
    }
    this->mass[0] = config->initial_feedback_test_star_mass;
    this->luminosity[0] = config->initial_feedback_test_luminosity;
  }

  return;
}

//------------------------------------------------------------------------------

void EnzoInitialFeedbackTest::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  p | num_particles;
  for (int dim = 0; dim < 3; dim++) p | position[dim];
  p | mass;

  return;
}

//----------------------------------------------------------------------

void EnzoInitialFeedbackTest::enforce_block
( Block * block, const Hierarchy * hierarchy) throw()
{

  ASSERT("EnzoInitialFeedbackTest","Block does not exist", block != NULL);

  if( !(block->is_leaf())) {
    block->initial_done();
    return;
  }

  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

  Field field = block->data()->field();

  // setup grackle if needed

  enzo_float * d  = (enzo_float *) field.values("density");
  enzo_float * ge = (enzo_float *) field.values("internal_energy");
  enzo_float * te = (enzo_float *) field.values("total_energy");

  enzo_float * d_HI = (enzo_float *) field.values("HI_density");
  enzo_float * d_HII = (enzo_float *) field.values("HII_density");
  enzo_float * d_HeI = (enzo_float *) field.values("HeI_density");
  enzo_float * d_HeII = (enzo_float *) field.values("HeII_density");
  enzo_float * d_HeIII = (enzo_float *) field.values("HeIII_density");
  enzo_float * d_electron = (enzo_float *) field.values("e_density");
  enzo_float * temperature = (enzo_float *) field.values("temperature");

  enzo_float * v3[3] = { (enzo_float *) field.values("velocity_x"),
                         (enzo_float *) field.values("velocity_y"),
                         (enzo_float *) field.values("velocity_z")};

  enzo_float * metal = (enzo_float *) field.values("metal_density");

  // Block size (excluding ghosts)
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  // Cell widths
  double xm,ym,zm;
  block->data()->lower(&xm,&ym,&zm);

  double xp,yp,zp;
  block->data()->upper(&xp,&yp,&zp);

  double hx,hy,hz;
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);

  // Ghost depths
  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;
  int ngz = nz + 2*gz;

  const enzo_float gamma = enzo::fluid_props()->gamma();
  const enzo_float mol_weight = enzo::fluid_props()->mol_weight();

  for (int iz = 0; iz < ngz; iz++){
    for (int iy = 0; iy < ngy; iy++){
      for (int ix = 0; ix < ngx; ix++){

         int i = INDEX(ix,iy,iz,ngx,ngy);

         // values are specified in CGS in the parameter file
         d[i]  = enzo_config->initial_feedback_test_density / enzo_units->density(); 
         d_HI[i]  = enzo_config->initial_feedback_test_HI_density / enzo_units->density();
         d_HII[i]  = enzo_config->initial_feedback_test_HII_density / enzo_units->density();
         d_HeI[i]  = enzo_config->initial_feedback_test_HeI_density / enzo_units->density();
         d_HeII[i]  = enzo_config->initial_feedback_test_HeII_density / enzo_units->density();
         d_HeIII[i]  = enzo_config->initial_feedback_test_HeIII_density / enzo_units->density();
         d_electron[i]  = enzo_config->initial_feedback_test_e_density / enzo_units->density();

         for (int dim = 0; dim < 3; dim++) v3[dim][i] = 0.0;

         ge[i] = (enzo_config->initial_feedback_test_temperature / mol_weight /
                  enzo_units->kelvin_per_energy_units() / (gamma - 1.0));

         for (int dim = 0; dim < 3; dim ++)
             te[i] = ge[i] + 0.5 * v3[dim][i] * v3[dim][i];

         metal[i] = enzo_config->initial_feedback_test_metal_fraction * d[i];

      }
    }
  }

  // drop in a particle

  ParticleDescr * particle_descr = cello::particle_descr();
  Particle particle              = block->data()->particle();

  int it = particle_descr->type_index("star");

  int ia_m = particle.attribute_index (it, "mass");
  int ia_x = particle.attribute_index (it, "x");
  int ia_y = particle.attribute_index (it, "y");
  int ia_z = particle.attribute_index (it, "z");
  int ia_vx = particle.attribute_index (it, "vx");
  int ia_vy = particle.attribute_index (it, "vy");
  int ia_vz = particle.attribute_index (it, "vz");
  int ia_cop  = particle.attribute_index (it, "is_copy");
  int ia_id   = particle.attribute_index (it, "id");
  
  int ia_L     = particle.has_attribute (it, "luminosity") ?
                 particle.attribute_index(it,"luminosity") : -1;

  int ia_to    = particle.has_attribute(it,"creation_time") ?
                 particle.attribute_index(it,"creation_time") : -1;

  int ia_l     = particle.has_attribute(it,"lifetime") ?
                 particle.attribute_index(it,"lifetime") : -1;

  int ia_metal = particle.has_attribute(it,"metal_fraction") ?
                 particle.attribute_index(it,"metal_fraction") : -1;

  int ib  = 0; // batch counter
  int ipp = 0; // particle counter

  // this will point to the particular value in the
  // particle attribute array
  enzo_float * pmass = 0;
  enzo_float * plum = 0;
  enzo_float * px   = 0;
  enzo_float * py   = 0;
  enzo_float * pz   = 0;
  enzo_float * pvx  = 0;
  enzo_float * pvy  = 0;
  enzo_float * pvz  = 0;
  enzo_float * pmetal = 0;
  enzo_float * plifetime = 0;
  enzo_float * pform     = 0;
  int64_t * is_copy = 0;
  int64_t * id = 0;


  // from global particle list, figure out which ones go on this block
  int * mask = new int [this->num_particles];
  int num_local_particles = 0;
  for (int i = 0; i < this->num_particles; i++){
    mask[i] =  block->check_position_in_block(this->position[0][i],
                                              this->position[1][i],
                                              this->position[2][i]);
    if(mask[i]) num_local_particles++;
  }


  int new_particle = particle.insert_particles(it, num_local_particles);
  particle.index(new_particle,&ib,&ipp);

  id   = (int64_t *) particle.attribute_array(it, ia_id, ib);
  pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);
  plum  = (enzo_float *) particle.attribute_array(it, ia_L, ib);
  px    = (enzo_float *) particle.attribute_array(it, ia_x, ib);
  py    = (enzo_float *) particle.attribute_array(it, ia_y, ib);
  pz    = (enzo_float *) particle.attribute_array(it, ia_z, ib);
  pvx   = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
  pvy   = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
  pvz   = (enzo_float *) particle.attribute_array(it, ia_vz, ib);
  pmetal      = (enzo_float *) particle.attribute_array(it, ia_metal, ib);
  plifetime  = (enzo_float *) particle.attribute_array(it, ia_l, ib);
  pform      = (enzo_float *) particle.attribute_array(it, ia_to, ib);
  is_copy   = (int64_t *) particle.attribute_array(it, ia_cop, ib);

  ipp = 0;
  for (int i = 0; i < this->num_particles; i++){
    if (mask[i]){

#ifdef DEBUG_INITIALFB
      CkPrintf("EnzoInitialFeedbackTest (%s): Mass MyPe id_counter istatic id %f %d %d %d %d\n", enzo_block->name().c_str(),
                                                                        this->mass[i],
                                                                        cello::index_static(), ParticleData::id_counter[cello::index_static()],
                                                                        CkNumPes(), CkMyPe() + (ParticleData::id_counter[cello::index_static()]) * CkNumPes());
#endif
      id[ipp] = CkMyPe() + (ParticleData::id_counter[cello::index_static()]++) * CkNumPes();
      pmass[ipp] = this->mass[i] * enzo_constants::mass_solar / enzo_units->mass();
      px[ipp]    = this->position[0][i];
      py[ipp]    = this->position[1][i];
      pz[ipp]    = this->position[2][i];
      plum[ipp]  = this->luminosity[i] * enzo_units->time(); // luminosity in photons/time
      pvx[ipp]   = 0.0;
      pvy[ipp]   = 0.0;
      pvz[ipp]   = 0.0;

      pmetal[ipp]    = 0.01;
      plifetime[ipp] = 1.00E9* enzo_constants::yr_s / enzo_units->time();
      pform[ipp]     = 1.0E-10 * enzo_constants::yr_s / enzo_units->time(); // really just needs to be non-zero

      is_copy[ipp] = 1;

      ipp++;
    }
  }

  block->initial_done();

  return;
}
