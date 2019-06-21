// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialFeedbackTest.cpp
/// @author   Andrew Emerick (emerick@astro.columbia.edu)
/// @date
/// @brief    [\ref Enzo] Initialization routine for Feedback test problem

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoInitialFeedbackTest::EnzoInitialFeedbackTest
(const EnzoConfig * config) throw ()
  : Initial(config->initial_cycle, config->initial_time)
{

  return;
}

//------------------------------------------------------------------------------

void EnzoInitialFeedbackTest::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  return;
}

void EnzoInitialFeedbackTest::enforce_block
( Block * block, const Hierarchy * hierarchy) throw()
{

  ASSERT("EnzoInitialFeedbackTest","Block does not exist", block != NULL);

  if( !(block->is_leaf())) return;

  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

  Field field = block->data()->field();

  // setup grackle if needed

  enzo_float * d  = (enzo_float *) field.values("density");
  enzo_float * ge = (enzo_float *) field.values("internal_energy");
  enzo_float * te = (enzo_float *) field.values("total_energy");

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

  for (int iz = 0; iz < ngz; iz++){
    for (int iy = 0; iy < ngy; iy++){
      for (int ix = 0; ix < ngx; ix++){

         int i = INDEX(ix,iy,iz,ngx,ngy);

         d[i]  = enzo_config->initial_feedback_test_density / enzo_units->density();

         for (int dim = 0; dim < 3; dim++) v3[dim][i] = 0.0;

         ge[i] = 1.0E4 / enzo_config->ppm_mol_weight / enzo_units->temperature() /
                         (enzo_config->field_gamma - 1.0);

         for (int dim = 0; dim < 3; dim ++)
             te[i] = ge[i] + 0.5 * v3[dim][i] * v3[dim][i];

         metal[i] = 0.01 * d[i]; // half solar

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

  int ia_to    = particle.is_attribute(it,"creation_time") ?
                 particle.attribute_index(it,"creation_time") : -1;

  int ia_l     = particle.is_attribute(it,"lifetime") ?
                 particle.attribute_index(it,"lifetime") : -1;

  int ia_metal = particle.is_attribute(it,"metal_fraction") ?
                 particle.attribute_index(it,"metal_fraction") : -1;

  int ib  = 0; // batch counter
  int ipp = 0; // particle counter

  // this will point to the particular value in the
  // particle attribute array
  enzo_float * pmass = 0;
  enzo_float * px   = 0;
  enzo_float * py   = 0;
  enzo_float * pz   = 0;
  enzo_float * pvx  = 0;
  enzo_float * pvy  = 0;
  enzo_float * pvz  = 0;
  enzo_float * pmetal = 0;
  enzo_float * plifetime = 0;
  enzo_float * pform     = 0;

  // just one particle for now

  if (   block->check_position_in_block(enzo_config->initial_feedback_test_position[0],
                                        enzo_config->initial_feedback_test_position[1],
                                        enzo_config->initial_feedback_test_position[2]) ){
    int new_particle = particle.insert_particles(it, 1);
    particle.index(new_particle,&ib,&ipp);

    pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);
    px    = (enzo_float *) particle.attribute_array(it, ia_x, ib);
    py    = (enzo_float *) particle.attribute_array(it, ia_y, ib);
    pz    = (enzo_float *) particle.attribute_array(it, ia_z, ib);
    pvx   = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
    pvy   = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
    pvz   = (enzo_float *) particle.attribute_array(it, ia_vz, ib);
    pmetal      = (enzo_float *) particle.attribute_array(it, ia_metal, ib);
    plifetime  = (enzo_float *) particle.attribute_array(it, ia_l, ib);
    pform      = (enzo_float *) particle.attribute_array(it, ia_to, ib);

    pmass[ipp] = enzo_config->initial_feedback_test_star_mass * cello::mass_solar / enzo_units->mass();
    px[ipp]    = enzo_config->initial_feedback_test_position[0];
    py[ipp]    = enzo_config->initial_feedback_test_position[1];
    pz[ipp]    = enzo_config->initial_feedback_test_position[2];
    pvx[ipp]   = 0.0;
    pvy[ipp]   = 0.0;
    pvz[ipp]   = 0.0;

    pmetal[ipp]    = 0.01;
    plifetime[ipp] = 1.00E9* cello::yr_s / enzo_units->time();
    pform[ipp]     = 1.0E-10 * cello::yr_s / enzo_units->time(); // really just needs to be non-zero


  }

  return;
}
