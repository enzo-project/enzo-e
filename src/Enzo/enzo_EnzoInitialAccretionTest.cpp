// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialShuCollapse.cpp
/// @author   Stefan Arridge (stefan.arridge@gmail.com)
/// @date     2022-03-09
/// @brief    Implementation of EnzoInitialAccretionTest, an initializer for an
///           accretion test problem, which puts an accreting sink particle with
///           a given initial position and velocity in a static medium.
///


#include "cello.hpp"
#include "enzo.hpp"

EnzoInitialAccretionTest::EnzoInitialAccretionTest
  (int cycle, double time,
   const double sink_position[3],
   const double sink_velocity[3],
   double sink_mass,
   double gas_density,
   double gas_pressure
   ) throw()
    : Initial(cycle,time),
      sink_mass_(sink_mass),
      gas_density_(gas_density),
      gas_pressure_(gas_pressure)
{
  sink_position_[0] = sink_position[0];
  sink_position_[1] = sink_position[1];
  sink_position_[2] = sink_position[2];
  sink_velocity_[0] = sink_velocity[0];
  sink_velocity_[1] = sink_velocity[1];
  sink_velocity_[2] = sink_velocity[2];

}

void EnzoInitialAccretionTest::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  PUParray(p,sink_position_,3);
  PUParray(p,sink_velocity_,3);
  p | sink_mass_;
  p | gas_density_;
  p | gas_pressure_;
}


void EnzoInitialAccretionTest::enforce_block
( Block * block, const Hierarchy * hierarchy ) throw()

{

  enzo::check_particle_attribute("sink","mass");

  // Check if accretion method is being used
    ASSERT("EnzoInitialAccretionTest",
	   "If accretion_test initializer is used, the accretion "
	   "method is required.",
	   enzo::problem()->method_exists("accretion"));

  // Check if mhd_vlct method is being used
  ASSERT("EnzoInitialAccretionTest",
	 "If accretion_test initializer is used, the mhd_vlct method is "
	 "required.",
         enzo::problem()->method_exists("mhd_vlct"));

  // Check that mhd_choice parameter is set to "no_bfield"
  ASSERT("EnzoInitialAccretionTest",
	 "Method:mhd_vlct:mhd_choice must be set to no_bfield",
         enzo::config()->method_vlct_mhd_choice == "no_bfield");

  // Check that riemann_solver parameter is set to "hllc"
  ASSERT("EnzoInitialAccretionTest",
	 "Method:mhd_vlct:mhd_choice must be set to hllc",
         enzo::config()->method_vlct_riemann_solver == "hllc");

  // Check if the initial density and pressure are at least as large as the
  // density and pressure floors set by mhd_vlct method
  ASSERT("EnzoInitialAccretionTest",
	 "Initial gas density must be at least as large as the density "
	 "floor set by the ppm method",
         gas_density_ >= enzo::config()->method_vlct_density_floor);

  ASSERT("EnzoInitialAccretionTest",
	 "Initial gas pressure must be at least as large as the pressure "
	 "floor set by the ppm method",
         gas_pressure_ >= enzo::config()->method_vlct_pressure_floor);

  // Check if we have periodic boundary conditions
  int px, py, pz;
  cello::hierarchy()->get_periodicity(&px,&py,&pz);

  ASSERT("EnzoInitialAccretionTest::EnzoInitialAccretionTest() ",
	 "accretion_test requires periodic boundary conditions.",
         px && py && pz);
  if (!block->is_leaf()) return;
  const EnzoConfig * enzo_config = enzo::config();
  ASSERT("EnzoInitialAccretionTest",
	 "Block does not exist",
	 block != NULL);
  
  // TODO: Check required fields?

  Field field = block->data()->field();

  // Get Field parameters
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);
  const int mx = nx + 2*gx;
  const int my = ny + 2*gy;
  const int mz = nz + 2*gz;
  const int m = mx*my*mz;

  // Block extents
  double bxm,bym,bzm;
  double bxp,byp,bzp;
  block->data()->lower(&bxm,&bym,&bzm);
  block->data()->upper(&bxp,&byp,&bzp);
  double hx,hy,hz;
  field.cell_width(bxm,bxp,&hx,
		   bym,byp,&hy,
		   bzm,bzp,&hz);

  // Get pointers to fields
  enzo_float *  d   = (enzo_float *) field.values ("density");
  enzo_float *  p   = (enzo_float *) field.values ("pressure");
  enzo_float * te   = (enzo_float *) field.values ("total_energy");
  enzo_float * ie   = (enzo_float *) field.values ("internal_energy");
  enzo_float * vx   = (enzo_float *) field.values ("velocity_x");
  enzo_float * vy   = (enzo_float *) field.values ("velocity_y");
  enzo_float * vz   = (enzo_float *) field.values ("velocity_z");

  enzo_float *  ds   = (enzo_float *) field.values ("density_source");
  enzo_float *  dsa  = (enzo_float *) field.values ("density_source_accumulate");
  enzo_float *  mxs  = (enzo_float *) field.values ("mom_dens_x_source");
  enzo_float *  mxsa = (enzo_float *) field.values ("mom_dens_x_source_accumulate");
  enzo_float *  mys  = (enzo_float *) field.values ("mom_dens_y_source");
  enzo_float *  mysa = (enzo_float *) field.values ("mom_dens_y_source_accumulate");
  enzo_float *  mzs  = (enzo_float *) field.values ("mom_dens_z_source");
  enzo_float *  mzsa = (enzo_float *) field.values ("mom_dens_z_source_accumulate");
  enzo_float *  tes  = (enzo_float *) field.values ("te_dens_source");
  enzo_float *  tesa = (enzo_float *) field.values ("te_dens_source_accumulate");

  // Initialise all fields to zero except for density and total_energy
  std::fill_n(p,m,0.0);
  std::fill_n(ie,m,0.0);
  std::fill_n(vx,m,0.0);
  std::fill_n(vy,m,0.0);
  std::fill_n(vz,m,0.0);
  std::fill_n(ds,m,0.0);
  std::fill_n(dsa,m,0.0);
  std::fill_n(mxs,m,0.0);
  std::fill_n(mxsa,m,0.0);
  std::fill_n(mys,m,0.0);
  std::fill_n(mysa,m,0.0);
  std::fill_n(mzs,m,0.0);
  std::fill_n(mzsa,m,0.0);
  std::fill_n(tes,m,0.0);
  std::fill_n(tesa,m,0.0);

  // Set density
  std::fill_n(d,m,gas_density_);

  const enzo_float gamma = EnzoBlock::Gamma[cello::index_static()];
  
  // Set total energy and internal energy
  const enzo_float te_value = gas_pressure_ / ((gamma - 1.0) * gas_density_);
  std::fill_n(te,m,te_value);

  // Create sink particle if its position is in the block
  if (block->check_position_in_block(sink_position_[0],
				     sink_position_[1],
				     sink_position_[2]))
    {
      ParticleDescr * particle_descr = cello::particle_descr();
      Particle particle              = block->data()->particle();

      // Attribute indices
      const int it   = particle.type_index("sink");
      const int ia_m = particle.attribute_index (it, "mass");
      const int ia_x = particle.attribute_index (it, "x");
      const int ia_y = particle.attribute_index (it, "y");
      const int ia_z = particle.attribute_index (it, "z");
      const int ia_vx = particle.attribute_index (it, "vx");
      const int ia_vy = particle.attribute_index (it, "vy");
      const int ia_vz = particle.attribute_index (it, "vz");
      int ia_copy  = particle.attribute_index (it, "is_copy");
      int ia_id   = particle.attribute_index (it, "id");

      /// Initialise pointers for particle attribute arrays
      enzo_float * pmass = 0;
      enzo_float * px   = 0;
      enzo_float * py   = 0;
      enzo_float * pz   = 0;
      enzo_float * pvx  = 0;
      enzo_float * pvy  = 0;
      enzo_float * pvz  = 0;
      int64_t * is_copy = 0;
      int64_t * id = 0;
     
      // insert particle
      int ib,ipp  = 0;
      const int new_particle = particle.insert_particles(it, 1);
      enzo::simulation()->data_insert_particles(1);
      particle.index(new_particle,&ib,&ipp);

      // Get pointers to particle attribute arrays
      pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);
      px    = (enzo_float *) particle.attribute_array(it, ia_x, ib);
      py    = (enzo_float *) particle.attribute_array(it, ia_y, ib);
      pz    = (enzo_float *) particle.attribute_array(it, ia_z, ib);
      pvx   = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
      pvy   = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
      pvz   = (enzo_float *) particle.attribute_array(it, ia_vz, ib);
      id   = (int64_t *) particle.attribute_array(it, ia_id, ib);
      is_copy   = (int64_t *) particle.attribute_array(it, ia_copy, ib);
      
      // Now assign values to attributes
      pmass[0] = sink_mass_;
      px[0] = sink_position_[0];
      py[0] = sink_position_[1];
      pz[0] = sink_position_[2];
      pvx[0] = sink_velocity_[0];
      pvy[0] = sink_velocity_[1];
      pvz[0] = sink_velocity_[2];
      id[0] = 1;
      is_copy[0] = 0;
      
    } // Is there a particle to place in this block?
  return;
}
