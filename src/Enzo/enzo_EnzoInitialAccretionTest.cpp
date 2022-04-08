// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialShuCollapse.cpp
/// @author   Stefan Arridge (stefan.arridge@gmail.com)
/// @date     2022-03-09
/// @brief    Implementation of EnzoInitialAccretionTest, an initializer for an
///           accretion test problem, which puts an accreting star particle with
///           a given initial position and velocity in a static medium.
///


#include "cello.hpp"
#include "enzo.hpp"

EnzoInitialAccretionTest::EnzoInitialAccretionTest
  (int cycle, double time,
   const double star_position[3],
   const double star_velocity[3],
   double star_mass,
   double gas_density,
   double gas_pressure
   ) throw()
    : Initial(cycle,time),
      star_mass_(star_mass),
      gas_density_(gas_density),
      gas_pressure_(gas_pressure)
{
  star_position_[0] = star_position[0];
  star_position_[1] = star_position[1];
  star_position_[2] = star_position[2];
  star_velocity_[0] = star_velocity[0];
  star_velocity_[1] = star_velocity[1];
  star_velocity_[2] = star_velocity[2];

}

void EnzoInitialAccretionTest::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  PUParray(p,star_position_,3);
  PUParray(p,star_velocity_,3);
  p | star_mass_;
  p | gas_density_;
  p | gas_pressure_;
}


void EnzoInitialAccretionTest::enforce_block
( Block * block, const Hierarchy * hierarchy ) throw()

{

  enzo::check_particle_attribute("star","mass");

  // Check if accretion_compute and accretion_remove_gas methods are being used,
  // and that accretion_compute precedes accretion_remove_gas
  ASSERT("EnzoInitialAccretionTest",
	 "If accretion_test initializer is used, the accretion_compute "
	 "and accretion_remove_gas methods are required, and "
	 "accretion_compute must precede accretion_remove_gas.",
         enzo::problem()->method_precedes("accretion_compute",
					  "accretion_remove_gas"));

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
  enzo_float *  d = (enzo_float *) field.values ("density");
  enzo_float *  da = (enzo_float *) field.values ("density_accreted");
  enzo_float *  daa = (enzo_float *) field.values ("density_accreted_accumulate");
  enzo_float *  p = (enzo_float *) field.values ("pressure");
  enzo_float * te = (enzo_float *) field.values ("total_energy");
  enzo_float * ie = (enzo_float *) field.values ("internal_energy");
  enzo_float * vx = (enzo_float *) field.values ("velocity_x");
  enzo_float * vy = (enzo_float *) field.values ("velocity_y");
  enzo_float * vz = (enzo_float *) field.values ("velocity_z");

  // Initialise all fields to zero except for density and total_energy
  std::fill_n(p,m,0.0);
  std::fill_n(ie,m,0.0);
  std::fill_n(vx,m,0.0);
  std::fill_n(vy,m,0.0);
  std::fill_n(vz,m,0.0);
  std::fill_n(da,m,0.0);
  std::fill_n(daa,m,0.0);

  // Set density
  std::fill_n(d,m,gas_density_);

  const enzo_float gamma = EnzoBlock::Gamma[cello::index_static()];
  
  // Set total energy and internal energy
  const enzo_float te_value = gas_pressure_ / ((gamma - 1.0) * gas_density_);
  std::fill_n(te,m,te_value);

  // Create star particle if its position is in the block
  if (block->check_position_in_block(star_position_[0],
				     star_position_[1],
				     star_position_[2]))
    {
      ParticleDescr * particle_descr = cello::particle_descr();
      Particle particle              = block->data()->particle();

      // Attribute indices
      const int it   = particle.type_index("star");
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
      pmass[0] = star_mass_;
      px[0] = star_position_[0];
      py[0] = star_position_[1];
      pz[0] = star_position_[2];
      pvx[0] = star_velocity_[0];
      pvy[0] = star_velocity_[1];
      pvz[0] = star_velocity_[2];
      id[0] = 1;
      is_copy[0] = 0;
      
    } // Is there a particle to place in this block?
  return;
}
