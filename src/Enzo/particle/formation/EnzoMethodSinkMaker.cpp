/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodSinkMaker.cpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @date       1 June 2022
/// @brief      Implements a method for forming sink particles based on
///             the method described in Krumholz et al 2004, ApJ, 611, 399 and
///             Federrath et al 2010, ApJ, 713, 269.

#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"
#include "Enzo/particle/particle.hpp"

#include <random>

//-------------------------------------------------------------------

EnzoMethodSinkMaker::EnzoMethodSinkMaker(ParameterGroup p) noexcept
  : Method(),
    jeans_length_resolution_cells_
      (p.value_float("jeans_length_resolution_cells", 4.0)),
    physical_density_threshold_cgs_
      (p.value_float("physical_density_threshold_cgs",1.0e-24)),
    check_density_maximum_(p.value_logical("check_density_maximum",true)),
    max_mass_fraction_(p.value_float("max_mass_fraction", 0.25)),
    min_sink_mass_solar_(p.value_float("min_sink_mass_solar",0.0)),
    max_offset_cell_fraction_(p.value_float("max_offset_cell_fraction",0.0)),
    offset_seed_shift_(0)
{
  int tmp_offset_seed_shift_input = p.value_integer("offset_seed_shift",0);
  ASSERT("EnzoMethodSinkMaker",
         "Method:sink_maker:offset_seed_shift must be >=0",
	 tmp_offset_seed_shift_input >= 0);
  offset_seed_shift_ = static_cast<uint64_t>(tmp_offset_seed_shift_input);

  // Check that we have 3 spatial dimensions.
  ASSERT("EnzoMethodSinkMaker::EnzoMethodSinkMaker()",
	 "Cannot use sink_maker method with rank < 3",
	 cello::rank() == 3);

  // Check if we are running in unigrid mode (will get rid of this in future)
  ASSERT("EnzoMethodSinkMaker::EnzoMethodSinkMaker()",
	 "EnzoMethodSinkMaker requires unigrid mode (Adapt : max_level = 0). "
	 "In future, we may put in a refinement condition that blocks containing "
	 "sink particles are at the highest refinement level, as well as a Jeans "
	 "length refinement condition.",
	 enzo::config()->mesh_max_level == 0);

  // Check the max_mass_fraction is between 0 and 1.
  ASSERT("EnzoMethodSinkMaker::EnzoMethodSinkMaker()",
	 "Method:sink_maker:max_mass_fraction must be between 0 and 1 inclusive.",
	  max_mass_fraction_ >= 0.0 && max_mass_fraction_ <= 1.0);

  // Check that the minimum sink mass is non-negative
  ASSERT("EnzoMethodSinkMaker::EnzoMethodSinkMaker()",
	 "Method:sink_maker:min_sink_mass_solar must be non-negative.",
	 min_sink_mass_solar_ >= 0.0);

  // Check that the Jeans length threshold is non-negative
  ASSERT("EnzoMethodSinkMaker::EnzoMethodSinkMaker()",
	 "Method:sink_maker:jeans_length_resolution_cells must be non-negative.",
	 jeans_length_resolution_cells_ >= 0);

  // Check that max_offset_cell_fraction_ is between 0.0 and 0.1 (inclusive)
  ASSERT("EnzoMethodSinkMaker::EnzoMethodSinkMaker()",
	 "Method:sink_maker:max_offset_cell_fraction must be between 0.0 and "
	 "0.1 (inclusive).",
	 (max_offset_cell_fraction_ >= 0.0) && (max_offset_cell_fraction_ <= 0.1));

  cello::simulation()->refresh_set_name(ir_post_,name());

  Refresh * refresh = cello::refresh(ir_post_);
  ParticleDescr * particle_descr = cello::particle_descr();
  refresh->add_all_fields();

}

//-------------------------------------------------------------------

void EnzoMethodSinkMaker::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);
  p | jeans_length_resolution_cells_;
  p | min_sink_mass_solar_;
  p | physical_density_threshold_cgs_;
  p | check_density_maximum_;
  p | max_mass_fraction_;
  p | min_sink_mass_solar_;
  p | max_offset_cell_fraction_;
  p | offset_seed_shift_;

  return;
}

//-------------------------------------------------------------------------------------

void EnzoMethodSinkMaker::compute ( Block *block) throw()
{
  if (cello::is_initial_cycle(InitCycleKind::fresh_or_noncharm_restart)) {
    do_checks_(block);
  }

  // Only call compute_ if block is on maximum refinement level.
  if (block->level() == enzo::config()->mesh_max_level){
    this->compute_(block);
  }
  block->compute_done();

  return;
}

//-------------------------------------------------------------------------------------

void EnzoMethodSinkMaker::compute_ ( Block *block) throw()

{
  // Get particle data
  Particle particle = block->data()->particle();

  const int it = particle.type_index ("sink");

  const int ia_m                = particle.attribute_index (it, "mass");
  const int ia_x                = particle.attribute_index (it, "x");
  const int ia_y                = particle.attribute_index (it, "y");
  const int ia_z                = particle.attribute_index (it, "z");
  const int ia_vx               = particle.attribute_index (it, "vx");
  const int ia_vy               = particle.attribute_index (it, "vy");
  const int ia_vz               = particle.attribute_index (it, "vz");
  const int ia_metal_fraction   = particle.attribute_index (it, "metal_fraction");
  const int ia_creation_time    = particle.attribute_index (it, "creation_time");
  const int ia_id               = particle.attribute_index (it, "id");
  const int ia_copy             = particle.attribute_index (it, "is_copy");

  // Attribrute stride lengths
  const int dm               = particle.stride(it, ia_m);
  const int dp               = particle.stride(it, ia_x);
  const int dv               = particle.stride(it, ia_vx);
  const int dmetal_fraction  =
    (ia_metal_fraction != -1) ? particle.stride(it, ia_metal_fraction) : 0;
  const int dcreation_time   = particle.stride(it,ia_creation_time);
  const int did              = particle.stride(it,ia_id);
  const int dcopy            = particle.stride(it, ia_copy);

  // Declare pointers for particle attribute arrays
  enzo_float * pmass;
  enzo_float * px  ;
  enzo_float * py  ;
  enzo_float * pz  ;
  enzo_float * pvx ;
  enzo_float * pvy ;
  enzo_float * pvz ;
  enzo_float * pmetal_fraction;
  enzo_float * pcreation_time;
  int64_t    * pid;
  int64_t    * is_copy;

  // Get field data
  Field field = block->data()->field();
  int gx,gy,gz;
  field.ghost_depth (0, &gx, &gy, &gz);
  int mx, my, mz;
  field.dimensions(0, &mx, &my, &mz);
  double hx, hy, hz;
  block->cell_width(&hx, &hy, &hz);
  const double cell_volume = hx * hy * hz;
  double xm, ym, zm;
  block->lower(&xm,&ym,&zm);

  // get pointers to field values
  enzo_float * density      = (enzo_float *) field.values("density");
  enzo_float * vx_gas       = (enzo_float *) field.values("velocity_x");
  enzo_float * vy_gas       = (enzo_float *) field.values("velocity_y");
  enzo_float * vz_gas       = (enzo_float *) field.values("velocity_z");

  enzo_float * metal_density = field.is_field("metal_density") ?
    (enzo_float *) field.values("metal_density") : nullptr;

  // Get the minimum sink mass in code units
  const double minimum_sink_mass =
    min_sink_mass_solar_ * enzo_constants::mass_solar / enzo::units()->mass();

  // Keep track of the number of sinks formed.
  int n_sinks_formed = 0;

  // This vector will store the indices of cells in which sinks are formed, and so will be
  // used later to set new values for density fields.
  std::vector<int> cells_forming_sinks;

  // This vector will be used to store the new densities of these cells
  std::vector<enzo_float> new_densities;

  // Get gravitational constant in code units
  const double const_G =
    enzo::grav_constant_cgs() * enzo::units()->density() *
    enzo::units()->time() * enzo::units()->time();
  
  // Get density threshold in code units for this cycle (value will change in
  // cosmological simultions.
  const double density_threshold =
    physical_density_threshold_cgs_ / enzo::units()->density();

  // Get global block index
  // ix/y/z_block is the x/y/z index of the block
  // nx/y/z_block is the number of blocks along the x/y/z axis
  int ix_block, iy_block, iz_block, nx_block, ny_block, nz_block;
  block->index_global(&ix_block, &iy_block, &iz_block,
		      &nx_block, &ny_block, &nz_block);

  const int n_blocks = nx_block * ny_block * nz_block;
  const int global_block_index = INDEX(ix_block,iy_block,iz_block,nx_block,ny_block);

  // Now get the total number of cells per block
  const int n_cells_per_block = mx * my * mz;

  // We want a unique ID for each cell for each cycle. This will be used to set sink particle
  // IDs and will be used for the seeds for the random number generator which generates the
  // initial particle positions.
  const uint64_t global_cell_index_start =
    (enzo::simulation()->cycle() * n_blocks + global_block_index) * n_cells_per_block;

  // Will be used to generate the random offsets
  std::uniform_real_distribution<double> distribution (0.0,1.0);
  std::mt19937_64 generator;

  for (int iz=gz; iz<mz - gz; iz++){
    for (int iy=gy; iy<my - gy; iy++){
      for (int ix=gx; ix<mx - gx; ix++){

	int block_cell_index = INDEX(ix,iy,iz,mx,my);

	// Check if density is greater than density threshold
	if (density[block_cell_index] <= density_threshold) continue;

	// Check if the mass of the potential sink particle is greater than the minimum sink
	// mass
	const enzo_float density_change =
	  std::min(density[block_cell_index] - density_threshold,
		   max_mass_fraction_ * density[block_cell_index]);

	const enzo_float sink_mass = density_change * cell_volume;

	if (sink_mass < minimum_sink_mass) continue;

	// Check whether Jeans length is resolved.
	if (!jeans_length_not_resolved_(block, block_cell_index, const_G)) continue;

	// Check whether flow is converging
	if (!flow_is_converging_(block, ix, iy, iz)) continue;

	// If check_density_maximum_ is true, check if this cell's density is a local maximum
	if (check_density_maximum_ && !density_is_local_maximum_(block,ix,iy,iz)) continue;

	// If we are here, this cell satisfies all the conditions for forming a sink particle.

	// Compute global cell index
	const uint64_t global_cell_index = global_cell_index_start + block_cell_index;

	// Append to the vectors
	cells_forming_sinks.push_back(block_cell_index);
	new_densities.push_back(density[block_cell_index] - density_change);

	// So. now create a sink particle
	n_sinks_formed++;

	// ip_block is the index of the particle in the block
	// ip_batch is the index if the particle in its batch
	// ibatch is the index of the batch
	int ip_block = particle.insert_particles(it, 1);
	int ip_batch;
	int ibatch;
	particle.index(ip_block, &ibatch, &ip_batch);

	// Set the mass of the sink particle to be sink_mass
	pmass = (enzo_float *) particle.attribute_array(it, ia_m, ibatch);
	pmass[ip_batch * dm] = sink_mass;

	// Get 3 seeds for the random number generator using the global cell index and
	// offset_seed_shift_
	uint64_t x_seed = offset_seed_shift_ + 3 * global_cell_index;
	uint64_t y_seed = offset_seed_shift_ + 3 * global_cell_index + 1;
	uint64_t z_seed = offset_seed_shift_ + 3 * global_cell_index + 2;

	// x_offset is in range [-hx*max_offset_cell_fraction_,hx*max_offset_cell_fraction_]
	generator.seed(x_seed);
	const double x_offset =
	  hx * max_offset_cell_fraction_ * (2.0 * distribution(generator) - 1.0);
	// y_offset is in range [-hy*max_offset_cell_fraction_,hy*max_offset_cell_fraction_]
	generator.seed(y_seed);
	const double y_offset =
	  hy * max_offset_cell_fraction_ * (2.0 * distribution(generator) - 1.0);
	// z_offset is in range [-hz*max_offset_cell_fraction_,hz*max_offset_cell_fraction_]
	generator.seed(z_seed);
	const double z_offset =
	  hz * max_offset_cell_fraction_ * (2.0 * distribution(generator) - 1.0);

	// Set particle position to be at centre of cell plus the offset
	px = (enzo_float *) particle.attribute_array(it, ia_x, ibatch);
	py = (enzo_float *) particle.attribute_array(it, ia_y, ibatch);
	pz = (enzo_float *) particle.attribute_array(it, ia_z, ibatch);

	px[ip_batch * dp] = xm + (ix - gx + 0.5) * hz + x_offset;
	py[ip_batch * dp] = ym + (iy - gy + 0.5) * hy + y_offset;
	pz[ip_batch * dp] = zm + (iz - gz + 0.5) * hz + z_offset;

	// Set particle velocity equal to gas velocity in cell
	pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ibatch);
	pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ibatch);
	pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ibatch);
	pvx[ip_batch * dv] = vx_gas[block_cell_index];
	pvy[ip_batch * dv] = vy_gas[block_cell_index];
	pvz[ip_batch * dv] = vz_gas[block_cell_index];

	// Set creation time equal to current time
	pcreation_time     =
	  (enzo_float *) particle.attribute_array(it, ia_creation_time, ibatch);
	pcreation_time[ip_batch * dcreation_time] = enzo::block(block)->time();

	// Set ID to be the global cell index
	pid = (int64_t * ) particle.attribute_array(it, ia_id, ibatch);
	pid[ip_batch * did] = global_cell_index;

	// If we are tracking metals, set metal fraction of sink particle
	if (metal_density){
	  pmetal_fraction  =
	    (enzo_float *) particle.attribute_array(it, ia_metal_fraction, ibatch);
	  pmetal_fraction[ip_batch * dmetal_fraction] =
	    metal_density[block_cell_index] / density[block_cell_index];
	}

	/* Specify that newly created particle is not a copy*/
	is_copy = (int64_t *) particle.attribute_array(it,ia_copy,ibatch);
	is_copy[ip_batch * dcopy] = 0;
      }
    }
  } // Loop over active cells

  // Add newly created sink particles to simulation.
  enzo::simulation()->data_insert_particles(n_sinks_formed);

  // Loop over the cells which formed sink particles and set new values for density fields.
  for (decltype(cells_forming_sinks.size()) i = 0; i < cells_forming_sinks.size(); i++){
    int ind = cells_forming_sinks[i];
    EnzoMethodStarMaker::rescale_densities(enzo::block(block),
					   ind,
					   new_densities[i] / density[ind]);
    density[ind] = new_densities[i];
  }

  return;
}

//-------------------------------------------------------------------------------------------

double EnzoMethodSinkMaker::timestep ( Block *block) const throw()
{
  return std::numeric_limits<double>::max();
}

//-------------------------------------------------------------------------------------------

bool EnzoMethodSinkMaker::jeans_length_not_resolved_(Block * block, int i,
						     double const_G) throw()
{
  Field field = block->data()->field();

  // First, let's get the square of the sound speed.
  // Do this by getting the specific internal energy in this cell, and then computing
  // the sound speed squared from that.

  // Get pointers to values of density field and specific total energy field
  enzo_float * density      = (enzo_float*) field.values("density");
  enzo_float * vx           = (enzo_float *) field.values("velocity_x");
  enzo_float * vy           = (enzo_float *) field.values("velocity_y");
  enzo_float * vz           = (enzo_float *) field.values("velocity_z");
  enzo_float * specific_te  = (enzo_float *) field.values("total_energy");

  // Check if dual-energy formalism is in use. If so, then this needs to use
  // then the simulation evolves a "specific internal energy" field.
  const bool dual_energy =
    ! enzo::fluid_props()->dual_energy_config().is_disabled();

  const enzo_float * specific_ie_field
    = dual_energy ? (enzo_float*) field.values("internal_energy") : nullptr;

  // Get pointers to magnetic field values if they exist.
  enzo_float * bx =
    field.is_field("bfield_x") ? (enzo_float*) field.values("bfield_x") : nullptr;
  enzo_float * by =
    field.is_field("bfield_y") ? (enzo_float*) field.values("bfield_y") : nullptr;
  enzo_float * bz =
    field.is_field("bfield_z") ? (enzo_float*) field.values("bfield_z") : nullptr;

  // Get the specific internal energy
  enzo_float specific_ie_cell;
  if (dual_energy) specific_ie_cell = specific_ie_field[i];
  else if (bx)
    specific_ie_cell = specific_te[i] -
      0.5 * ( vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i] ) -
      0.5 * ( bx[i] * bx[i] + by[i] * by[i] + bz[i] * bz[i] ) / density[i];
  else specific_ie_cell = specific_te[i] -
	 0.5 * ( vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i] );

  // Now compute the square of the sound speed
  const double gamma = enzo::fluid_props()->gamma();
  const double c_s_2 =  gamma * (gamma - 1.0) * specific_ie_cell;

  // Now compute the Jeans length
  const enzo_float jeans_length = sqrt(cello::pi * c_s_2 / (const_G * density[i]));

  // Get the maximum cell width
  double hx, hy, hz;
  block->cell_width(&hx, &hy, &hz);
  const double max_cell_width = std::max(hx,std::max(hy,hz));

  return (jeans_length < jeans_length_resolution_cells_ * max_cell_width);
}

//-----------------------------------------------------------------------------------------

bool EnzoMethodSinkMaker::flow_is_converging_(Block * block,
					     int ix, int iy, int iz) throw()
{

  // Get field data
  Field field = block->data()->field();
  int mx, my, mz;
  field.dimensions(0, &mx, &my, &mz);
  double hx, hy, hz;
  block->cell_width(&hx, &hy, &hz);
  double xm, ym, zm;
  block->lower(&xm,&ym,&zm);

  // Get pointers to velocity fields
  enzo_float * vx           = (enzo_float *) field.values("velocity_x");
  enzo_float * vy           = (enzo_float *) field.values("velocity_y");
  enzo_float * vz           = (enzo_float *) field.values("velocity_z");

  // a_{ij} = 0.5*(dv_i/dx_j + dv_j/dx_i)
  // First, compute the diagonal terms and the trace (sum of diagonal terms).
  // Since trace is also equal to the sum of the eigenvalues, if the trace is non-negative,
  // at least one of the eigenvalues is non-negative.
  const double a_11 = (vx[INDEX(ix+1,iy,iz,mx,my)] - vx[INDEX(ix-1,iy,iz,mx,my)]) / (2.0 * hx);
  const double a_22 = (vy[INDEX(ix,iy+1,iz,mx,my)] - vy[INDEX(ix,iy-1,iz,mx,my)]) / (2.0 * hy);
  const double a_33 = (vz[INDEX(ix,iy,iz+1,mx,my)] - vz[INDEX(ix,iy,iz-1,mx,my)]) / (2.0 * hz);

  if (a_11 + a_22 + a_33 >= 0) return false;

  // Compute other terms in the tensor / matrix
  const double a_12 =
    0.5 * ((vx[INDEX(ix,iy+1,iz,mx,my)] - vx[INDEX(ix,iy-1,iz,mx,my)]) / (2.0 * hy)
	   + (vy[INDEX(ix+1,iy,iz,mx,my)] - vy[INDEX(ix-1,iy,iz,mx,my)]) /  (2.0 * hx));

  const double a_13 =
    0.5 * ((vx[INDEX(ix,iy,iz+1,mx,my)] - vx[INDEX(ix,iy,iz-1,mx,my)]) / (2.0 * hz)
	   + (vz[INDEX(ix+1,iy,iz,mx,my)] - vz[INDEX(ix-1,iy,iz,mx,my)]) / (2.0 * hx));

  const double a_23 =
    0.5 * ((vy[INDEX(ix,iy,iz+1,mx,my)] - vy[INDEX(ix,iy,iz-1,mx,my)]) / (2.0 * hz)
	   + (vz[INDEX(ix,iy+1,iz,mx,my)] - vz[INDEX(ix,iy-1,iz,mx,my)]) / (2.0 * hy));

  // Set the coefficients of the cubic equation which gives the
  // eigenvalues, i.e. lambda^3 + A*lambda^2 + B*lambda + C = 0

  const double A = -1.0 * (a_11 + a_22 + a_33);
  const double B = -1.0 * (  a_11 * a_22 + a_11 * a_33
			   + a_22 * a_33 + a_12 * a_12
			   + a_23 * a_23 + a_13 * a_13 );
  const double C =
    a_11 * a_23 * a_23 +
    a_22 * a_13 * a_13 +
    a_33 * a_12 * a_12 -
    a_11 * a_22 * a_33;

  // Equation can be transformed to the form t^3 + beta*t + gamma = 0
  // where t = lambda - alpha, with alpha, beta, gamma defined as
  // follows

  const double alpha = -3.0 * A;
  const double beta  = B - A * A / 3.0;
  const double gamma = C + 2.0 * A * A * A / 27.0 - A * B / 3.0;

  /// Can transform this to a trigonometric equation by taking
  /// t = 2*sqrt(-beta/3) * cos(theta). See
  /// https://en.wikipedia.org/wiki/Cubic_equation for derivation

  /// We check if any of the eigenvalues are non-negative, if so, return false
  for (int i = 0; i < 3; i++){
    const double lambda_k = alpha + 2.0 * sqrt(-1.0 * beta / 3.0) *
      cos( acos(1.5 * gamma / beta * sqrt(-3.0 / beta)) / 3.0 + i * 2.0 * cello::pi / 3.0);
    if (lambda_k >= 0) return false;
  }

  // If we are here, all the eigenvalues are negative, so return true

  return true;
}

//---------------------------------------------------------------------------------------

bool EnzoMethodSinkMaker::density_is_local_maximum_(Block * block,
						    int ix, int iy, int iz) throw()
{

  // Get field dimensions
  Field field = block->data()->field();
  int mx, my, mz;
  field.dimensions (0, &mx, &my, &mz);
  int gx,gy,gz;
  field.ghost_depth (0, &gx, &gy, &gz);

  // Index of the central cell
  const int central_cell_index = INDEX(ix,iy,iz,mx,my);

  // Get pointer to density field
  enzo_float * density = (enzo_float *) field.values("density");

  // Loop over neighboring cells.
  // If any of them have a density larger than the central cell, return false
  for (int jz = iz - 1; jz <= iz + 1; jz++){
    for (int jy = iy - 1; jy <= iy + 1; jy++){
      for (int jx = ix - 1; jx <= ix + 1; jx++){

	const int j = INDEX(jx,jy,jz,mx,my);

	if (density[j] > density[central_cell_index])  return false;
      }
    }
  }

  // If we get here then this cell must be a local maximum.
  return true;

}

//-------------------------------------------------------------------------------------------


void EnzoMethodSinkMaker::do_checks_(Block *block) throw()
{
  // Check sink particle attributes
  cello::particle_descr()->check_particle_attribute("sink","mass");
  cello::particle_descr()->check_particle_attribute("sink","x");
  cello::particle_descr()->check_particle_attribute("sink","y");
  cello::particle_descr()->check_particle_attribute("sink","z");
  cello::particle_descr()->check_particle_attribute("sink","vx");
  cello::particle_descr()->check_particle_attribute("sink","vy");
  cello::particle_descr()->check_particle_attribute("sink","vz");
  cello::particle_descr()->check_particle_attribute("sink","creation_time");
  cello::particle_descr()->check_particle_attribute("sink","id");
  cello::particle_descr()->check_particle_attribute("sink","is_copy");

  // If there is a "metal_density" field, check if sink particles have a
  // "metal_fraction" attribute.
  Field field = block->data()->field();
  if (field.is_field("metal_density"))
    cello::particle_descr()->check_particle_attribute("sink","metal_fraction");

  // Check if merge_sinks method precedes pm_update method
  ASSERT("EnzoMethodSinkMaker",
	 "sink_maker must precede pm_update",
	 enzo::problem()->method_precedes("sink_maker",
					  "pm_update"));

  // Check if either PPM or VL+CT method is being used.
  ASSERT("EnzoMethodSinkMaker",
	 "accretion requires ppm or vlct methods.",
	 enzo::problem()->method_exists("mhd_vlct") ||
	 enzo::problem()->method_exists("ppm"));

  // Get density floor set by the hydro method.
  //
  // TODO: remove the use of density_dbl_prec. This is a temporary workaround
  const double density_floor =
    enzo::fluid_props()->fluid_floor_config().density_dbl_prec();

  // Get density threshold in code units.
  const double density_threshold =
    physical_density_threshold_cgs_ / enzo::units()->density();

  // In a cosmological simulation, the density unit is the mean matter density
  // of the universe, which decreases with time, which means that the value of
  // a fixed physical density quantity will increase with time. So if the density
  // threshold is above the density floor at the start of the simulation, it is
  // guaranteed to be above the density floor at all times.
  ASSERT2("EnzoMethodSinkMaker",
	  "Density threshold must be at least as large as the density "
	  "floor set by the VL+CT method. At start of simulation, "
	  "density threshold is %g and density floor is %g "
	  "(code density units).",
	  density_threshold,
	  density_floor,
	  density_threshold >= density_floor);

  return;
}
