/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodSinkMaker.cpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @date       1 June 2022
/// @brief      Implements a method for forming sink particles based on
///             the method described in Krumholz+ 2004 and Federrath+ 2010.

#include "cello.hpp"
#include "enzo.hpp"

//-------------------------------------------------------------------

EnzoMethodSinkMaker::EnzoMethodSinkMaker
(double min_control_volume_cells,
 double max_control_volume_cells,
 double jeans_length_resolution_cells,
 double density_threshold,
 double max_mass_fraction,
 double min_sink_mass_solar)
  : Method(),
    min_control_volume_cells_(min_control_volume_cells),
    max_control_volume_cells_(max_control_volume_cells),
    jeans_length_resolution_cells_(jeans_length_resolution_cells),
    density_threshold_(density_threshold),
    max_mass_fraction_(max_mass_fraction),
    min_sink_mass_solar_(min_sink_mass_solar)
{
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

  // Check that min_control_volume_cells_ is larger than 1 and
  // max_control_volume_cells_ is less than the minimum ghost depth.
  const int * ghost_depth = enzo::config()->field_ghost_depth;
  const int min_ghost_depth = std::min(ghost_depth[0],
				       std::min(ghost_depth[1],ghost_depth[2]));

  ASSERT("EnzoMethodSinkMaker::EnzoMethodSinkMaker()",
	 "Method:sink_maker:min_control_volume_cells must be at least 1.",
	 (min_control_volume_cells_ >= 1));

  ASSERT("EnzoMethodSinkMaker::EnzoMethodSinkMaker()",
	 "Method:sink_maker:max_control_volume_cells must be no greater than the "
	 "minimum ghost depth (4 by default).",
	 (max_control_volume_cells_ <= min_ghost_depth));

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

  p | min_control_volume_cells_;
  p | max_control_volume_cells_;
  p | jeans_length_resolution_cells_;
  p | min_sink_mass_solar_;

  return;
}

//-------------------------------------------------------------------------------------

void EnzoMethodSinkMaker::compute ( Block *block) throw()
{
  if (enzo::simulation()->cycle() == enzo::config()->initial_cycle)
    do_checks_(block);

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
  const int ia_id               = particle.attribute_index (it, "ID");
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
    min_sink_mass_solar_ * cello::mass_solar / enzo::units()->mass();

  // Keep track of the number of sinks formed.
  int n_sinks_formed = 0;

  // This vector will store the indices of cells in which sinks are formed, and so will be
  // used later to set new values for density fields.
  std::vector<int> cells_forming_sinks;

  // This vector will be used to store the new densities of these cells
  std::vector<enzo_float> new_densities;

  // Get gravitational constant in code units
  const double const_G =
    cello::grav_constant * enzo::units()->mass() * enzo::units()->time() * enzo::units()->time()
    / (enzo::units()->length() * enzo::units()->length() * enzo::units()->length());

  for (int iz=gz; iz<mz - gz; iz++){
    for (int iy=gy; iy<my - gy; iy++){
      for (int ix=gx; ix<mx - gx; ix++){

	int cell_index = INDEX(ix,iy,iz,mx,my);

	// Check if density is greater than density_threshold_
	if (density[cell_index] <= density_threshold_) continue;

	// Check if the mass of the potential sink particle is greater than the minimum sink
	// mass
	const enzo_float density_change =
	  std::min(density[cell_index] - density_threshold_,
		   max_mass_fraction_ * density[cell_index]);

	const enzo_float sink_mass = density_change * cell_volume;

	if (sink_mass < minimum_sink_mass) continue;

	// Check whether Jeans length is resolved.
	enzo_float jeans_length;
	if (!jeans_length_not_resolved_(block, cell_index, const_G, &jeans_length)) continue;

	// Check whether flow is converging
	if (!flow_is_converging_(block, ix, iy, iz)) continue;

	// Check if this cell's density is the maximum density within the control volume
	if (!density_is_local_maximum_(block,jeans_length,ix,iy,iz)) continue;

	// If we are here, this cell satisfies all the conditions for forming a sink particle.

	// Append to the vectors
	cells_forming_sinks.push_back(cell_index);
	new_densities.push_back(density[cell_index] - density_change);

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

	// Set particle position to be at centre of cell
	px = (enzo_float *) particle.attribute_array(it, ia_x, ibatch);
	py = (enzo_float *) particle.attribute_array(it, ia_y, ibatch);
	pz = (enzo_float *) particle.attribute_array(it, ia_z, ibatch);

	px[ip_batch * dp] = xm + (ix - gx + 0.5) * hz;
	py[ip_batch * dp] = ym + (iy - gy + 0.5) * hy;
	pz[ip_batch * dp] = zm + (iz - gz + 0.5) * hz;

	// Set particle velocity equal to gas velocity in cell
	pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ibatch);
	pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ibatch);
	pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ibatch);
	pvx[ip_batch * dv] = vx_gas[cell_index];
	pvy[ip_batch * dv] = vy_gas[cell_index];
	pvz[ip_batch * dv] = vz_gas[cell_index];

	// Set creation time equal to current time
	pcreation_time     =
	  (enzo_float *) particle.attribute_array(it, ia_creation_time, ibatch);
	pcreation_time[ip_batch * dcreation_time] = enzo::block(block)->time();

	// Set ID.
	// Copyied from EnzoMethodStarMakerStochasticSF
	pid = (int64_t * ) particle.attribute_array(it, ia_id, ibatch);
	pid[ip_batch * did] =
	  CkMyPe() + (ParticleData::id_counter[cello::index_static()]++) * CkNumPes();

	// If we are tracking metals, set metal fraction of sink particle
	if (metal_density){
	  pmetal_fraction  =
	    (enzo_float *) particle.attribute_array(it, ia_metal_fraction, ibatch);
	  pmetal_fraction[ip_batch * dmetal_fraction] =
	    metal_density[cell_index] / density[cell_index];
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
    int cell_index = cells_forming_sinks[i];
    EnzoMethodStarMaker::rescale_densities(enzo::block(block),
					   cell_index,
					   new_densities[i] / density[cell_index]);
    density[cell_index] = new_densities[i];
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
						     double const_G,
						     enzo_float* jeans_length) throw()
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

  // Not sure if this is the best way to do this
  // EnzoConfig has three different "dual energy" attributes.
  // One of them is `method_hydro_dual_energy` but it seems that EnzoMethodHydro is not
  // actually used.
  // Assume that if either of the other two are true, there an "specific internal energy" field.
  bool dual_energy =
    enzo::config()->ppm_dual_energy ||  enzo::config()->method_vlct_dual_energy;
  const enzo_float * specific_ie_field
    = dual_energy ? (enzo_float*) field.values("internal_energy") : nullptr;

  // Get pointers to magnetic field values if they exist.
  enzo_float * bx =
    field.is_field("bfield_x") ? (enzo_float*) field.values("b_field_x") : nullptr;
  enzo_float * by =
    field.is_field("bfield_y") ? (enzo_float*) field.values("b_field_y") : nullptr;
  enzo_float * bz =
    field.is_field("bfield_z") ? (enzo_float*) field.values("b_field_z") : nullptr;

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
  const double gamma = enzo::config()->field_gamma;
  const double c_s_2 =  gamma * (gamma - 1.0) * specific_ie_cell;

  // Now compute the Jeans length
  *jeans_length = sqrt(cello::pi * c_s_2 / (const_G * density[i]));

  // Get the maximum cell width
  double hx, hy, hz;
  block->cell_width(&hx, &hy, &hz);
  const double max_cell_width = std::max(hx,std::max(hy,hz));

  return (*jeans_length < jeans_length_resolution_cells_ * max_cell_width);
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

bool EnzoMethodSinkMaker::density_is_local_maximum_(Block * block, enzo_float jeans_length,
						    int ix, int iy, int iz) throw()
{

  // Get field data
  Field field = block->data()->field();
  int mx, my, mz;
  field.dimensions (0, &mx, &my, &mz);
  double hx, hy, hz;
  block->cell_width(&hx, &hy, &hz);
  int gx,gy,gz;
  field.ghost_depth (0, &gx, &gy, &gz);

  // Get maximum and minimum cell widths
  const double min_cell_width = std::min(hx,std::min(hy,hz));
  const double max_cell_width = std::max(hx,std::max(hy,hz));

  // Compute the control volume radius
  const double control_volume_radius =
    std::min(max_control_volume_cells_* min_cell_width,
	     std::max(jeans_length, min_control_volume_cells_ * max_cell_width));

  // Get the "bounding indices" of the region containing the control volume
  const int min_ind_x =  ceil(ix - control_volume_radius / hx) + gx;
  const int min_ind_y =  ceil(iy - control_volume_radius / hy) + gy;
  const int min_ind_z =  ceil(iz - control_volume_radius / hz) + gz;
  const int max_ind_x = floor(ix + control_volume_radius / hx) + gx;
  const int max_ind_y = floor(iy + control_volume_radius / hy) + gy;
  const int max_ind_z = floor(iz + control_volume_radius / hz) + gz;

  // Index of the central cell
  const int central_cell_index = INDEX(ix,iy,iz,mx,my);

  // Get pointer to density field
  enzo_float * density = (enzo_float *) field.values("density");

  // Loop over cells contained within the "bounding indices".
  // For each cell, check if it is within the control volume radius. If so, check if
  // the density is greater than the density of the central cell.
  // If it is, return false
  for (int jz = min_ind_z; jz <= max_ind_z; jz++){
    for (int jy = min_ind_y; jy <= max_ind_y; jy++){
      for (int jx = min_ind_x; jx <= max_ind_x; jx++){

	const int j = INDEX(jx,jy,jz,mx,my);

	const double distance2 =
	  hx * hx * (jx - ix) * (jx - ix) +
	  hy * hy * (jy - iy) * (jy - iy) +
	  hz * hz * (jz - iz) * (jz - iz);

	if ((distance2 <= control_volume_radius * control_volume_radius)
	    && (density[j] > density[central_cell_index])) return false;
      }
    }
  }

  // If we get here then this cell must have largest density within the control volume,
  // so return true.
  return true;

}

//-------------------------------------------------------------------------------------------


void EnzoMethodSinkMaker::do_checks_(const Block *block) throw()
{
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

    // Check if density threshold is at least as large as the density floor
    // set by the hydro method.
    const double density_floor = enzo::problem()->method_exists("mhd_vlct") ?
      enzo::config()->method_vlct_density_floor :
      enzo::config()->ppm_density_floor ;

    ASSERT("EnzoMethodSinkMaker",
	   "Density threshold must be at least as large as the density "
	   "floor set by the hydro method",
	   density_threshold_ >= density_floor);

    // The minimum control volume radius is min_control_volume_cells_ times the maximum
    // cell width. The maximum control volume radius is max_control_volume_cells_ times
    // the minimum cell width. To ensure that the maximum control volume radius is
    // greater than or equal to the minimum control volume radius, we check that the ratio of the
    // maximum cell width to the minimum cell width is less than or equal to the ratio of
    // max_control_volume_cells_ to min_control_volume_cells_.
    // Note that the ratio of max cell width to min cell width is the same at all
    // refinement levels.

    // Get the cell widths
    double hx, hy, hz;
    block->cell_width(&hx, &hy, &hz);
    const double min_cell_width = std::min(hx,std::min(hy,hz));
    const double max_cell_width = std::max(hx,std::max(hy,hz));

    const double cell_width_ratio = max_cell_width / min_cell_width;
    const double control_volume_ratio = max_control_volume_cells_ / min_control_volume_cells_;

    ASSERT2("EnzoMethodSinkMaker",
	    "Ratio of max cell width to min cell width must be less than or equal to the ratio "
	    "of Method:sink_maker:max_control_volume_cells to "
	    "Method:sink_maker:min_control_volume_cells. Currently, the former ratio "
	    "is %g and the latter ratio is %g \n",
	    cell_width_ratio, control_volume_ratio,
	    cell_width_ratio <= control_volume_ratio);

    return;
}
