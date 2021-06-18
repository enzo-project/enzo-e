// See LICENSE_CELLO file for license and copyright information
/// @file     enzo_EnzoInitialCloud.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon May 27 2019
/// @brief    [\ref Enzo] Implementation of EnzoInitialCloud for initializing a
///           a spherical cloud in a cloudroutine for spherical cloud in a wind

#include <random>
#include <limits>

#include "enzo.hpp"
#include "charm_enzo.hpp"
#include "cello.hpp"

//----------------------------------------------------------------------

void prep_subcell_offsets_(double* subcell_offsets, double cell_width,
			   int subsample_n)
{
  int num_subcells_axis = (int)std::pow(2,subsample_n);
  double cur_frac = 1./std::pow(2,subsample_n+1);
  subcell_offsets[0] = cur_frac*cell_width;
  for (int i = 1; i <num_subcells_axis; i++){
    cur_frac += 1./std::pow(2,subsample_n);
    subcell_offsets[i] = cur_frac*cell_width;
  }
}

class SphereRegion
{
public:
  SphereRegion(Data *data, int subsample_n,
	       double center_x, double center_y, double center_z,
	       double radius)
    : subsample_n_(subsample_n),
      num_subcells_axis_((int)std::pow(2,subsample_n)),
      num_subsampled_cells_(std::pow(std::pow(2,subsample_n),3)),
      center_x_(center_x),
      center_y_(center_y),
      center_z_(center_z),
      radius_(radius),
      sqr_radius_(radius*radius),
      cell_xf_(nullptr),
      cell_yf_(nullptr),
      cell_zf_(nullptr),
      subcell_xoffsets_(nullptr),
      subcell_yoffsets_(nullptr),
      subcell_zoffsets_(nullptr)
  {
    ASSERT("SphereRegion", "subsample_num must be >= 0", subsample_n>=0);
    ASSERT("SphereRegion", "radius must be >0", radius>0);

    double hx, hy, hz;
    data->field_cell_width(&hx,&hy,&hz);

    ASSERT("EnzoSphereRegion",
	   "the cell width along each axis is assumed to be the same",
	   (fabs((hy - hx)/hx) <= INIT_CLOUD_TOLERANCE) &&
           (fabs((hz - hx)/hx) <= INIT_CLOUD_TOLERANCE));

    Field field = data->field();
    int nx,ny,nz;
    field.size(&nx,&ny,&nz);
    int gx,gy,gz;
    field.ghost_depth(field.field_id("density"),&gx,&gy,&gz);
    int ndx, ndy, ndz;
    ndx = nx + 2*gx;
    ndy = ny + 2*gy;
    ndz = nz + 2*gz;

    cell_xf_ = new double[ndx + 1];
    cell_yf_ = new double[ndy + 1];
    cell_zf_ = new double[ndz + 1];

    data->field_cell_faces(cell_xf_,cell_yf_,cell_zf_,gx,gy,gz,1,1,1);

    subcell_xoffsets_ = new double[ndx];
    subcell_yoffsets_ = new double[ndy];
    subcell_zoffsets_ = new double[ndz];

    prep_subcell_offsets_(subcell_xoffsets_, hx, subsample_n_);
    prep_subcell_offsets_(subcell_yoffsets_, hy, subsample_n_);
    prep_subcell_offsets_(subcell_zoffsets_, hz, subsample_n_);
  }

  ~SphereRegion(){
    delete[] cell_xf_;
    delete[] cell_yf_;
    delete[] cell_zf_;

    delete[] subcell_xoffsets_;
    delete[] subcell_yoffsets_;
    delete[] subcell_zoffsets_;
  }


  /// returns whether a point lies within the sphere
  bool check_point(double x, double y, double z)
  {
    return (x*x + y*y + z*z) <= sqr_radius_;
  }

  /// returns the fraction of the cell at (iz, iy, ix) that is enclosed by the
  /// sphere 
  double cell_fraction_enclosed(int iz, int iy, int ix){
    // Get the x,y,z on the lower left and upper right corners of the cell
    double left_x = cell_xf_[ix];  double right_x = cell_xf_[ix+1];
    double left_y = cell_yf_[iy];  double right_y = cell_yf_[iy+1];
    double left_z = cell_zf_[iz];  double right_z = cell_zf_[iz+1];

    // Find the closest and furthest points on the cell from the sphere center
    double closest_x, closest_y, closest_z;
    double furthest_x, furthest_y, furthest_z;
    find_closest_furthest_points_(left_x, right_x, center_x_, &closest_x,
				  &furthest_x);
    find_closest_furthest_points_(left_y, right_y, center_y_, &closest_y,
				  &furthest_y);
    find_closest_furthest_points_(left_z, right_z, center_z_, &closest_z,
				  &furthest_z);
    
    
    if (check_point(furthest_x, furthest_y, furthest_z)){
      // if the furthest point on the cell is enclosed in the sphere, then the
      // entire cell is enclosed
      return 1.0;
    } else if (!check_point(closest_x, closest_y, closest_z)){
      // if the closest point on the cell is not enclosed then, then the none of
      // cell is enclosed
      return 0.0;
    }

    // now let's determine the fraction of the cell that is enclosed

    int n_enclosed = 0;
    for (int sub_iz=0; sub_iz<num_subcells_axis_; sub_iz++){
      double sub_zc = left_z+subcell_zoffsets_[sub_iz];

      for (int sub_iy=0; sub_iy<num_subcells_axis_; sub_iy++){
	double sub_yc = left_y + subcell_yoffsets_[sub_iy];

	for (int sub_ix=0; sub_ix<num_subcells_axis_; sub_ix++){
	  double sub_xc = left_x+subcell_xoffsets_[sub_ix];

	  if (check_point(sub_xc, sub_yc, sub_zc)){
	    // if the subcell center is enclosed, increment n_enclosed
	    n_enclosed++;
	  }
	}
      }
    }

    return ((double)n_enclosed)/num_subsampled_cells_;

  }

private:

  /// Finds the closest and nearest coordinate values in a cell
  /// (This get's called separately for x, y, and z)
  void find_closest_furthest_points_(double cell_min_pos, double cell_max_pos,
				     double sphere_center,
				     double *nearest_point,
				     double *furthest_point)
  {
    if (sphere_center <= cell_min_pos){
      *nearest_point = cell_min_pos;
      *furthest_point = cell_max_pos;
    } else if (sphere_center >= cell_max_pos) {
      *nearest_point = cell_max_pos;
      *furthest_point = cell_min_pos;
    } else {
      *nearest_point = sphere_center;
      if ((sphere_center - cell_min_pos) > (cell_max_pos - sphere_center)){
	*furthest_point = cell_min_pos;
      } else {
	*furthest_point = cell_max_pos;
      }
    }
  }

public:

  /// subsample_n
  int subsample_n_;
  /// number of subsampled_cells per axis (2^subsample_n_)
  int num_subcells_axis_;
  /// total number of subsampled_cells
  double num_subsampled_cells_;

  /// center of spher
  double center_x_;
  double center_y_;
  double center_z_;
  /// sphere radius
  double radius_;
  double sqr_radius_;

  /// x, y, and z cell face positions. Along axis i, if there are mi cells
  /// (including ghost cells), then the corresponding array has mi+1 elements
  double *cell_xf_;
  double *cell_yf_;
  double *cell_zf_;

  /// offsets from the left  cell face that indicates the center of a given
  /// sub-sampled cell
  double *subcell_xoffsets_;
  double *subcell_yoffsets_;
  double *subcell_zoffsets_;

};

//----------------------------------------------------------------------

class TruncatedGaussianNoise
{
public:

  TruncatedGaussianNoise(Block* block, unsigned int default_seed, double mean,
			 double standard_deviation, double trunc_deviations)
  {
    // first compute the random seed for the current block
    int ix, iy, iz, nx, ny, nz;
    block->index_global(&ix, &iy, &iz, &nx, &ny, &nz);
    unsigned int index = (unsigned int)(ix + nx*(iy + ny*iz));

    unsigned int wraparound_diff =
      std::numeric_limits<unsigned int>::max() - default_seed;
    if (index > wraparound_diff){
      gen_.seed(index - wraparound_diff - 1);
    } else {
      gen_.seed(default_seed + index);
    }

    ASSERT("TruncatedGaussianNoise", "trunc_deviations must be >= zero",
	   trunc_deviations>=0);
    truncate_ = (trunc_deviations != 0.);
    upper_lim_ = mean + trunc_deviations * standard_deviation;
    lower_lim_ = mean - trunc_deviations * standard_deviation;
    dist_ = std::normal_distribution<double>{mean,standard_deviation};
  }

  double operator()(){
    if (truncate_){
      double out = dist_(gen_);
      while (out > upper_lim_ || out < lower_lim_){
	out = dist_(gen_);
      }
      return out;
    } else {
      return dist_(gen_);
    }
  }

private:

  std::mt19937 gen_;
  std::normal_distribution<double> dist_;
  bool truncate_;
  double upper_lim_;
  double lower_lim_;
};

//----------------------------------------------------------------------

// Helper function that checks the assumption that the bfield is constant in
// the active_zone
void check_uniform_bfield_(EFlt3DArray bfield,
			   const std::string field_name,
			   const int gx, const int gy, const int gz)
{
  EFlt3DArray temp = bfield.subarray(CSlice(gz,-gz), CSlice(gy,-gy),
				     CSlice(gx,-gx));
  const enzo_float val = temp(0,0,0);

  for (int iz = 0; iz<temp.shape(0); iz++){
    for (int iy = 0; iy<temp.shape(1); iy++){
      for (int ix = 0; ix<temp.shape(2); ix++){

	if (temp(iz,iy,ix) != val){
	  ERROR1("EnzoInitialCloud",
		 ("Currently %s must have a constant value throughout the "
		  "entire active zone"),
		 field_name.c_str());
	}

      }
    }
  }
}

//----------------------------------------------------------------------

class MHDHandler
{
  // at the moment this only supports uniform bfields
  //
  // this currently supports 2 means of initialization:
  //   1. This handler can initialize the bfields uniform bfields
  //   2. (Deprecated) bfields can be pre-initialized. This is the historic
  //      approach - but doing this usually involves invoking InitialValue
  //      which is not thread-safe. We may choose to drop this.
  //
  // Note: Testing indicates that these approaches may cause slightly different
  // values to be initialized for non-trivial bfield values (non-zero and
  // not-easily rounded). It appears that values initialized by this class
  // have the exact expected value whereas the values initialized by the value
  // initializer are slightly offset (there could potentially be a cast to
  // float causing round-off errors)

public:
  MHDHandler(Block *block, bool has_bfield, bool has_interface_bfield,
             bool initialize_uniform_bfield, double uniform_bfield[3])
    : has_bfield_(has_bfield),
      has_interface_bfield_(has_interface_bfield),
      initialize_uniform_bfield_(initialize_uniform_bfield)
  {
    for (int i = 0; i < 3; i++){
      uniform_bfield_[i] = (enzo_float) uniform_bfield[i];
    }
    EnzoFieldArrayFactory array_factory(block);
    if (has_bfield_){
      bfield_x = array_factory.from_name("bfield_x");
      bfield_y = array_factory.from_name("bfield_y");
      bfield_z = array_factory.from_name("bfield_z");
    }
    if (has_interface_bfield_){
      bfieldi_x = array_factory.from_name("bfieldi_x");
      bfieldi_y = array_factory.from_name("bfieldi_y");
      bfieldi_z = array_factory.from_name("bfieldi_z");
    }

    if (initialize_uniform_bfield_){
      ASSERT("MHDHandler::MHDHandler",
             ("to initialize uniform bfields, the \"bfield_x\", \"bfield_y\", "
              "and \"bfield_z\" fields must be defined"), has_bfield_);
      magnetic_edens_wind_ = 0.5 * (uniform_bfield_[0] * uniform_bfield_[0] +
                                    uniform_bfield_[1] * uniform_bfield_[1] +
                                    uniform_bfield_[2] * uniform_bfield_[2]);
      this->initialize_bfields();
    } else if (has_bfield_) {
      // compute the specific magnetic energy of the wind in the lower left
      // corner of the active zone
      Field field = block->data()->field();
      int gx,gy,gz;
      field.ghost_depth(field.field_id("bfield_x"),&gx,&gy,&gz);

      // Currently bfields are assumed to be constant throughout the active zone
      // We check this assumption below:
      check_uniform_bfield_(bfield_x, "bfield_x", gx, gy, gz);
      check_uniform_bfield_(bfield_y, "bfield_y", gx, gy, gz);
      check_uniform_bfield_(bfield_z, "bfield_z", gx, gy, gz);

      // now compute specific magnetic energy of the wind
      magnetic_edens_wind_ = 0.5 * (bfield_x(gz,gy,gx) * bfield_x(gz,gy,gx) +
                                    bfield_y(gz,gy,gx) * bfield_y(gz,gy,gx) +
                                    bfield_z(gz,gy,gx) * bfield_z(gz,gy,gx));
    } else {
      magnetic_edens_wind_ = 0;
    }
  }

  enzo_float magnetic_edens_wind() const { return magnetic_edens_wind_; }

  enzo_float magnetic_edens(int iz, int iy, int ix){
    if (!has_bfield_) { return 0;}
    return 0.5*(bfield_x(iz,iy,ix)*bfield_x(iz,iy,ix)+
                bfield_y(iz,iy,ix)*bfield_y(iz,iy,ix)+
                bfield_z(iz,iy,ix)*bfield_z(iz,iy,ix));
  }

  void initialize_bfields()
  {
    // in the future, this may need to be refactorred
    if (initialize_uniform_bfield_) {
      if (has_bfield_){ // always true
        initialize_uniform_bfield_helper_(bfield_x, bfield_y, bfield_z);
      }
      if (has_interface_bfield_) {
        initialize_uniform_bfield_helper_(bfieldi_x, bfieldi_y, bfieldi_z);
      }
    }

  }

  static MHDHandler construct_MHDHandler(Block *block,
                                         bool initialize_uniform_bfield,
                                         double uniform_bfield[3])
  {
    Field field = block->data()->field();
    bool has_bfield = (field.is_field("bfield_x") ||
                       field.is_field("bfield_y") ||
                       field.is_field("bfield_z"));
    bool has_interface_bfield = (field.is_field("bfieldi_x") ||
                                 field.is_field("bfieldi_y") ||
                                 field.is_field("bfieldi_z"));
    return MHDHandler(block, has_bfield, has_interface_bfield,
                      initialize_uniform_bfield, uniform_bfield);
  }

private:

  void initialize_uniform_bfield_helper_(EFlt3DArray &bx, EFlt3DArray &by,
                                         EFlt3DArray &bz)
  {
    // Handle them separately because we might be handling face-centered fields
    for (int iz = 0; iz<bx.shape(0); iz++){
      for (int iy = 0; iy<bx.shape(1); iy++){
        for (int ix = 0; ix<bx.shape(2); ix++){
          bx(iz,iy,ix) = uniform_bfield_[0];
        }
      }
    }

    for (int iz = 0; iz<by.shape(0); iz++){
      for (int iy = 0; iy<by.shape(1); iy++){
        for (int ix = 0; ix<by.shape(2); ix++){
          by(iz,iy,ix) = uniform_bfield_[1];
        }
      }
    }

    for (int iz = 0; iz<bz.shape(0); iz++){
      for (int iy = 0; iy<bz.shape(1); iy++){
        for (int ix = 0; ix<bz.shape(2); ix++){
          bz(iz,iy,ix) = uniform_bfield_[2];
        }
      }
    }
  }

private: // atributes
  const bool has_bfield_;
  const bool has_interface_bfield_;
  const bool initialize_uniform_bfield_;
  enzo_float uniform_bfield_[3];

  enzo_float magnetic_edens_wind_;
  EFlt3DArray bfield_x, bfield_y, bfield_z;
  EFlt3DArray bfieldi_x, bfieldi_y, bfieldi_z;
};

//----------------------------------------------------------------------

void EnzoInitialCloud::enforce_block
(Block * block, const Hierarchy * hierarchy) throw()
{
  EnzoFieldArrayFactory array_factory(block);
  EFlt3DArray density = array_factory.from_name("density");
  EFlt3DArray velocity_x = array_factory.from_name("velocity_x");
  EFlt3DArray velocity_y = array_factory.from_name("velocity_y");
  EFlt3DArray velocity_z = array_factory.from_name("velocity_z");
  EFlt3DArray total_energy = array_factory.from_name("total_energy");

  // Currently assume pressure equilibrium and pre-initialized uniform B-field
  SphereRegion sph(block->data(), subsample_n_, cloud_center_x_,
		   cloud_center_y_, cloud_center_z_, cloud_radius_);

  // Passively advected scalar "cloud_dye" denotes the material originally
  // located in the cloud
  FieldDescr * field_descr = cello::field_descr();
  const bool use_cloud_dye
    = (field_descr->is_field("cloud_dye") &&
       field_descr->groups()->is_in("cloud_dye", "color"));
  EFlt3DArray cloud_dye_density;
  if (use_cloud_dye){
    cloud_dye_density = array_factory.from_name("cloud_dye");
  }

  const bool set_metal_density
    = (field_descr->is_field("metal_density") &&
       field_descr->groups()->is_in("metal_density","color"));
  if (metal_mass_frac_>0.){
    ASSERT("EnzoInitialCloud::enforce_block",
	   "metal_denisty field must be in the \"color\" group since a metal "
	   "mass fraction was specified.", set_metal_density);
  }
  EFlt3DArray metal_density;
  if (set_metal_density){
    metal_density = array_factory.from_name("metal_density");
  }

  Field field = block->data()->field();
  
  // Handle magnetic fields (mhd indicates whether the fields are present)
  MHDHandler mhd_handler = MHDHandler::construct_MHDHandler
    (block, initialize_uniform_bfield_, uniform_bfield_);

  // Handle internal_energy & compute eint_density
  // dual_energy indicates if there is an internal_energy to initialize
  const bool dual_energy = field.is_field("internal_energy");
  EFlt3DArray internal_energy; // internal_energy field
  double eint_density; // internal energy density (constant over full domain).
                       // used to initialize total energy in cells overlapping
                       // with the cloud

  TruncatedGaussianNoise random_factor_gen(block, perturb_seed_, 1.,
					   perturb_stddev_, truncate_dev_);
  const bool use_random_factor = (perturb_stddev_ != 0);

  if (dual_energy){
    internal_energy = array_factory.from_name("internal_energy");
    double temp = (eint_wind_ + 0.5*velocity_wind_*velocity_wind_ +
		   mhd_handler.magnetic_edens_wind() / density_wind_);
    ASSERT2("EnzoInitialCloud::enforce_block",
	    ("Relative error of the wind's specific etot computed from "
	     "specified eint, velocity, density & preinitialized B-fields, "
	     "w.r.t. the specified etot, has a magnitude of %e (exceeding %e)"),
	    fabs((temp - etot_wind_)/etot_wind_), INIT_CLOUD_TOLERANCE,
	    fabs((temp - etot_wind_)/etot_wind_) <= INIT_CLOUD_TOLERANCE);
    eint_density = eint_wind_ * density_wind_;
  } else {
    ASSERT("EnzoInitialCloud::enforce_block",
	   ("The internal energy of the wind has been specified, but there "
	    "is no \"internal_energy\" field."),
	   eint_wind_ == 0);
    eint_density = ((etot_wind_ - 0.5 * velocity_wind_ * velocity_wind_)
		    * density_wind_ - mhd_handler.magnetic_edens_wind());
  }

  ASSERT("EnzoInitialCloud::enforce_block",
	 "Internal Energy Density must be positive", eint_density > 0);

  for (int iz = 0; iz<density.shape(0); iz++){
    for (int iy = 0; iy<density.shape(1); iy++){
      for (int ix = 0; ix<density.shape(2); ix++){
	velocity_y(iz,iy,ix) = 0.;
	velocity_z(iz,iy,ix) = 0.;

	// In the case of overlap, we use a volume weighted average density
	// For other primitive quantities, like velocity_x (and eventually
	// specific total energy), we use mass weighted average quantities
	//
	// let f = fraction_enclosed  -> V_cloud = f * V_cell
	// V_wind = (1-f) * V_cell
	//
	// rho_cell = f*rho_cl + (1-f)*rho_wind
	// cell_mass:  M_cell = V_cell * rho_cell
	// mass_from cloud:
	//   M_cloud = V_cloud * rho_cloud = f * V_cell * rho_cloud
	// mass from wind:
	//   M_wind = V_wind * rho_wind = (1-f) * V_cell * rho_wind
	//
	// cloud_mass_weight = M_cloud/M_cell
	//   = f * rho_cloud / rho_cell
	// cloud_mass_weight = M_wind/M_cell
	//   = (1-f) * rho_wind / rho_cell

	double frac_enclosed = sph.cell_fraction_enclosed(iz, iy, ix);

	double random_factor = 1.;
	if (use_random_factor && frac_enclosed > 0){
	  random_factor = random_factor_gen();
	}

	double avg_density = (frac_enclosed * density_cloud_ * random_factor +
			      (1. - frac_enclosed) * density_wind_);

	density(iz,iy,ix) = avg_density;
	if (use_cloud_dye){
	  cloud_dye_density(iz,iy,ix) = (frac_enclosed * density_cloud_ *
					 random_factor);
	}
	if (set_metal_density){
	  metal_density(iz,iy,ix) = metal_mass_frac_ * avg_density;
	}

	//cloud_mass_weight = frac_enclosed * density_cloud_ / avg_density;

	double wind_to_average_ratio = density_wind_ / avg_density;
	double wind_mass_weight = (1. - frac_enclosed) * wind_to_average_ratio;
	velocity_x(iz,iy,ix) = (wind_mass_weight * velocity_wind_);

	if (dual_energy) {
	  internal_energy(iz,iy,ix) = eint_wind_ * wind_to_average_ratio;
	}

	if (frac_enclosed == 0){
	  total_energy(iz,iy,ix) = etot_wind_;
	} else {
	  double magnetic_edens = mhd_handler.magnetic_edens(iz,iy,ix);

	  total_energy(iz,iy,ix)
	    = ((eint_density + magnetic_edens) / avg_density +
	       0.5 * velocity_x(iz,iy,ix) * velocity_x(iz,iy,ix));
	}
      }
    }
  }
}
