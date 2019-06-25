// See LICENSE_CELLO file for license and copyright information
/// @file     enzo_EnzoInitialCloud.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon May 27 2019
/// @brief    [\ref Enzo] Implementation of EnzoInitialCloud for initializing a
///           a spherical cloud in a cloudroutine for spherical cloud in a wind

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
      cell_xf_(NULL),
      cell_yf_(NULL),
      cell_zf_(NULL),
      subcell_xoffsets_(NULL),
      subcell_yoffsets_(NULL),
      subcell_zoffsets_(NULL)
  {
    ASSERT("SphereRegion", "subsample_num must be >= 0", subsample_n>=0);
    ASSERT("SphereRegion", "radius must be >0", radius>0);

    double hx, hy, hz;
    data->field_cell_width(&hx,&hy,&hz);
    ASSERT("EnzoSphereRegion",
	   "the cell width along each axis is assumed to be the same",
	   hx == hy && hy == hz);

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
       field_descr->groups()->is_in("cloud_dye", "colour"));
  EFlt3DArray cloud_dye_density;
  if (use_cloud_dye){
    cloud_dye_density = array_factory.from_name("cloud_dye");
  }
  
  // Handle magnetic fields
  bool mhd; // Whether or not we use MHD
  EFlt3DArray bfield_x, bfield_y, bfield_z; // Cell-centered bfields
  enzo_float magnetic_edens_wind; // magnetic energy density of the wind

  Field field = block->data()->field();
  if (field.is_field("bfield_x")){
    mhd = true;
    // initialize bfields
    bfield_x = array_factory.from_name("bfield_x");
    bfield_y = array_factory.from_name("bfield_y");
    bfield_z = array_factory.from_name("bfield_z");

    // compute the specific magnetic energy of the wind in the lower left
    // corner of the active zone
    int gx,gy,gz;
    field.ghost_depth(field.field_id("bfield_x"),&gx,&gy,&gz);

    // Currently bfields are assumed to be constant throughout the active zone
    // We check this assumption below:
    check_uniform_bfield_(bfield_x, "bfield_x", gx, gy, gz);
    check_uniform_bfield_(bfield_y, "bfield_y", gx, gy, gz);
    check_uniform_bfield_(bfield_z, "bfield_z", gx, gy, gz);

    // now compute specific magnetic energy of the wind
    double frac_enclosed = sph.cell_fraction_enclosed(gz, gy, gx);
    ASSERT("EnzoInitialCloud",
	   ("The lower left corner of the active zone is assumed to lie "
	    "outside the cloud."), frac_enclosed == 0);
    magnetic_edens_wind = 0.5*(bfield_x(gz,gy,gx)*bfield_x(gz,gy,gx)+
			       bfield_y(gz,gy,gx)*bfield_y(gz,gy,gx)+
			       bfield_z(gz,gy,gx)*bfield_z(gz,gy,gx));
  } else {
    mhd = false;
    // So that we don't get yelled at by the compiler, set the bfield variables
    // to array [[[0]]] and use them to compute magnetic_wind of 0
    bfield_x = EFlt3DArray(1,1,1);
    bfield_y = EFlt3DArray(1,1,1);
    bfield_z = EFlt3DArray(1,1,1);
    
    magnetic_edens_wind = 0.5*(bfield_x(0,0,0)*bfield_x(0,0,0)+
			       bfield_y(0,0,0)*bfield_y(0,0,0)+
			       bfield_z(0,0,0)*bfield_z(0,0,0));
  }

  double eint_density = ((etot_wind_ - 0.5*velocity_wind_*velocity_wind_)
			 * density_wind_ - magnetic_edens_wind);

  ASSERT("EnzoInitialCloud", "Internal Energy Density must be positive",
	 eint_density > 0);


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

	double wind_mass_weight;
	double frac_enclosed = sph.cell_fraction_enclosed(iz, iy, ix);

	double avg_density = (frac_enclosed * density_cloud_ +
			      (1. - frac_enclosed) * density_wind_);
	density(iz,iy,ix) = avg_density;
	if (use_cloud_dye){
	  cloud_dye_density(iz,iy,ix) = frac_enclosed * density_cloud_;
	}

	//cloud_mass_weight = frac_enclosed * density_cloud_ / avg_density;
	wind_mass_weight = (1. - frac_enclosed) * density_wind_ / avg_density;

	velocity_x(iz,iy,ix) = (wind_mass_weight * velocity_wind_);

	if (frac_enclosed == 0){
	  total_energy(iz,iy,ix) = etot_wind_;
	} else {

	  double magnetic_edens = 0;
	  if (mhd){
	    magnetic_edens = 0.5*(bfield_x(iz,iy,ix)*bfield_x(iz,iy,ix)+
				  bfield_y(iz,iy,ix)*bfield_y(iz,iy,ix)+
				  bfield_z(iz,iy,ix)*bfield_z(iz,iy,ix));
	  }

	  total_energy(iz,iy,ix)
	    = ((eint_density + magnetic_edens) / avg_density +
	       0.5 * velocity_x(iz,iy,ix) * velocity_x(iz,iy,ix));
	}
      }
    }
  }
}
