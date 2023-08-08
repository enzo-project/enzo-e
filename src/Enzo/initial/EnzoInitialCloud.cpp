// See LICENSE_CELLO file for license and copyright information
/// @file     enzo_EnzoInitialCloud.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon May 27 2019
/// @brief    [\ref Enzo] Implementation of EnzoInitialCloud for initializing a
///           a spherical cloud in a cloudroutine for spherical cloud in a wind

#include <random>
#include <limits>
#include <array>

#include "enzo.hpp"
#include "charm_enzo.hpp"
#include "cello.hpp"


// return random value drawn from a uniform distribution on the unit interval
double uniform_dist_transform_(std::minstd_rand &generator,
                               bool include_zero, bool include_one){

  static_assert((generator.max() <= UINT32_MAX) && (generator.min() == 1), 
                "Unexpected PRNG Property"); // sanity-check!

  // cast to double since they can perfectly represent all values of uint32_t
  double raw = static_cast<double>(generator());

  double range;
  if (include_zero && include_one){
    range = static_cast<double>(generator.max()) - 1.;
    raw--;
  } else if (include_zero) {
    range = static_cast<double>(generator.max());
    raw--;
  } else if (include_one) {
    range = static_cast<double>(generator.max());

  } else {
    range = static_cast<double>(generator.max()) + 1.;

  }
  return raw / range;
}


// draw samples from a prng and transform it to a gaussian distribution
// For simplicity, using the standard Box-muller (truncates at ~6.6*sigma)
//
// std::normal_distribution isn't portable across library versions & platforms
std::pair<double, double> normal_dist_transform_(std::minstd_rand &generator){

  double x1 = uniform_dist_transform_(generator, false, false);
  double x2 = uniform_dist_transform_(generator, false, false);

  double coef = std::sqrt(-2.*std::log(x1));
  return std::pair<double,double>(coef * std::cos(2. * cello::pi * x2),
                                  coef * std::sin(2. * cello::pi * x2));
}

std::array<double,3> sample_sphere_points_(std::minstd_rand &generator){
  // draw a sample from uniform distribution on the surface of unit sphere
  // and return corresponding unit vector  

  double x,y,z;
  const int MAX_ITER = 10000000; 
  for (int i =0; i < MAX_ITER; i++){

    std::pair<double,double> pair1 = normal_dist_transform_(generator);
    std::pair<double,double> pair2 = normal_dist_transform_(generator);

    x = pair1.first;
    y = pair1.second;
    z = pair2.first;

    if ((x != 0.) || (y != 0.) || (z!=0.)){
      break;
    } else if (i+1 == MAX_ITER){
      abort(); // almost certainly won't occur unless there's a bug
    }
  }

  double magnitude = std::sqrt(x*x + y*y + z*z);
  return {x/magnitude,y/magnitude,z/magnitude};
}

class WavePerturbation {
  /// This density perturbations modelled by summing a series of inclined
  /// planar cosine waves

public:
  WavePerturbation(std::size_t Nwaves, uint32_t seed, double amplitude,
                   double min_lambda, double max_lambda) noexcept
  {
    Nwaves_ = Nwaves;
    amplitude_ = amplitude;
    if (Nwaves_ > 0 && amplitude_ > 0.){
      std::minstd_rand generator(seed);

      // These are sanity-checks for our assumptions for our about the
      // generator
      static_assert((generator.max() <= UINT32_MAX) &&(generator.min() == 1), 
                    "Unexpected PRNG Property"); // sanity-check!

      for (std::size_t i=0; i < Nwaves; i++){
        // draw wavelength from [min_lambda, max_lambda]
        double lambda = min_lambda +
          (max_lambda - min_lambda) * uniform_dist_transform_(generator,
                                                              true, true);
        std::array<double,3> khat = sample_sphere_points_(generator);
        kx_.push_back(khat[0] * 2 * cello::pi / lambda);
        ky_.push_back(khat[1] * 2 * cello::pi / lambda);
        kz_.push_back(khat[2] * 2 * cello::pi / lambda);

        // draw phase from [0,pi). Don't need to include values between pi and
        // 2*pi since k_x, k_y, and k_z can each be positive or negative
        phi_.push_back(cello::pi * uniform_dist_transform_(generator, true,
                                                           false) );
      }
    }
  }

  void print_summary() const noexcept
  {
    if (Nwaves_ == 0){
      cello::monitor()->print ("EnzoInitialCloud", "Not perturbing cloud.");
    } else {
      for (std::size_t i=0; i < Nwaves_; i++){
        cello::monitor()->print
          ("EnzoInitialCloud",
           "Perturb %2u/%d: kx=%+1.8e ky=%+1.8e kz=%+1.8e phi=%1.8e",
           i,Nwaves_,kx_[i],ky_[i],kz_[i],phi_[i]);
      }
    }
  }


  double operator() (double xc, double yc, double zc,
                     double hx, double hy, double hz) const noexcept
  {
    // The total instantaneous perturbation at location (x,y,z) is given by:
    //    delta = amplitude * Sum( cos(kx_i * x + ky_i*y + kz_i*z + phi_i)),
    // in which the sum runs from i=0 to i = N_waves_-1
    //
    // The volume averaged perturbation, for a cell at center (xc,yc,zc) withs
    // widths (hx,hy,hz) is given by
    //     alpha * Sum( ci * cos(kx_i*xc + ky_i*yc + kz_i*zc + phi_i) ),
    // where:
    //   - the sum runs from i=0 to i = N_waves_-1
    //   - alpha = (8 * amplitue / (hx*hy*hz) )
    //   - ci = sin(kx_i*hx/2)*sin(ky_i*hy/2)*sin(kz_i*hz/2) / (kx_i*ky_i*kz_i)

    double total = 0.;
    double alpha = 8. * amplitude_ / (hx * hy * hz);
    for (std::size_t i = 0; i < Nwaves_; i++){
      double ci = ( std::sin(kx_[i] * hx * 0.5) *
                    std::sin(ky_[i] * hy * 0.5) *
                    std::sin(kz_[i] * hz * 0.5) ) / (kx_[i] * ky_[i] * kz_[i]);

      total += ci * std::cos(kx_[i]*xc + ky_[i]*yc + kz_[i]*zc + phi_[i]);
    }
    return alpha * total;
  }

private:
  /// number of waves
  std::size_t Nwaves_;

  /// Wavelength amplitude
  double amplitude_;

  /// vector holding sets of wavenumber components
  std::vector<double> kx_;
  std::vector<double> ky_;
  std::vector<double> kz_;

  /// vector holding the different phases of the different waves
  std::vector<double> phi_;

};

//----------------------------------------------------------------------

enum class CellIntersection{ enclosed_cell, partial_overlap, no_overlap };

//----------------------------------------------------------------------

class SphereRegion
{
public:
  SphereRegion(double center_x, double center_y, double center_z, double radius)
  {
    center_[0] = center_x;
    center_[1] = center_y;
    center_[2] = center_z;
    ASSERT("SphereRegion", "radius must be >0", radius>0);
    sqr_radius_ = radius*radius;
  }

  /// returns whether a point lies within the sphere
  inline bool check_point(double x, double y, double z) const noexcept
  {
    double dx = x - center_[0];
    double dy = y - center_[1];
    double dz = z - center_[2];
    return (dx*dx + dy*dy + dz*dz) <= sqr_radius_;
  }

  /// determine whether a cell intersects with the sphere
  /// past the (x,y,z) coordinates of the lower left and upper right corners
  CellIntersection check_intersect(std::array<double,3> left,
                                   std::array<double,3> right) const noexcept 
  {
    std::array<double, 3> nearest;    // nearest point in cell to center_
    std::array<double, 3> furthest;   // furthest point in cell from center_

    for (std::size_t i = 0; i<3; i++){
      if (center_[i] <= left[i]){
        nearest[i]  = left[i];
        furthest[i] = right[i];
      } else if (center_[i] >= right[i]) {
        nearest[i] = right[i];
        furthest[i] = left[i];
      } else {
        nearest[i] = center_[i];
        if ((center_[i] - left[i]) > (right[i] - center_[i])){
          furthest[i] = left[i];
        } else {
          furthest[i] = right[i];
        }
      }
    }

    if (check_point(furthest[0], furthest[1], furthest[2])){
      return CellIntersection::enclosed_cell;
    } else if (!check_point(nearest[0], nearest[1], nearest[2])){
      return CellIntersection::no_overlap;
    } else {
      return CellIntersection::partial_overlap;
    }
  }

private:

  /// center x,y, and z of sphere
  std::array<double,3> center_;
  /// squared sphere radius
  double sqr_radius_;
};

void prep_subcell_offsets_(std::vector<double> &subcell_offsets,
                           double cell_width, int subsample_n)
{
  int num_subcells_axis = (int)std::pow(2,subsample_n);
  double cur_frac = 1./std::pow(2,subsample_n+1);
  subcell_offsets[0] = cur_frac*cell_width;
  for (int i = 1; i <num_subcells_axis; i++){
    cur_frac += 1./std::pow(2,subsample_n);
    subcell_offsets[i] = cur_frac*cell_width;
  }
}

class CloudInitHelper {

public:

  CloudInitHelper(Block *block, int subsample_n, SphereRegion region,
                  WavePerturbation perturb_generator)
    : region_(region),
      perturb_generator_(perturb_generator),
      subsample_n_(subsample_n),
      num_subcells_axis_((int)std::pow(2,subsample_n)),
      num_subsampled_cells_(std::pow(std::pow(2,subsample_n),3)),
      cell_widths_(),
      subcell_widths_(),
      cell_xf_(),
      cell_yf_(),
      cell_zf_(),
      subcell_xoffsets_(),
      subcell_yoffsets_(),
      subcell_zoffsets_()
  {
    ASSERT("SphereRegion", "subsample_num must be >= 0", subsample_n>=0);

    Data* data = block->data();

    double hx, hy, hz;
    data->field_cell_width(&hx,&hy,&hz);

    ASSERT("EnzoSphereRegion",
	   "the cell width along each axis is assumed to be the same",
	   (fabs((hy - hx)/hx) <= INIT_CLOUD_TOLERANCE) &&
           (fabs((hz - hx)/hx) <= INIT_CLOUD_TOLERANCE));

    double n_subcells_axis = static_cast<double>(num_subcells_axis_);
    cell_widths_ = {hx, hy, hz};
    subcell_widths_ =
      {hx/n_subcells_axis, hy/n_subcells_axis, hz/n_subcells_axis};

    Field field = data->field();
    int nx,ny,nz;
    field.size(&nx,&ny,&nz);
    int gx,gy,gz;
    field.ghost_depth(field.field_id("density"),&gx,&gy,&gz);
    int ndx, ndy, ndz;
    ndx = nx + 2*gx;
    ndy = ny + 2*gy;
    ndz = nz + 2*gz;

    cell_xf_ = std::vector<double>(ndx + 1);
    cell_yf_ = std::vector<double>(ndy + 1);
    cell_zf_ = std::vector<double>(ndz + 1);

    data->field_cell_faces(cell_xf_.data(), cell_yf_.data(), cell_zf_.data(),
                           gx,gy,gz,1,1,1);

    subcell_xoffsets_ = std::vector<double>(num_subcells_axis_);
    subcell_yoffsets_ = std::vector<double>(num_subcells_axis_);
    subcell_zoffsets_ = std::vector<double>(num_subcells_axis_);

    prep_subcell_offsets_(subcell_xoffsets_, hx, subsample_n_);
    prep_subcell_offsets_(subcell_yoffsets_, hy, subsample_n_);
    prep_subcell_offsets_(subcell_zoffsets_, hz, subsample_n_);

    if (block->index().is_root()) { perturb_generator_.print_summary(); }
  }

  void query_cell(int iz, int iy, int ix, double &frac_enclosed,
                  double &perturbation) const noexcept
  {
    // Get the x,y,z on the lower left and upper right corners of the cell
    std::array<double,3> left, right;
    left[0] = cell_xf_[ix];    right[0] = cell_xf_[ix+1];
    left[1] = cell_yf_[iy];    right[1] = cell_yf_[iy+1];
    left[2] = cell_zf_[iz];    right[2] = cell_zf_[iz+1];

    switch (region_.check_intersect(left,right)){
    case CellIntersection::enclosed_cell :
      {
        frac_enclosed = 1.0;
        perturbation = perturb_generator_
          (0.5*(left[0] + right[0]), 0.5*(left[1] + right[1]),
           0.5*(left[2] + right[2]), cell_widths_[0], cell_widths_[1],
           cell_widths_[2]);
        return;
      }
    case CellIntersection::no_overlap :
      {
        frac_enclosed = 0.0;
        perturbation = 0.0;
        return;
      }
    case CellIntersection::partial_overlap :
      {
        // determine the fraction of the cell that is enclosed and calculate
        // the average density perturbation of that region
        double perturb_sum = 0.;

        int n_enclosed = 0;
        for (int sub_iz=0; sub_iz<num_subcells_axis_; sub_iz++){
          double sub_zc = left[2] + subcell_zoffsets_[sub_iz];

          for (int sub_iy=0; sub_iy<num_subcells_axis_; sub_iy++){
            double sub_yc = left[1] + subcell_yoffsets_[sub_iy];

            for (int sub_ix=0; sub_ix<num_subcells_axis_; sub_ix++){
              double sub_xc = left[0] + subcell_xoffsets_[sub_ix];

              if (region_.check_point(sub_xc, sub_yc, sub_zc)){
                // if the subcell center is enclosed, increment n_enclosed
                n_enclosed++;

                perturb_sum += perturb_generator_(sub_xc, sub_yc, sub_zc,
                                                  subcell_widths_[0],
                                                  subcell_widths_[1],
                                                  subcell_widths_[2]);
              }
            }
          }
        }

        frac_enclosed = static_cast<double>(n_enclosed)/num_subsampled_cells_;
        if (n_enclosed > 0){
          perturbation = perturb_sum / static_cast<double>(n_enclosed);
        } else { // I don't think that this should ever be evaluated
          perturbation = 0.;
        }
      }

    }
  }

  
private:

  const SphereRegion region_;
  const WavePerturbation perturb_generator_;

  /// subsample_n
  int subsample_n_;
  /// number of subsampled_cells per axis (2^subsample_n_)
  int num_subcells_axis_;
  /// total number of subsampled_cells
  double num_subsampled_cells_;

  /// x, y, and z cell widths
  std::array<double,3> cell_widths_;
  /// x, y, and z subcell widths
  std::array<double,3> subcell_widths_;

  /// x, y, and z cell face positions. Along axis i, if there are mi cells
  /// (including ghost cells), then the corresponding array has mi+1 elements
  std::vector<double> cell_xf_;
  std::vector<double> cell_yf_;
  std::vector<double> cell_zf_;

  /// offsets from the left  cell face that indicates the center of a given
  /// sub-sampled cell
  std::vector<double> subcell_xoffsets_;
  std::vector<double> subcell_yoffsets_;
  std::vector<double> subcell_zoffsets_;
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
    Field field = block->data()->field();
    if (has_bfield_){
      bfield_x = field.view<enzo_float>("bfield_x");
      bfield_y = field.view<enzo_float>("bfield_y");
      bfield_z = field.view<enzo_float>("bfield_z");
    }
    if (has_interface_bfield_){
      bfieldi_x = field.view<enzo_float>("bfieldi_x");
      bfieldi_y = field.view<enzo_float>("bfieldi_y");
      bfieldi_z = field.view<enzo_float>("bfieldi_z");
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
  Field field = block->data()->field();

  EFlt3DArray density = field.view<enzo_float>("density");
  EFlt3DArray velocity_x = field.view<enzo_float>("velocity_x");
  EFlt3DArray velocity_y = field.view<enzo_float>("velocity_y");
  EFlt3DArray velocity_z = field.view<enzo_float>("velocity_z");
  EFlt3DArray total_energy = field.view<enzo_float>("total_energy");

  // Currently assume pressure equilibrium and pre-initialized uniform B-field
  CloudInitHelper init_helper(block, subsample_n_,
                              SphereRegion(cloud_center_x_, cloud_center_y_,
                                           cloud_center_z_, cloud_radius_),
                              WavePerturbation(perturb_Nwaves_, perturb_seed_,
                                               perturb_amplitude_,
                                               perturb_min_wavelength_,
                                               perturb_max_wavelength_)
                              );

  // Passively advected scalar "cloud_dye" denotes the material originally
  // located in the cloud
  FieldDescr * field_descr = cello::field_descr();
  const bool use_cloud_dye
    = (field_descr->is_field("cloud_dye") &&
       field_descr->groups()->is_in("cloud_dye", "color"));
  EFlt3DArray cloud_dye_density;
  if (use_cloud_dye){
    cloud_dye_density = field.view<enzo_float>("cloud_dye");
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
    metal_density = field.view<enzo_float>("metal_density");
  }
  
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

  if (dual_energy){
    internal_energy = field.view<enzo_float>("internal_energy");
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

	double frac_enclosed, perturbation;
        init_helper.query_cell(iz, iy, ix, frac_enclosed, perturbation);
        perturbation += 1.;

	double avg_density = (frac_enclosed * density_cloud_ * perturbation +
			      (1. - frac_enclosed) * density_wind_);

	density(iz,iy,ix) = avg_density;
	if (use_cloud_dye){
	  cloud_dye_density(iz,iy,ix) = (frac_enclosed * density_cloud_ *
                                         perturbation);
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
  block->initial_done();
}
