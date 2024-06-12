// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialInclinedWave.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2019-04-27
/// @brief    Implementation of EnzoInitialInclinedWave for initializing
///           inclined linear MHD waves and inclined circularly polarized
///           alfven waves detailed in Gardiner & Stone (2008). These are used
///           to test the VL+CT MHD integrator. An inclined sound wave can also
///           be initialized (its initial conditions were inspired by the
///           Athena test suite)

// Using the coordinate systems given by eqn 68 of Gardiner & Stone (2008)
// (except that index numbers start at 0)
// The system of the mesh are x,y,z
// The system of the linear wave is x0, x1, x2
// The linear waves are defined along x0
// {{x0},   {{ cos(a), 0, sin(a)},   {{ cos(b), sin(b), 0},   {{x},
//  {x1}, =  {     0,  1,      0}, .  {-sin(b), cos(b), 0}, .  {y},
//  {x2}}    {-sin(a), 0, cos(a)}}    {      0,      0, 1}}    {z}}
//                    R2                          R1
//  R1, the second matix, rotates a 3D vector about it's third dimension by
//  angle -b. In this case, if z-axis points out of the page, then rotate
//  axes clockwise. The rotated y-axis now points along the x1 dimension
//
//  R2, the first matrix, rotates a 3D vector about it's second axis by angle
//  -a. In this case, if the axis along the second dimension (x1) points out of
//  the page, then rotate the other axes clockwise. The rotated x- and z- axes
//  now point along the x0- and x2- axes
//
//  Note that they give transformation of x0,x1,x2 -> x,y,z while the
//  implementation primarily use x,y,z -> x0,x1,x2 (since the the rotation is
//  easier to visualize)
//
//
//  In general, this implementation is more complicated than it needs to be.
//  It could definitely be simplified (although the Rotation class probably
//  shouldn't be touched since it will be reused)

#include <sstream> // std::stringstream

#include "Enzo/initial/initial.hpp"
#include "Enzo/enzo.hpp"
#include "Cello/cello.hpp"


//----------------------------------------------------------------------

// This will be reused to implement Cosmic Rays
// When we do that this should probably be made into a class template
// so that the values don't necessarily have to be cast to doubles
// (We uses doubles by default for it's current use because the convention
//  is to perform all operations on position using Double precision)
class Rotation {
public:
  Rotation(double a, double b)
  {
    matrix_[0][0] = std::cos(a)*std::cos(b);
    matrix_[0][1] = std::cos(a)*std::sin(b);
    matrix_[0][2] = std::sin(a);
    matrix_[1][0] = -1.*std::sin(b);
    matrix_[1][1] = std::cos(b);
    matrix_[1][2] = 0.;
    matrix_[2][0] = -1.*std::sin(a)*std::cos(b);
    matrix_[2][1] = -1.*std::sin(a)*std::sin(b);
    matrix_[2][2] = std::cos(a);
  }

  Rotation()
  {
    for (int j = 0; j<3; j++){
      for (int i = 0; i<3; i++){
        matrix_[j][i] = (j == i); // initialize identity matrix
      }
    }
  }

  // rotates vector {v0, v1, v2} to {rot0, rot1, rot2}
  void rot(double v0, double v1, double v2,
	   double &rot0, double &rot1, double &rot2) const noexcept
  {
    rot0 = matrix_[0][0]*v0 + matrix_[0][1]*v1 + matrix_[0][2]*v2;
    rot1 = matrix_[1][0]*v0 + matrix_[1][1]*v1 + matrix_[1][2]*v2;
    rot2 = matrix_[2][0]*v0 + matrix_[2][1]*v1 + matrix_[2][2]*v2;
  }

  // rotates vector {rot0, rot1, rot2} to {v0, v1, v2}
  void inv_rot(double rot0, double rot1, double rot2,
	       double &v0, double &v1, double &v2) const noexcept
  {
    v0 = matrix_[0][0]*rot0 + matrix_[1][0]*rot1 + matrix_[2][0]*rot2;
    v1 = matrix_[0][1]*rot0 + matrix_[1][1]*rot1 + matrix_[2][1]*rot2;
    v2 = matrix_[0][2]*rot0 + matrix_[1][2]*rot1 + matrix_[2][2]*rot2;
  }

  // For debugging
  void print_matrix()
  {
    for (int j = 0; j<3; j++){
      if (j == 0){
	CkPrintf("{");
      } else {
	CkPrintf(" ");
      }
      CkPrintf(" %.5g, %.5g, %.5g}",
	       matrix_[j][0], matrix_[j][1], matrix_[j][2]);
      if (j == 2){
	CkPrintf("}\n");
      } else {
	CkPrintf("\n");
      }
    }
  }

private:
  double matrix_[3][3];
};

//----------------------------------------------------------------------

// Use initializer functors to set up the grid
//   - for scalars: density, total energy density
//   - for vectors: momentum, magnetic vector potential
// These functors are most frequently most easily defined along axes that are
// rotated with respect to the grid. RotatedScalarInit and
// RotatedVectorInit are defined wrap around such initializers
// There are large similarities between Scalar and Vector Inits, there
// may be a way to consolidate

class ScalarInit{
public:
  virtual ~ScalarInit()
  {}

  virtual double operator() (double x0, double x1, double x2)=0;
};

//----------------------------------------------------------------------

class LinearScalarInit : public ScalarInit{

public:
  LinearScalarInit(double background, double eigenvalue,
		   double amplitude, double lambda,
                   bool use_cosine = true)
    : background_(background),
      eigenvalue_(eigenvalue),
      amplitude_(amplitude),
      lambda_(lambda),
      use_cosine_(use_cosine)
  {}

  double operator() (double x0, double x1, double x2){
    double tmp = x0*2.*cello::pi/lambda_;
    return (background_ + amplitude_ * eigenvalue_
	    * (use_cosine_*std::cos(tmp) + (1 - use_cosine_) * std::sin(tmp)));
  }

protected:
  // value of the background state
  double background_;
  // value of the eigenvalue for a given mode
  double eigenvalue_;
  // perturbation amplitude
  double amplitude_;
  // wavelength
  double lambda_;
  // indicates whether cosine or sine is used
  bool use_cosine_;
};

//----------------------------------------------------------------------

template<class F>
class WrapperScalarInit : public ScalarInit{

public:
  WrapperScalarInit(F functor)
    : functor_(functor)
  {}

  double operator() (double x0, double x1, double x2){
    return functor_(x0, x1, x2);
  }

protected:
  F functor_;
};

//----------------------------------------------------------------------

class RotatedScalarInit : public ScalarInit {
public:
  RotatedScalarInit(double a, double b, ScalarInit *inner)
    : ScalarInit()
  {
    rot_ = new Rotation(a,b);
    inner_ = inner;
  }

  ~RotatedScalarInit()
  {
    delete inner_;
    delete rot_;
  }

  double operator() (double x, double y, double z)
  {
    double x0,x1,x2;
    rot_->rot(x,y,z,x0,x1,x2);
    return (*inner_)(x0,x1,x2);
  }

private:
  Rotation *rot_;
  ScalarInit *inner_;
};

//----------------------------------------------------------------------

class VectorInit {
public:
  virtual ~VectorInit()
  {}

  virtual void operator() (double x0, double x1, double x2,
			   double &v0, double &v1, double &v2)=0;
};

//----------------------------------------------------------------------

class LinearVectorInit : public VectorInit {
public:
  LinearVectorInit(double back0, double back1, double back2,
		   double ev0, double ev1, double ev2,
		   double amplitude, double lambda,
                   bool use_cosine = true)
    : VectorInit(),
      back0_(back0), back1_(back1), back2_(back2),
      ev0_(ev0), ev1_(ev1), ev2_(ev2),
      amplitude_(amplitude),
      lambda_(lambda),
      use_cosine_(use_cosine)
  {}

  void operator() (double x0, double x1, double x2,
		   double &v0, double &v1, double &v2)
  {
    double tmp = x0*2.*cello::pi/lambda_;
    double trig_term = (use_cosine_ * std::cos(tmp) +
                        (1 - use_cosine_) * std::sin(tmp));
    v0 = back0_ + amplitude_*ev0_*trig_term;
    v1 = back1_ + amplitude_*ev1_*trig_term;
    v2 = back2_ + amplitude_*ev2_*trig_term;
  }

protected:
  // value of the background state
  double back0_, back1_, back2_;
  // value of the eigenvalue for a given mode
  double ev0_, ev1_, ev2_;
  // perturbation amplitude
  double amplitude_;
  // wavelength
  double lambda_;
  // whether to use cosine or sine
  bool use_cosine_;
};

//----------------------------------------------------------------------

class RotatedVectorInit : public VectorInit {
public:
  RotatedVectorInit(double a, double b, VectorInit *inner)
    : VectorInit()
  {
    rot_ = new Rotation(a,b);
    inner_ = inner;
  }

  ~RotatedVectorInit()
  {
    delete rot_;
    delete inner_;
  }

  void operator() (double x, double y, double z,
		   double &v0, double &v1, double &v2)
  {
    double x0,x1,x2, rot0, rot1, rot2;
    rot_->rot(x,y,z,x0,x1,x2);
    (*inner_)(x0,x1,x2,rot0,rot1,rot2);
    rot_->inv_rot(rot0,rot1,rot2,v0,v1,v2);
  }

private:
  Rotation *rot_;
  VectorInit *inner_;
};

//----------------------------------------------------------------------

class LinearVectorPotentialInit : public VectorInit {
public:
  // The absence of ev_B0 is intentional
  LinearVectorPotentialInit(double back_B0, double back_B1,
			    double back_B2, double ev_B1,
			    double ev_B2, double amplitude,
			    double lambda)
    :VectorInit(),
     back_B0_(back_B0), back_B1_(back_B1), back_B2_(back_B2),
     ev_B1_(ev_B1), ev_B2_(ev_B2),
     amplitude_(amplitude),
     lambda_(lambda)
  {}

  void operator() (double x0, double x1, double x2,
		   double &v0, double &v1, double &v2)
  {
    v0 = (x2 * amplitude_ * ev_B1_ * std::cos(2. * cello::pi * x0/ lambda_) -
	  x1 * amplitude_ * ev_B2_ * std::cos(2. * cello::pi * x0/ lambda_));
    v1 = back_B2_ * x0; // For Gardiner & Stone (2008) setup, should be 0
    v2 = back_B0_ * x1 - back_B1_ * x0;
  }

protected:
  // value of the magnetic field background state
  double back_B0_, back_B1_, back_B2_;
  // magnetic field eigenvalue for a given mode (its always 0, along dim 0)
  double ev_B1_, ev_B2_;
  // perturbation amplitude
  double amplitude_;
  // wavelength
  double lambda_;
};

//----------------------------------------------------------------------

template<class F>
class WrapperVectorInit : public VectorInit{

public:
  WrapperVectorInit(F functor)
    : functor_(functor)
  {}

  void operator() (double x0, double x1, double x2,
                   double &v0, double &v1, double &v2)
  { functor_(x0,x1,x2, v0,v1,v2); }

protected:
  F functor_;
};

//----------------------------------------------------------------------

template<class F>
static WrapperVectorInit<F>* alloc_vector_init_wrapper_(F &functor)
{ return new WrapperVectorInit<F>(functor); }

//----------------------------------------------------------------------

HydroInitPack::HydroInitPack(ScalarInit* density_init,
                             VectorInit* velocity_init,
                             ScalarInit* pressure_init,
                             InitializerForm initializer_form)
  : density_init(density_init),
    velocity_init(velocity_init),
    pressure_init(pressure_init),
    initializer_form(initializer_form)
{}

//----------------------------------------------------------------------

class MeshPos {
  // TODO: Remove this class. We can accomplish the same thing with the
  //       Data::field_cells method

public:
  MeshPos(Block *block)
  {
    Field field  = block->data()->field();

    // lower left edge of active region
    block->data()->lower(&x_left_edge_,&y_left_edge_,&z_left_edge_);

    // cell widths
    block->data()->field_cell_width(&dx_,&dy_,&dz_);

    // ghost depth 
    field.ghost_depth(0,&gx_,&gy_,&gz_);
  }

  // compute positions at cell-centers
  double x(int k, int j, int i){
    return x_left_edge_ + dx_* (0.5+(double)(i-gx_));}
  double y(int k, int j, int i){
    return y_left_edge_ + dy_* (0.5+(double)(j-gy_));}
  double z(int k, int j, int i){
    return z_left_edge_ + dz_* (0.5+(double)(k-gz_));}

  // computes the x, y, or z position for the cell corner at the
  // (i-1/2, j,k),  (i,j-1/2,k) or (i,j, k-1/2)
  double x_face(int k, int j, int i){
    return x_left_edge_ + dx_* (double)(i-gx_);}
  double y_face(int k, int j, int i){
    return y_left_edge_ + dy_* (double)(j-gy_);}
  double z_face(int k, int j, int i){
    return z_left_edge_ + dz_* (double)(k-gz_);}

  double dx() {return dx_;}
  double dy() {return dy_;}
  double dz() {return dz_;}
  
private:
  // the starting position of the edge of the active grid
  double x_left_edge_, y_left_edge_, z_left_edge_;
  // ghost depth (and the index at the start of the active region)
  int gx_, gy_, gz_;
  double dx_, dy_, dz_;
};

//----------------------------------------------------------------------

// Components of the magnetic field from the vector potential
void setup_bfield(Block * block, VectorInit *a, MeshPos &pos,
		  int mx, int my, int mz)
{
  Field field = block->data()->field();
  ASSERT("setup_bfield", ("Can only currently handle initialization of "
			  "B-fields if interface values are used"),
	 field.is_field("bfieldi_x") &&
	 field.is_field("bfieldi_y") &&
	 field.is_field("bfieldi_z"));

  EFlt3DArray bfieldi_x = field.view<enzo_float>("bfieldi_x");
  EFlt3DArray bfieldi_y = field.view<enzo_float>("bfieldi_y");
  EFlt3DArray bfieldi_z = field.view<enzo_float>("bfieldi_z");

  if (a == NULL){
    for (int iz=0; iz<mz+1; iz++){
      for (int iy=0; iy<my+1; iy++){
	for (int ix=0; ix<mx+1; ix++){
	  if ((iz != mz) && (iy != my)) { bfieldi_x(iz,iy,ix) = 0; }
	  if ((iz != mz) && (ix != mx)) { bfieldi_y(iz,iy,ix) = 0; }
	  if ((iy != my) && (ix != mx)) { bfieldi_z(iz,iy,ix) = 0; }
	}
      }
    }
    EnzoInitialBCenter::initialize_bfield_center(block);
    return;
  }


  // allocate corner-centered arrays for the magnetic vector potentials
  // Ax, Ay, and Az are always cell-centered along the x, y, and z dimensions,
  // respectively
  CelloView<double,3> Ax(mz+1,my+1,mx);
  CelloView<double,3> Ay(mz+1,my,mx+1);
  CelloView<double,3> Az(mz,my+1,mx+1);

  // Compute the Magnetic Vector potential at all points on the grid
  for (int iz=0; iz<mz+1; iz++){
    for (int iy=0; iy<my+1; iy++){
      for (int ix=0; ix<mx+1; ix++){

	double temp_ax, temp_ay, temp_az;

	if (ix != mx){
	  (*a)(pos.x(iz,iy,ix),  pos.y_face(iz,iy,ix),  pos.z_face(iz,iy,ix),
	       Ax(iz,iy,ix), temp_ay, temp_az);
	}
	if (iy != my){
	  (*a)(pos.x_face(iz,iy,ix),  pos.y(iz,iy,ix),  pos.z_face(iz,iy,ix),
	       temp_ay, Ay(iz,iy,ix), temp_az);
	}
	if (iz != mz){
	  (*a)(pos.x_face(iz,iy,ix),  pos.y_face(iz,iy,ix),  pos.z(iz,iy,ix),
	       temp_ax, temp_ay, Az(iz,iy,ix));
	}
      }
    }
  }

  // Compute the Interface B-fields
  EnzoInitialBCenter::initialize_bfield_interface(block, Ax, Ay, Az);

  // Compute the Cell-Centered B-fields
  EnzoInitialBCenter::initialize_bfield_center(block);

}

//----------------------------------------------------------------------

void setup_eint_(Block *block)
{
  // because cell-centered bfields for CT are computed in terms of enzo_float,
  // this calculation is handled in terms of enzo_float
  // This operation could be split off and placed in a separate initializer
  Field field = block->data()->field();

  EFlt3DArray density, etot, eint;
  density = field.view<enzo_float>("density");
  etot = field.view<enzo_float>("total_energy");
  eint = field.view<enzo_float>("internal_energy");

  EFlt3DArray velocity_x, velocity_y, velocity_z;
  velocity_x = field.view<enzo_float>("velocity_x");
  velocity_y = field.view<enzo_float>("velocity_y");
  velocity_z = field.view<enzo_float>("velocity_z");

  const bool mhd = field.is_field("bfield_x");
  EFlt3DArray bfield_x, bfield_y, bfield_z;
  if (mhd){
    bfield_x = field.view<enzo_float>("bfield_x");
    bfield_y = field.view<enzo_float>("bfield_y");
    bfield_z = field.view<enzo_float>("bfield_z");
  }

  for (int iz=0; iz<density.shape(0); iz++){
    for (int iy=0; iy<density.shape(1); iy++){
      for (int ix=0; ix<density.shape(2); ix++){
	enzo_float kinetic, magnetic;
	kinetic = 0.5 * (velocity_x(iz,iy,ix) * velocity_x(iz,iy,ix) +
		         velocity_y(iz,iy,ix) * velocity_y(iz,iy,ix) +
		         velocity_z(iz,iy,ix) * velocity_z(iz,iy,ix));
	if (mhd){
	  magnetic =
	    0.5 * (bfield_x(iz,iy,ix)*bfield_x(iz,iy,ix) +
		   bfield_y(iz,iy,ix)*bfield_y(iz,iy,ix) +
		   bfield_z(iz,iy,ix)*bfield_z(iz,iy,ix)) / density(iz,iy,ix);
	} else {
	  magnetic = 0.;
	}
	eint(iz,iy,ix) = etot(iz,iy,ix) - kinetic - magnetic;
      }
    }
  }
}

//----------------------------------------------------------------------

static void setup_fluid_(Block *block, HydroInitPack hydro_init_pack,
                         MeshPos &pos, int mx, int my, int mz, double gamma)
{
  Field field = block->data()->field();
  EFlt3DArray density = field.view<enzo_float>("density");
  EFlt3DArray specific_total_energy = field.view<enzo_float>("total_energy");

  EFlt3DArray velocity_x = field.view<enzo_float>("velocity_x");
  EFlt3DArray velocity_y = field.view<enzo_float>("velocity_y");
  EFlt3DArray velocity_z = field.view<enzo_float>("velocity_z");

  ScalarInit* density_init = hydro_init_pack.density_init.get();
  if (hydro_init_pack.initializer_form == InitializerForm::conserved){

    VectorInit* momentum_init = hydro_init_pack.velocity_init.get();
    ScalarInit* total_energy_density_init = hydro_init_pack.pressure_init.get();

    for (int iz=0; iz<mz; iz++){
      for (int iy=0; iy<my; iy++){
        for (int ix=0; ix<mx; ix++){
          double x,y,z;
          double rho, px, py, pz, etot_dens;
          x = pos.x(iz,iy,ix); y = pos.y(iz,iy,ix); z = pos.z(iz,iy,ix);
          rho = (*density_init)(x,y,z);
          density(iz,iy,ix) = (enzo_float)rho;
          (*momentum_init)(x, y, z, px, py, pz);
          velocity_x(iz,iy,ix) = (enzo_float)(px/rho);
          velocity_y(iz,iy,ix) = (enzo_float)(py/rho);
          velocity_z(iz,iy,ix) = (enzo_float)(pz/rho);

          etot_dens = (*total_energy_density_init)(x,y,z);
          specific_total_energy(iz,iy,ix) = (enzo_float)(etot_dens/rho);
        }
      }
    }

    // If present, setup internal energy (to test dual energy formalism):
    if (field.is_field("internal_energy")) { setup_eint_(block); }

  } else { // hydro_init_pack.initializer_form == InitializerForm::primitive

    VectorInit* velocity_init = hydro_init_pack.velocity_init.get();
    ScalarInit* pressure_init = hydro_init_pack.pressure_init.get();

    const bool dual_energy = field.is_field("internal_energy");
    EFlt3DArray specific_internal_energy;
    if (dual_energy){
      specific_internal_energy = field.view<enzo_float>("internal_energy");
    }

    const bool mag = field.is_field("bfield_x");
    using RdOnlyEFlt3DArray = CelloView<const enzo_float, 3>;
    const RdOnlyEFlt3DArray bfield_x = (mag) ?
      RdOnlyEFlt3DArray(field.view<enzo_float>("bfield_x")) :
      RdOnlyEFlt3DArray();
    const RdOnlyEFlt3DArray bfield_y = (mag) ?
      RdOnlyEFlt3DArray(field.view<enzo_float>("bfield_y")) :
      RdOnlyEFlt3DArray();
    const RdOnlyEFlt3DArray bfield_z = (mag) ?
      RdOnlyEFlt3DArray(field.view<enzo_float>("bfield_z")) :
      RdOnlyEFlt3DArray();

    const double inv_gm1 = 1.0 / (gamma - 1.0);

    for (int iz=0; iz<mz; iz++){
      for (int iy=0; iy<my; iy++){
        for (int ix=0; ix<mx; ix++){
          double x = pos.x(iz,iy,ix);
          double y = pos.y(iz,iy,ix);
          double z = pos.z(iz,iy,ix);

          double rho = (*density_init)(x,y,z);
          density(iz,iy,ix) = (enzo_float)rho;

          double vx, vy, vz;
          (*velocity_init)(x, y, z, vx, vy, vz);
          velocity_x(iz,iy,ix) = (enzo_float)(vx);
          velocity_y(iz,iy,ix) = (enzo_float)(vy);
          velocity_z(iz,iy,ix) = (enzo_float)(vz);

          
          double inv_rho = 1.0 / rho;
          double eint = inv_gm1 * (*pressure_init)(x,y,z) * inv_rho;
          if (dual_energy) {
            specific_internal_energy(iz,iy,ix) = (enzo_float)eint;
          }

          double non_thermal_e = 0.5 * ((vx*vx) + (vy*vy) + (vz*vz));
          if (mag){
            double b2 = ((bfield_x(iz,iy,ix) * bfield_x(iz,iy,ix)) +
                         (bfield_y(iz,iy,ix) * bfield_y(iz,iy,ix)) +
                         (bfield_z(iz,iy,ix) * bfield_z(iz,iy,ix)));
            non_thermal_e += 0.5 * b2 * inv_rho;
          }
          specific_total_energy(iz,iy,ix) = (enzo_float)(eint + non_thermal_e);
          if (dual_energy) { specific_internal_energy(iz,iy,ix) = eint;}
        }
      }
    }
  }
}

//----------------------------------------------------------------------

EnzoInitialInclinedWave::EnzoInitialInclinedWave(int cycle, double time,
                                                 ParameterGroup p) throw()
  : Initial (cycle,time),
    alpha_(p.value_float("alpha",0.0)),
    beta_(p.value_float("beta",0.0)),
    gamma_(enzo::fluid_props()->gamma()), // TODO: drop this as a parameter
    amplitude_(p.value_float("amplitude",1.e-6)),
    lambda_(p.value_float("lambda",1.0)),
    parallel_vel_(p.value_float("parallel_vel",
                                std::numeric_limits<double>::min())),
    pos_vel_(p.value_logical("positive_vel",true)),
    wave_type_("")
{
  const std::string wave_type_param_name = p.full_name("wave_type");
  if (p.param("wave_type") == nullptr){
    ERROR1("EnzoInitialInclinedWave::EnzoInitialInclinedWave",
           "%s must be specified", wave_type_param_name.c_str());
  } else {
    wave_type_ = p.value_string("wave_type", "");
  }

  std::vector<std::string> mhd_waves = mhd_waves_();
  std::vector<std::string> hd_waves  =  hd_waves_();

  bool is_MHD_wave = (std::find(mhd_waves.begin(), mhd_waves.end(), wave_type_)
                      != mhd_waves.end());

  if (is_MHD_wave | (wave_type_ == "jeans")){
    ASSERT1("EnzoInitialInclinedWave",
	    "parallel_vel isn't currently supported for wave_type: \"%s\"",
	    wave_type_.c_str(), !specified_parallel_vel_());
  } else if (std::find(hd_waves.begin(), hd_waves.end(), wave_type_)
	     == hd_waves.end()) {
    // wave_type_ isn't a known type. Raise error with list of known
    // wave_types (follows https://stackoverflow.com/a/1430774):

    std::stringstream ss;
    for (std::vector<std::string>::size_type i = 0; i < mhd_waves.size(); i++){
      if (i != 0) { ss << ", "; }
      ss << "'" << mhd_waves[i] << "'";
    }
    for (std::vector<std::string>::size_type i = 0; i < hd_waves.size(); i++){
      ss << ", '" << hd_waves[i] << "'";
    }
    ss << ", 'jeans'";

    std::string s = ss.str();
    ERROR3("EnzoInitialInclinedWave",
           "%s has invalid value, \"%s\". Must be %s",
	   wave_type_param_name.c_str(), wave_type_.c_str(), s.c_str());
  }
}

//----------------------------------------------------------------------

std::vector<std::string> EnzoInitialInclinedWave::mhd_waves_() const throw()
{
  std::vector<std::string> v = {"fast", "alfven", "slow", "mhd_entropy",
				"circ_alfven"};
  return v;
}

//----------------------------------------------------------------------

std::vector<std::string> EnzoInitialInclinedWave::hd_waves_() const throw()
{
  std::vector<std::string> v = {"sound",
				"hd_entropy",            // perturb v0
				"hd_transv_entropy_v1",  // perturb v1
				"hd_transv_entropy_v2"}; // perturb v2
  return v;
}

//----------------------------------------------------------------------

void EnzoInitialInclinedWave::enforce_block(Block * block,
                                            const Hierarchy * hierarchy) throw()
{
  // Set up the test problem
  // Only currently works on unigrid and only currently supports hydro methods
  // and VLCT (PPML initial conditions are much more complicated)
  VectorInit *a_init = nullptr;

  Field field = block->data()->field();
  const bool mhd = field.is_field("bfield_x");

  std::vector<std::string> hd_waves = hd_waves_();

  const bool is_root_block = block->index().is_root();

  HydroInitPack hydro_inits;

  if (wave_type_ == "jeans"){
    // linear jeans wave
    hydro_inits = prepare_jeans_initializers_(is_root_block);
  } else if (std::find(hd_waves.begin(), hd_waves.end(), wave_type_) !=
             hd_waves.end()) {
    hydro_inits = prepare_HD_initializers_(is_root_block);
  } else {
    ASSERT("EnzoInitialInclinedWave::enforce_block",
	   "A MHD wave requires fields to store magnetic field values.", mhd);
    hydro_inits = prepare_MHD_initializers_(&a_init);
  }


  if (block->is_leaf()){
    // load the dimensions and initialize object to compute positions
    int mx,my,mz;
    const int id = field.field_id ("density");
    field.dimensions (id,&mx,&my,&mz);
    MeshPos pos(block);

    if (mhd) {
      // Initialize bfields. If the wave is purely hydrodynamical, they will be
      // all be set to 0.
      setup_bfield(block, a_init, pos, mx, my, mz);
    }

    // it's important that this happens after setting up the magnetic fields
    setup_fluid_(block, hydro_inits, pos, mx, my, mz, gamma_);
  }

  if (a_init != nullptr) { delete a_init; }

  block->initial_done();
}

//----------------------------------------------------------------------

HydroInitPack EnzoInitialInclinedWave::prepare_jeans_initializers_
(bool is_root_block) const noexcept
{
  ASSERT("EnzoInitialInclinedWave::prepare_jeans_initializers_",
         "The gravity method must be used with the jeans wave.",
         enzo::problem()->method_exists("gravity"));
  ASSERT("EnzoInitialInclinedWave::prepare_jeans_initializers_",
         "There currently isn't support for specifying a velocity with a "
         "Jeans wave wave",
         pos_vel_);

  // this was originally based on the description in section 4.2 of Mullen,
  // Hanawa, & Gammie (2021). Unfortunately, we encountered problems with
  // implementing those equations (we have subsequently confirmed that the
  // problem was a transcription error related to initializing pressure).
  // The current implementation uses the same initial conditions as Athena and
  // (earlier versions of) Athena++
  double lambda = lambda_;
  double wave_number = (2*cello::pi/lambda);

  // compute the gravitational constant in code units
  const double G_code = enzo::grav_constant_codeU();

  // background properties:
  // -> density = 1.0
  // -> all velocity components = 0
  // -> pressure = 1/gamma
  // This means cs_0^2 = 1

  // lambda_J = sqrt(pi * cs_0^2 / (G_code * rho_bkg))
  //   -> based on the initial conditions: cs_0^2 -> 1
  //                                       rho_bkg -> 1
  const double lambda_J_sq = cello::pi / G_code;

  const bool stable_wave = lambda_J_sq >= (lambda * lambda);

  // The dispersion relation is given by:
  //   angular_freq ^ 2 = (wave_number*cs_0)^2 - 4 * pi * G_code;
  //                    = wave_number^2 * cs_0^2 * (1 - (lambda / lambda_J)^2)
  // -> lambda < lambda_J, the Jean waves is stable
  // -> lambda > lambda_J, the Jeans wave is unstable
  //
  // to initialize the problem, we need: sqrt(|angular_freq^2| / wave_number^2)
  // -> this is essentially just just the magnitude of the phase velocity
  //    (recall: v_phase = angular_freq / wave_number)
  // -> v_phase = sqrt(|cs_0^2 * (1 - (lambda / lambda_J)^2)|). Recall that
  //    our choice of initial condtions give: cs_0^2 = 1

  double v_phase = std::sqrt(std::fabs(1 - (lambda*lambda) / lambda_J_sq));
  if (is_root_block) {
    double angular_freq_sq =
      wave_number*wave_number * (1.0 - (lambda*lambda)/ lambda_J_sq);
    const char* adjective = (stable_wave) ? "stable" : "unstable";
    cello::monitor()->print
      ("Initial", ("Initializing %s Jeans wave: lambda = %e, lambda_J = %e, "
                   "angular_freq^2 = %e"),
       adjective, lambda_, std::sqrt(lambda_J_sq), angular_freq_sq);
  }

  // the modern version of Athena++ uses the following initializer. When the
  // Jeans wave test was originally introduced to Athena++ (for use with their
  // Multigrid gravity solver), the amplitude^2 term was not present. The
  // amplitude^2 term also was not present in the C-version of Athena. The
  // amplitude^2 term only got introduced to Athena++ in a commit related to
  // introducing a FFT gravity solver.
  //
  //  double amplitude = amplitude_;
  //  Rotation rot(alpha_, beta_);
  //  auto functor = [=](double x, double y, double z) -> double
  //    {
  //      double x0, x1, x2;
  //      rot.rot(x,y,z,x0,x1,x2);
  //      double density_bkg = 1.0;
  //      double tmp = x0*wave_number;
  //      return density_bkg*(1.0 + amplitude * std::sin(tmp) +
  //                          amplitude * amplitude * std::sin(x*wave_number));
  //    };
  //  ScalarInit* density_init = new WrapperScalarInit<decltype(functor)>
  //    (functor);


  // our choice of density initializer, is the same, but it omits the
  // amplitude^2 term
  double density_bkg = 1.0;
  ScalarInit* density_init = new RotatedScalarInit
    (alpha_, beta_, new LinearScalarInit(density_bkg, 1.0,
                                         amplitude_, lambda_, false));


  // coefficient for momentum component along axis 0
  double mom0_coef = (stable_wave) ? 0.0 : density_bkg * v_phase;
  VectorInit* momentum_init = new RotatedVectorInit
    (alpha_, beta_, new LinearVectorInit(0.0, 0.0, 0.0,
                                         mom0_coef, 0.0, 0.0,
                                         amplitude_, lambda));

  double p_bkg = 1.0 / gamma_;
  double eint_dens_bkg = p_bkg / (gamma_ - 1.0);

  ScalarInit* etotdens_init = new RotatedScalarInit
    (alpha_, beta_, new LinearScalarInit(eint_dens_bkg, gamma_,
                                         amplitude_, lambda_, false));
                                         
  return HydroInitPack(density_init, momentum_init, etotdens_init,
                       InitializerForm::conserved);
}

//----------------------------------------------------------------------

HydroInitPack EnzoInitialInclinedWave::prepare_MHD_initializers_
(VectorInit **a_init) const noexcept
{
  double lambda = lambda_;

  if (wave_type_ == "circ_alfven"){

    ASSERT("EnzoInitialInclinedWave::prepare_MHD_initializers_",
           "There currently isn't support for initializing circularly "
           "polarized Alfven wave with a negative velocity.",
           pos_vel_);

    // Construct vector potential initializer. This was defined such that
    // B0 = 1.0, B1 = 0.1 * sin(2*pi*x0/lambda), B2 = 0.1 * cos(2*pi*x0/lambda)
    auto a_func_ = [=](double x0, double x1, double x2,
                       double &A0, double &A1, double &A2)
      {
        A0 = (x2 * 0.1 * std::sin(2. * cello::pi * x0/ lambda_) -
              x1 * 0.1 * std::cos(2. * cello::pi * x0/ lambda_));
        A1 = 0.0;
        A2 = x1;
      };
    *a_init = new RotatedVectorInit(alpha_,beta_,
                                    alloc_vector_init_wrapper_(a_func_));

    // construct hydro initializers
    //
    // it's important that we define the hydro initializers in primitive form
    // - doing so let's us force pressure to be 0.1 and then dynamically
    //   compute the total energy based on the already initialized Bfield
    //   values at the cell-centers
    // - the alternative would be to compute the total energy using the
    //   analytic formulas used to initialize the problem. To initialize the
    //   cell-centered Bfield, we perform finite differences of the Vector
    //   potential. With that in mind, rounding error may make the
    //   cell-centered Bfield slightly different than what we expect. This
    //   isn't ideal since Gardiner & Stone (2008) explicitly states that the
    //   truncation error of B_perp**2/P is important

    // density is 1.0 everywhere
    ScalarInit* rho_init = new LinearScalarInit(1.0,0.0,0.0,lambda);
    // pressure is 0.1 everywhere
    ScalarInit* pressure_init = new LinearScalarInit(0.1,0.0,0.0,lambda);

    auto vel_func_ = [=](double x0, double x1, double x2,
                         double &vel_0, double &vel_1, double &vel_2)
      {
        vel_0 = 0.0;
        vel_1 = 0.1 * std::sin(2. * cello::pi * x0/ lambda);
        vel_2 = 0.1 * std::cos(2. * cello::pi * x0/ lambda);
      };
    VectorInit* velocity_init = new RotatedVectorInit
      (alpha_,beta_, alloc_vector_init_wrapper_(vel_func_));

    return HydroInitPack(rho_init, velocity_init, pressure_init,
                         InitializerForm::primitive);
  } else {
    // wsign indicates direction of propogation.
    double wsign = 1.;
    if (!pos_vel_){
      wsign = -1;
    }
    // initialize the background values:
    // density = 1, pressure = 1/gamma, mom1 = 0, mom2 = 0, B0 = 1, B1=1.5, B2=0
    // For entropy wave, mom0 = 1. Otherwise, mom0 = 0.
    double density_back = 1;
    double mom0_back = 0;
    double mom1_back = 0;
    double mom2_back = 0;
    double b0_back = 1.;
    double b1_back = 1.5;
    double b2_back = 0.0;
    // etot = pressure/(gamma-1)+0.5*rho*v^2+.5*B^2
    double etot_back = (1./gamma_)/(gamma_-1.)+1.625;
    if (wave_type_ == "mhd_entropy"){
      mom0_back = wsign;
      etot_back += 0.5;
    }

    // Get the eigenvalues for density, velocity, and etot
    double density_ev, mom0_ev, mom1_ev, mom2_ev, etot_ev, b1_ev, b2_ev;
    // intentionally omit b0_ev

    if (wave_type_ == "fast"){
      double coef = 0.5/std::sqrt(5.);
      density_ev = 2. *coef;
      mom0_ev = wsign*4. * coef;
      mom1_ev = -1.*wsign*2. * coef;
      mom2_ev = 0;
      etot_ev = 9. * coef;
      b1_ev = 4. * coef;
      b2_ev = 0;
    } else if (wave_type_ == "alfven") {
      density_ev = 0;
      mom0_ev = 0;
      mom1_ev = 0;
      mom2_ev = -1.*wsign*1;
      etot_ev = 0;
      b1_ev = 0;
      b2_ev = 1.;
    } else if (wave_type_ == "slow") {
      double coef = 0.5/std::sqrt(5.);
      density_ev = 4. *coef;
      mom0_ev = wsign*2. * coef;
      mom1_ev = wsign*4. * coef;
      mom2_ev = 0;
      etot_ev = 3. * coef;
      b1_ev = -2. * coef;
      b2_ev = 0;
    } else {
      // wave_type_ == "mhd_entropy"
      double coef = 0.5;
      density_ev = 2. *coef;
      mom0_ev = 2. * coef * wsign;
      mom1_ev = 0;
      mom2_ev = 0;
      etot_ev = 1. * coef;
      b1_ev = 0;
      b2_ev = 0;
    }
 
    double amplitude = amplitude_;

    // construct vector potential initializer
    *a_init = new RotatedVectorInit
      (alpha_, beta_, new LinearVectorPotentialInit(b0_back, b1_back, b2_back,
                                                    b1_ev, b2_ev, amplitude,
                                                    lambda));
    // construct the hydro initializers:
    return build_linear_HD_inits_(density_back, density_ev, etot_back, etot_ev,
				  mom0_back, mom1_back, mom2_back,
				  mom0_ev, mom1_ev, mom2_ev);
  }
}

//----------------------------------------------------------------------

HydroInitPack EnzoInitialInclinedWave::prepare_HD_initializers_
(bool is_root_block) const noexcept
{

  // wsign indicates direction of propogation.
  double wsign = 1.;
  if (!pos_vel_){
    wsign = -1;
  }

  // the values for initializing hydrodynamical equations comes from columns of
  // the matrix in B3 of Stone+08. We assume the background values of
  // density=1, pressure = 1/gamma
  //
  // Because the calculation of the perturbations is MUCH simpler for HD than
  // for MHD, we allow the user to specify custom velocities. For completeness
  // and easy comparision to the MHD wave initialization, we specify the values
  // the eigenvalues, dU = (rho,rho*v0, rho*v1, rho*v2, etot_dens), under the
  // assumptions that v1=v2=0 and gamma = 5/3 below:
  //   - for sound wave we assume bkg v0 = 0. To propagate at +-cs (or +-1)
  //     In this case, dU = (1,+-1,0,0,1.5)
  //   - For all 3 entropy waves, bkg v0 = +- 1. The requirements are:
  //       - hd_entropy needs dU = (1,+-1,0,0,0.5) where the sign of the second
  //         element corresponds to the sign of the bkg v0
  //       - hd_transv_entropy_v1: (0,0,1,0,0)
  //       - hd_transv_entropy_v2: (0,0,0,1,0)

  // determine the background velocity values:
  double v0_back, v1_back, v2_back;
  v0_back = 0;            v1_back = 0;            v2_back = 0;

  if (specified_parallel_vel_()){
    v0_back = parallel_vel_;
  } else if (wave_type_ != "sound"){
    // For all types of entropy waves, if parallel_vel_ is not specified ,
    // pos_vel_ is used to initialize the wave to 1 or -1. But if parallel_vel_
    // is specified, then pos_vel_ is ignored
    v0_back = wsign;
  }

  // Setup the conserved background values
  double squared_v_back = v0_back*v0_back + v1_back*v1_back + v2_back*v2_back;
  double density_back, mom0_back, mom1_back, mom2_back, etot_back;
  density_back = 1;
  mom0_back = v0_back;    mom1_back = v1_back;    mom2_back = v2_back;
  etot_back = ((1. / gamma_) / (gamma_ - 1.) + 0.5 * squared_v_back);

  // now compute eigenvalues:
  double density_ev, mom0_ev, mom1_ev, mom2_ev, etot_ev;

  if (wave_type_ == "sound") {
    // under assumption that density_back = 1 and pressure_back = 1/gamma:
    double h_back = 1/(gamma_ - 1.) + 0.5 * squared_v_back; // enthalpy
    double signed_cs = wsign * 1;

    density_ev = 1;
    mom0_ev = v0_back + signed_cs;
    mom1_ev = v1_back;
    mom2_ev = v2_back;
    etot_ev = h_back + v0_back * signed_cs;
  } else if (wave_type_ == "hd_entropy") {
    // entropy wave where v0 is perturbed
    // Equivalent to "mhd_entropy" without bfield
    density_ev = 1;
    mom0_ev = v0_back;
    mom1_ev = v1_back;
    mom2_ev = v2_back;
    etot_ev = 0.5*squared_v_back;
  } else if (wave_type_ == "hd_transv_entropy_v1") {
    // entropy wave where v1 is perturbed
    density_ev = 0;
    mom0_ev = 0;
    mom1_ev = 1;
    mom2_ev = 0;
    etot_ev = v1_back;
  } else { // (wave_type_ == "hd_transv_entropy_v2")
    // entropy wave where v2 is perturbed
    density_ev = 0;
    mom0_ev = 0;
    mom1_ev = 0;
    mom2_ev = 1;
    etot_ev = v2_back;
  }

  if (is_root_block) {
    cello::monitor()->print("Initial", "Initializing \"%s\" hydro wave.",
                            wave_type_.c_str());
    cello::monitor()->print("Initial",
                            "U = [rho, rho*v0, rho*v1, rho*v2, etot_dens]");
    cello::monitor()->print("Initial", "bkg vals: [%lf, %lf, %lf, %lf, %lf]",
                            density_back, mom0_back, mom1_back, mom2_back,
                            etot_back);
    cello::monitor()->print("Initial", "eigenvals: [%lf, %lf, %lf, %lf, %lf]",
                            density_ev, mom0_ev, mom1_ev, mom2_ev, etot_ev);
  }

  // Now allocate the actual Initializers
  return build_linear_HD_inits_(density_back, density_ev, etot_back, etot_ev,
				mom0_back, mom1_back, mom2_back,
				mom0_ev, mom1_ev, mom2_ev);
}

//----------------------------------------------------------------------

HydroInitPack EnzoInitialInclinedWave::build_linear_HD_inits_
(double density_back, double density_ev,
 double etot_back, double etot_ev,
 double mom0_back, double mom1_back, double mom2_back,
 double mom0_ev, double mom1_ev, double mom2_ev) const
{
  ScalarInit* density_init = new RotatedScalarInit
    (alpha_, beta_, new LinearScalarInit(density_back, density_ev,
                                         amplitude_, lambda_));

  ScalarInit* etot_dens_init = new RotatedScalarInit
    (alpha_, beta_, new LinearScalarInit(etot_back, etot_ev,
                                         amplitude_, lambda_));

  VectorInit* momentum_init = new RotatedVectorInit
    (alpha_, beta_, new LinearVectorInit(mom0_back, mom1_back, mom2_back,
					 mom0_ev, mom1_ev, mom2_ev,
					 amplitude_, lambda_));

  return HydroInitPack(density_init, momentum_init, etot_dens_init,
                       InitializerForm::conserved);
}
