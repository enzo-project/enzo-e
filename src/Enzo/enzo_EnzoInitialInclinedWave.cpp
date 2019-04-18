// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialInclinedWave.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2019-04-27
/// @brief    Implementation of EnzoInitialInclinedWave for initializing
///           inclined linear MHD waves and inclined circularly polarized
///           alfven waves detailed in Gardiner & Stone (2008). These are used
///           to test the VL+CT MHD integrator.

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

#include "enzo.hpp"
#include "charm_enzo.hpp"
#include "cello.hpp"

//----------------------------------------------------------------------

// This will be reused to implement Cosmic Rays
class Rotation {
public:
  Rotation(enzo_float a, enzo_float b)
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

  // rotates vector {v0, v1, v2} to {rot0, rot1, rot2}
  void rot(enzo_float v0, enzo_float v1, enzo_float v2,
	      enzo_float &rot0, enzo_float &rot1, enzo_float &rot2)
  {
    rot0 = matrix_[0][0]*v0 + matrix_[0][1]*v1 + matrix_[0][2]*v2;
    rot1 = matrix_[1][0]*v0 + matrix_[1][1]*v1 + matrix_[1][2]*v2;
    rot2 = matrix_[2][0]*v0 + matrix_[2][1]*v1 + matrix_[2][2]*v2;
  }

  // rotates vector {rot0, rot1, rot2} to {v0, v1, v2}
  void inv_rot(enzo_float rot0, enzo_float rot1, enzo_float rot2,
	       enzo_float &v0, enzo_float &v1, enzo_float &v2)
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

  virtual enzo_float operator() (enzo_float x0, enzo_float x1, enzo_float x2)=0;
};

//----------------------------------------------------------------------

class LinearScalarInit : public ScalarInit{

public:
  LinearScalarInit(enzo_float background, enzo_float eigenvalue,
		   enzo_float amplitude, enzo_float lambda)
    : background_(background),
      eigenvalue_(eigenvalue),
      amplitude_(amplitude),
      lambda_(lambda)
  {}

  enzo_float operator() (enzo_float x0, enzo_float x1, enzo_float x2){
    return (background_ + amplitude_ * eigenvalue_
	    * std::cos(x0*2.*cello::pi/lambda_));
  }

protected:
  // value of the background state
  enzo_float background_;
  // value of the eigenvalue for a given mode
  enzo_float eigenvalue_;
  // perturbation amplitude
  enzo_float amplitude_;
  // wavelength
  enzo_float lambda_;
};

//----------------------------------------------------------------------

class RotatedScalarInit : public ScalarInit {
public:
  RotatedScalarInit(enzo_float a, enzo_float b, ScalarInit *inner)
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

  enzo_float operator() (enzo_float x, enzo_float y, enzo_float z)
  {
    enzo_float x0,x1,x2;
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

  virtual void operator() (enzo_float x0, enzo_float x1, enzo_float x2,
			   enzo_float &v0, enzo_float &v1, enzo_float &v2)=0;
};

//----------------------------------------------------------------------

class LinearVectorInit : public VectorInit {
public:
  LinearVectorInit(enzo_float back0, enzo_float back1, enzo_float back2,
		   enzo_float ev0, enzo_float ev1, enzo_float ev2,
		   enzo_float amplitude, enzo_float lambda)
    : VectorInit(),
      back0_(back0), back1_(back1), back2_(back2),
      ev0_(ev0), ev1_(ev1), ev2_(ev2),
      amplitude_(amplitude),
      lambda_(lambda)
  {}

  void operator() (enzo_float x0, enzo_float x1, enzo_float x2,
		   enzo_float &v0, enzo_float &v1, enzo_float &v2)
  {
    v0 = (back0_ + amplitude_ * ev0_* std::cos(x0*2.*cello::pi/lambda_));
    v1 = (back1_ + amplitude_ * ev1_* std::cos(x0*2.*cello::pi/lambda_));
    v2 = (back2_ + amplitude_ * ev2_* std::cos(x0*2.*cello::pi/lambda_));
  }

protected:
  // value of the background state
  enzo_float back0_, back1_, back2_;
  // value of the eigenvalue for a given mode
  enzo_float ev0_, ev1_, ev2_;
  // perturbation amplitude
  enzo_float amplitude_;
  // wavelength
  enzo_float lambda_;
};

//----------------------------------------------------------------------

class RotatedVectorInit : public VectorInit {
public:
  RotatedVectorInit(enzo_float a, enzo_float b, VectorInit *inner)
    : VectorInit()
  {
    rot_ = new Rotation(a,b);
    inner_ = inner;
  }

  ~RotatedVectorInit()
  {
    delete inner_;
    delete rot_;
  }

  void operator() (enzo_float x, enzo_float y, enzo_float z,
		   enzo_float &v0, enzo_float &v1, enzo_float &v2)
  {
    enzo_float x0,x1,x2, rot0, rot1, rot2;
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
  LinearVectorPotentialInit(enzo_float back_B0, enzo_float back_B1,
			    enzo_float back_B2, enzo_float ev_B1,
			    enzo_float ev_B2, enzo_float amplitude,
			    enzo_float lambda)
    :VectorInit(),
     back_B0_(back_B0), back_B1_(back_B1), back_B2_(back_B2),
     ev_B1_(ev_B1), ev_B2_(ev_B2),
     amplitude_(amplitude),
     lambda_(lambda)
  {}

  void operator() (enzo_float x0, enzo_float x1, enzo_float x2,
		   enzo_float &v0, enzo_float &v1, enzo_float &v2)
  {
    v0 = (x2 * amplitude_ * ev_B1_ * std::cos(2. * cello::pi * x0/ lambda_) -
	  x1 * amplitude_ * ev_B2_ * std::cos(2. * cello::pi * x0/ lambda_));
    v1 = back_B2_ * x0; // For Gardiner & Stone (2008) setup, should be 0
    v2 = back_B0_ * x1 - back_B1_ * x0;
  }

protected:
  // value of the magnetic field background state
  enzo_float back_B0_, back_B1_, back_B2_;
  // magnetic field eigenvalue for a given mode (its always 0, along dim 0)
  enzo_float ev_B1_, ev_B2_;
  // perturbation amplitude
  enzo_float amplitude_;
  // wavelength
  enzo_float lambda_;
};

//----------------------------------------------------------------------

class CircAlfvenMomentumInit : public VectorInit {
public:
  // The absence of ev_B0 is intentional
  CircAlfvenMomentumInit(enzo_float p0,
			 enzo_float lambda)
    :VectorInit(),
     p0_(p0),
     lambda_(lambda)
  {}

  void operator() (enzo_float x0, enzo_float x1, enzo_float x2,
		   enzo_float &v0, enzo_float &v1, enzo_float &v2)
  {
    v0 = p0_;
    v1 = 0.1 * std::sin(2. * cello::pi * x0/ lambda_);
    v2 = 0.1 * std::cos(2. * cello::pi * x0/ lambda_);
  }

protected:
  // momentum along p0
  enzo_float p0_;
  // wavelength
  enzo_float lambda_;
};

//----------------------------------------------------------------------

class CircAlfvenVectorPotentialInit : public VectorInit {
public:
  // The absence of ev_B0 is intentional
  CircAlfvenVectorPotentialInit(enzo_float lambda)
    :VectorInit(),
     lambda_(lambda)
  {}

  void operator() (enzo_float x0, enzo_float x1, enzo_float x2,
		   enzo_float &v0, enzo_float &v1, enzo_float &v2)
  {
    v0 = (x2 * 0.1 * std::sin(2. * cello::pi * x0/ lambda_) -
	  x1 * 0.1 * std::cos(2. * cello::pi * x0/ lambda_));
    v1 = 0.0;
    v2 = x1;
  }

protected:
  // wavelength
  enzo_float lambda_;
};

//----------------------------------------------------------------------

class MeshPos {

public:
  MeshPos(Block *block)
  {
    Field field  = block->data()->field();

    // lower left edge of active region
    double xmb,ymb,zmb;
    block->data()->lower(&xmb,&ymb,&zmb);

    x_left_edge_ = (enzo_float)xmb;
    y_left_edge_ = (enzo_float)ymb;
    z_left_edge_ = (enzo_float)zmb;

    // cell widths
    double xpb,ypb,zpb,hx,hy,hz;
    block->data()->upper(&xpb,&ypb,&zpb);
    field.cell_width(xmb,xpb,&hx,
		     ymb,ypb,&hy,
		     zmb,zpb,&hz);


    dx_ = (enzo_float) hx;
    dy_ = (enzo_float) hy;
    dz_ = (enzo_float) hz;

    // ghost depth 
    field.ghost_depth(0,&gx_,&gy_,&gz_);
  }

  // compute positions at cell-centers
  enzo_float x(int k, int j, int i){
    return x_left_edge_ + dx_* (0.5+(enzo_float)(i-gx_));}
  enzo_float y(int k, int j, int i){
    return y_left_edge_ + dy_* (0.5+(enzo_float)(j-gy_));}
  enzo_float z(int k, int j, int i){
    return z_left_edge_ + dz_* (0.5+(enzo_float)(k-gz_));}

  // computes the x, y, or z position for the cell corner at the
  // (i-1/2, j,k),  (i,j-1/2,k) or (i,j, k-1/2)
  enzo_float x_face(int k, int j, int i){
    return x_left_edge_ + dx_* (enzo_float)(i-gx_);}
  enzo_float y_face(int k, int j, int i){
    return y_left_edge_ + dy_* (enzo_float)(j-gy_);}
  enzo_float z_face(int k, int j, int i){
    return z_left_edge_ + dz_* (enzo_float)(k-gz_);}

  enzo_float dx() {return dx_;}
  enzo_float dy() {return dy_;}
  enzo_float dz() {return dz_;}
  
private:
  // the starting position of the edge of the active grid
  enzo_float x_left_edge_, y_left_edge_, z_left_edge_;
  // ghost depth (and the index at the start of the active region)
  int gx_, gy_, gz_;
  enzo_float dx_, dy_, dz_;
};

//----------------------------------------------------------------------

void bfieldi_helper_(EnzoArray<enzo_float,3> &bfield,
		     EnzoArray<enzo_float,3> &Aj,
		     EnzoArray<enzo_float,3> &Ak, int dim, enzo_float dj,
		     enzo_float dk)
{

  // Aj is centered on edges of i and k dimension but cell-centered along j
  // Ak is centered on edges of i and j dimension but cell-centered along k
  // Aj_right(k,j,i) = Aj(k+1/2,j,i-1/2)    Aj_left(k,j,i) = Aj(k-1/2,j,i-1/2)
  // Ak_right(k,j,i) = Ak(k,j+1/2,i-1/2)    Ak_left(k,j,i) = Ak(k,j-1/2,i-1/2)

  // get dimensions of interface field
  int fc_mx = bfield.shape(2);
  int fc_my = bfield.shape(1);
  int fc_mz = bfield.shape(0);

  EnzoArray<enzo_float,3> Ak_left, Ak_right,Aj_right, Aj_left;

  if (dim == 0){
    // Aj_right = Ay(iz+1/2,iy,ix-1/2), Ak_right = Az(iz, iy+1/2,ix-1/2)
    Aj_right = Aj.subarray(ESlice(1, fc_mz+1), ESlice(0, fc_my),
			   ESlice(0, fc_mx));
    Ak_right = Ak.subarray(ESlice(0, fc_mz), ESlice(1, fc_my+1),
			   ESlice(0, fc_mx));
  } else if (dim == 1){
    // Aj_right = Az(iz,iy-1/2,ix+1/2), Ak_right = Ax(iz+1/2,iy-1/2,ix)
    Aj_right = Aj.subarray(ESlice(0, fc_mz), ESlice(0, fc_my),
			   ESlice(1, fc_mx+1));
    Ak_right = Ak.subarray(ESlice(1, fc_mz+1), ESlice(0, fc_my),
			   ESlice(0, fc_mx));
  } else {
    // Aj_right = Ax(iz-1/2,iy+1/2,ix), Ak_right = Ay(iz-1/2,iy,ix+1/2)
    Aj_right = Aj.subarray(ESlice(0, fc_mz), ESlice(1, fc_my+1),
			   ESlice(0, fc_mx));
    Ak_right = Ak.subarray(ESlice(0, fc_mz), ESlice(0, fc_my),
			   ESlice(1, fc_mx+1));
  }

  Aj_left = Aj.subarray(ESlice(0, fc_mz), ESlice(0, fc_my), ESlice(0, fc_mx));
  Ak_left = Ak.subarray(ESlice(0, fc_mz), ESlice(0, fc_my), ESlice(0, fc_mx));

  for (int iz=0;iz<fc_mz;iz++){
    for (int iy=0; iy<fc_my; iy++){
      for (int ix=0; ix<fc_mx; ix++){
	// Bi(k,j,i-1/2) =
	//    ( Ak(    k, j+1/2, i-1/2) - Ak(    k, j-1/2, i-1/2) )/dj -
	//    ( Aj(k+1/2,     j, i-1/2) - Aj(k+1/2,     j, i-1/2) )/dk
	bfield(iz,iy,ix) = ((Ak_right(iz,iy,ix) - Ak_left(iz,iy,ix))/dj -
			    (Aj_right(iz,iy,ix) - Aj_left(iz,iy,ix))/dk);
	
      }
    }
  }
}

//----------------------------------------------------------------------

// Components of the magnetic field from the vector potential
// Bx(k, j, i-1/2) =
//    ( Az(    k, j+1/2, i-1/2) - Az(    k, j-1/2, i-1/2) )/dy -
//    ( Ay(k+1/2,     j, i-1/2) - Ay(k+1/2,     j, i-1/2) )/dz
// By(k, j-1/2, i) =
//    ( Ax(k+1/2, j-1/2,     i) - Ax(k+1/2, j-1/2,     i) )/dz -
//    ( Az(    k, j-1/2, i+1/2) - Az(    k, j-1/2, i-1/2) )/dx
// Bz(k-1/2, j, i) =
//    ( Ay(k-1/2,     j, i+1/2) - Ay(k-1/2,     j, i-1/2) )/dx -
//    ( Ax(k-1/2, j+1/2,     i) - Ax(k-1/2, j-1/2,     i) )/dy
void setup_bfield(Block * block, VectorInit *a, MeshPos &pos,
		  int mx, int my, int mz)
{
  EnzoArray<enzo_float,3> bfieldi_x, bfieldi_y, bfieldi_z;
  EnzoFieldArrayFactory array_factory(block);
  bfieldi_x = array_factory.from_name("bfieldi_x");
  bfieldi_y = array_factory.from_name("bfieldi_y");
  bfieldi_z = array_factory.from_name("bfieldi_z");

  // cell widths
  enzo_float dx = pos.dx();
  enzo_float dy = pos.dy();
  enzo_float dz = pos.dz();

  // allocate corner-centered arrays for the magnetic vector potentials
  // Ax, Ay, and Az are always cell-centered along the x, y, and z dimensions,
  // respectively
  EnzoArray<enzo_float,3> Ax(mz+1,my+1,mx);
  EnzoArray<enzo_float,3> Ay(mz+1,my,mx+1);
  EnzoArray<enzo_float,3> Az(mz,my+1,mx+1);

  // Compute the Magnetic Vector potential at all points on the grid
  for (int iz=0; iz<mz+1; iz++){
    for (int iy=0; iy<my+1; iy++){
      for (int ix=0; ix<mx+1; ix++){

	enzo_float temp_ax, temp_ay, temp_az;

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
  bfieldi_helper_(bfieldi_x, Ay, Az, 0, dy,dz);
  bfieldi_helper_(bfieldi_y, Az, Ax, 1, dz,dx);
  bfieldi_helper_(bfieldi_z, Ax, Ay, 2, dx,dy);

  // Compute the Cell-Centered B-fields
  EnzoInitialBCenter::initialize_bfield_center(block);

}

//----------------------------------------------------------------------

void setup_fluid(Block *block, ScalarInit *density_init,
		 ScalarInit *total_energy_init, 
		 VectorInit *momentum_init,
		 MeshPos &pos, int mx, int my, int mz, enzo_float gamma)
{
  EnzoArray<enzo_float,3> density, pressure;
  EnzoFieldArrayFactory array_factory(block);
  density = array_factory.from_name("density");
  pressure = array_factory.from_name("pressure");

  EnzoArray<enzo_float,3> velocity_x, velocity_y, velocity_z;
  velocity_x = array_factory.from_name("velocity_x");
  velocity_y = array_factory.from_name("velocity_y");
  velocity_z = array_factory.from_name("velocity_z");

  // Required for computing Pressure
  EnzoArray<enzo_float,3> bfieldc_x, bfieldc_y, bfieldc_z;
  bfieldc_x = array_factory.from_name("bfield_x");
  bfieldc_y = array_factory.from_name("bfield_y");
  bfieldc_z = array_factory.from_name("bfield_z");

  for (int iz=0; iz<mz; iz++){
    for (int iy=0; iy<my; iy++){
      for (int ix=0; ix<mx; ix++){
	enzo_float x,y,z;
	enzo_float rho, px, py, pz, etot, p2, b2;
	x = pos.x(iz,iy,ix); y = pos.y(iz,iy,ix); z = pos.z(iz,iy,ix);
	rho = (*density_init)(x,y,z);
	density(iz,iy,ix) = rho;
	(*momentum_init)(x, y, z, px, py, pz);
	velocity_x(iz,iy,ix) = px/rho;
	velocity_y(iz,iy,ix) = py/rho;
	velocity_z(iz,iy,ix) = pz/rho;

	etot = (*total_energy_init)(x,y,z);
	p2 = (px*px + py*py + pz*pz);
	b2 = (bfieldc_x(iz,iy,ix)*bfieldc_x(iz,iy,ix) +
	      bfieldc_y(iz,iy,ix)*bfieldc_y(iz,iy,ix) +
	      bfieldc_z(iz,iy,ix)*bfieldc_z(iz,iy,ix));

	pressure(iz,iy,ix) = (gamma-1.)*(etot - 0.5*(p2/density(iz,iy,ix) +
						     b2));
      }
      fflush(stdout);
    }
  }
}

//----------------------------------------------------------------------

void setup_circ_polarized_alfven(Block *block, ScalarInit *density_init, 
				 VectorInit *momentum_init,
				 MeshPos &pos, int mx, int my, int mz,
				 enzo_float gamma)
{
  // this function directly sets pressure = 0.1 by hand
  // Gardiner & Stone (2008) explicitly as states that the truncation error of
  // B_perp**2/P is important
  EnzoArray<enzo_float,3> density, pressure;
  EnzoFieldArrayFactory array_factory(block);
  density = array_factory.from_name("density");
  pressure = array_factory.from_name("pressure");

  EnzoArray<enzo_float,3> velocity_x, velocity_y, velocity_z;
  velocity_x = array_factory.from_name("velocity_x");
  velocity_y = array_factory.from_name("velocity_y");
  velocity_z = array_factory.from_name("velocity_z");

  for (int iz=0; iz<mz; iz++){
    for (int iy=0; iy<my; iy++){
      for (int ix=0; ix<mx; ix++){
	enzo_float x,y,z;
	enzo_float rho, px, py, pz;
	x = pos.x(iz,iy,ix); y = pos.y(iz,iy,ix); z = pos.z(iz,iy,ix);
	rho = (*density_init)(x,y,z);
	density(iz,iy,ix) = rho;
	(*momentum_init)(x, y, z, px, py, pz);
	velocity_x(iz,iy,ix) = px/rho;
	velocity_y(iz,iy,ix) = py/rho;
	velocity_z(iz,iy,ix) = pz/rho;

	pressure(iz,iy,ix) = 0.1;
      }
      fflush(stdout);
    }
  }
}


//----------------------------------------------------------------------

EnzoInitialInclinedWave::EnzoInitialInclinedWave(int cycle, double time,
					     double alpha, double beta,
					     double gamma, double amplitude,
					     double lambda, bool pos_vel,
					     std::string wave_type) throw()
  : Initial (cycle,time),
    alpha_(alpha),
    beta_(beta),
    gamma_(gamma),
    amplitude_(amplitude),
    lambda_(lambda),
    pos_vel_(pos_vel),
    wave_type_(wave_type)
{
  ASSERT("EnzoInitialInclinedWave",
	 ("Invalid wave_type specified (must be 'fast', 'slow', 'alfven', "
	  "'entropy', or 'circ_alfven'"),
	 (wave_type_ == "fast") || (wave_type_ == "alfven") ||
	 (wave_type_ == "slow") || (wave_type_ == "entropy") ||
	 (wave_type_ == "circ_alfven"));
}

//----------------------------------------------------------------------

void EnzoInitialInclinedWave::enforce_block(Block * block,
					  const Hierarchy * hierarchy) throw()
{
  // Set up the test problem
  // This will only work for the VLCT method for an adiabatic fluid
  // does NOT support AMR
  // (PPML initial conditions are much more complicated)

  Field field  = block->data()->field();
  ASSERT("EnzoInitialInclinedWave",
	 "Insufficient number of fields",
	 field.field_count() >= 8);
  
  ScalarInit *density_init = NULL;
  ScalarInit *total_energy_init = NULL;
  VectorInit *momentum_init = NULL;
  VectorInit *a_init = NULL;

  prepare_initializers_(&density_init, &total_energy_init, &momentum_init,
			&a_init);

  // Try to load the dimensions, again
  int mx,my,mz;
  const int id = field.field_id ("density");
  field.dimensions (id,&mx,&my,&mz);
  
  // Initialize object to compute positions
  MeshPos pos(block);

  setup_bfield(block, a_init, pos, mx, my, mz);
  if (wave_type_ != "circ_alfven"){
    setup_fluid(block, density_init, total_energy_init, momentum_init,
		pos, mx, my, mz, gamma_);
    
  } else {
    setup_circ_polarized_alfven(block, density_init, momentum_init,
				pos, mx, my, mz, gamma_);
  }

  delete density_init;
  delete momentum_init;
  delete total_energy_init;
  delete a_init;
}

//----------------------------------------------------------------------

void EnzoInitialInclinedWave::prepare_initializers_(ScalarInit **density_init,
						  ScalarInit **etot_init,
						  VectorInit **momentum_init,
						  VectorInit **a_init)
{
  enzo_float lambda = lambda_;
  

  if (wave_type_ == "circ_alfven"){
    *density_init = new LinearScalarInit(1.0,0.0,0.0,lambda);
    // we don't use etot_init for circularly polarized alfven waves
    *etot_init = new LinearScalarInit(0.66,0.0,0.0,lambda);
    *momentum_init = new RotatedVectorInit(alpha_,beta_,
					   new CircAlfvenMomentumInit(0.0,
								      lambda));
    *a_init = new RotatedVectorInit(alpha_,beta_,
				    new CircAlfvenVectorPotentialInit(lambda));

  } else {
    // wsign indicates direction of propogation.
    enzo_float wsign = 1.;
    if (!pos_vel_){
      wsign = -1;
    }
    // initialize the background values:
    // density = 1, pressure = 1/gamma, mom1 = 0, mom2 = 0, B0 = 1, B1=1.5, B2=0
    // For entropy wave, mom0 = 1. Otherwise, mom0 = 0.
    enzo_float density_back = 1;
    enzo_float mom0_back = 0;
    enzo_float mom1_back = 0;
    enzo_float mom2_back = 0;
    enzo_float b0_back = 1.;
    enzo_float b1_back = 1.5;
    enzo_float b2_back = 0.0;
    // etot = pressure/(gamma-1)+0.5*rho*v^2+.5*B^2
    enzo_float etot_back = (1./gamma_)/(gamma_-1.)+1.625;
    if (wave_type_ == "entropy"){
      mom0_back = wsign;
      b0_back*=wsign;
      b1_back*=wsign;
      b2_back*=wsign;
      etot_back += 0.5;
    }


    // Get the eigenvalues for density, velocity, and etot
    enzo_float density_ev, mom0_ev, mom1_ev, mom2_ev, etot_ev, b1_ev, b2_ev;
    // intentionally omit b0_ev

    if (wave_type_ == "fast"){
      enzo_float coef = 0.5/std::sqrt(5.);
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
      enzo_float coef = 0.5/std::sqrt(5.);
      density_ev = 4. *coef;
      mom0_ev = wsign*2. * coef;
      mom1_ev = wsign*4. * coef;
      mom2_ev = 0;
      etot_ev = 3. * coef;
      b1_ev = -2. * coef;
      b2_ev = 0;
    } else {
      // (wave_type_ == "entropy")
      enzo_float coef = 0.5;
      density_ev = 2. *coef;
      mom0_ev = 2. * coef;
      mom1_ev = 0;
      mom2_ev = 0;
      etot_ev = 1. * coef;
      b1_ev = 0;
      b2_ev = 0;
    }
  
    enzo_float amplitude = amplitude_;
    // Now allocate the actual Initializers
  
    *density_init = new RotatedScalarInit(alpha_, beta_,
					  new LinearScalarInit(density_back,
							       density_ev,
							       amplitude,
							       lambda));
      
    *etot_init = new RotatedScalarInit(alpha_, beta_,
				       new LinearScalarInit(etot_back,
							    etot_ev,
							    amplitude,
							    lambda));

    *momentum_init = new RotatedVectorInit(alpha_, beta_,
					   new LinearVectorInit(mom0_back,
								mom1_back,
								mom2_back,
								mom0_ev,
								mom1_ev,
								mom2_ev,
								amplitude,
								lambda));
    *a_init = new RotatedVectorInit(alpha_, beta_,
				    new LinearVectorPotentialInit(b0_back,
								  b1_back,
								  b2_back,
								  b1_ev,
								  b2_ev,
								  amplitude,
								  lambda));
  }
}
