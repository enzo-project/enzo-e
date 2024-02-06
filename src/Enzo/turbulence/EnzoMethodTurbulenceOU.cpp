// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodTurbulenceOU.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @author   Alexei Kritsuk (akritsuk@ucsd.edu)
/// @date     2021-09-01
/// @brief    Implements the EnzoMethodTurbulenceOU class

#include "cello.hpp"
#include "enzo.hpp"
#include "charm_enzo.hpp"

// #define DEBUG_FIELDS

#ifdef DEBUG_FIELDS
#   define CHECK_FIELD(VALUES,NAME)             \
  ASSERT1("CHECK_FIELD",                        \
          "Field %s must be defined",           \
          NAME,                                 \
          (VALUES != nullptr));

#   define FIELD_STATS(NAME,VALUES,mx,my,mz,gx,gy,gz)           \
  {                                                             \
   double avg=0.0, max=-1.0, min=1e9;                           \
   int count=0;                                                 \
   for (int iz=gz; iz<mz-gz; iz++) {                            \
     for (int iy=gy; iy<my-gy; iy++) {                          \
       for (int ix=gx; ix<mx-gx; ix++) {                        \
         const int i=ix+mx*(iy+my*iz);                          \
         avg += VALUES[i];                                      \
         max = std::max(max,VALUES[i]);                         \
         min = std::min(min,VALUES[i]);                         \
         count++;                                               \
       }                                                        \
     }                                                          \
   }                                                            \
   avg /= count;                                                \
   CkPrintf ("FIELD_STATS %s  %g %g %g\n",NAME,min,avg,max);    \
  }
#else
#   define CHECK_FIELD(VALUES,NAME) /* ... */
#   define FIELD_STATS(NAME,VALUES,mx,my,mz,gx,gy,gz) /* ... */
#endif

int EnzoMethodTurbulenceOU::iupdate_phases_ = 1;

//----------------------------------------------------------------------
EnzoMethodTurbulenceOU::EnzoMethodTurbulenceOU
(double gamma,
 const double domain_lower[3],
 const double domain_upper[3],
 bool apply_cooling,
 bool apply_forcing,
 bool apply_injection_rate,
 int cooling_term,
 double hc_alpha,
 double hc_sigma,
 double injection_rate,
 double kfi,
 double kfa,
 double mach,
 int olap,
 bool read_sol,
 double sol_weight,
 double totemp,
 bool update_solution)
  : Method(),
    gamma_(gamma),
    apply_cooling_(apply_cooling),
    apply_forcing_(apply_forcing),
    apply_injection_rate_(apply_injection_rate),
    cooling_term_(cooling_term),
    hc_alpha_(hc_alpha),
    hc_sigma_(hc_sigma),
    injection_rate_(injection_rate),
    kfi_(kfi),
    kfa_(kfa),
    mach_(mach),
    olap_(0),
    read_sol_(read_sol),
    sol_weight_(sol_weight),
    totemp_(totemp),
    update_solution_(update_solution)
{
  // Initialize default Refresh object

  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);

  cello::define_field("jacobian");
  cello::define_field("resid_density");
  cello::define_field("resid_energy");
  cello::define_field("resid_total_energy");
  cello::define_field("resid_velocity_x");
  cello::define_field("resid_velocity_y");
  cello::define_field("resid_velocity_z");
  cello::define_field("acceleration_x");
  cello::define_field("acceleration_y");
  cello::define_field("acceleration_z");
  cello::define_field("energy");

  //  refresh->add_all_fields();

  refresh->add_field ("acceleration_x");
  refresh->add_field ("acceleration_y");
  refresh->add_field ("acceleration_z");
  refresh->add_field ("bfieldx");
  refresh->add_field ("bfieldx_rx");
  refresh->add_field ("bfieldx_ry");
  refresh->add_field ("bfieldx_rz");
  refresh->add_field ("bfieldy");
  refresh->add_field ("bfieldy_rx");
  refresh->add_field ("bfieldy_ry");
  refresh->add_field ("bfieldy_rz");
  refresh->add_field ("bfieldz");
  refresh->add_field ("bfieldz_rx");
  refresh->add_field ("bfieldz_ry");
  refresh->add_field ("bfieldz_rz");
  refresh->add_field ("density");
  refresh->add_field ("velocity_x");
  refresh->add_field ("velocity_y");
  refresh->add_field ("velocity_z");

  // Call fortran initializer

  double domain_size[3] =
    { domain_upper[0]-domain_lower[0],
      domain_upper[1]-domain_lower[1],
      domain_upper[2]-domain_lower[2] };

  int is_root = (CkMyPe() == 0) ? 1 : 0;
  int rank = cello::rank();
  int iapply_injection_rate = apply_injection_rate_ ? 1 : 0;
  int iread_sol = read_sol_ ? 1 : 0;

  double weight_norm = 0;
  FORTRAN_NAME(cello_init_turbulence_ou)
    (&is_root,
     &rank,
     domain_size,
     &gamma_,
     &iapply_injection_rate,
     &cooling_term_,
     &injection_rate_,
     &kfi_,
     &kfa_,
     &mach_,
     &iread_sol,
     &sol_weight_,
     &weight_norm);
  static bool first_call = true;
  if (first_call && CkMyPe() == 0) {
    CkPrintf ("WeightNorm = %24.18g\n",weight_norm);
    first_call = false;
  }

}

//----------------------------------------------------------------------

void EnzoSimulation::get_turbou_state()
{
  // Create scalar variables for storing state for checkpoint/restart

  // Allocate state arrays if needed
  int n_size_double, n_size_int;
  FORTRAN_NAME(cello_turbou_state_size)(&n_size_double,&n_size_int);
  if (turbou_real_state_.size() < n_size_double) {
    turbou_real_state_.resize(n_size_double);
  }
  if (turbou_int_state_.size() < n_size_int) {
    turbou_int_state_.resize(n_size_int);
  }

  // Save state to arrays
  FORTRAN_NAME(cello_get_turbou_state)
    (turbou_real_state_.data(), turbou_int_state_.data());
}

//----------------------------------------------------------------------

void EnzoSimulation::put_turbou_state()
{
  // Create scalar variables for storing state for checkpoint/restart

  // Save arrays to state
  if (turbou_real_state_.size() > 0) {
    FORTRAN_NAME(cello_put_turbou_state)
      (turbou_real_state_.data(), turbou_int_state_.data());
  }
}

//----------------------------------------------------------------------

void EnzoMethodTurbulenceOU::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);
  p | gamma_;
  p | apply_cooling_;
  p | apply_forcing_;
  p | apply_injection_rate_;
  p | cooling_term_;
  p | hc_alpha_;
  p | hc_sigma_;
  p | injection_rate_;
  p | kfi_;
  p | kfa_;
  p | mach_;
  p | olap_;
  p | read_sol_;
  p | sol_weight_;
  p | totemp_;
  p | update_solution_;
}

//----------------------------------------------------------------------

void EnzoMethodTurbulenceOU::compute ( Block * block) throw()
{
  long double r_gvld1[5] = {4.0, 0.0, 0.0, 0.0, 0.0};

  if (block->is_leaf()) {

    // Restore Fortran phase arrays if restarting
    if (enzo::config()->initial_restart &&
        EnzoMethodTurbulenceOU::iupdate_phases_ == 1) {
      enzo::simulation()->put_turbou_state();
    }

    Field field = block->data()->field();
    int mx,my,mz;
    int nx,ny,nz;
    field.dimensions(0,&mx,&my,&mz);
    field.size(&nx,&ny,&nz);
    double * w;
    double * jac;
    int iupdate_sol = update_solution_ ? 1 : 0;
    int iapply_cooling = apply_cooling_ ? 1 : 0;
    int iapply_forcing = apply_forcing_ ? 1 : 0;
    int iapply_injection_rate = apply_injection_rate_ ? 1 : 0;
    int cello_apply_injection_rate;
    double r_gv[4] = {0.0,0.0,0.0,0.0};

    double time = cello::simulation()->time();
    double dt   = cello::simulation()->dt();

    double * field_density = (double *)field.values("density");
    CHECK_FIELD(field_density,"density");

    double xm,ym,zm;
    double xp,yp,zp;
    block->lower(&xm,&ym,&zm);
    block->upper(&xp,&yp,&zp);
    double hx,hy,hz;
    field.cell_width (xm,xp,&hx, ym,yp,&hy, zm,zp,&hz);
    int gx,gy,gz;
    field.ghost_depth(0,&gx,&gy,&gz);

    // Restore work array
    double * field_work_1 =  (double *)field.values("work_1");
    double * field_work_2 =  (double *)field.values("work_2");
    double * field_work_3 =  (double *)field.values("work_3");
    CHECK_FIELD(field_work_1,"work_1");
    CHECK_FIELD(field_work_2,"work_2");
    CHECK_FIELD(field_work_3,"work_3");
    const int m = mx*my*mz;
    double * array_work = new double [3*m];
    std::copy_n (field_work_1,m,array_work+0*m);
    std::copy_n (field_work_2,m,array_work+1*m);
    std::copy_n (field_work_3,m,array_work+2*m);

    double * grid = new double [3*m];
    for (int iz=0; iz<mz; iz++) {
      double z = zm+(iz-gz+0.5)*hz;
      for (int iy=0; iy<my; iy++) {
        double y = ym+(iy-gy+0.5)*hy;
        for (int ix=0; ix<mx; ix++) {
          double x = xm+(ix-gx+0.5)*hx; // ix=4 -> xm + 0.5*hx
          int i=3*(ix+mx*(iy+my*iz));
          grid[i]=x;
          grid[i+1]=y;
          grid[i+2]=z;
        }
      }
    }

    FIELD_STATS("density force start",field_density,mx,my,mz,gx,gy,gz);
    // index of first non-ghost value
    const int i0 = gx + mx*(gy + my*gz);

    FORTRAN_NAME(turbforceou)
      (&mx, &my, &mz,
       &nx, &ny, &nz,
       field_density+i0,
       grid+3*i0,
       array_work+i0, &time, &dt,
       &iupdate_sol,
       &iapply_cooling,
       &iapply_forcing,
       &iapply_injection_rate,
       &EnzoMethodTurbulenceOU::iupdate_phases_,
       &cooling_term_,
       &gamma_,
       &injection_rate_,
       &olap_,
       r_gv);

    EnzoMethodTurbulenceOU::iupdate_phases_ = 0;

    FIELD_STATS("density force stop ",field_density,mx,my,mz,gx,gy,gz);

    for (int i=0; i<4; i++) r_gvld1[i+1] = r_gv[i];

    // Save work array
    std::copy_n (array_work+0*m,m,field_work_1);
    std::copy_n (array_work+1*m,m,field_work_2);
    std::copy_n (array_work+2*m,m,field_work_3);
    FIELD_STATS("work_1 force",field_work_1,mx,my,mz,gx,gy,gz);
    FIELD_STATS("work_2 force",field_work_2,mx,my,mz,gx,gy,gz);
    FIELD_STATS("work_3 force",field_work_3,mx,my,mz,gx,gy,gz);
    delete [] array_work;
    delete [] grid;

  }

  CkCallback callback = CkCallback
    (CkIndex_EnzoBlock::r_method_turbulence_ou_shift(nullptr),
     enzo::block_array());
  int n = 4;

  block->contribute((n+1)*sizeof(long double), r_gvld1,
                    sum_long_double_n_type, callback);

}

//----------------------------------------------------------------------

void EnzoBlock::r_method_turbulence_ou_shift(CkReductionMsg *msg)
{
  EnzoMethodTurbulenceOU * method = static_cast<EnzoMethodTurbulenceOU*> (this->method());
  method->compute_shift(this, msg);
}

//----------------------------------------------------------------------

void EnzoMethodTurbulenceOU::compute_shift
(EnzoBlock * enzo_block, CkReductionMsg * msg)
{
  long double * data = (long double *) msg->getData();
  long double r_gvld0[4];
  int id=0;
  int n = data[id++];
  ASSERT1 ("EnzoMethodTurbulenceOU::compute_shift()",
           "Expected length 4 but actual length %d",
           n,(n==4));
  for (int i=0; i<n; i++) r_gvld0[i] = data[id++];
  delete msg;

  const int num_reduce = 4;
  long double r_avld1[num_reduce+1];
  std::fill_n(r_avld1,num_reduce+1,0.0);
  r_avld1[0] = num_reduce;

  if (enzo_block->is_leaf()) {

    double * w;
    int iupdate_sol = update_solution_ ? 1 : 0;
    int iapply_cooling = apply_cooling_ ? 1 : 0;
    int iapply_forcing = apply_forcing_ ? 1 : 0;
    int iapply_injection_rate = apply_injection_rate_ ? 1 : 0;
    int cello_apply_injection_rate;
    double r_gv[4] = {r_gvld0[0],r_gvld0[1],r_gvld0[2],r_gvld0[3]};
    double r_av[num_reduce];
    std::fill_n(r_av,num_reduce,0.0);

    Field field = enzo_block->data()->field();
    double * field_density     = (double *)field.values("density");
    double * field_momentum_x  = (double *)field.values("velocity_x");
    double * field_momentum_y  = (double *)field.values("velocity_y");
    double * field_momentum_z  = (double *)field.values("velocity_z");
    double * field_jacobian    = (double *)field.values("jacobian");
    CHECK_FIELD(field_density,"density");
    CHECK_FIELD(field_momentum_x,"velocity_x");
    CHECK_FIELD(field_momentum_y,"velocity_y");
    CHECK_FIELD(field_momentum_z,"velocity_z");
    CHECK_FIELD(field_jacobian,"jacobian");

    // convert to conservative form
    int mx,my,mz;
    int nx,ny,nz;
    field.dimensions(0,&mx,&my,&mz);
    field.size(&nx,&ny,&nz);
    for (int iz=0; iz<mz; iz++) {
      for (int iy=0; iy<my; iy++) {
        for (int ix=0; ix<mx; ix++) {
          int i=ix+mx*(iy+my*iz);
          double density = field_density[i];
          field_momentum_x[i] *= density;
          field_momentum_y[i] *= density;
          field_momentum_z[i] *= density;
        }
      }
    }

    // Restore work array
    double * field_work_1 = (double *)field.values("work_1");
    double * field_work_2 = (double *)field.values("work_2");
    double * field_work_3 = (double *)field.values("work_3");
    CHECK_FIELD(field_work_1,"work_1");
    CHECK_FIELD(field_work_2,"work_2");
    CHECK_FIELD(field_work_3,"work_3");
    const int m = mx*my*mz;
    double * array_work = new double [3*m];
    std::copy_n (field_work_1,m,array_work+0*m);
    std::copy_n (field_work_2,m,array_work+1*m);
    std::copy_n (field_work_3,m,array_work+2*m);

    int gx,gy,gz;
    field.ghost_depth(0,&gx,&gy,&gz);
    FIELD_STATS("density shift start",field_density,mx,my,mz,gx,gy,gz);
    FIELD_STATS("momentum_x shift start",field_momentum_x,mx,my,mz,gx,gy,gz);
    FIELD_STATS("momentum_y shift start",field_momentum_y,mx,my,mz,gx,gy,gz);
    FIELD_STATS("momentum_z shift start",field_momentum_z,mx,my,mz,gx,gy,gz);
    FIELD_STATS("work_1 shift start",field_work_1,mx,my,mz,gx,gy,gz);
    FIELD_STATS("work_2 shift start",field_work_2,mx,my,mz,gx,gy,gz);
    FIELD_STATS("work_3 shift start",field_work_3,mx,my,mz,gx,gy,gz);
    // index of first non-ghost value
    const int i0 = gx + mx*(gy + my*gz);
    FORTRAN_NAME(turbforceshift)
      (&mx, &my, &mz,
       &nx, &ny, &nz,
       field_density+i0,
       field_momentum_x+i0,
       field_momentum_y+i0,
       field_momentum_z+i0,
       field_jacobian+i0,
       array_work+i0,
       &iupdate_sol,
       &iapply_injection_rate,
       &olap_,
       &injection_rate_,
       r_gv,
       r_av);

    FIELD_STATS("density shift stop ",field_density,mx,my,mz,gx,gy,gz);
    FIELD_STATS("momentum_x shift stop ",field_momentum_x,mx,my,mz,gx,gy,gz);
    FIELD_STATS("momentum_y shift stop ",field_momentum_y,mx,my,mz,gx,gy,gz);
    FIELD_STATS("momentum_z shift stop ",field_momentum_z,mx,my,mz,gx,gy,gz);

    for (int i=0; i<2; i++) r_avld1[i+1] = r_av[i];

    // revert to code units
    for (int iz=0; iz<mz; iz++) {
      for (int iy=0; iy<my; iy++) {
        for (int ix=0; ix<mx; ix++) {
          int i=ix+mx*(iy+my*iz);
          double density_inverse = 1.0/field_density[i];
          field_momentum_x[i] *= density_inverse;
          field_momentum_y[i] *= density_inverse;
          field_momentum_z[i] *= density_inverse;
        }
      }
    }

    // Save work array
    std::copy_n (array_work+0*m,m,field_work_1);
    std::copy_n (array_work+1*m,m,field_work_2);
    std::copy_n (array_work+2*m,m,field_work_3);
    FIELD_STATS("work_1 shift stop ",field_work_1,mx,my,mz,gx,gy,gz);
    FIELD_STATS("work_2 shift stop ",field_work_2,mx,my,mz,gx,gy,gz);
    FIELD_STATS("work_3 shift stop ",field_work_3,mx,my,mz,gx,gy,gz);
    delete [] array_work;

  }

  Field field = enzo_block->data()->field();
  const bool include_ppml_divb = (field.values("bfieldx") != nullptr);
  if (include_ppml_divb && enzo_block->is_leaf()) {
    enzo_float * bfield[3] =
      { (enzo_float *) field.values("bfieldx"),
        (enzo_float *) field.values("bfieldy"),
        (enzo_float *) field.values("bfieldz") };

    enzo_float * bfield_rx[3] =
      { (enzo_float *) field.values("bfieldx_rx"),
        (enzo_float *) field.values("bfieldy_rx"),
        (enzo_float *) field.values("bfieldz_rx") };

    enzo_float * bfield_ry[3] =
      { (enzo_float *) field.values("bfieldx_ry"),
        (enzo_float *) field.values("bfieldy_ry"),
        (enzo_float *) field.values("bfieldz_ry") };

    enzo_float * bfield_rz[3] =
      { (enzo_float *) field.values("bfieldx_rz"),
        (enzo_float *) field.values("bfieldy_rz"),
        (enzo_float *) field.values("bfieldz_rz") };

    CHECK_FIELD(bfield[0],"bfieldx");
    CHECK_FIELD(bfield[1],"bfieldy");
    CHECK_FIELD(bfield[2],"bfieldz");
    CHECK_FIELD(bfield_rx[0],"bfieldx_rx");
    CHECK_FIELD(bfield_ry[1],"bfieldx_ry");
    CHECK_FIELD(bfield_rz[2],"bfieldx_rz");
    CHECK_FIELD(bfield_rx[0],"bfieldy_rx");
    CHECK_FIELD(bfield_ry[1],"bfieldy_ry");
    CHECK_FIELD(bfield_rz[2],"bfieldy_rz");
    CHECK_FIELD(bfield_rx[0],"bfieldz_rx");
    CHECK_FIELD(bfield_ry[1],"bfieldz_ry");
    CHECK_FIELD(bfield_rz[2],"bfieldz_rz");

    int mx,my,mz;
    int nx,ny,nz;
    int gx,gy,gz;
    field.size(&nx,&ny,&nz);
    field.ghost_depth(0,&gx,&gy,&gz);
    field.dimensions(0,&mx,&my,&mz);
    const int dx = 1;
    const int dy = mx;
    const int dz = mx*my;
    r_avld1[3] = 0.0;
    r_avld1[4] = 0.0;
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
          int i = (ix+gx) + mx*((iy+gy) + my*(iz+gz));
          ++r_avld1[4];
          r_avld1[3] += fabs(bfield_rx[0][i+dx] - bfield_rx[0][i-dx] +
                             bfield_rx[0][i+dx+dy] - bfield_rx[0][i-dx+dy] +
                             bfield_rx[0][i+dx+dz] - bfield_rx[0][i-dx+dz] +
                             bfield_rx[0][i+dx+dy+dz] - bfield_rx[0][i-dx+dy+dz]
                             +
                             bfield_ry[1][i+dy] - bfield_ry[1][i-dy] +
                             bfield_ry[1][i+dy+dx] - bfield_ry[1][i-dy+dx] +
                             bfield_ry[1][i+dy+dz] - bfield_ry[1][i-dy+dz] +
                             bfield_ry[1][i+dy+dx+dz] - bfield_ry[1][i-dy+dx+dz]
                             +
                             bfield_rz[2][i+dz] - bfield_rz[2][i-dz] +
                             bfield_rz[2][i+dz+dx] - bfield_rz[2][i-dz+dx] +
                             bfield_rz[2][i+dz+dy] - bfield_rz[2][i-dz+dy] +
                             bfield_rz[2][i+dz+dx+dy] - bfield_rz[2][i-dz+dx+dy]);
        }
      }
    }
  } else {
    r_avld1[3] = r_avld1[4] = 0.0;
  }

  CkCallback callback = CkCallback
    (CkIndex_EnzoBlock::r_method_turbulence_ou_update(nullptr),
     enzo::block_array());
  enzo_block->contribute((num_reduce+1)*sizeof(long double), r_avld1,
                         sum_long_double_n_type, callback);
}

void EnzoBlock::r_method_turbulence_ou_update(CkReductionMsg *msg)
{
  EnzoMethodTurbulenceOU * method = static_cast<EnzoMethodTurbulenceOU*> (this->method());
  method->compute_update(this,msg);
}

//----------------------------------------------------------------------

void EnzoMethodTurbulenceOU::compute_update
(EnzoBlock * enzo_block, CkReductionMsg *msg)
{
  long double * data = (long double *) msg->getData();
  const int num_reduce = 4;
  long double r_avld0[num_reduce+1];
  int id=0;
  int n = data[id++];
  ASSERT2 ("EnzoMethodTurbulenceOU::compute_shift()",
           "Expected length %d but actual length %d",
           num_reduce,n,(n==num_reduce));
  for (int i=0; i<n; i++) r_avld0[i] = data[id++];
  delete msg;

  double dt   = cello::simulation()->dt();
  int iupdate_sol = update_solution_ ? 1 : 0;
  int iapply_cooling = apply_cooling_ ? 1 : 0;
  int iapply_forcing = apply_forcing_ ? 1 : 0;
  int iapply_injection_rate = apply_injection_rate_ ? 1 : 0;
  int cello_apply_injection_rate;
  double r_av[2] = {r_avld0[0],r_avld0[1]};

  Field field = enzo_block->data()->field();
  double * field_density     = (double *)field.values("density");
  double * field_momentum_x  = (double *)field.values("velocity_x");
  double * field_momentum_y  = (double *)field.values("velocity_y");
  double * field_momentum_z  = (double *)field.values("velocity_z");
  double * field_energy_total      = (double *)field.values("total_energy");
  double * resid_density     = (double *)field.values("resid_density");
  double * resid_momentum_x  = (double *)field.values("resid_velocity_x");
  double * resid_momentum_y  = (double *)field.values("resid_velocity_y");
  double * resid_momentum_z  = (double *)field.values("resid_velocity_z");
  double * resid_energy_total      = (double *)field.values("resid_total_energy");
  double * field_temperature = (double *)field.values("temperature");
  CHECK_FIELD(field_density,"density");
  CHECK_FIELD(field_momentum_x,"velocity_x");
  CHECK_FIELD(field_momentum_y,"velocity_y");
  CHECK_FIELD(field_momentum_z,"velocity_z");
  CHECK_FIELD(field_energy_total,"total_energy");
  CHECK_FIELD(resid_density,"resid_density");
  CHECK_FIELD(resid_momentum_x,"resid_velocity_x");
  CHECK_FIELD(resid_momentum_y,"resid_velocity_y");
  CHECK_FIELD(resid_momentum_z,"resid_velocity_z");
  CHECK_FIELD(resid_energy_total,"resid_energy_total");
  CHECK_FIELD(field_temperature,"temperature");

  const bool using_ppml = (field.values("bfieldx") != nullptr);
  if (using_ppml && enzo_block->index().is_root()) {
    Monitor * monitor = cello::monitor();
    monitor->print ("Method","<|div(b)|>  %.17Lg",
		    r_avld0[2]/8.0/r_avld0[3]);
  }
  int mx,my,mz;
  int nx,ny,nz;
  field.dimensions(0,&mx,&my,&mz);
  field.size(&nx,&ny,&nz);
  const int m = mx*my*mz;

  double * turbAcc = new double [4*m];
  double * field_acceleration_x = (double *)field.values("acceleration_x");
  double * field_acceleration_y = (double *)field.values("acceleration_y");
  double * field_acceleration_z = (double *)field.values("acceleration_z");
  double * field_energy = (double *)field.values("energy");
  CHECK_FIELD(field_acceleration_x,"acceleration_x");
  CHECK_FIELD(field_acceleration_y,"acceleration_y");
  CHECK_FIELD(field_acceleration_z,"acceleration_z");
  CHECK_FIELD(field_energy,"energy");

  std::fill_n(resid_density,m,0.0);
  std::fill_n(resid_momentum_x,m,0.0);
  std::fill_n(resid_momentum_y,m,0.0);
  std::fill_n(resid_momentum_z,m,0.0);
  std::fill_n(resid_energy_total,m,0.0);

  // convert to conservative form
  for (int iz=0; iz<mz; iz++) {
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {
        int i=ix+mx*(iy+my*iz);
        double density = field_density[i];
        field_momentum_x[i] *= density;
        field_momentum_y[i] *= density;
        field_momentum_z[i] *= density;
        field_energy_total[i] *= density;
      }
    }
  }
  // Update off-centered velocities (as momenta) if they exist
  double * field_px[3] = {nullptr, nullptr, nullptr};
  double * field_py[3] = {nullptr, nullptr, nullptr};
  double * field_pz[3] = {nullptr, nullptr, nullptr};
  if (using_ppml) {
    field_px[0] = (double *)field.values("velox_rx");
    field_px[1] = (double *)field.values("velox_ry");
    field_px[2] = (double *)field.values("velox_rz");
    field_py[0] = (double *)field.values("veloy_rx");
    field_py[1] = (double *)field.values("veloy_ry");
    field_py[2] = (double *)field.values("veloy_rz");
    field_pz[0] = (double *)field.values("veloz_rx");
    field_pz[1] = (double *)field.values("veloz_ry");
    field_pz[2] = (double *)field.values("veloz_rz");
    for (int iz=0; iz<mz; iz++) {
      for (int iy=0; iy<my; iy++) {
        for (int ix=0; ix<mx; ix++) {
          int i=ix+mx*(iy+my*iz);
          double density = field_density[i];
          for (int k=0; k<3; k++) {
            field_px[k][i] *= density;
            field_py[k][i] *= density;
            field_pz[k][i] *= density;
          }
        }
      }
    }
  }

  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);
  const int i0 = gx + mx*(gy + my*gz);

  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
        int i=ix+mx*(iy+my*iz);
        turbAcc[0+4*i] = field_acceleration_x[i+i0];
        turbAcc[1+4*i] = field_acceleration_y[i+i0];
        turbAcc[2+4*i] = field_acceleration_z[i+i0];
        turbAcc[3+4*i] = field_energy[i+i0];
      }
    }
  }
  // Restore work array
  double * field_work_1 =  (double *)field.values("work_1");
  double * field_work_2 =  (double *)field.values("work_2");
  double * field_work_3 =  (double *)field.values("work_3");
  CHECK_FIELD(field_work_1,"work_1");
  CHECK_FIELD(field_work_2,"work_2");
  CHECK_FIELD(field_work_3,"work_3");
  double * array_work = new double [3*m];
  std::copy_n (field_work_1,m,array_work+0*m);
  std::copy_n (field_work_2,m,array_work+1*m);
  std::copy_n (field_work_3,m,array_work+2*m);

  FIELD_STATS("density update start",field_density,mx,my,mz,gx,gy,gz);
  FIELD_STATS("momentum_x update start",field_momentum_x,mx,my,mz,gx,gy,gz);
  FIELD_STATS("momentum_y update start",field_momentum_y,mx,my,mz,gx,gy,gz);
  FIELD_STATS("momentum_z update start",field_momentum_z,mx,my,mz,gx,gy,gz);
  FIELD_STATS("energy_total update start",field_energy_total,mx,my,mz,gx,gy,gz);
  FIELD_STATS("work_1 update start",field_work_1,mx,my,mz,gx,gy,gz);
  FIELD_STATS("work_2 update start",field_work_2,mx,my,mz,gx,gy,gz);
  FIELD_STATS("work_3 update start",field_work_3,mx,my,mz,gx,gy,gz);
  // index of first non-ghost value

  int have_faces = using_ppml ? 1 : 0;
  FORTRAN_NAME(turbforceupdate)
    (&mx, &my, &mz,
     &nx, &ny, &nz,
     field_density+i0,
     field_momentum_x+i0,
     field_momentum_y+i0,
     field_momentum_z+i0,
     &have_faces,
     field_px[0]+i0,
     field_py[0]+i0,
     field_pz[0]+i0,
     field_px[1]+i0,
     field_py[1]+i0,
     field_pz[1]+i0,
     field_px[2]+i0,
     field_py[2]+i0,
     field_pz[2]+i0,
     field_energy_total+i0,
     resid_density+i0,
     resid_momentum_x+i0,
     resid_momentum_y+i0,
     resid_momentum_z+i0,
     resid_energy_total+i0,
     field_temperature+i0,
     array_work+i0, &dt,
     turbAcc,
     &iupdate_sol,
     &iapply_injection_rate,
     &injection_rate_,
     &cooling_term_,
     &iapply_cooling,
     &gamma_,
     &hc_alpha_,
     &hc_sigma_,
     &totemp_,
     r_av );
  FIELD_STATS("density update stop ",field_density,mx,my,mz,gx,gy,gz);
  FIELD_STATS("momentum_x update stop ",field_momentum_x,mx,my,mz,gx,gy,gz);
  FIELD_STATS("momentum_y update stop ",field_momentum_y,mx,my,mz,gx,gy,gz);
  FIELD_STATS("momentum_z update stop ",field_momentum_z,mx,my,mz,gx,gy,gz);
  FIELD_STATS("energy_total update stop ",field_energy_total,mx,my,mz,gx,gy,gz);

  // revert to code units
  for (int iz=0; iz<mz; iz++) {
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {
        int i=ix+mx*(iy+my*iz);
        double density_inverse = 1.0/field_density[i];
        field_momentum_x[i] *= density_inverse;
        field_momentum_y[i] *= density_inverse;
        field_momentum_z[i] *= density_inverse;
        field_energy_total[i]     *= density_inverse;
      }
    }
  }
  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
        int i=ix+mx*(iy+my*iz);
        field_acceleration_x[i+i0] = turbAcc[0+4*i];
        field_acceleration_y[i+i0] = turbAcc[1+4*i];
        field_acceleration_z[i+i0] = turbAcc[2+4*i];
        field_energy[i+i0] = turbAcc[3+4*i];
      }
    }
  }

  if (using_ppml) {
    for (int iz=0; iz<mz; iz++) {
      for (int iy=0; iy<my; iy++) {
        for (int ix=0; ix<mx; ix++) {
          int i=ix+mx*(iy+my*iz);
          double density_inverse = 1.0/field_density[i];
          for (int k=0; k<3; k++) {
            field_px[k][i] *= density_inverse;
            field_py[k][i] *= density_inverse;
            field_pz[k][i] *= density_inverse;
          }
        }
      }
    }
  }
  // Save work array
  std::copy_n (array_work+0*m,m,field_work_1);
  std::copy_n (array_work+1*m,m,field_work_2);
  std::copy_n (array_work+2*m,m,field_work_3);
  FIELD_STATS("work_1 update stop ",field_work_1,mx,my,mz,gx,gy,gz);
  FIELD_STATS("work_2 update stop ",field_work_2,mx,my,mz,gx,gy,gz);
  FIELD_STATS("work_3 update stop ",field_work_3,mx,my,mz,gx,gy,gz);
  delete [] array_work;
  delete [] turbAcc;

  EnzoMethodTurbulenceOU::iupdate_phases_ = 1;

  // Save Fortran phase and random number state for checkpointing
  enzo::simulation()->get_turbou_state();

  compute_reductions_(enzo_block);
}

//----------------------------------------------------------------------

void EnzoMethodTurbulenceOU::compute_reductions_(EnzoBlock * enzo_block)
{
  Field field = enzo_block->data()->field();

  const EnzoConfig * enzo_config = enzo::config();

  EnzoComputeTemperature compute_temperature
    (enzo::fluid_props(), enzo_config->physics_cosmology);

  compute_temperature.compute(enzo_block);

  enzo_float * d  = (enzo_float *) field.values("density");
  enzo_float * vx = (enzo_float *) field.values("velox");
  enzo_float * vy = (enzo_float *) field.values("veloy");
  enzo_float * vz = (enzo_float *) field.values("veloz");
  enzo_float * ax = (enzo_float *) field.values("acceleration_x");
  enzo_float * ay = (enzo_float *) field.values("acceleration_y");
  enzo_float * az = (enzo_float *) field.values("acceleration_z");
  enzo_float * t  = (enzo_float *) field.values("temperature");
  enzo_float * bx = (enzo_float *) field.values("bfieldx");
  enzo_float * by = (enzo_float *) field.values("bfieldy");
  enzo_float * bz = (enzo_float *) field.values("bfieldz");
  enzo_float * e  = (enzo_float *) field.values("energy"); // turbAcc(4)
  enzo_float * p  = (enzo_float *) field.values("pressure");

  int nx,ny,nz;
  int gx,gy,gz;
  field.size(&nx,&ny,&nz);
  field.ghost_depth(0,&gx,&gy,&gz);
  const int mx = nx + 2*gx;
  const int my = ny + 2*gy;
  const int mz = nz + 2*gz;
  const int dx = 1;
  const int dy = mx;
  const int dz = mx*my;

  // Allocate zero field if rank < 3 to remove if's from loops
  const int rank = cello::rank();
  enzo_float * zero = rank <3 ? new enzo_float[mx*my*mz] : nullptr;
  if (rank < 2) {
    vy = zero;
    ay = zero;
    by = zero;
  }
  if (rank < 3) {
    vz = zero;
    az = zero;
    bz = zero;
  }

  // Allocate and initialize reduction array
  const int n = 30;
  long double g[n+1];
  std::fill_n(g,n+1,0.0);
  g[0] = n;

  const double c2 = 29979245800.0 * 29979245800.0;
  if (enzo_block->is_leaf()) {

    double xm,ym,zm;
    double xp,yp,zp;
    enzo_block->lower(&xm,&ym,&zm);
    enzo_block->upper(&xp,&yp,&zp);
    double hx,hy,hz;
    field.cell_width (xm,xp,&hx, ym,yp,&hy, zm,zp,&hz);
    double hxi=1.0/hx;
    double hyi=(cello::rank() >= 2) ? 1.0/hy : 0.0;
    double hzi=(cello::rank() >= 3) ? 1.0/hz : 0.0;

    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {

	  int i = (ix+gx) + mx*((iy+gy) + my*(iz+gz));

          g[ 1] += d[i];
          g[ 2] += log(d[i]);
          g[ 3] += d[i]*d[i];
          g[ 4] += d[i]*log(d[i]);
          g[ 5] += vx[i];
          g[ 6] += vy[i];
          g[ 7] += vz[i];
          const double v2 = vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
          g[ 8] += d[i]*v2;
          g[ 9] += v2;
          const double di = 1.0/d[i];
          g[10] += vx[i]*di;
          g[11] += vy[i]*di;
          g[12] += vz[i]*di;
          // compute partial derivatives
          const double dvxdx = 0.5*(vx[i+dx] - vx[i-dx])*hxi;
          const double dvxdy = 0.5*(vx[i+dy] - vx[i-dy])*hyi;
          const double dvxdz = 0.5*(vx[i+dz] - vx[i-dz])*hzi;
          const double dvydx = 0.5*(vy[i+dx] - vy[i-dx])*hxi;
          const double dvydy = 0.5*(vy[i+dy] - vy[i-dy])*hyi;
          const double dvydz = 0.5*(vy[i+dz] - vy[i-dz])*hzi;
          const double dvzdx = 0.5*(vz[i+dx] - vz[i-dx])*hxi;
          const double dvzdy = 0.5*(vz[i+dy] - vz[i-dy])*hyi;
          const double dvzdz = 0.5*(vz[i+dz] - vz[i-dz])*hzi;
          // compute div velocity
          const double vdiv = (dvxdx + dvydy + dvzdz);
          g[13] += vdiv;
          g[14] += vdiv*vdiv;
          // compute velocity curl components
          const double cvx = dvzdy - dvydz;
          const double cvy = dvxdz - dvzdx;
          const double cvz = dvydx - dvxdy;
          g[15] += cvx*cvx+cvy*cvy+cvz*cvz;
          // g[16] += p[i];
          // g[17] += p[i]/d[i]
          g[18] += bx[i];
          g[19] += by[i];
          g[20] += bz[i];
          const double b2 = bx[i]*bx[i] + by[i]*by[i] + bz[i]*bz[i];
          g[21] += b2;
          g[22] += b2*di;
          // compute bfield curl components
          const double dbxdx = 0.5*(bx[i+dx] - bx[i-dx])*hxi;
          const double dbydy = 0.5*(by[i+dy] - by[i-dy])*hyi;
          const double dbzdz = 0.5*(bz[i+dz] - bz[i-dz])*hzi;
          g[23] += abs(dbxdx) + abs(dbydy) + abs(dbzdz);
          // or is it g[23] += abs(dbxdx+dbydy+dbzdz) ?
          g[24] += e[i];
          g[25] += v2/c2;
          g[26] += p[i]*vdiv;
          g[27] += p[i]*vdiv*vdiv;
          g[28] += hx*hy*hz;
          g[29] += dbxdx + dbydy + dbzdz;
	}
      }
    }
    // Scale by volume
    const double v = hx*hy*hz;
    for (int i=1; i<=30; i++) g[i] *= v;
  }

  delete [] zero;

  CkCallback callback (CkIndex_EnzoBlock::r_method_turbulence_ou_end(NULL),
		       enzo_block->proxy_array());
  enzo_block->contribute((n+1)*sizeof(long double),g,sum_long_double_n_type,callback);
}

//----------------------------------------------------------------------

void EnzoBlock::r_method_turbulence_ou_end(CkReductionMsg * msg)
{
  performance_start_(perf_compute,__FILE__,__LINE__);
  static_cast<EnzoMethodTurbulenceOU*> (method())->compute_reduce (this,msg);
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoMethodTurbulenceOU::compute_reduce
(EnzoBlock * enzo_block,CkReductionMsg * msg)
{
  long double * g = (long double *)msg->getData();

  Data * data = enzo_block->data();
  Field field = data->field();


  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  int n = nx*ny*nz;

  double dt = ((Block*)enzo_block)->dt();

  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);
  const int rank = ((mz == 1) ? ((my == 1) ? 1 : 2) : 3);

  double xdm,ydm,zdm;
  data->lower(&xdm,&ydm,&zdm);
  double xdp,ydp,zdp;
  data->upper(&xdp,&ydp,&zdp);

  if (enzo_block->index().is_root()) {

    Monitor * m = cello::monitor();
#define FMT "%.17Lg"
#define NAME "Turbulence"
    
    int i=1;
    m->print (NAME,"stat %02d rho " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d log(rho) " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d rho^2 " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d rho*log(rho) " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d rho*u " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d rho*v " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d rho*w " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d rho*(u^2 + v^2 + w^2) " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d u^2 + v^2 + w^2 " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d u " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d v " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d w " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d div u " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d (div u)^2 " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d Omega^2 " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d p " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d p/rho " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d Bx " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d By " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d Bz " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d Bx^2 + By^2 + Bz^2 " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d (Bx^2 + By^2 + Bz^2)/rho " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d |div b| " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d E_turb " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d (u^2 + v^2 + w^2)/c^2 " FMT, i,sqrt(g[i])); i++;
    m->print (NAME,"stat %02d p*(div u) " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d (p*(div u))^2 " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d total volume " FMT, i,g[i]); i++;
    m->print (NAME,"stat %02d div b " FMT, i,g[i]); i++;
  }
  enzo_block->compute_done();
}


