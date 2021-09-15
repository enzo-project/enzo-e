// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodTurbulenceOU.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @author   Alexei Kritsuk (akritsuk@ucsd.edu)
/// @date     2021-09-01
/// @brief    Implements the EnzoMethodTurbulenceOU class

#include "cello.hpp"
#include "enzo.hpp"
#include "charm_enzo.hpp"

#define CALL_FORTRAN
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

  refresh->add_all_fields();

  // Call fortran initializer
  
  double domain_size[3] =
    { domain_upper[0]-domain_lower[0],
      domain_upper[1]-domain_lower[1],
      domain_upper[2]-domain_lower[2] };

  int is_root = (CkMyPe() == 0) ? 1 : 0;
  int rank = cello::rank();
  int iapply_injection_rate = apply_injection_rate_ ? 1 : 0;
  int iread_sol = read_sol_ ? 1 : 0;
#ifdef CALL_FORTRAN    
  FORTRAN_NAME(cello_init_turbulence_ou)
    (&is_root,                  // (*)
     &rank,                     // (*)
     domain_size,               // (*)
     &gamma_,                   // (*)
     &iapply_injection_rate,    // (*)
     &cooling_term_,            // (*)
     &injection_rate_,          // (*)
     &kfi_,                     // (*)
     &kfa_,                     // (*)
     &mach_,                    // (*)
     &iread_sol,                // (*)
     &sol_weight_);             // (*)
#endif  
}

//----------------------------------------------------------------------

void EnzoMethodTurbulenceOU::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

}

//----------------------------------------------------------------------

void EnzoMethodTurbulenceOU::compute ( Block * block) throw()
{
  long double r_gvld1[5] = {4.0, 0.0, 0.0, 0.0, 0.0};
  if (block->is_leaf()) {
#ifdef CALL_FORTRAN    
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
    CHECK_FIELD(field_work_1,"work_1");
    double * field_work_2 =  (double *)field.values("work_2");
    CHECK_FIELD(field_work_2,"work_2");
    double * field_work_3 =  (double *)field.values("work_3");
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
       r_gv
       );
    EnzoMethodTurbulenceOU::iupdate_phases_ = 0;

    FIELD_STATS("density force stop ",field_density,mx,my,mz,gx,gy,gz);
#endif    

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
  
  long double r_avld1[2+1] = {2.0, 0.0,0.0};

  if (enzo_block->is_leaf()) {
  
#ifdef CALL_FORTRAN    
    double * w;
    int iupdate_sol = update_solution_ ? 1 : 0;
    int iapply_cooling = apply_cooling_ ? 1 : 0;
    int iapply_forcing = apply_forcing_ ? 1 : 0;
    int iapply_injection_rate = apply_injection_rate_ ? 1 : 0;
    int cello_apply_injection_rate;
    double r_gv[4] = {r_gvld0[0],r_gvld0[1],r_gvld0[2],r_gvld0[3]};
    double r_av[2] = {0.0,0.0};

    Field field = enzo_block->data()->field();
    double * field_density     = (double *)field.values("density");
    CHECK_FIELD(field_density,"density");
    double * field_momentum_x  = (double *)field.values("velocity_x");
    CHECK_FIELD(field_momentum_x,"velocity_x");
    double * field_momentum_y  = (double *)field.values("velocity_y");
    CHECK_FIELD(field_momentum_y,"velocity_y");
    double * field_momentum_z  = (double *)field.values("velocity_z");
    CHECK_FIELD(field_momentum_z,"velocity_z");
    double * field_jacobian    = (double *)field.values("jacobian");
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
    CHECK_FIELD(field_work_1,"work_1");
    double * field_work_2 = (double *)field.values("work_2");
    CHECK_FIELD(field_work_2,"work_2");
    double * field_work_3 = (double *)field.values("work_3");
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

    for (int i=0; i<4+1; i++) r_avld1[i+1] = r_av[i];

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
#endif  
  }
  
  CkCallback callback = CkCallback
    (CkIndex_EnzoBlock::r_method_turbulence_ou_update(nullptr), 
     enzo::block_array());
  n = 2;
  enzo_block->contribute((n+1)*sizeof(long double), r_avld1, 
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
  long double r_avld0[2];
  int id=0;
  int n = data[id++];
  ASSERT1 ("EnzoMethodTurbulenceOU::compute_shift()",
           "Expected length 2 but actual length %d",
           n,(n==2));
  for (int i=0; i<n; i++) r_avld0[i] = data[id++];
  delete msg;

#ifdef CALL_FORTRAN    
  double dt   = cello::simulation()->dt();
  int iupdate_sol = update_solution_ ? 1 : 0;
  int iapply_cooling = apply_cooling_ ? 1 : 0;
  int iapply_forcing = apply_forcing_ ? 1 : 0;
  int iapply_injection_rate = apply_injection_rate_ ? 1 : 0;
  int cello_apply_injection_rate;
  double r_av[2] = {r_avld0[0],r_avld0[1]};

  Field field = enzo_block->data()->field();
  double * field_density     = (double *)field.values("density");
  CHECK_FIELD(field_density,"density");
  double * field_momentum_x  = (double *)field.values("velocity_x");
  CHECK_FIELD(field_momentum_x,"velocity_x");
  double * field_momentum_y  = (double *)field.values("velocity_y");
  CHECK_FIELD(field_momentum_y,"velocity_y");
  double * field_momentum_z  = (double *)field.values("velocity_z");
  CHECK_FIELD(field_momentum_z,"velocity_z");
  double * field_energy      = (double *)field.values("total_energy");
  CHECK_FIELD(field_energy,"total_energy");
  double * resid_density     = (double *)field.values("resid_density");
  CHECK_FIELD(resid_density,"resid_density");
  double * resid_momentum_x  = (double *)field.values("resid_velocity_x");
  CHECK_FIELD(resid_momentum_x,"resid_velocity_x");
  double * resid_momentum_y  = (double *)field.values("resid_velocity_y");
  CHECK_FIELD(resid_momentum_y,"resid_velocity_y");
  double * resid_momentum_z  = (double *)field.values("resid_velocity_z");
  CHECK_FIELD(resid_momentum_z,"resid_velocity_z");
  double * resid_energy      = (double *)field.values("resid_total_energy");
  CHECK_FIELD(resid_energy,"resid_energy");
  double * field_temperature = (double *)field.values("temperature");
  CHECK_FIELD(field_temperature,"temperature");

  int mx,my,mz;
  int nx,ny,nz;
  field.dimensions(0,&mx,&my,&mz);
  field.size(&nx,&ny,&nz);
  const int m = mx*my*mz;

  double * turbAcc = new double [4*m];
  double * field_work_4 = (double *)field.values("work_4");
  CHECK_FIELD(field_work_4,"work_4");
  double * field_work_5 = (double *)field.values("work_5");
  CHECK_FIELD(field_work_5,"work_5");
  double * field_work_6 = (double *)field.values("work_6");
  CHECK_FIELD(field_work_6,"work_6");
  double * field_work_7 = (double *)field.values("work_7");
  CHECK_FIELD(field_work_7,"work_7");
  std::copy_n (field_work_4,m,turbAcc+0*m);
  std::copy_n (field_work_5,m,turbAcc+1*m);
  std::copy_n (field_work_6,m,turbAcc+2*m);
  std::copy_n (field_work_7,m,turbAcc+3*m);
  
  std::fill_n(resid_density,m,0.0);
  std::fill_n(resid_momentum_x,m,0.0);
  std::fill_n(resid_momentum_y,m,0.0);
  std::fill_n(resid_momentum_z,m,0.0);
  std::fill_n(resid_energy,m,0.0);

        // convert to conservative form
  for (int iz=0; iz<mz; iz++) {
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {
        int i=ix+mx*(iy+my*iz);
        double density = field_density[i];
        field_momentum_x[i] *= density;
        field_momentum_y[i] *= density;
        field_momentum_z[i] *= density;
        field_energy[i]     *= density;
      }
    }
  }

  // Restore work array
  double * field_work_1 =  (double *)field.values("work_1");
  CHECK_FIELD(field_work_1,"work_1");
  double * field_work_2 =  (double *)field.values("work_2");
  CHECK_FIELD(field_work_2,"work_2");
  double * field_work_3 =  (double *)field.values("work_3");
  CHECK_FIELD(field_work_3,"work_3");
  double * array_work = new double [3*m];
  std::copy_n (field_work_1,m,array_work+0*m);
  std::copy_n (field_work_2,m,array_work+1*m);
  std::copy_n (field_work_3,m,array_work+2*m);

    int gx,gy,gz;
    field.ghost_depth(0,&gx,&gy,&gz);
    FIELD_STATS("density update start",field_density,mx,my,mz,gx,gy,gz);
    FIELD_STATS("momentum_x update start",field_momentum_x,mx,my,mz,gx,gy,gz);
    FIELD_STATS("momentum_y update start",field_momentum_y,mx,my,mz,gx,gy,gz);
    FIELD_STATS("momentum_z update start",field_momentum_z,mx,my,mz,gx,gy,gz);
    FIELD_STATS("energy update start",field_energy,mx,my,mz,gx,gy,gz);
    FIELD_STATS("work_1 update start",field_work_1,mx,my,mz,gx,gy,gz);
    FIELD_STATS("work_2 update start",field_work_2,mx,my,mz,gx,gy,gz);
    FIELD_STATS("work_3 update start",field_work_3,mx,my,mz,gx,gy,gz);
    // index of first non-ghost value
    const int i0 = gx + mx*(gy + my*gz);
  FORTRAN_NAME(turbforceupdate)
    (&mx, &my, &mz,
     &nx, &ny, &nz,
     field_density+i0,
     field_momentum_x+i0,
     field_momentum_y+i0,
     field_momentum_z+i0,
     field_energy+i0,
     resid_density+i0,
     resid_momentum_x+i0,
     resid_momentum_y+i0,
     resid_momentum_z+i0,
     resid_energy+i0,
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
    FIELD_STATS("energy update stop ",field_energy,mx,my,mz,gx,gy,gz);
#endif

    std::copy_n (turbAcc+0*m,m,field_work_4);
    std::copy_n (turbAcc+1*m,m,field_work_5);
    std::copy_n (turbAcc+2*m,m,field_work_6);
    std::copy_n (turbAcc+3*m,m,field_work_7);

  // revert to code units
  for (int iz=0; iz<mz; iz++) {
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {
        int i=ix+mx*(iy+my*iz);
        double density_inverse = 1.0/field_density[i];
        field_momentum_x[i] *= density_inverse;
        field_momentum_y[i] *= density_inverse;
        field_momentum_z[i] *= density_inverse;
        field_energy[i]     *= density_inverse;
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

  enzo_block->compute_done();
}

