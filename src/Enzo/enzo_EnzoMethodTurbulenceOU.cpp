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
    olap_(olap),
    read_sol_(read_sol),
    sol_weight_(sol_weight),
    totemp_(totemp),
    update_solution_(update_solution)
{
  // Initialize default Refresh object

  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);

  refresh->add_field("acceleration_x");
  refresh->add_field("acceleration_y");
  refresh->add_field("acceleration_z");

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
  CkPrintf ("DEBUG_TURBULENCE_OU compute() enter %s\n",block->name().c_str());
  fflush(stdout);
  long double r_gv128[5] = {4.0, 0.0, 0.0, 0.0, 0.0};
  if (block->is_leaf()) {
#ifdef CALL_FORTRAN    
    int mx,my,mz;
    block->data()->field().dimensions(0,&mx,&my,&mz);
    double * w;
    double * jac;
    double * turbAcc;
    double dt;
    int iupdate_sol = update_solution_ ? 1 : 0;
    int iapply_cooling = apply_cooling_ ? 1 : 0;
    int iapply_forcing = apply_forcing_ ? 1 : 0;
    CkPrintf ("iapply_forcing %d\n",iapply_forcing);
    int iapply_injection_rate = apply_injection_rate_ ? 1 : 0;
    int cello_apply_injection_rate;
    double cello_injection_rate;
    double r_gv[4] = {0.0,0.0,0.0,0.0};

    Field field = block->data()->field();
    double * field_density = cello::field(block,"density");
   
    double xm,ym,zm;
    double xp,yp,zp;
    block->lower(&xm,&ym,&zm);
    block->upper(&xp,&yp,&zp);
    double hx,hy,hz;
    field.cell_width
      (xm,xp,&hx, ym,yp,&hy, zm,zp,&hz);
    int gx,gy,gz;
    field.ghost_depth(0,&gx,&gy,&gz);

    // Restore work array
    double * field_work_1 =  cello::field(block,"work_1");
    double * field_work_2 =  cello::field(block,"work_2");
    double * field_work_3 =  cello::field(block,"work_3");
    const int m = mx*my*mz;
    double * array_work = new double [3*m];
    std::copy_n (field_work_1,m,array_work+0*m);
    std::copy_n (field_work_2,m,array_work+1*m);
    std::copy_n (field_work_3,m,array_work+2*m);

    double * grid = new double [3*m];
    for (int iz=0; iz<mz; iz++) {
      double z = (iz-gz+0.5)*hz;
      for (int iy=0; iy<my; iy++) {
        double y = (iy-gy+0.5)*hy;
        for (int ix=0; ix<mx; ix++) {
          double x = (ix-gx+0.5)*hx;
          int i=3*(ix+mx*(iy+my*iz));
          grid[i]=x;
          grid[i+1]=y;
          grid[i+2]=z;
        }
      }
    }
    
    CkPrintf ("array_work forcing in %g %g %g\n",array_work[0],array_work[m],array_work[2*m]);

    FORTRAN_NAME(turbforceou)
      (&mx, &my, &mz,          // (*)
       field_density,          // (*)
       grid,                   // (*)
       array_work,  &dt,       // ( )
       turbAcc, &iupdate_sol,  // ( )
       &iapply_cooling,        // ( )
       &iapply_forcing,        // ( )
       &iapply_injection_rate, // ( )
       &cooling_term_,         // ( )
       &gamma_,                // ( )
       &injection_rate_,       // ( )
       &olap_,                 // ( )
       r_gv                    // ( )
       );

    CkPrintf ("array_work forcing out %g %g %g\n",array_work[0],array_work[m],array_work[2*m]);
#endif    

    for (int i=0; i<4; i++) r_gv128[i+1] = r_gv[i];

    // Save work array
    std::copy_n (array_work+0*m,m,field_work_1);
    std::copy_n (array_work+1*m,m,field_work_2);
    std::copy_n (array_work+2*m,m,field_work_3);
    delete [] array_work;
    delete [] grid;

  }
  
  CkCallback callback = CkCallback
    (CkIndex_EnzoBlock::r_method_turbulence_ou_shift(nullptr), 
     enzo::block_array());
  int n = 4;
  CkPrintf ("DEBUG_TURBULENCE_OU compute() contribute %s [%Lg: %Lg %Lg %Lg %Lg]\n",
            block->name().c_str(),r_gv128[0],r_gv128[1],r_gv128[2],r_gv128[3],r_gv128[4]);
  fflush(stdout);
  
  block->contribute((n+1)*sizeof(long double), r_gv128, 
                    sum_long_double_n_type, callback);
  
}

//----------------------------------------------------------------------

void EnzoBlock::r_method_turbulence_ou_shift(CkReductionMsg *msg)
{
  CkPrintf ("DEBUG_TURBULENCE_OU shift() (block) %s\n",this->name().c_str());
  fflush(stdout);
  EnzoMethodTurbulenceOU * method = static_cast<EnzoMethodTurbulenceOU*> (this->method());
  method->compute_shift(this, msg);
}

//----------------------------------------------------------------------

void EnzoMethodTurbulenceOU::compute_shift
(EnzoBlock * enzo_block, CkReductionMsg * msg)
{
  CkPrintf ("DEBUG_TURBULENCE_OU shift() (method) %s\n",enzo_block->name().c_str());
  fflush(stdout);
  long double * data = (long double *) msg->getData();
  long double r_gv128[4];
  int id=0;
  int n = data[id++];
  ASSERT1 ("EnzoMethodTurbulenceOU::compute_shift()",
           "Expected length 4 but actual length %d",
           n,(n==4));
  for (int i=0; i<n; i++) r_gv128[i] = data[id++];
  delete msg;
  
  CkPrintf ("DEBUG_TURBULENCE_OU shift() reduced %s [%Lg %Lg %Lg %Lg]\n",
            enzo_block->name().c_str(),r_gv128[0],r_gv128[1],r_gv128[2],r_gv128[3]);

  long double r_av128[2+1] = {2.0, 0.0,0.0};

  if (enzo_block->is_leaf()) {
  
#ifdef CALL_FORTRAN    
    double * w;
    double * turbAcc;
    double dt;
    int iupdate_sol = update_solution_ ? 1 : 0;
    int iapply_cooling = apply_cooling_ ? 1 : 0;
    int iapply_forcing = apply_forcing_ ? 1 : 0;
    int iapply_injection_rate = apply_injection_rate_ ? 1 : 0;
    int cello_apply_injection_rate;
    double cello_injection_rate;
    double r_gv[4] = {r_gv128[1],r_gv128[2],r_gv128[3],r_gv128[4]};
    CkPrintf ("r_gv = %g %g %g %g\n",r_gv[0],r_gv[1],r_gv[2],r_gv[3]);
    double r_av[2] = {0.0,0.0};

    double * field_density     = cello::field(enzo_block,"density");
    double * field_momentum_x  = cello::field(enzo_block,"velocity_x");
    double * field_momentum_y  = cello::field(enzo_block,"velocity_y");
    double * field_momentum_z  = cello::field(enzo_block,"velocity_z");
    double * field_jacobian      = cello::field(enzo_block,"total_jacobian");
    
    // convert to conservative form
    int mx,my,mz;
    enzo_block->data()->field().dimensions(0,&mx,&my,&mz);
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
    double * field_work_1 = cello::field(enzo_block,"work_1");
    double * field_work_2 = cello::field(enzo_block,"work_2");
    double * field_work_3 = cello::field(enzo_block,"work_3");
    const int m = mx*my*mz;
    double * array_work = new double [3*m];
    std::copy_n (field_work_1,m,array_work+0*m);
    std::copy_n (field_work_2,m,array_work+1*m);
    std::copy_n (field_work_3,m,array_work+2*m);

  CkPrintf ("array_work shift in %g %g %g\n",array_work[0],array_work[m],array_work[2*m]);

  FORTRAN_NAME(turbforceshift)
      (&mx, &my, &mz,       // ( )
       field_density,            // ( )
       field_momentum_x,         // ( )
       field_momentum_y,         // ( )
       field_momentum_z,         // ( )
       field_jacobian,           // (*)
       array_work,               // ( )
       turbAcc,                  // ( )
       &iupdate_sol,             // ( )
       &iapply_injection_rate,   // ( )
       &olap_,                   // ( )
       &cello_injection_rate,    // ( )
       r_gv,                     // ( )
       r_av);                    // ( )

  CkPrintf ("array_work shift out %g %g %g\n",array_work[0],array_work[m],array_work[2*m]);
    for (int i=0; i<4+1; i++) r_av128[i+1] = r_av[i];

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
    delete [] array_work;
#endif  
  }
  
  CkCallback callback = CkCallback
    (CkIndex_EnzoBlock::r_method_turbulence_ou_update(nullptr), 
     enzo::block_array());
  n = 2;
  CkPrintf ("DEBUG_TURBULENCE_OU shift() contribute %s [%Lg: %Lg %Lg]\n",
            enzo_block->name().c_str(),r_av128[0],r_av128[1],r_av128[2]);
  enzo_block->contribute((n+1)*sizeof(long double), r_av128, 
                         sum_long_double_n_type, callback);
}

void EnzoBlock::r_method_turbulence_ou_update(CkReductionMsg *msg)
{
  CkPrintf ("DEBUG_TURBULENCE_OU update() (block) %s\n",this->name().c_str());
  fflush(stdout);
  EnzoMethodTurbulenceOU * method = static_cast<EnzoMethodTurbulenceOU*> (this->method());
  method->compute_update(this,msg);
}

void EnzoMethodTurbulenceOU::compute_update
(EnzoBlock * enzo_block, CkReductionMsg *msg)
{
  CkPrintf ("DEBUG_TURBULENCE_OU update() (method) %s\n",enzo_block->name().c_str());
  fflush(stdout);
  long double * data = (long double *) msg->getData();
  long double r_av128[2];
  int id=0;
  int n = data[id++];
  ASSERT1 ("EnzoMethodTurbulenceOU::compute_shift()",
           "Expected length 2 but actual length %d",
           n,(n==2));
  for (int i=0; i<n; i++) r_av128[i] = data[id++];
  delete msg;

  CkPrintf ("DEBUG_TURBULENCE_OU update() reduced %s [%Lg %Lg]\n",
            enzo_block->name().c_str(),r_av128[0],r_av128[1]);
  
#ifdef CALL_FORTRAN    
  int ni;
  int nj;
  int nk;
  double * w;
  double * turbAcc;
  double dt;
  double * res;
  int iupdate_sol = update_solution_ ? 1 : 0;
  int iapply_cooling = apply_cooling_ ? 1 : 0;
  int iapply_forcing = apply_forcing_ ? 1 : 0;
  int iapply_injection_rate = apply_injection_rate_ ? 1 : 0;
  int cello_apply_injection_rate;
  double cello_injection_rate;
  double r_av[4] = {0.0,0.0,0.0,0.0};

  double * field_density     = cello::field(enzo_block,"density");
  double * field_momentum_x  = cello::field(enzo_block,"velocity_x");
  double * field_momentum_y  = cello::field(enzo_block,"velocity_y");
  double * field_momentum_z  = cello::field(enzo_block,"velocity_z");
  double * field_energy      = cello::field(enzo_block,"total_energy");
  double * resid_density     = cello::field(enzo_block,"resid_density_r");
  double * resid_momentum_x  = cello::field(enzo_block,"resid_velocity_x");
  double * resid_momentum_y  = cello::field(enzo_block,"resid_velocity_y");
  double * resid_momentum_z  = cello::field(enzo_block,"resid_velocity_z");
  double * resid_energy      = cello::field(enzo_block,"total_energy");
  double * field_temperature = cello::field(enzo_block,"temperature");

  int mx,my,mz;
  enzo_block->data()->field().dimensions(0,&mx,&my,&mz);
  const int m = mx*my*mz;

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
  double * field_work_1 =  cello::field(enzo_block,"work_1");
  double * field_work_2 =  cello::field(enzo_block,"work_2");
  double * field_work_3 =  cello::field(enzo_block,"work_3");
  double * array_work = new double [3*m];
  std::copy_n (field_work_1,m,array_work+0*m);
  std::copy_n (field_work_2,m,array_work+1*m);
  std::copy_n (field_work_3,m,array_work+2*m);

  CkPrintf ("array_work update in %g %g %g\n",array_work[0],array_work[m],array_work[2*m]);
  FORTRAN_NAME(turbforceupdate)
    (&mx, &my, &mz,       // ( )
     field_density,            // ( )
     field_momentum_x,         // ( )
     field_momentum_y,         // ( )
     field_momentum_z,         // ( )
     field_energy,             // ( )
     resid_density,            // ( )
     resid_momentum_x,         // ( )
     resid_momentum_y,         // ( )
     resid_momentum_z,         // ( )
     resid_energy,             // ( )
     field_temperature,        // ( )
     array_work, &dt, res,     // ( )
     turbAcc,                  // ( )
     &iupdate_sol,             // ( )
     &iapply_injection_rate,   // ( )
     &cello_injection_rate,    // ( )
     &cooling_term_,           // ( )
     &iapply_cooling,          // ( )
     &gamma_,                  // ( )
     &hc_alpha_,               // ( )
     &hc_sigma_,               // ( )
     &totemp_,                 // ( )
     r_av );                   // ( )
  CkPrintf ("array_work update out %g %g %g\n",array_work[0],array_work[m],array_work[2*m]);
#endif

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
  delete [] array_work;

  CkPrintf ("DEBUG_TURBULENCE_OU update() exit %s\n",enzo_block->name().c_str());
  enzo_block->compute_done();
}

