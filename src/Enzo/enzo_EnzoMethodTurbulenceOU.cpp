// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodTurbulenceOU.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @author   Alexei Kritsuk (akritsuk@ucsd.edu)
/// @date     2021-09-01
/// @brief    Implements the EnzoMethodTurbulenceOU class

#include "cello.hpp"
#include "enzo.hpp"
#include "charm_enzo.hpp"

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
  FORTRAN_NAME(cello_init_turbulence_ou)
    (&is_root,
     &rank,
     domain_size,
     &gamma,
     &iapply_injection_rate,
     &cooling_term,
     &injection_rate,
     &kfi,
     &kfa,
     &mach,
     &iread_sol,
     &sol_weight);
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
    int nc;
    int ni;
    int nj;
    int nk;
    double * w;
    double * grid;
    double * jac;
    double * temperature;
    double * turbAcc;
    double * wk;
    double dt;
    double * res;
    int iupdate_sol = update_solution_ ? 1 : 0;
    int iapply_cooling = apply_cooling_ ? 1 : 0;
    int iapply_forcing = apply_forcing_ ? 1 : 0;
    int iapply_injection_rate = apply_injection_rate_ ? 1 : 0;
    int cello_apply_injection_rate;
    double cello_gamma;
    double cello_hc_alpha;
    double cello_hc_sigma;
    double cello_injection_rate;
    double cello_totemp;
    double r_gv[4];

#ifdef CALL_FORTRAN    
    FORTRAN_NAME(turbforceou)
      (& nc, & ni, & nj, & nk,
       w,  grid,  jac,  temperature,
       wk,  & dt,  res,
       turbAcc, & iupdate_sol,
       & iapply_cooling,
       & iapply_forcing,
       & iapply_injection_rate,
       & cooling_term_,
       &gamma_,
       &injection_rate_,
       &olap_,
       r_gv
       );
#endif    
    for (int i=0; i<4+1; i++) r_gv128[i+1] = r_gv[i];

    r_gv128[1] = 1.0;
    r_gv128[2] = 2.0;
    r_gv128[3] = 3.0;
    r_gv128[4] = 4.0;
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

  long double r_av128[2+1] = {2.0,0.0,0.0};

  if (enzo_block->is_leaf()) {
    double r_av[2] = {0.0,0.0};;
    // call turbulence_shift

    for (int i=0; i<2; i++) r_av128[i+1] = r_av[i];
    r_av128[1] = 1.0;
    r_av128[2] = 2.0;
    
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
  
  CkPrintf ("DEBUG_TURBULENCE_OU update() exit %s\n",enzo_block->name().c_str());
  enzo_block->compute_done();
}

