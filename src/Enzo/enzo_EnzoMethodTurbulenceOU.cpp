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
     &hc_alpha, 
     &hc_sigma,
     &injection_rate,
     &kfi,
     &kfa,
     &mach,
     &iread_sol,
     &sol_weight,
     &totemp);
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
    double * dt;
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
    INCOMPLETE ("MethodTurbulenceOU::compute()");
    FORTRAN_NAME(turbforceou)
      (& nc, & ni, & nj, & nk,
       w,  grid,  jac,  temperature,
       wk,  dt,  res,
       turbAcc, & iupdate_sol,
       & iapply_cooling,
       & iapply_forcing,
       & iapply_injection_rate,
       & cooling_term_,
       &gamma_,
       &hc_alpha_,
       &hc_sigma_,
       &injection_rate_,
       &olap_,
       &totemp_
       );
      }
  block->compute_done();
  
}

//----------------------------------------------------------------------
