// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoPhysicsCosmology.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-05-23
/// @brief    Porting cosmology from Enzo

#include "enzo.hpp"

enzo_float EnzoPhysicsCosmology::time_from_redshift (enzo_float redshift) const
{
  enzo_float eta, time_hubble_0 = -1.0;
 
  /* Find Omega due to curvature. */
 
  enzo_float omega_curvature_now = 1 - omega_matter_now_ - omega_lambda_now_;
 
  /* 1) For a flat universe with omega_matter_now_ = 1, things are easy. */
 
  if (omega_matter_now_ == 1 && omega_lambda_now_ == 0)
    time_hubble_0 = 2.0/3.0/pow(1.0+redshift, 1.5);

  /* 2) For omega_matter_now_ < 1 and omega_lambda_now_ == 0 see
        Peebles 1993, eq. 13-3, 13-10. */
 
  if (omega_matter_now_ < 1 && omega_lambda_now_ == 0) {
    eta = acosh(1 + 2*(1-omega_matter_now_)/omega_matter_now_/(1+redshift));
    time_hubble_0 = omega_matter_now_/(2*pow(1.0-omega_matter_now_, 1.5))*
		  (sinh(eta) - eta);
  }
 
  /* 3) For omega_matter_now_ > 1 && omega_lambda_now_ == 0, use sin/cos. */
 
  if (omega_matter_now_ > 1 && omega_lambda_now_ == 0) {
    eta = acos(1.0 - 2.0*(1.0-omega_matter_now_)/
	       omega_matter_now_/(1.0+redshift));
    time_hubble_0 = omega_matter_now_/(2.0*pow(1.0-omega_matter_now_, 1.5))*
                  (eta - sin(eta));
  }
 
  /* 4) For flat universe, with non-zero omega_lambda_now_, see eq. 13-20. */
 
  if (fabs(omega_curvature_now) < 1.0e-3 && omega_lambda_now_ != 0)
    time_hubble_0 = 2.0/3.0/sqrt(1.0-omega_matter_now_)*
		    asinh(sqrt((1.0-omega_matter_now_)/omega_matter_now_)/
		           pow(1.0+redshift,1.5)             );
 
  /* Now convert from Time * H0 to code units (see also CosmologyGetUnits). */
 
  enzo_float time_units = 2.52e17/sqrt(omega_matter_now_)/
    hubble_constant_now_/pow(1.0 + initial_redshift_, 1.5);
 
  return time_hubble_0 / (hubble_constant_now_*3.24e-18) /
    time_units;
}

//======================================================================

void EnzoPhysicsCosmology::compute_expansion_factor
(enzo_float *a, enzo_float *dadt, enzo_float time) const
{
 
  /* Error check. */

  ASSERT ("EnzoPhysicsCosmology::compute_expansion_factor",
	  "Initial time in code units is 0",
	  (initial_time_in_code_units() != 0) );
 
  /* Find Omega due to curvature. */
 
  enzo_float omega_curvature_now_ = 1 - omega_matter_now_ - omega_lambda_now_;
 
  /* Convert the time from code units to Time * H0 (c.f. CosmologyGetUnits). */
 
  enzo_float time_units = 2.52e17/sqrt(omega_matter_now_)/hubble_constant_now_/
                    pow(1 + initial_redshift_,enzo_float(1.5));

  enzo_float time_hubble_0 = time * time_units * (hubble_constant_now_*3.24e-18);
 
  /* 1) For a flat universe with omega_matter_now_ = 1, it's easy. */
 
  if (fabs(omega_matter_now_-1) < OMEGA_TOLERANCE &&
      omega_lambda_now_ < OMEGA_TOLERANCE) {
    *a    = pow(time/initial_time_in_code_units(), enzo_float(2.0/3.0));
  }
 
  /* 2) For omega_matter_now_ < 1 and omega_lambda_now_ == 0 see
        Peebles 1993, eq. 13-3, 13-10.
	Actually, this is a little tricky since we must solve an equation
	of the form eta - sinh(eta) + x = 0..*/
 
  if (omega_matter_now_ < 1 && omega_lambda_now_ < OMEGA_TOLERANCE) {
 
    enzo_float eta, eta_old, x;
    int i;

    x = 2*time_hubble_0*pow(1.0 - omega_matter_now_, 1.5) / omega_matter_now_;
 
    /* Compute eta in a three step process, first from a third-order
       Taylor expansion of the formula above, then use that in a fifth-order
       approximation.  Then finally, iterate on the formula itself, solving for
       eta.  This works well because parts 1 & 2 are an excellent approximation
       when x is small and part 3 converges quickly when x is large. */
 
    eta = pow(6*x, enzo_float(1.0/3.0));                     // part 1
    eta = pow(120*x/(20+eta*eta), enzo_float(1.0/3.0));      // part 2
    for (i = 0; i < 40; i++) {                          // part 3
      eta_old = eta;
      eta = asinh(eta + x);
      if (fabs(eta-eta_old) < ETA_TOLERANCE) break;
    }
    if (i == 40) {
      fprintf(stderr, "Case 2 -- no convergence after %" ISYM " iterations.\n", i);
      return;
    }
 
    /* Now use eta to compute the expansion factor (eq. 13-10, part 2). */
 
    *a = omega_matter_now_/(2*(1 - omega_matter_now_))*(cosh(eta) - 1);
    *a *= (1 + initial_redshift_);    // to convert to code units, divide by [a]
  }
 
  /* 3) For omega_matter_now_ > 1 && omega_lambda_now_ == 0, use sin/cos.
        Easy, but skip it for now. */
 
  if (omega_matter_now_ > 1 && omega_lambda_now_ < OMEGA_TOLERANCE) {
  }
 
  /* 4) For flat universe, with non-zero omega_lambda_now_, see eq. 13-20. */
 
  if (fabs(omega_curvature_now_) < OMEGA_TOLERANCE &&
      omega_lambda_now_ > OMEGA_TOLERANCE) {
    *a = pow(enzo_float(omega_matter_now_/(1 - omega_matter_now_)), enzo_float(1.0/3.0)) *
         pow(enzo_float(sinh(1.5 * sqrt(1.0 - omega_matter_now_)*time_hubble_0)),
	     enzo_float(2.0/3.0));
    *a *= (1 + initial_redshift_);    // to convert to code units, divide by [a]
  }
 
  /* Compute the derivative of the expansion factor (Peebles93, eq. 13.3). */
 
  enzo_float temp = (*a)/(1 + initial_redshift_);
  *dadt = sqrt( 2.0/(3.0*omega_matter_now_*(*a)) *
	       (omega_matter_now_ + omega_curvature_now_*temp +
		omega_lambda_now_*temp*temp*temp));
 
  /* Someday, we'll implement the general case... */
 
}

//---------------------------------------------------------------------- 
 
void EnzoPhysicsCosmology::compute_expansion_timestep
(enzo_float *dt_expansion, enzo_float time) const
{
 
  /* Error check. */
 
  ASSERT ("EnzoPhysicsCosmology::compute_expansion_timestep",
	  "Initial time in code units is 0",
	  initial_time_in_code_units() != 0);
 
  /* Compute the expansion factors. */
 
  enzo_float a, dadt;
  compute_expansion_factor(&a, &dadt, time);
 
  /* Compute the maximum allwed timestep given the maximum allowed
     expansion factor. */
 
  *dt_expansion = max_expansion_rate_*a/dadt;
}

//----------------------------------------------------------------------

