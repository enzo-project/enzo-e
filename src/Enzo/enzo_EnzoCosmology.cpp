// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoCosmology.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-05-23
/// @brief    Porting cosmology from Enzo

#include "enzo.hpp"

double EnzoCosmology::time_from_redshift_ (double redshift) const
{
 
  /* Error check [INCOMPLETE]. */
 
  double eta, time_hubble_0 = -1.0;
 
  /* Find Omega due to curvature. */
 
  double omega_curvature_now = 1 - omega_matter_now_ - omega_lambda_now_;
 
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
 
  double time_units = 2.52e17/sqrt(omega_matter_now_)/
    hubble_constant_now_/pow(1.0 + initial_redshift_, 1.5);
 
  return time_hubble_0 / (hubble_constant_now_*3.24e-18) /
    time_units;
}

