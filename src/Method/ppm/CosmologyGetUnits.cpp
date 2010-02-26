/** 
 *********************************************************************
 *
 * @file      CosmologyGetUnits.cpp
 * @brief     Compute and return the cosmology units
 * @author    Greg Bryan
 * @date      April, 1995
 * @ingroup   Enzo
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * DESCRIPTION:
 *
 *     Returns the cosmology units:
 *     
 *              time:        utim = 1 / sqrt(4 * \pi * G * \rho_0 * (1+zri)^3)
 *              density:     urho = \rho_0 * (1+z)^3
 *              length:      uxyz = (1 Mpc) * box / h  (1+z)
 *              velocity:    uvel = uaye * uxyz / utim  (since u = a * dx/dt)
 *         (*)  temperature: utem = m_H * \mu / k * uvel**2
 *              a(t):        uaye = 1 / (1 + zri)
 *
 *                where:
 *                  box     - size of simulation box in Mpc/h
 *                  zri     - initial redshift (start of simulation)
 *                  \rho_0  = 3*\Omega_0*H_0^2/(8*\pi*G)
 *                  Omega_0 - the fraction of non-relativistic matter at z=0
 *     
 *                Note that two definitions are dependent on redshift (urho
 *                  and uxyz) so make sure to call this routine immediately
 *                  before writing.
 *
 *                * - the utem given below assumes that \mu = 1, so you must
 *                    multiply the resulting temperature field by \mu.
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */
 
#include "cello_hydro.h"

/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time)
{
 
  /* From the time, compute the current redshift. */
 
  FLOAT a, dadt;
  if (CosmologyComputeExpansionFactor(Time, &a, &dadt) == FAIL) {
    fprintf(stderr, "Error in ComputeExpansionFactor.\n");
    return FAIL;
  }
 
  /* Compute the current redshift (remember a(init) = 1). */
 
  FLOAT CurrentRedshift = (1 + InitialRedshift)/a - 1;
 
  /* Determine the units. */
 
  *DensityUnits     = 1.88e-29*OmegaMatterNow*pow(HubbleConstantNow,2)*
                      pow(1 + CurrentRedshift,3);
 
  *LengthUnits      = 3.086e24*ComovingBoxSize/HubbleConstantNow/
                      (1 + CurrentRedshift);
 
  *TemperatureUnits = 1.88e6*pow(ComovingBoxSize,2)*OmegaMatterNow*
                      (1 + InitialRedshift);
 
  *TimeUnits        = 2.52e17/sqrt(OmegaMatterNow)/HubbleConstantNow/
                      pow(1 + InitialRedshift,FLOAT(1.5));
 
  *VelocityUnits    = 1.225e7*ComovingBoxSize*sqrt(OmegaMatterNow)*
                      sqrt(1 + InitialRedshift);
 
  return SUCCESS;
}
