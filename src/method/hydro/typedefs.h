#ifndef __typedefs_h_
#define __typedefs_h_
/***********************************************************************
/
/  MISCELANEOUS TYPEDEFS AND ENUMERATIONS
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include "CoolData.h"
#include "RateData.h"
#include "RadiationFieldData.h"

/* These are the different types of baryon fields. */

#ifdef SMALL_INTS
typedef int field_type;
typedef int boundary_type;
typedef int gravity_boundary_type;
typedef int interpolation_type;
typedef int hydro_method;
#endif

#ifdef LARGE_INTS
typedef long_int field_type;
typedef long_int boundary_type;
typedef long_int gravity_boundary_type;
typedef long_int interpolation_type;
typedef long_int hydro_method;
#endif

const field_type 
  Density         = 0,
  TotalEnergy     = 1,
  InternalEnergy  = 2,
  Pressure        = 3,
  Velocity1       = 4,
  Velocity2       = 5,
  Velocity3       = 6,
  ElectronDensity = 7,
  HIDensity       = 8,
  HIIDensity      = 9,
  HeIDensity      = 10,
  HeIIDensity     = 11,
  HeIIIDensity    = 12,
  HMDensity       = 13,
  H2IDensity      = 14,
  H2IIDensity     = 15,
  DIDensity       = 16,
  DIIDensity      = 17,
  HDIDensity      = 18,
  Metallicity     = 19,
  ExtraType0      = 20,
  ExtraType1      = 21,
  GravPotential   = 22,
  Acceleration0   = 23,
  Acceleration1   = 24,
  Acceleration2   = 25,
  RadiationFreq0  = 26,
  RadiationFreq1  = 27,
  RadiationFreq2  = 28,
  RadiationFreq3  = 29,
  RadiationFreq4  = 30,
  RadiationFreq5  = 31,
  FieldUndefined  = 32,

/* these pseudo-fields are used to access grid data 
   the "g" prefix is to avoid namespace conflict */

  gParticlePosition     = 100,
  gParticleVelocity     = 101,
  gParticleMass         = 102,
  gParticleAcceleration = 103,
  gParticleNumber       = 104,
  gParticleType         = 105,
  gParticleAttribute    = 106,
  gPotentialField       = 107,
  gAccelerationField    = 108,
  gGravitatingMassField = 109,
  gFlaggingField        = 110,
  gVelocity             = 111;

  
/*
enum field_type {Density, TotalEnergy, InternalEnergy, Pressure,
		 Velocity1, Velocity2, Velocity3, 
		 ElectronDensity, HIDensity, HIIDensity,  HeIDensity, 
		 HeIIDensity, HeIIIDensity, HMDensity, H2IDensity, 
		 H2IIDensity, DIDensity, DIIDensity, HDIDensity,
                 Metallicity, ExtraType0, ExtraType1, GravPotential,
		 Acceleration0, Acceleration1,Acceleration2,
		 RadiationFreq0, RadiationFreq1, RadiationFreq2, 
		 RadiationFreq3, RadiationFreq4, RadiationFreq5,
		 FieldUndefined};
*/

#define FieldTypeIsDensity(A) (((A) >= TotalEnergy && (A) <= Velocity3) ? FALSE : TRUE)

/* These are the different types of fluid boundary conditions. */

const boundary_type
  reflecting        = 0,
  outflow           = 1,
  inflow            = 2,
  periodic          = 3,
  BoundaryUndefined = 4;

// enum boundary_type {reflecting, outflow, inflow, periodic, BoundaryUndefined};

/* These are the different types of gravity boundary conditions. */

const gravity_boundary_type
  TopGridPeriodic  = 0,
  TopGridIsolated  = 1,
  SubGridIsolated  = 2,
  GravityUndefined = 3;

// enum gravity_boundary_type {TopGridPeriodic, TopGridIsolated, 
// 				    SubGridIsolated, GravityUndefined};

/* Interpolation types. */

const interpolation_type
  ThirdOrderA            = 0,
  SecondOrderA           = 1,
  SecondOrderB           = 2,
  SecondOrderC           = 3,
  FirstOrderA            = 4,
  InterpolationUndefined = 5;


// enum interpolation_type {ThirdOrderA, SecondOrderA, SecondOrderB, SecondOrderC,
// 			 FirstOrderA, InterpolationUndefined};

/* Hydrodynamics methods. */

const hydro_method
  PPM_DirectEuler      = 0,
  PPM_LagrangeRemap    = 1,
  Zeus_Hydro           = 2,
  HydroMethodUndefined = 3;

// enum hydro_method {PPM_DirectEuler, PPM_LagrangeRemap, Zeus_Hydro};

/* Define a float/int union. */

union float_int {
  long_int ival;
  float fval;
  FLOAT FVAL;
};

#endif
