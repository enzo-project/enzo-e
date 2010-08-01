#ifndef CELLO_HYDRO_H
#define CELLO_HYDRO_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defines.h"
#include "monitor.hpp"

//----------------------------------------------------------------------
// DEFINES
//----------------------------------------------------------------------

#define WARNING(MESSAGE) printf ("%s:%d WARNING: %s\n",__FILE__,__LINE__,MESSAGE);

#define ENZO_FAIL                           0 // Error handling
#define ENZO_SUCCESS                        1 // Error handling

#define FALSE                               0 // Needed for fortran
#define TRUE                                1 // Needed for fortran

#define ENZO_FLOAT                          double // Scalar
#define ENZO_FLOAT_UNDEFINED              -99999.0 // use NaN: CosmologyComputeExpansionFactor()
#define ISYM                              "d" // Scalar

#define MAX_DIMENSION                       3 // for array declarations and loops in SolveHydro
#define MAX_NUMBER_OF_BARYON_FIELDS         28 // for array declarations and loops in SolveHydro

#define SIGN(A)   ((A) >  0  ?  1  : -1 )     // upper-case inline function
#define NINT(A)   ((int) ((A) + 0.5*SIGN(A)) ) // rename to round(), upper-case inline function

// defined in cello.h
// #define MIN(A,B)  ((A) < (B) ? (A) : (B))      // upper-case inline function
// #define MAX(A,B)  ((A) > (B) ? (A) : (B))      // upper-case inline function

// #define FORTRAN_NAME(NAME) NAME
#define FORTRAN_NAME(NAME) NAME##_

//----------------------------------------------------------------------

const int Density         = 0;  // Field identifiers: use Field's instead
const int TotalEnergy     = 1;
const int InternalEnergy  = 2;
const int Velocity1       = 4;
const int Velocity2       = 5;
const int Velocity3       = 6;
const int ElectronDensity = 7;
const int HIDensity       = 8;
const int HIIDensity      = 9;
const int HeIDensity      = 10;
const int HeIIDensity     = 11;
const int HeIIIDensity    = 12;
const int HMDensity       = 13;
const int H2IDensity      = 14;
const int H2IIDensity     = 15;
const int DIDensity       = 16;
const int DIIDensity      = 17;
const int HDIDensity      = 18;
const int Metallicity     = 19;

//----------------------------------------------------------------------
// TYPEDEFS
//----------------------------------------------------------------------

typedef int            Eint32;     // c_message only
typedef long long      long_int;   // use long long
typedef long long int  Elong_int;  // use long long
typedef long long unsigned  global_index; // 

#endif /* CELLO_HYDRO_H */

