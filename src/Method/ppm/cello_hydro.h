#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defines.h"
#include "error.hpp"

//----------------------------------------------------------------------
// DEFINES
//----------------------------------------------------------------------

#define WARNING(MESSAGE) printf ("%s:%d WARNING: %s\n",__FILE__,__LINE__,MESSAGE);

#define FAIL                                0 // Error handling
#define SUCCESS                             1 // Error handling

#define FALSE                               0 // Needed for fortran
#define TRUE                                1 // Needed for fortran

#define FLOAT                          double // Scalar
#define FLOAT_UNDEFINED              -99999.0 // use NaN: CosmologyComputeExpansionFactor()
#define ISYM                              "d" // Scalar

#define MAX_DIMENSION                       3 // for array declarations and loops in SolveHydro
#define MAX_NUMBER_OF_BARYON_FIELDS         5 // for array declarations and loops in SolveHydro

#define sign(A)   ((A) >  0  ?  1  : -1 )     // upper-case inline function
#define nint(A)   ((int) ((A) + 0.5*sign(A)) ) // rename to round(), upper-case inline function
#define min(A,B)  ((A) < (B) ? (A) : (B))      // upper-case inline function
#define max(A,B)  ((A) > (B) ? (A) : (B))      // upper-case inline function

// #define FORTRAN_NAME(NAME) NAME
#define FORTRAN_NAME(NAME) NAME##_


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

//----------------------------------------------------------------------
// EXTERNS
//----------------------------------------------------------------------

// Cosmology

extern int    ComovingCoordinates;
extern int    UseMinimumPressureSupport;
extern float  MinimumPressureSupportParameter;
extern float  ComovingBoxSize;
extern float  HubbleConstantNow;
extern float  OmegaLambdaNow;
extern float  OmegaMatterNow;
extern float  MaxExpansionRate;

// Chemistry

extern int    MultiSpecies;

// Gravity

extern int    GravityOn;
extern float *AccelerationField[MAX_DIMENSION];

// Physics

extern int    PressureFree;
extern float  Gamma;
extern float  GravitationalConstant;

// Problem-specific

extern int    ProblemType;

// Method PPM

extern int    PPMFlatteningParameter;
extern int    PPMDiffusionParameter;
extern int    PPMSteepeningParameter;

// Parallel

extern int    ProcessorNumber;

// Numerics

extern int    DualEnergyFormalism;
extern float  DualEnergyFormalismEta1;
extern float  DualEnergyFormalismEta2;
extern float  pressure_floor;
extern float  density_floor;
extern float  number_density_floor;
extern float  temperature_floor;

// Control

extern float  time_stop;   // stop after this simulation time
extern int    cycle_stop;  // stop after this number of cycles

extern float  CourantSafetyNumber;
extern FLOAT  InitialRedshift;
extern FLOAT  InitialTimeInCodeUnits;
extern FLOAT  Time;
extern FLOAT  OldTime;

// Domain

extern FLOAT  DomainLeftEdge [MAX_DIMENSION];
extern FLOAT  DomainRightEdge[MAX_DIMENSION];

// Grid

extern int    GridRank;
extern int    GridDimension[MAX_DIMENSION];
extern int    GridStartIndex[MAX_DIMENSION];
extern int    GridEndIndex[MAX_DIMENSION];
extern FLOAT  GridLeftEdge[MAX_DIMENSION];
extern FLOAT *CellWidth[MAX_DIMENSION];

// Fields

extern int field_density;
extern int field_total_energy;
extern int field_internal_energy;
extern int field_velocity_x;
extern int field_velocity_y;
extern int field_velocity_z;
extern int field_color;

extern int    NumberOfBaryonFields;
extern float *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS];
extern float *OldBaryonField[MAX_NUMBER_OF_BARYON_FIELDS];
extern int    FieldType[MAX_NUMBER_OF_BARYON_FIELDS];

// Boundary

enum bc_type 
  {
    bc_unknown    = 0, 
    bc_reflecting = 1, 
    bc_outflow    = 2, 
    bc_inflow     = 3, 
    bc_periodic   = 4
  };

extern int      BoundaryRank;
extern int      BoundaryDimension[MAX_DIMENSION];
extern int      BoundaryFieldType[MAX_NUMBER_OF_BARYON_FIELDS];
extern bc_type *BoundaryType[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2];
extern float   *BoundaryValue[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2];  

//----------------------------------------------------------------------
// FUNCTION PROTOTYPES
//----------------------------------------------------------------------

int ComputeGammaField(float *GammaField);

int ComputePressure(FLOAT time, float *pressure);

int ComputePressureDualEnergyFormalism(FLOAT time, float *pressure);

int ComputeTemperatureField(float *temperature);

int IdentifyPhysicalQuantities(int &DensNum, int &GENum,   int &Vel1Num, 
			       int &Vel2Num, int &Vel3Num, int &TENum);

int IdentifySpeciesFields(int &DeNum, int &HINum, int &HIINum, 
			  int &HeINum, int &HeIINum, int &HeIIINum,
			  int &HMNum, int &H2INum, int &H2IINum,
			  int &DINum, int &DIINum, int &HDINum);

int SetMinimumSupport(float &MinimumSupportEnergyCoefficient);

int SolveHydroEquations(int cycle, float dt);
int SolveMHDEquations(int cycle, float dt);

float ComputeTimeStep();
int SetExternalBoundaryValues();
int SetExternalBoundary(int FieldRank, 
			     int GridDims[],
			     int GridOffset[],
			     int StartIndex[], 
			     int EndIndex[],
			     float *Field, 
			     int FieldType);



void initialize_hydro ();

// PPM

void initialize_implosion (int size, int cycles);
void initialize_implosion3 (int size, int cycles);
void initialize_image ();
void initialize_color ();

// PPML

void initialize_ppml (int size, int cycles);

float sum_field (int field);
float print_field (int field);
void data_dump (const char * file_root, int cycle);
