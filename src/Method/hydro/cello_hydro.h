#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fortran.def"

//----------------------------------------------------------------------
// DEFINES
//----------------------------------------------------------------------

#define WARNING(MESSAGE) printf ("%s:%d WARNING: %s\n",__FILE__,__LINE__,MESSAGE);

#define FAIL                                0 // Error handling
#define FALSE                               0 // obsolete: use false
#define FLOAT                          double // Scalar
#define FLOAT_UNDEFINED              -99999.0 // use NaN: CosmologyComputeExpansionFactor()
#define ISYM                              "d" // Scalar
#define MAX_DIMENSION                       3 // for array declarations and loops in SolveHydro
#define MAX_NUMBER_OF_BARYON_FIELDS        20 // for array declarations and loops in SolveHydro
#define SUCCESS                             1 // Error handling
#define TRUE                                1 // obsolete: use true

#define sign(A)   ((A) >  0  ?  1  : -1 )     // upper-case inline function
#define nint(A)   ((int) ((A) + 0.5*sign(A)) ) // rename to round(), upper-case inline function
#define min(A,B)  ((A) < (B) ? (A) : (B))      // upper-case inline function
#define max(A,B)  ((A) > (B) ? (A) : (B))      // upper-case inline function

// #define FORTRAN_NAME(NAME) NAME
#define FORTRAN_NAME(NAME) NAME##_

const int PPM_DirectEuler = 0;  // Select PPM method: not necessary
const int Zeus_Hydro      = 2;  // Select Zeus_Hydro: not necessary

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

extern int   ComovingCoordinates;
extern int   DualEnergyFormalism;
extern int   MultiSpecies;
extern int   GravityOn;
extern int   PressureFree;
extern int   ProblemType;
extern int   UseMinimumPressureSupport;
extern float ComovingBoxSize;
extern float DualEnergyFormalismEta1;
extern float DualEnergyFormalismEta2;
extern float Gamma;
extern float GravitationalConstant;
extern float HubbleConstantNow;
extern float MinimumPressureSupportParameter;
extern float OmegaLambdaNow;
extern float OmegaMatterNow;
extern float pressure_floor;
extern float density_floor;
extern float number_density_floor;
extern float temperature_floor;
extern float CourantSafetyNumber;
extern float MaxExpansionRate;
extern FLOAT DomainLeftEdge[MAX_DIMENSION];
extern FLOAT InitialRedshift;
extern FLOAT InitialTimeInCodeUnits;
//----------------------------------------------------------------------
// CLASSES
//----------------------------------------------------------------------

extern int    GridRank;
extern int    GridDimension[MAX_DIMENSION];
extern int    GridStartIndex[MAX_DIMENSION];
extern int    GridEndIndex[MAX_DIMENSION];
extern FLOAT  GridLeftEdge[MAX_DIMENSION];

struct fluxes;
class grid
{
  float  dtFixed;                        // current (fixed) timestep
  FLOAT  Time;                           // current problem time
  FLOAT  OldTime;                        // time corresponding to OldBaryonField

  int    NumberOfBaryonFields;                        // active baryon fields
  float *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS];    // pointers to arrays
  float *OldBaryonField[MAX_NUMBER_OF_BARYON_FIELDS]; // pointers to old arrays
  int    FieldType[MAX_NUMBER_OF_BARYON_FIELDS];
  FLOAT *CellWidth[MAX_DIMENSION];

  int    PPMFlatteningParameter;
  int    PPMDiffusionParameter;
  int    PPMSteepeningParameter;
  float *AccelerationField[MAX_DIMENSION]; // cell cntr acceleration at n+1/2
  int    ProcessorNumber;

 public:

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
  int SolveHydroEquations(int CycleNumber, int NumberOfSubgrids, 
			  fluxes *SubgridFluxes[], int level);
  float ComputeTimeStep();

};

//----------------------------------------------------------------------
// STRUCTS
//----------------------------------------------------------------------

struct fluxes
{
  // first index selects axis orthogonal to face
  // second index selects coordinate of corner point
  // equivalently:
  //    Point lower [axis][face] ("Start")
  //    Point upper [axis][face] ("End")
  // where face is 0 or 1 for Left or Right

  global_index LeftFluxStartGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
  global_index LeftFluxEndGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
  global_index RightFluxStartGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
  global_index RightFluxEndGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];

  // Array of fluxes for each baryon field and each face
  // equivalently:
  //   float * fluxes [field][axis][face]
  float       *LeftFluxes[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION];
  float       *RightFluxes[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION];
};


//----------------------------------------------------------------------
// FUNCTIONS
//----------------------------------------------------------------------

void initialize_hydro ();

