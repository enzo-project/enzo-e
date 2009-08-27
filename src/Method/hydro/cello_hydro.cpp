#include "cello_hydro.h"

// Global parameters

int ComovingCoordinates;
int DualEnergyFormalism;
int MultiSpecies;
int GravityOn;
int PressureFree;
int ProblemType;
int UseMinimumPressureSupport;
float ComovingBoxSize;
float DualEnergyFormalismEta1;
float DualEnergyFormalismEta2;
float Gamma;
float GravitationalConstant;
float HubbleConstantNow;
float MinimumPressureSupportParameter;
float OmegaMatterNow;
float OmegaLambdaNow;
float pressure_floor;
float density_floor;
float number_density_floor;
float temperature_floor;
float CourantSafetyNumber;
float MaxExpansionRate;
FLOAT DomainLeftEdge[MAX_DIMENSION];
FLOAT InitialRedshift;
FLOAT InitialTimeInCodeUnits;

// Grid variables

int GridRank;
int GridDimension[MAX_DIMENSION];   // total dimensions of all grids
int GridStartIndex[MAX_DIMENSION];  // starting index of the active region
int GridEndIndex[MAX_DIMENSION];    // stoping index of the active region
FLOAT GridLeftEdge[MAX_DIMENSION];    // starting pos (active problem space)

// class grid
// {

//   float  dtFixed;                        // current (fixed) timestep
//   FLOAT  Time;                           // current problem time
//   FLOAT  OldTime;                        // time corresponding to OldBaryonField

//   int    NumberOfBaryonFields;                        // active baryon fields
//   float *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS];    // pointers to arrays
//   float *OldBaryonField[MAX_NUMBER_OF_BARYON_FIELDS]; // pointers to old arrays
//   int    FieldType[MAX_NUMBER_OF_BARYON_FIELDS];
//   FLOAT *CellWidth[MAX_DIMENSION];

//   int    PPMFlatteningParameter;
//   int    PPMDiffusionParameter;
//   int    PPMSteepeningParameter;
//   float *AccelerationField[MAX_DIMENSION]; // cell cntr acceleration at n+1/2
//   int    ProcessorNumber;
// };



