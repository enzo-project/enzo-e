// $Id$
// See LICENSE_ENZO file for license and copyright information

#include "cello_hydro.h"

// Cosmology

int   ComovingCoordinates;
int   UseMinimumPressureSupport;
float MinimumPressureSupportParameter;
float ComovingBoxSize;
float HubbleConstantNow;
float OmegaMatterNow;
float OmegaLambdaNow;
float MaxExpansionRate;

// Chemistry

int MultiSpecies;

// Gravity

int GravityOn;
float *AccelerationField[MAX_DIMENSION]; // cell cntr acceleration at n+1/2

// Physics

int PressureFree;
float Gamma;
float GravitationalConstant;

// Problem-specific

int ProblemType;

// Method PPM

int    PPMFlatteningParameter;
int    PPMDiffusionParameter;
int    PPMSteepeningParameter;

// Parallel

int    ProcessorNumber;

// Numerics

int   DualEnergyFormalism;
float DualEnergyFormalismEta1;
float DualEnergyFormalismEta2;
float pressure_floor;
float density_floor;
float number_density_floor;
float temperature_floor;

// Control

float  time_stop;   // stop after this simulation time
int    cycle_stop;  // stop after this number of cycles

float  CourantSafetyNumber;
ENZO_FLOAT  InitialRedshift;
ENZO_FLOAT  InitialTimeInCodeUnits;
ENZO_FLOAT  Time;
ENZO_FLOAT  OldTime;

// Domain

ENZO_FLOAT DomainLeftEdge [MAX_DIMENSION];
ENZO_FLOAT DomainRightEdge[MAX_DIMENSION];

// Grid

int field_density;
int field_total_energy;
int field_internal_energy;
int field_velocity_x;
int field_velocity_y;
int field_velocity_z;
int field_color;
int field_magnetic_x;
int field_magnetic_y;
int field_magnetic_z;

int field_density_xp;
int field_velocity_x_xp;
int field_velocity_y_xp;
int field_velocity_z_xp;
int field_magnetic_x_xp;
int field_magnetic_y_xp;
int field_magnetic_z_xp;

int field_density_yp;
int field_velocity_x_yp;
int field_velocity_y_yp;
int field_velocity_z_yp;
int field_magnetic_x_yp;
int field_magnetic_y_yp;
int field_magnetic_z_yp;

int field_density_zp;
int field_velocity_x_zp;
int field_velocity_y_zp;
int field_velocity_z_zp;
int field_magnetic_x_zp;
int field_magnetic_y_zp;
int field_magnetic_z_zp;


int    GridRank;
int    GridDimension[MAX_DIMENSION];   // total dimensions of all grids
int    GridStartIndex[MAX_DIMENSION];  // starting index of the active region
int    GridEndIndex[MAX_DIMENSION];    // stoping index of the active region
ENZO_FLOAT  GridLeftEdge[MAX_DIMENSION];    // starting pos (active problem space)
ENZO_FLOAT *CellWidth[MAX_DIMENSION];

// Fields

int    NumberOfBaryonFields;                        // active baryon fields
float *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS];    // pointers to arrays
float *OldBaryonField[MAX_NUMBER_OF_BARYON_FIELDS]; // pointers to old arrays
int    FieldType[MAX_NUMBER_OF_BARYON_FIELDS];

// Boundary

int      BoundaryRank;
int      BoundaryDimension[MAX_DIMENSION];
int      BoundaryFieldType[MAX_NUMBER_OF_BARYON_FIELDS];
bc_type *BoundaryType[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2];
float   *BoundaryValue[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2];  




