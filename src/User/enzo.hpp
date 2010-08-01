// $Id: Enzo.hpp 1394 2010-04-22 20:52:54Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef ENZO_HPP
#define ENZO_HPP

/// @file     Enzo.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Brief description of file Enzo.hpp

#include "global.hpp"
#include "cello_hydro.h"

class Enzo {

  /// @class    Enzo
  /// @ingroup  User
  /// @brief    Object for storing cross-cutting objects, including monitor, error, parameters

public: // interface

  /// Constructor
  Enzo() throw()
  {
  }

  /// Destructor
  ~Enzo() throw()
  {
  }


  float ComputeTimeStep();
  float sum_field (int field);
  int ComputeGammaField(float *GammaField);
  int ComputePressureDualEnergyFormalism(ENZO_FLOAT time, float *pressure);
  int ComputePressure(ENZO_FLOAT time, float *pressure);
  int ComputeTemperatureField(float *temperature);
  int CosmologyComputeExpansionFactor(ENZO_FLOAT time, ENZO_FLOAT *a, ENZO_FLOAT *dadt);
  int CosmologyComputeExpansionTimestep(ENZO_FLOAT time, float *dtExpansion);
  int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,			      float *TemperatureUnits, float *TimeUnits,			      float *VelocityUnits, ENZO_FLOAT Time);
  int FindField(int field, int farray[], int numfields);
  int IdentifyPhysicalQuantities(int &DensNum, int &GENum, int &Vel1Num,				 int &Vel2Num, int &Vel3Num, int &TENum);
  int IdentifySpeciesFields(int &DeNum, int &HINum, int &HIINum,			    int &HeINum, int &HeIINum, int &HeIIINum,			    int &HMNum, int &H2INum, int &H2IINum,			    int &DINum, int &DIINum, int &HDINum);
    int SetExternalBoundaryValues();
  int SetMinimumSupport(float &MinimumSupportEnergyCoefficient);
  int SolveHydroEquations ( int CycleNumber,  float dt);
  void print_field (int field);

private: // prohibit copy constructor

  /// Copy constructor
  Enzo(const Enzo & enzo) throw()
  {
  }

private: // prohibit assignment

  /// Assignment operator
  Enzo & operator= (const Enzo & enzo) throw();

public: // public attributes (!!)

  
  int   ComovingCoordinates;
//   int   UseMinimumPressureSupport;
//   float MinimumPressureSupportParameter;
//   float ComovingBoxSize;
//   float HubbleConstantNow;
//   float OmegaMatterNow;
//   float OmegaLambdaNow;
//   float MaxExpansionRate;

//   // Chemistry

//   int MultiSpecies;

//   // Gravity

//   int GravityOn;
//   float *AccelerationField[MAX_DIMENSION]; // cell cntr acceleration at n+1/2

//   // Physics

//   int PressureFree;
//   float Gamma;
//   float GravitationalConstant;

//   // Problem-specific

//   int ProblemType;

//   // Method PPM

//   int    PPMFlatteningParameter;
//   int    PPMDiffusionParameter;
//   int    PPMSteepeningParameter;

//   // Parallel

//   int    ProcessorNumber;

//   // Numerics

//   int   DualEnergyFormalism;
//   float DualEnergyFormalismEta1;
//   float DualEnergyFormalismEta2;
//   float pressure_floor;
//   float density_floor;
//   float number_density_floor;
//   float temperature_floor;

//   float  CourantSafetyNumber;
//   ENZO_FLOAT  InitialRedshift;
//   ENZO_FLOAT  InitialTimeInCodeUnits;
//   ENZO_FLOAT  Time;
//   ENZO_FLOAT  OldTime;

//   // Domain

//   ENZO_FLOAT DomainLeftEdge [MAX_DIMENSION];
//   ENZO_FLOAT DomainRightEdge[MAX_DIMENSION];

//   // Grid

//   int field_density;
//   int field_total_energy;
//   int field_internal_energy;
//   int field_velocity_x;
//   int field_velocity_y;
//   int field_velocity_z;
//   int field_color;
//   int field_magnetic_x;
//   int field_magnetic_y;
//   int field_magnetic_z;

//   int field_density_xp;
//   int field_velocity_x_xp;
//   int field_velocity_y_xp;
//   int field_velocity_z_xp;
//   int field_magnetic_x_xp;
//   int field_magnetic_y_xp;
//   int field_magnetic_z_xp;

//   int field_density_yp;
//   int field_velocity_x_yp;
//   int field_velocity_y_yp;
//   int field_velocity_z_yp;
//   int field_magnetic_x_yp;
//   int field_magnetic_y_yp;
//   int field_magnetic_z_yp;

//   int field_density_zp;
//   int field_velocity_x_zp;
//   int field_velocity_y_zp;
//   int field_velocity_z_zp;
//   int field_magnetic_x_zp;
//   int field_magnetic_y_zp;
//   int field_magnetic_z_zp;


//   int    GridRank;
//   int    GridDimension[MAX_DIMENSION];   // total dimensions of all grids
//   int    GridStartIndex[MAX_DIMENSION];  // starting index of the active region
//   int    GridEndIndex[MAX_DIMENSION];    // stoping index of the active region
//   ENZO_FLOAT  GridLeftEdge[MAX_DIMENSION];    // starting pos (active problem space)
//   ENZO_FLOAT *CellWidth[MAX_DIMENSION];
//   int    ghost_depth[MAX_DIMENSION];

//   // Fields

//   int    NumberOfBaryonFields;                        // active baryon fields
//   float *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS];    // pointers to arrays
//   float *OldBaryonField[MAX_NUMBER_OF_BARYON_FIELDS]; // pointers to old arrays
//   int    FieldType[MAX_NUMBER_OF_BARYON_FIELDS];

//   // Boundary

//   int      BoundaryRank;
//   int      BoundaryDimension[MAX_DIMENSION];
//   int      BoundaryFieldType[MAX_NUMBER_OF_BARYON_FIELDS];
//   bc_type *BoundaryType[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2];
//   float   *BoundaryValue[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2];  

//   // problem

//   int   CycleNumber;
//   float dt;

};

#endif /* ENZO_HPP */

