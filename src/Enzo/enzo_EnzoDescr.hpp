// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef ENZO_ENZO_DESCR_HPP
#define ENZO_ENZO_DESCR_HPP

/// @file     enzo_EnzoDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @todo     Comment Enzo functions / attributes
/// @todo     Introduce EnzoBlock for block-dependent attributes
/// @brief    [\ref Enzo] Interface for the Enzo descriptor class

enum bc_enum 
{ // explicitly enumerated to match what Enzo expects
  bc_unknown    = 0, 
  bc_reflecting = 1, 
  bc_outflow    = 2, 
  bc_inflow     = 3, 
  bc_periodic   = 4 
};

struct fluxes
{
  long_int LeftFluxStartGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
  long_int LeftFluxEndGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
  long_int RightFluxStartGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
  long_int RightFluxEndGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
  enzo_float *LeftFluxes[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION];
  enzo_float *RightFluxes[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION];
};

class EnzoDescr {

  /// @class    EnzoDescr
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Object for storing cross-cutting objects,
  /// including monitor, error, parameters

public: // interface

  /// Constructor
  EnzoDescr() throw();

  /// Destructor
  ~EnzoDescr() throw();

  /// Initialize EnzoDescr values from parameter file
  void initialize(Parameters * parameters) throw();

  /// Write attributes, e.g. to stdout for debugging
  void write(FILE *fp) throw ();

  //----------------------------------------------------------------------
  // Original Enzo functions
  //----------------------------------------------------------------------

  enzo_float ComputeTimeStep();
  enzo_float sum_field (int field);
  int ComputeGammaField(enzo_float *GammaField);
  int ComputePressureDualEnergyFormalism(enzo_float time, enzo_float *pressure);
  int ComputePressure(enzo_float time, enzo_float *pressure);
  int ComputeTemperatureField(enzo_float *temperature);
  int CosmologyComputeExpansionFactor(enzo_float time, enzo_float *a, enzo_float *dadt);
  int CosmologyComputeExpansionTimestep(enzo_float time, enzo_float *dtExpansion);
  int CosmologyGetUnits(enzo_float *DensityUnits, enzo_float *LengthUnits, enzo_float *TemperatureUnits, enzo_float *TimeUnits, enzo_float *VelocityUnits, enzo_float Time);
  int FindField(int field, int farray[], int numfields);
  int IdentifyPhysicalQuantities(int &DensNum, int &GENum, int &Vel1Num, int &Vel2Num, int &Vel3Num, int &TENum);
  int IdentifySpeciesFields(int &DeNum, int &HINum, int &HIINum, int &HeINum, int &HeIINum, int &HeIIINum, int &HMNum, int &H2INum, int &H2IINum, int &DINum, int &DIINum, int &HDINum);
  int SetExternalBoundaryValues();
  int SetMinimumSupport(enzo_float &MinimumSupportEnergyCoefficient);
  int SolveHydroEquations (DataBlock * data_block, int CycleNumber, enzo_float dt);
  void print_field (int field);
  int SetExternalBoundary(int FieldRank, int GridDims[], int GridOffset[], int StartIndex[], int EndIndex[], enzo_float *Field, int FieldType);
  void image_dump(const char * file_root, int cycle, double lower, double upper, Monitor * monitor);

  void initialize_hydro ();
  void initialize_image ();
  void initialize_implosion3 (int size_param);
  void initialize_implosion (int size_param);
  void initialize_ppml_implosion3 (int size_param);

  int SolveMHDEquations(DataBlock * data_block, int cycle, enzo_float dt);
  void initialize_ppml (int size_param);

private: // prohibit copy constructor

  /// Copy constructor
  EnzoDescr(const EnzoDescr & enzo) throw()
  {
  }

private: // prohibit assignment

  /// Assignment operator
  EnzoDescr & operator= (const EnzoDescr & enzo) throw();

public: // PUBLIC ATTRIBUTES (!!)

 
  int ComovingCoordinates;
  int UseMinimumPressureSupport;
  enzo_float MinimumPressureSupportParameter;
  enzo_float ComovingBoxSize;
  enzo_float HubbleConstantNow;
  enzo_float OmegaMatterNow;
  enzo_float OmegaLambdaNow;
  enzo_float MaxExpansionRate;

  // Chemistry

  int MultiSpecies;

  // Gravity

  int GravityOn;
  enzo_float *AccelerationField[MAX_DIMENSION]; // cell cntr acceleration at n+1/2

  // Physics

  int PressureFree;
  enzo_float Gamma;
  enzo_float GravitationalConstant;

  // Problem-specific

  int ProblemType;

  // Method PPM

  int PPMFlatteningParameter;
  int PPMDiffusionParameter;
  int PPMSteepeningParameter;

  // Parallel

  int ProcessorNumber;

  // Numerics

  int DualEnergyFormalism;
  enzo_float DualEnergyFormalismEta1;
  enzo_float DualEnergyFormalismEta2;

  enzo_float pressure_floor;
  enzo_float density_floor;
  enzo_float number_density_floor;
  enzo_float temperature_floor;

  enzo_float CourantSafetyNumber;
  enzo_float InitialRedshift;
  enzo_float InitialTimeInCodeUnits;
  enzo_float Time;
  enzo_float OldTime;

  // Domain

  enzo_float DomainLeftEdge [MAX_DIMENSION];
  enzo_float DomainRightEdge[MAX_DIMENSION];

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


  int GridRank;
  int GridDimension[MAX_DIMENSION]; // total dimensions of all grids
  int GridStartIndex[MAX_DIMENSION]; // starting index of the active region
  int GridEndIndex[MAX_DIMENSION]; // stoping index of the active region
  enzo_float GridLeftEdge[MAX_DIMENSION]; // starting pos (active problem space)
  enzo_float CellWidth[MAX_DIMENSION];
  int ghost_depth[MAX_DIMENSION];

  // Fields

  int NumberOfBaryonFields;      // active baryon fields
  enzo_float *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS]; // pointers to arrays
  enzo_float *OldBaryonField[MAX_NUMBER_OF_BARYON_FIELDS]; // pointers to old arrays
  int FieldType[MAX_NUMBER_OF_BARYON_FIELDS];

  // Boundary

  int  BoundaryRank;
  int  BoundaryDimension[MAX_DIMENSION];
  int  BoundaryFieldType[MAX_NUMBER_OF_BARYON_FIELDS];
  bc_enum *BoundaryType[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2];
  enzo_float *BoundaryValue[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2]; 

  // problem

  int CycleNumber;
  enzo_float dt;

  // Fluxes
  fluxes ** SubgridFluxes;

};


#endif /* ENZO_ENZO_DESCR_HPP */

