// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar  3 23:01:51 PST 2011
/// @brief    [\ref Enzo] Declaration of the EnzoBlock class

#ifndef ENZO_ENZO_BLOCK_HPP
#define ENZO_ENZO_BLOCK_HPP

//----------------------------------------------------------------------

class EnzoBlock : public Block

{

  /// @class    EnzoBlock
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] An EnzoBlock is a Block with Enzo data

  friend class IoEnzoBlock;

  friend class EnzoSimulationCharm;
  friend class EnzoSimulationMpi;
  friend class EnzoTimestep;
  friend class EnzoTimestepPpm;
  friend class EnzoTimestepPpml;
  friend class EnzoInitialImplosion2;

  //----------------------------------------------------------------------
  // functions

  static void initialize (Parameters * parameters, FieldDescr *);

  //----------------------------------------------------------------------
  // variables

  /// Boundary

  static int  BoundaryRank;
  static int  BoundaryDimension[MAX_DIMENSION];
  static int  BoundaryFieldType[MAX_NUMBER_OF_BARYON_FIELDS];
  static bc_enum *BoundaryType[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2];
  static enzo_float *BoundaryValue[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2]; 

  /// Cosmology

  static int ComovingCoordinates;
  static int UseMinimumPressureSupport;
  static enzo_float MinimumPressureSupportParameter;
  static enzo_float ComovingBoxSize;
  static enzo_float HubbleConstantNow;
  static enzo_float OmegaMatterNow;
  static enzo_float OmegaLambdaNow;
  static enzo_float MaxExpansionRate;

  // Chemistry

  static int MultiSpecies;

  // Gravity

  static int GravityOn;

  // Physics

  static int PressureFree;
  static enzo_float Gamma;
  static enzo_float GravitationalConstant;

  // Problem-specific

  static int ProblemType;

  // Method PPM

  static int PPMFlatteningParameter;
  static int PPMDiffusionParameter;
  static int PPMSteepeningParameter;

  // Parallel

  //  static int ProcessorNumber;

  // Numerics

  static int DualEnergyFormalism;
  static enzo_float DualEnergyFormalismEta1;
  static enzo_float DualEnergyFormalismEta2;

  static enzo_float pressure_floor;
  static enzo_float density_floor;
  static enzo_float number_density_floor;
  static enzo_float temperature_floor;

  static enzo_float CourantSafetyNumber;
  static enzo_float InitialRedshift;
  static enzo_float InitialTimeInCodeUnits;

  // Domain

  static enzo_float DomainLeftEdge [MAX_DIMENSION];
  static enzo_float DomainRightEdge[MAX_DIMENSION];

  // PPM

  static int field_density;
  static int field_total_energy;
  static int field_internal_energy;
  static int field_velocity_x;
  static int field_velocity_y;
  static int field_velocity_z;

  static int field_color;

  // PPM

  static int field_velox;
  static int field_veloy;
  static int field_veloz;
  static int field_bfieldx;
  static int field_bfieldy;
  static int field_bfieldz;

  static int field_dens_rx;
  static int field_velox_rx;
  static int field_veloy_rx;
  static int field_veloz_rx;
  static int field_bfieldx_rx;
  static int field_bfieldy_rx;
  static int field_bfieldz_rx;

  static int field_dens_ry;
  static int field_velox_ry;
  static int field_veloy_ry;
  static int field_veloz_ry;
  static int field_bfieldx_ry;
  static int field_bfieldy_ry;
  static int field_bfieldz_ry;

  static int field_dens_rz;
  static int field_velox_rz;
  static int field_veloy_rz;
  static int field_veloz_rz;
  static int field_bfieldx_rz;
  static int field_bfieldy_rz;
  static int field_bfieldz_rz;


  static int GridRank;

  static int ghost_depth[MAX_DIMENSION];

  // Fields

  static int NumberOfBaryonFields;      // active baryon fields

  static int FieldType[MAX_NUMBER_OF_BARYON_FIELDS];



public: // interface

  /// Initialize the EnzoBlock chare array
  EnzoBlock
  (
   int ix, int iy, int iz,
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   double xm, double ym, double zm,
   double hx, double hy, double hz,
#ifdef CONFIG_USE_CHARM
   CkChareID proxy_patch,
#endif
   int patch_id,
   int patch_rank,
   int num_field_blocks) throw();

#ifdef CONFIG_USE_CHARM
  /// Initialize a migrated Block
  EnzoBlock (CkMigrateMessage *m) {};

  /// Initialize the EnzoBlock chare array
  EnzoBlock
  (
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   double xm, double ym, double zm,
   double hx, double hy, double hz,
   CkChareID proxy_patch,
   int patch_id,
   int patch_rank,
   int num_field_blocks) throw();

  /// Pack / unpack the Block in a CHARM++ program
  void pup(PUP::er &p)
  { ERROR("EnoBlock::pup","Not Implemented"); };

#endif

  /// Destructor
  virtual ~EnzoBlock() throw();


  /// Write attributes, e.g. to stdout for debugging
  void write(FILE *fp=stdout) throw ();

  //----------------------------------------------------------------------
  // Enzo attribute access functions
  //----------------------------------------------------------------------

  /// When Enzo accesses Time, refresh Cello time_ to Time_
  double Time() { Time_ = time_; return Time_; };

  //----------------------------------------------------------------------
  // Original Enzo functions
  //----------------------------------------------------------------------

  //  enzo_float ComputeTimeStep();

  /// Compute the ratio of specific heats
  int ComputeGammaField(enzo_float *GammaField);

  /// Compute the pressure field at the given time) - dual energy
  int ComputePressureDualEnergyFormalism
  ( enzo_float time, enzo_float *pressure );

  /// Compute the pressure field at the given time
  int ComputePressure(enzo_float time, enzo_float *pressure);

  /// Compute the temperature field
  int ComputeTemperatureField(enzo_float *temperature);

  /// Computes the expansion factors (a & dadt) at the requested time
  int CosmologyComputeExpansionFactor
  (enzo_float time, enzo_float *a, enzo_float *dadt);

  /// Computes the maximum allowed expansion timestep at given time
  int CosmologyComputeExpansionTimestep
  (enzo_float time, enzo_float *dtExpansion);

  /// Compute and return the cosmology units
  int CosmologyGetUnits
  ( enzo_float *DensityUnits, 
    enzo_float *LengthUnits, 
    enzo_float *TemperatureUnits, 
    enzo_float *TimeUnits, 
    enzo_float *VelocityUnits, 
    enzo_float Time);

  /// Find field type field in array field_type, returning the index
  /// into the field array or -1 if it is not there
  int FindField(int field, int farray[], int numfields);

  /// Identify certain commonly used variables from the list
  int IdentifyPhysicalQuantities
  ( int &DensNum, 
    int &GENum, 
    int &Vel1Num, 
    int &Vel2Num, 
    int &Vel3Num, 
    int &TENum);

  /// Identify the multi-species fields
  int IdentifySpeciesFields
  ( int &DeNum, 
    int &HINum, 
    int &HIINum, 
    int &HeINum, 
    int &HeIINum, 
    int &HeIIINum, 
    int &HMNum, 
    int &H2INum, 
    int &H2IINum, 
    int &DINum, 
    int &DIINum, 
    int &HDINum);

  /// Copy the current baryon fields to the old baryon fields
  int SetExternalBoundaryValues();

  /// Set the energy to provide minimal pressure support
  int SetMinimumSupport(enzo_float &MinimumSupportEnergyCoefficient);

  /// Solve the hydro equations, saving subgrid fluxes
  int SolveHydroEquations ( int CycleNumber, enzo_float time, enzo_float dt);

  /// Set external boundary values
  int SetExternalBoundary
  ( int FieldRank, 
    int GridDims[], 
    int GridOffset[], 
    int StartIndex[], 
    int EndIndex[], 
    enzo_float *Field, 
    int FieldType);

  /// Solve the mhd equations (with ppml), saving subgrid fluxes
  int SolveMHDEquations(FieldDescr *,  enzo_float dt);

public: // functions (TEMPORARILY PUBLIC)

  /// Set Block's cycle
  virtual void set_cycle (int cycle) throw();

  /// Set Block's time
  virtual void set_time (double time) throw();

  /// Set Block's dt
  virtual void set_dt (double dt) throw();

  /// Initialize Block
  virtual void initialize () throw();

public: // attributes (YIKES!)

  enzo_float Time_;

  int CycleNumber;

  enzo_float OldTime;
  enzo_float dt;

  /// cell cntr acceleration at n+1/2
  enzo_float *AccelerationField[MAX_DIMENSION]; 

  /// Fluxes
  fluxes ** SubgridFluxes;

  /// starting pos (active problem space)
  enzo_float GridLeftEdge[MAX_DIMENSION]; 
  /// total dimensions of all grids
  int GridDimension[MAX_DIMENSION]; 
  /// starting index of the active region
  int GridStartIndex[MAX_DIMENSION]; 
  /// stoping index of the active region
  int GridEndIndex[MAX_DIMENSION]; 
  enzo_float CellWidth[MAX_DIMENSION];
  /// pointers to arrays
  enzo_float *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS]; 
  /// pointers to old arrays
  enzo_float *OldBaryonField[MAX_NUMBER_OF_BARYON_FIELDS]; 

};

#endif /* ENZO_ENZO_BLOCK_HPP */

