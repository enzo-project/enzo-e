// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar  3 23:01:51 PST 2011
/// @brief    [\ref Enzo] Declaration of the EnzoBlock class

#ifndef ENZO_ENZO_BLOCK_HPP
#define ENZO_ENZO_BLOCK_HPP

enum field_type {
  field_bfieldx,
  field_bfieldx_rx,
  field_bfieldx_ry,
  field_bfieldx_rz,
  field_bfieldy,
  field_bfieldy_rx,
  field_bfieldy_ry,
  field_bfieldy_rz,
  field_bfieldz,
  field_bfieldz_rx,
  field_bfieldz_ry,
  field_bfieldz_rz,
  field_color,
  field_density,
  field_dens_rx,
  field_dens_ry,
  field_dens_rz,
  field_internal_energy,
  field_total_energy,
  field_velocity_x,
  field_velocity_y,
  field_velocity_z,
  field_velox,
  field_velox_rx,
  field_velox_ry,
  field_velox_rz,
  field_veloy,
  field_veloy_rx,
  field_veloy_ry,
  field_veloy_rz,
  field_veloz,
  field_veloz_rx,
  field_veloz_ry,
  field_veloz_rz
};

//----------------------------------------------------------------------

class EnzoBlock : public CommBlock

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
  friend class EnzoInitialSedovArray3;

  //----------------------------------------------------------------------
  // functions

  static void initialize (EnzoConfig * enzo_config, FieldDescr *);

  //----------------------------------------------------------------------
  // variables

public:
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

  static int field_index_[NUM_FIELDS];

  static int GridRank;

  static int ghost_depth[MAX_DIMENSION];

  // Fields

  static int NumberOfBaryonFields;      // active baryon fields

  static int FieldType[MAX_NUMBER_OF_BARYON_FIELDS];

public: // interface

#if defined CONFIG_USE_CHARM && ! defined PREPARE_AMR

  /// Initialize the EnzoBlock chare array
  EnzoBlock
  (
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   double xm, double ym, double zm,
   double hx, double hy, double hz,
   int num_field_blocks) throw();

#else /* defined(CONFIG_USE_CHARM) && ! defined (PREPARE_AMR) */

  /// Initialize the EnzoBlock chare array
  EnzoBlock
  (
   int ix, int iy, int iz,
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   double xm, double ym, double zm,
   double hx, double hy, double hz,
   int num_field_blocks) throw();


#endif /* defined(CONFIG_USE_CHARM) && ! defined (PREPARE_AMR) */

#ifdef CONFIG_USE_CHARM
  /// Initialize a migrated EnzoBlock
  EnzoBlock (CkMigrateMessage *m) 
    : CommBlock (m)
  {
    TRACE("CkMigrateMessage");
    //    initialize();
  };

#endif  /* CONFIG_USE_CHARM */

#ifdef CONFIG_USE_CHARM

  /// Pack / unpack the EnzoBlock in a CHARM++ program
  void pup(PUP::er &p)
  { 

    TRACEPUP;
    TRACE ("BEGIN EnzoBlock::pup()");

    CommBlock::pup(p);


    p | Time_;
    p | CycleNumber;
    p | OldTime;
    p | dt;

    WARNING("EnzoBlock::pup()", "skipping AccelerationField_ (not used)");
    WARNING("EnzoBlock::pup()", "skipping SubgridFluxes (not used)");

    PUParray(p,GridLeftEdge,MAX_DIMENSION); 
    PUParray(p,GridDimension,MAX_DIMENSION); 
    PUParray(p,GridStartIndex,MAX_DIMENSION); 
    PUParray(p,GridEndIndex,MAX_DIMENSION); 
    PUParray(p,CellWidth,MAX_DIMENSION);

    if (p.isUnpacking()) {
      for (int field = 0; field < EnzoBlock::NumberOfBaryonFields; field++) {
	BaryonField[field] = (enzo_float *)block_->field_block(0)->field_values(field);
      }
    }

    WARNING("EnzoBlock::pup()", "skipping OldBaryonField[] [not used]");

    TRACE ("END EnzoBlock::pup()");

  };

#endif /*CONFIG_USE_CHARM */

  /// Implementation of initialization in constructors
  void initialize_enzo_();

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

  /// Return the Cello FieldBlock index for the given field type
  int index (enum field_type type) const
  { return field_index_[type];}

  /// Set EnzoBlock's cycle
  virtual void set_cycle (int cycle) throw();

  /// Set EnzoBlock's time
  virtual void set_time (double time) throw();

  /// Set EnzoBlock's dt
  virtual void set_dt (double dt) throw();

  /// Initialize EnzoBlock
  virtual void initialize () throw();

private: // attributes

  enzo_float Time_;

public: // attributes (YIKES!)

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

