// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar  3 23:01:51 PST 2011
/// @brief    [\ref Enzo] Declaration of the EnzoBlock class

#ifndef ENZO_ENZO_BLOCK_HPP
#define ENZO_ENZO_BLOCK_HPP

class CProxy_EnzoBlock;

//----------------------------------------------------------------------

#include "enzo.decl.h"

class EnzoBlock : public Block

{

  /// @class    EnzoBlock
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] An EnzoBlock is a Block with Enzo data

  friend class IoEnzoBlock;

  friend class EnzoSimulation;
  friend class EnzoTimestep;
  friend class EnzoTimestepPpm;
  friend class EnzoTimestepPpml;
  friend class EnzoInitialGrackleTest;
  friend class EnzoInitialImplosion2;
  friend class EnzoInitialSedovArray2;

  //----------------------------------------------------------------------
  // functions

  static void initialize (EnzoConfig * enzo_config, FieldDescr * field_descr);

  //----------------------------------------------------------------------
  // variables

public:

  /// Cosmology

  static int UseMinimumPressureSupport;
  static enzo_float MinimumPressureSupportParameter;
  static enzo_float ComovingBoxSize;
  static enzo_float HubbleConstantNow;
  static enzo_float OmegaMatterNow;
  static enzo_float OmegaLambdaNow;
  static enzo_float MaxExpansionRate;

  // Chemistry

  static int MultiSpecies;

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

  static int GridRank;

  static int ghost_depth[MAX_DIMENSION];

  // Fields

  static int NumberOfBaryonFields;      // active baryon fields

public: // interface

  /// Initialize the EnzoBlock chare array
  EnzoBlock
  (
   Index index,
   int nx, int ny, int nz,
   int num_field_blocks,
   int count_adapt,
   int cycle, double time, double dt,
   int narray, char * array, int op_array,
   int num_face_level, int * face_level,
   bool testing=false) throw();

  /// Initialize an empty EnzoBlock
  EnzoBlock()  { };

  /// Initialize a migrated EnzoBlock
  EnzoBlock (CkMigrateMessage *m) 
    : Block (m)
  {
    TRACE("CkMigrateMessage");
    //    initialize();
  };

  /// Pack / unpack the EnzoBlock in a CHARM++ program
  void pup(PUP::er &p);

  /// Implementation of initialization in constructors
  void initialize_enzo_();

  /// Destructor
  virtual ~EnzoBlock() throw();

  //--------------------------------------------------
  // Charm++ virtual 
  //--------------------------------------------------

  virtual const CProxy_Block proxy_array() const 
  { return thisProxy; }

  virtual const CProxyElement_Block proxy_element() const 
  { return thisProxy[thisIndex]; }

  /// Write attributes, e.g. to stdout for debugging
  void write(FILE *fp=stdout) throw ();

  //----------------------------------------------------------------------
  // Original Enzo functions
  //----------------------------------------------------------------------

  //  enzo_float ComputeTimeStep();

  /// Compute the ratio of specific heats
  int ComputeGammaField(enzo_float *GammaField,
			int comoving_coordinates);

  /// Compute the pressure field at the given time) - dual energy
  int ComputePressureDualEnergyFormalism
  ( enzo_float time, 
    enzo_float *pressure,
    int comoving_coordinates);

  /// Compute the pressure field at the given time
  int ComputePressure(enzo_float time, enzo_float *pressure, 
		      int comoving_coordinates);

  /// Compute the temperature field
  int ComputeTemperatureField (enzo_float *temperature, 
			       int comoving_coordinates);

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

  /// Set the energy to provide minimal pressure support
  int SetMinimumSupport(enzo_float &MinimumSupportEnergyCoefficient,
			int comoving_coordinates);

  /// Solve the hydro equations using PPM
  int SolveHydroEquations ( enzo_float time, 
			    enzo_float dt,
			    int comoving_coordinates);

  /// Solve the hydro equations using Enzo 3.0 PPM
  int SolveHydroEquations3 ( enzo_float time, enzo_float dt);

  /// Solve the mhd equations (with ppml), saving subgrid fluxes
  int SolveMHDEquations(const FieldDescr *,  enzo_float dt);

  /// Set EnzoBlock's dt
  virtual void set_dt (double dt) throw();

  /// Set EnzoBlock's stopping criteria
  virtual void set_stop (bool stop) throw();

  /// Initialize EnzoBlock
  virtual void initialize () throw();

public: /// entry methods

  /// Compute sum, min, and max of g values for EnzoMethodTurbulence
  void method_turbulence_begin();

  /// Perform the necessary reductions
  CkReductionMsg * r_method_turbulence(int n, CkReductionMsg ** msgs);

  /// Compute sum, min, and max of g values for EnzoMethodTurbulence
  void p_method_turbulence_end(CkReductionMsg *msg);

  /// EnzoMethodGravityCg entry method: DOT ==> refresh P
  template <class T>
  void r_cg_loop_0a (CkReductionMsg * msg) ;  

  /// EnzoMethodGravityCg entry method: ==> refresh P
  template <class T>
  void r_cg_loop_0b (CkReductionMsg * msg) ;  

  /// EnzoMethodGravityCg entry method: DOT(R,R) after shift
  template <class T>
  void r_cg_shift_1 (CkReductionMsg * msg) ;

  /// EnzoMethodGravityCg entry method: refresh P for MATVEC
  template <class T>
  void r_cg_loop_1 (CkReductionMsg * msg) ;

  /// EnzoMethodGravityCg entry method: DOT(P,AP)
  template <class T>
  void r_cg_loop_3 (CkReductionMsg * msg) ;

  /// EnzoMethodGravityCg entry method: DOT(R,R)
  template <class T>
  void r_cg_loop_5 (CkReductionMsg * msg) ;

  /// EnzoMethodGravityCg entry method: 
  /// perform the necessary reductions for shift
  CkReductionMsg * r_method_gravity_cg(int n, CkReductionMsg ** msgs);

  /// EnzoMethodGravityBiCGStab entry method: SUM(B) and COUNT(B)
  template <class T>
  void r_gravity_bicgstab_start_1(CkReductionMsg* msg);  

  /// EnzoMethodGravityBiCGStab entry method: DOT(R,R)
  template <class T>
  void r_gravity_bicgstab_start_3(CkReductionMsg* msg);  

  /// EnzoMethodGravityBiCGStab entry method: refresh P
  void p_gravity_bicgstab_loop_1();  

  /// EnzoMethodGravityBiCGStab entry method: refresh Y
  void p_gravity_bicgstab_loop_3();

  /// EnzoMethodGravityBiCGStab entry method: DOT(V,R0)
  template <class T>
  void r_gravity_bicgstab_loop_5(CkReductionMsg* msg);

  /// EnzoMethodGravityBiCGStab entry method: refresh Q
  void p_gravity_bicgstab_loop_7();

  /// EnzoMethodGravityBiCGStab entry method: refresh Y
  void p_gravity_bicgstab_loop_9();

  /// EnzoMethodGravityBiCGStab entry method: DOT(U,U) and DOT(U,Q)
  template <class T>
  void r_gravity_bicgstab_loop_11(CkReductionMsg* msg);

  /// EnzoMethodGravityBiCGStab entry method: DOT(R,R) and DOT(R,R0)
  template <class T>
  void r_gravity_bicgstab_loop_13(CkReductionMsg* msg);

  /// EnzoMethodGravityBiCGStab entry method: ITER++
  template <class T>
  void r_gravity_bicgstab_loop_15(CkReductionMsg* msg);

  /// EnzoMethodGravityBiCGStab entry method: 
  /// perform the necessary reductions for shift
  CkReductionMsg* r_method_gravity_bicgstab(int n, CkReductionMsg** msgs);

  void p_enzo_matvec()
  {      enzo_matvec_(); }
  void r_enzo_matvec(CkReductionMsg * msg)
  {      enzo_matvec_(); delete msg; }

  /// EnzoMethodSolverMlat entry method: receive face data for refresh
  void p_mg_receive_face
  (int n, char buffer[],  int type_refresh, 
   int if3[3], int ic3[3], int count = 0);

  /// EnzoMethodSolverMg0
  template <class T>
  void p_mg0_pre_smooth(CkReductionMsg * msg);
  template <class T>
  void p_mg0_restrict_send(CkReductionMsg * msg);
  template <class T>
  void p_mg0_restrict_recv(FieldMsg * msg);
  template <class T>
  void p_mg0_prolong_recv(FieldMsg * msg);
  template <class T>
  void p_mg0_post_smooth(CkReductionMsg * msg);

  void mg_sync_reset()             { mg_sync_.reset(); }
  void mg_sync_set_stop(int value) { mg_sync_.set_stop(value); }
  bool mg_sync_next()         { return mg_sync_.next(); };
  int  mg_sync_value()        { return mg_sync_.value(); };
  int  mg_sync_stop()         { return mg_sync_.stop(); };

  void mg_iter_clear() { mg_iter_ = 0; }
  void mg_iter_increment() { ++mg_iter_; }
  int mg_iter() const {return mg_iter_; }

protected: // functions

  void enzo_matvec_() ;
  void gravity_bicgstab_matvec_1_();
  void gravity_bicgstab_matvec_2_();

protected: // attributes

  
  // MG SOLVER
  Sync mg_sync_;

  // MG iteration count
  int mg_iter_;

public: // attributes (YIKES!)

  union {
    enzo_float dt;
    enzo_float dtFixed;
  };

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

  /// Data for turbulence reductions
  double method_turbulence_data [MAX_TURBULENCE_ARRAY];

};

#endif /* ENZO_ENZO_BLOCK_HPP */

