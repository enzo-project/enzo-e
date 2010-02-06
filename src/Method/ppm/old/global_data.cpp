/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @bug       
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */
/***********************************************************************
/
/  GLOBAL DATA DECLARATIONS
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified:   Robert Harkness
/  date:       February, 2004
/
/  PURPOSE:
/    This is the global data, which should be held to a minimum.  Any changes
/    in this file require changes in: WriteGlobalData,
/    ReadGlobalData and InitializeNew.  
/    This file is dual-purposed:
/        1) read with    DEFINE_STORAGE defined for the (single) definition
/        2) read without DEFINE_STORAGE defined for external linkage
/
************************************************************************/

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"

/* debugging, extraction flags */

int debug;
int extract;

/* Problem: 00 = None                    01 = ShockTube
            02 = WavePool                03 = ShockPool  
	    04 = Double-Mach reflection  05 = ShockInABox
	    06 = Implosion               07 = SedovBlast
	    08 = KelvinHelmholtz instability
            09 = Noh Problem
	    20 = 1D Zeldovich Pancake    21 = 1D pressureless collapse
	    22 = Adiabatic expansion     23 = TestGravity
            24 = Spherical infall        25 = TestGravitySphere
	    26 = GravityEquilibriumTest  27 = CollapseTest
	    28 = TestGravityMotion
	    30 = Cosmology simulation
	    50 = ThermalInstability simulation
	    51 = ThermalPancake test
	    60 = TurbulenceSimulation
	                                                                  */
int ProblemType;

/* Hydrodynamics method:
       0 - PPM_DE      1 - PPM_LR (not working)    2 - ZEUS        */

hydro_method HydroMethod;

/* Large and small numbers (i.e. compared to any real quantity).  This may
   be machine and problem dependent. */

float huge_number, tiny_number;

/* Gamma: Ideal gas law constant. */

float Gamma;

/* Flag indicating if the gas is pressureless. */

int PressureFree;

/* Factor to refine by */

int RefineBy;

/* Maximum refinement level (0 = topgrid). */

int MaximumRefinementLevel;
int MaximumGravityRefinementLevel;
int MaximumParticleRefinementLevel;

/* Cell Flagging method:  0 = None
                          1 = FlagCellsToBeRefinedBySlope
			  2 = FlagCellsToBeRefinedByMass (baryon only)
			  3 = FlagCellsToBeRefinedByShocks
			  4 = FlagCellsToBeRefinedByMass (particles only)
	     (disabled)	  5 = FlagCellsToBeRefinedByOverdensity (baryon only)
			  6 = FlagCellsToBeRefinedByJeansLength
                          7 = FlagCellsToBeRefinedByCoolingTime
                          8 = FlagCellsToBeRefinedByMustRefineParticles
                          9 = FlagCellsToBeRefinedByShear
 */

int CellFlaggingMethod[MAX_FLAGGING_METHODS];

/* Flag indicating if the flux correction should be applied. */

int FluxCorrection;

/* This specifies the interpolation method (see typedefs.h). */

interpolation_type InterpolationMethod;
int ConservativeInterpolation;

/* This is the minimum efficiency of combined grid needs to achieve in
   order to be considered better than the two grids from which it formed. */

float MinimumEfficiency;

/* This is the minimum allowable edge size for a new subgrid (>=4) */

int MinimumSubgridEdge;

/* This is the maximum allowable size for a new subgrid (>=2000) */

int MaximumSubgridSize;

/* The number of zones that will be refined around each flagged zone. */

int NumberOfBufferZones;

/* The left and right boundaries of the entire computational domain. */

FLOAT DomainLeftEdge[MAX_DIMENSION], DomainRightEdge[MAX_DIMENSION];

/* Velocity of entire computational domain. */

float GridVelocity[MAX_DIMENSION];

/* HDF names for labels and scales. */

char *DimUnits[MAX_DIMENSION], *DimLabels[MAX_DIMENSION];
char *DataLabel[MAX_NUMBER_OF_BARYON_FIELDS];
char *DataUnits[MAX_NUMBER_OF_BARYON_FIELDS];

/* Region in which refinement is allowed (in problem space). */

FLOAT RefineRegionLeftEdge[MAX_DIMENSION], 
             RefineRegionRightEdge[MAX_DIMENSION];

/* Uniform gravity: on/off flag, direction, and strength. */

int UniformGravity, UniformGravityDirection;
float UniformGravityConstant;

/* point source gravity: on/off flag position, and strength. */

int PointSourceGravity;
FLOAT PointSourceGravityPosition[MAX_DIMENSION];
float PointSourceGravityConstant;
float PointSourceGravityCoreRadius;

#ifdef ISO_GRAV
/* AverageDensity used for ComovingGravitySourceTerm and 
   isolating gravitational potential solve */
float AverageDensity;
#endif

/* SelfGravity (TRUE or FALSE) */

int SelfGravity;

/* CopyGravPotential (TRUE or FALSE) */

int CopyGravPotential;

/* Flag indicating whether or not to use the baryon self-gravity approximation
   (subgrid cells influence are approximated by their projection to the
   current grid). */

int BaryonSelfGravityApproximation;

/* Coefficient in front of source term in Poisson's equations.
   (i.e. Del^phi = GravitationConstant * density, usually 4*Pi*G). */

float GravitationalConstant;

/* S2 Particle size in top grid cell units (usually around 3).  The S2
   particle is S(r) = A*(a/2-r) (if r < a/2, 0 otherwise).  The constant
   A depends on the dimension: 1D) 4/a^2,  2D) 24/(Pi*a^3)  3D) 48/(Pi*a^3). */

float S2ParticleSize;

/* Gravity resolution factor is a float indicating the comparative resolution
   of the gravitational computation compared to the grid (1-2 or so). */

float GravityResolution;

/* Flag to indicate if gravitational potential field should be computed
   and stored. */

int ComputePotential;

/* Flag to indicate output for gravitational potential field. */

int WritePotential;

/* Maximum number of GreensFunctions that will be stored in any time.
   This number must be less than MAX_NUMBER_OF_GREENS_FUNCTIONS. */

int GreensFunctionMaxNumber;

/* Maximum number of words associated with GreensFunction storage
   (Not currently implemented). */

int GreensFunctionMaxSize;

/* Dual energy formalism (TRUE or FALSE). */

int DualEnergyFormalism;

/* Two parameters for the dual energy formalism. */

float DualEnergyFormalismEta1;
float DualEnergyFormalismEta2;

/* This is the particle equivalent of the Courant factor.  It is the maximum
   number of cells a particle is allowed to travel in a single timestep. */

float ParticleCourantSafetyNumber;

/* Radiative cooling on/off flag and associated data. */

int RadiativeCooling;
CoolDataType CoolData;

/* Coupled Radiation-Hydrodynamics on/off flag */

int RadiationHydrodynamics;

/* Gadget Equilibrium cooling on/off flag */

int GadgetEquilibriumCooling;

/* Random Forcing on/off flag and associated data. */ //AK

int     RandomForcing;
FLOAT   RandomForcingEdot;
FLOAT   RandomForcingMachNumber;
fpos_t  BaryonFileNamePosition;

/* Multi-species rate equation flag and associated data. */

int MultiSpecies;
RateDataType RateData;

/* Multi-element metallicity field flag and count. */

int MultiMetals;

/* Type of radiation field. 
   0 - none,                    1 - Haardt & Madau alpha=-1.5
   2 - H&M alpha = -1.8       
   10 - homogenous internal radiation field (a la Renyue's work) */

int RadiationFieldType;
int AdjustUVBackground; 
float SetUVBAmplitude;
float SetHeIIHeatingScale;
RadiationFieldDataType RadiationData;
int RadiationFieldLevelRecompute;

/* ZEUS Hydro artificial viscosity parameters (C1, C2 of Stone & Norman). */

float ZEUSLinearArtificialViscosity;
float ZEUSQuadraticArtificialViscosity;

/* Parameters for MinimumPressureSupport. */

int UseMinimumPressureSupport;
float MinimumPressureSupportParameter;

/* Parameters for statically refined regions. */

FLOAT StaticRefineRegionLeftEdge[MAX_STATIC_REGIONS][MAX_DIMENSION];
FLOAT StaticRefineRegionRightEdge[MAX_STATIC_REGIONS][MAX_DIMENSION];
int   StaticRefineRegionLevel[MAX_STATIC_REGIONS];

/* Processor identifier for this thread/processor */

int MyProcessorNumber;
int NumberOfProcessors;
float CommunicationTime;
int CommunicationDirection;

/* Parameter to indicate if top grid should do parallel IO
   (currently only works for ProblemType == 30). */

int ParallelRootGridIO;
int ParallelParticleIO;
int Unigrid;
int CubeDumpEnabled;
int PartitionNestedGrids;
int ExtractFieldsOnly;
int First_Pass;

/************************************************/
/* Global data for specific problems or methods */
/************************************************/

/* For CellFlaggingMethod = 1,
   The minimum relative slope (da/dx over a) required for refinement. */

float MinimumSlopeForRefinement;

/* For CellFlaggingMethod = 2,
   The minimum refined mass for the ByMass refining scheme
   (Usually, user sets OverDensity and code sets MinimumMass but this can be
    overridden by directely setting MinimumMass). 
   The LevelExponent is used to change the minimum mass with level,
   the formula is MinimumMassForRefinement*pow(RefineBy, level*LevelExponent)*/

float MinimumOverDensityForRefinement[MAX_FLAGGING_METHODS];
float MinimumMassForRefinement[MAX_FLAGGING_METHODS];
float MinimumMassForRefinementLevelExponent[MAX_FLAGGING_METHODS];

/* For CellFlaggingMethod = 3,
   The minimum pressure jump required to be a shock.
   The minimum internal/total energy ratio for a shock. */

float MinimumPressureJumpForRefinement, MinimumEnergyRatioForRefinement;

/* For CellFlaggingMethod = 6,
   The number of cells by which the Jeans length should be resolved. */

float RefineByJeansLengthSafetyFactor;

/* For CellFlaggingMethod = 8,
   The level to which the must refine particles apply */

int   MustRefineParticlesRefineToLevel;

/* For CellFlaggingMethod = 9,   
   The minimum shear (roughly, dv accross two zones) required for 
   refinement.    */

float MinimumShearForRefinement;

/* Noh problem switch: Upper-Right quadrant or full domain */

int NohProblemFullBox;

/* A boolean flag indicating if we are using coordinate comoving with the
   expansion of the universe. */

int   ComovingCoordinates;

/* A flag indicating if we are using star particles. */

int   StarParticleCreation;
int   StarParticleFeedback;
int   NumberOfParticleAttributes;

/* Parameters governing certain time or redshift-dependent actions. */

int   TimeActionType[MAX_TIME_ACTIONS];
FLOAT TimeActionTime[MAX_TIME_ACTIONS];
FLOAT TimeActionRedshift[MAX_TIME_ACTIONS];
float TimeActionParameter[MAX_TIME_ACTIONS];

/* Parameters for direct unigrid dumps of entire top grid */

char  *CubeDumps[MAX_CUBE_DUMPS];

/* Parameters governing whether tracer particles are on or off. */

int   TracerParticleOn;
FLOAT TracerParticleCreationSpacing;
FLOAT TracerParticleCreationLeftEdge[MAX_DIMENSION];
FLOAT TracerParticleCreationRightEdge[MAX_DIMENSION];

int   ParticleTypeInFile;

int   ExternalBoundaryIO;
int   ExternalBoundaryTypeIO;
int   ExternalBoundaryValueIO;
int   ExternalBoundaryField;
int   SimpleConstantBoundary;

Eint64 TaskMemory[MAX_NUMBER_OF_TASKS];
int    TaskMap[MAX_NUMBER_OF_TASKS];

/* Zhiling Lan's modified code */

#ifdef MPI_INSTRUMENTATION
double GlobalCommunication;
double RecvComm;
double WaitComm;
double timer[40];
int counter[40];
char name[20];
FILE *filePtr;
char tracename[20];
FILE *tracePtr;
int traceMPI;
char memtracename[20];
FILE *memtracePtr;
int traceMEM;
double starttime, endtime;
double Start_Wall_Time, End_Wall_Time, WallTime;
int flagging_count, in_count, out_count, moving_count;
float flagging_pct, moving_pct;
#endif /* MPI_INSTRUMENTATION */

/* Storage Resource Broker */

char *SRBprefix;
char *SRBcwd;

/* New Movie Data */

int MovieDataField[MAX_MOVIE_FIELDS];
int MovieSkipTimestep;
char *NewMovieName;
int NewMovieDumpNumber;
int NewMovieEntries;
long *MovieEntriesPP;
int MaxMovieFilenum;
int NewMovieParticleOn;

/* ran1 initialization flag for star_maker5 */

int ran1_init;
