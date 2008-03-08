/***********************************************************************
/
/  GRID CLASS
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified:   Many times by AK, DC, RH, JB, DR...
/
/  PURPOSE:
/
************************************************************************/

#ifndef GRID_DEFINED__
#define GRID_DEFINED__

#include "ProtoSubgrid.h"
#include "ListOfParticles.h"
#include "region.h"
#include "FastSiblingLocator.h"

#ifdef FLUX_FIX
#include "TopGridData.h"
#endif

#ifdef DC_OPT_SIBSUB
struct HierarchyEntry;
#endif

#ifdef ANALYSIS_TOOLS
#   include "AnalyzeClusters.h"
#endif

//extern int CommunicationDirection;

struct LevelHierarchyEntry;

#include "EnzoArray.h"

class grid
{
 private:
//
//  General grid class data
//
  int GridRank;                        // number of dimensions
  int GridDimension[MAX_DIMENSION];    // total dimensions of all grids
  int GridStartIndex[MAX_DIMENSION];   // starting index of the active region
                                       //   (zero based)
  int GridEndIndex[MAX_DIMENSION];     // stoping index of the active region
                                       //   (zero based)
  FLOAT GridLeftEdge[MAX_DIMENSION];   // starting pos (active problem space)
  FLOAT GridRightEdge[MAX_DIMENSION];  // ending pos (active problem space)
  float dtFixed;                       // current (fixed) timestep
  FLOAT Time;                          // current problem time
  FLOAT OldTime;                       // time corresponding to OldBaryonField
  int   SubgridsAreStatic;             // 
//
//  Baryon grid data
//
  int    NumberOfBaryonFields;                        // active baryon fields
  float *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS];    // pointers to arrays
  float *OldBaryonField[MAX_NUMBER_OF_BARYON_FIELDS]; // pointers to old arrays
  float *RandomForcingField[MAX_DIMENSION];           // pointers to arrays //AK
  int    FieldType[MAX_NUMBER_OF_BARYON_FIELDS];
  FLOAT *CellLeftEdge[MAX_DIMENSION];
  FLOAT *CellWidth[MAX_DIMENSION];
  fluxes *BoundaryFluxes;

  float  CourantSafetyNumber;                       // Hydro parameter
  int    PPMFlatteningParameter;                    // PPM parameter
  int    PPMDiffusionParameter;                     // PPM parameter
  int    PPMSteepeningParameter;                    // PPM parameter
//
//  Particle data
//
  int    NumberOfParticles;
  FLOAT *ParticlePosition[MAX_DIMENSION];  // pointers to position arrays
  float *ParticleVelocity[MAX_DIMENSION];  // pointers to velocity arrays
  float *ParticleAcceleration[MAX_DIMENSION+1];  // 
  float *ParticleMass;                     // pointer to mass array
  int   *ParticleNumber;                   // unique identifier
  int   *ParticleType;                     // type of particle
  float *ParticleAttribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
//
//  Gravity data
// 
  float *PotentialField;
  float *AccelerationField[MAX_DIMENSION]; // cell cntr acceleration at n+1/2
  float *GravitatingMassField;
  FLOAT  GravitatingMassFieldLeftEdge[MAX_DIMENSION];
  int    GravitatingMassFieldDimension[MAX_DIMENSION];
  FLOAT  GravitatingMassFieldCellSize;     // all dimensions must be the same
  float *GravitatingMassFieldParticles;     // for particles only
  FLOAT  GravitatingMassFieldParticlesLeftEdge[MAX_DIMENSION];
  FLOAT  GravitatingMassFieldParticlesCellSize;
  int    GravitatingMassFieldParticlesDimension[MAX_DIMENSION];
  gravity_boundary_type GravityBoundaryType;
  float  PotentialSum;
  int    GravityGhostZones;
//
//  Top grid parallelism
//
  int ProcLayout[MAX_DIMENSION];
  int ProcLocation[MAX_DIMENSION];
  int ProcNeighbors[MAX_DIMENSION][2];
//
//  Rebuild Hierarchy Temporaries
//
  int *FlaggingField;              // Boolean flagging field (for refinement)
  float *MassFlaggingField;        // Used by mass flagging criterion
//
//  Parallel Information
//
  int ProcessorNumber;
//
// Movie Data Format
//
  int TimestepsSinceCreation; 	// Not really since creation anymore... 
  				// resets everytime the grid outputs
//
// Friends
//
  friend int ExternalBoundary::Prepare(grid *TopGrid);
  friend int ProtoSubgrid::CopyFlaggedZonesFromGrid(grid *Grid);

 public:

// -------------------------------------------------------------------------
//  Main hydro/AMR functions
//

/* Grid constructor (Set all data to null/default state). */

   grid();

/* Grid deconstructor (free up memory usage) */

   ~grid();

/* Read grid data from a file (returns: success/failure) */

  int ReadGrid(FILE *main_file_pointer, int GridID);

/* Get field or particle data based on name or integer 
   defined in typedefs.h. Details are in Grid_CreateFieldArray.C. */

  EnzoArrayInt *CreateFieldArrayInt(field_type field);
  EnzoArrayInt *CreateFieldArrayInt(char *field_name);
  
  EnzoArrayFloat *CreateFieldArrayFloat(field_type field);
  EnzoArrayFloat *CreateFieldArrayFloat(char *field_name);
  
  EnzoArrayFLOAT *CreateFieldArrayFLOAT(field_type field);
  EnzoArrayFLOAT *CreateFieldArrayFLOAT(char *field_name);

/* Write unigrid cubes to a file (returns: success/failure) */

   int WriteCube(char *base_name, int grid_id, int TGdims[]);

/* Write grid data to a file (returns: success/failure) */

   int WriteGrid(FILE *main_file_pointer, char *base_name, int grid_id);

/* Write grid data to separate files (returns: success/failure) */

   int WriteGridX(FILE *main_file_pointer, char *base_name, int grid_id);

/* Write grid memory data to map file */

   int WriteMemoryMap(FILE *file_pointer, char *base_name, int grid_id);

/* Write grid hierarchy only */

   int WriteStuff(FILE *main_file_pointer, char *base_name, int grid_id);

/* Interpolate to specified time and write unigrid cube data to a file
   (returns: success/failure). */

   int WriteCubeInterpolate(FLOAT WriteTime, char *base_name, int grid_id, int TGdims[]);

/* Interpolate to specified time and write grid data to a file
   (returns: success/failure). */

   int WriteGridInterpolate(FLOAT WriteTime, FILE *main_file_pointer, 
			    char *base_name, int grid_id);

/* Compute the timestep constraint for this grid
    (for steps #3 and #4) */

   float ComputeTimeStep();

/* Set the timestep in this grid to the timestep in the argument
    (for step #3) */

   void SetTimeStep(float dt) {dtFixed = dt;};

/* Check timestep (dtFixed) against argument (return fail if dtFixed > dt).
    (for step #4) */

   int CheckTimeStep(float dt) {return ((dtFixed > dt) ? FAIL : SUCCESS);};

/* Return time, timestep */

   FLOAT ReturnTime() {return Time;};
   float ReturnTimeStep() {return dtFixed;};

/* Baryons: Interpolate (parental) grid in argument to current grid.
            (returns success or fail).
    (for step #16) */

   int InterpolateBoundaryFromParent(grid *ParentGrid);

/* Baryons: Copy current solution to Old solution (returns success/fail)
    (for step #16) */

   int CopyBaryonFieldToOldBaryonField();

/* Copy potential field to baryon potential for output purposes. */

   int CopyPotentialToBaryonField();

/* Baryons: Update boundary according to the external boundary values
    (for step #16) */

   int SetExternalBoundaryValues(ExternalBoundary *Exterior);

/* Baryons: solve hydro equations in this grid (returns: the fluxes of the
           subgrids in the argument).  Returns SUCCESS or FAIL.
    (for step #16) */

   int SolveHydroEquations(int CycleNumber, int NumberOfSubgrids, 
			   fluxes *SubgridFluxes[], int level);

/* Baryons: return pointer to the BoundaryFluxes of this grid */

   int ReturnFluxDims(fluxes &f, int RefinementFactors[]);

/* Baryons: prepare and clear the accumulated boundary fluxes for this grid.
    (for step #16) */

   void PrepareBoundaryFluxes();
   void ClearBoundaryFluxes();

/* Baryons: projected solution in current grid to the grid in the 
           argument which must have a lower resolution (i.e. downsample 
           the current grid to the appropriate level).
    (for step #18) */

   int ProjectSolutionToParentGrid(grid &ParentGrid);

/* Baryons: return boundary fluxes from this grid.  Downsample them to
           the refinement factors specified in the argument.
	   Returns FAIL or SUCCESS.
    (for step #19) */

   int GetProjectedBoundaryFluxes(grid *ParentGrid, fluxes &ProjectedFluxes);

/* Return the refinement factors as compared to the grid in the argument
   (integer version) (for step #19) */

   void ComputeRefinementFactors(grid *SubGrid, int RefinementFactors[]) {
     int dim;
     for (dim = 0; dim < GridRank; dim++) RefinementFactors[dim] = 
	 int( CellWidth[dim][0] / SubGrid->CellWidth[dim][0] + 0.5);
     for (dim = GridRank; dim < MAX_DIMENSION; dim++)
       RefinementFactors[dim] = 1;
   };

/* Return the refinement factors as compared to the grid in the argument
   (float version) (for step #19) */

   void ComputeRefinementFactorsFloat(grid *SubGrid, float Factors[]) {
     int dim;
     for (dim = 0; dim < GridRank; dim++) Factors[dim] = 
       (*CellWidth[dim]) / (*(SubGrid->CellWidth[dim]));;
     for (dim = GridRank; dim < MAX_DIMENSION; dim++)
       Factors[dim] = 1.0;
   };

/* Baryons: Search for redundant overlap between two sets of fluxes (other
            and refined).  If found, set the refined fluxes equal to the
	    initial fluxes so there will be no double corrections.(step #19) */

   void CorrectRedundantFluxes(fluxes *OtherFluxes, fluxes *InitialFluxes, 
                               fluxes *RefinedFluxes);

/* Baryons: correct for better flux estimates produced by subgrids
           (i.e given the initial flux estimates and the subgrid flux 
	   estimates, correct the grid to account for the subgrid 
	   flux estimates).  Returns SUCCESS or FAIL.
    (for step #19) */

#ifdef FLUX_FIX
   int CorrectForRefinedFluxes(fluxes *InitialFluxes, fluxes *RefinedFluxes,
			       fluxes *BoundaryFluxesThisTimeStep,
			       int SUBlingGrid,
			       TopGridData *MetaData);
#else
   int CorrectForRefinedFluxes(fluxes *InitialFluxes, fluxes *RefinedFluxes,
			       fluxes *BoundaryFluxesThisTimeStep);
#endif

/* Baryons: add the fluxes pointed to by the argument to the boundary fluxes
            of this grid (sort of for step #16).  Note that the two fluxes
	    must have the same size. */

   int AddToBoundaryFluxes(fluxes *BoundaryFluxesToBeAdded);

/* set new time (time += dt)
    (step #21) */

   void SetTimeNextTimestep() {Time += dtFixed;};

/* set time of this grid (used in setup) */

   void SetTime(FLOAT NewTime) {Time = NewTime;};

/* set hydro parameters (used in setup) */

   void SetHydroParameters(float co, int p1, int p2, int p3) 
     {
       CourantSafetyNumber    = co;
       PPMFlatteningParameter = p1;
       PPMDiffusionParameter  = p2;
       PPMSteepeningParameter = p3;
     }

/* Baryons: compute the pressure at the requested time. */

   int ComputePressure(FLOAT time, float *pressure);

/* Baryons: compute the pressure at the requested time using the dual energy
            formalism. */

   int ComputePressureDualEnergyFormalism(FLOAT time, float *pressure);

/* Baryons: compute the temperature. */

   int ComputeTemperatureField(float *temperature);

/* Baryons: compute the temperature at the requested time using
   Gadget equilibrium cooling. */

   int GadgetComputeTemperature(FLOAT time, float *temperature);

/* Baryons: compute the temperatre at the requested time using the dual energy
   formalism when using Gadget equilibrium cooling. */

   int GadgetComputeTemperatureDEF(FLOAT time, float *temperature);

/* Baryons: compute X-ray emissivity in specified band. */

   int ComputeXrayEmissivity(float *temperature,
			     float *xray_emissivity, float keV1, float keV2,
			     char *XrayFileName);

/* Baryons: compute number density of ionized elements (just O7 and O8). */

   int ComputeElementalDensity(float *temperature, float *elemental_density,
			       int Type);

/* Baryons: compute the ratio of specific heats. */

   int ComputeGammaField(float *GammaField);

/* Baryons: compute the cooling time. */

   int ComputeCoolingTime(float *cooling_time);

/* Baryons & DualEnergyFormalism: Restore consistency between total and
                                  internal energy fields. */

   int RestoreEnergyConsistency(int Region);

/* Returns some grid info. */

   int ReturnGridInfo(int *Rank, int Dims[], FLOAT Left[], FLOAT Right[]);

/* Subtracts kinetic component from total energy. */

   int ConvertTotalEnergyToGasEnergy();

/* Sets the energy to provide Jean's level support (Zeus: returns coeff). */
   
   int SetMinimumSupport(float &MinimumSupportEnergyCoefficient);

/* Debugging support. */

   int DebugCheck(char *message = "Debug");

// -------------------------------------------------------------------------
// Functions used for analysis
//

/* Calculate the angular momentum of a grid (given center). */

   int CalculateAngularMomentum(FLOAT Center[], float AngularMomentum[],
				float MeanVelocity[], float DMVelocity[],
				FLOAT CenterOfMass[], FLOAT DMCofM[]);

/* Find and track density peaks. */

   int AnalyzeTrackPeaks(int level, int ReportLevel);

/* Project some of the fields to a plane. */

   int ProjectToPlane(FLOAT ProjectedFieldLeftEdge[], 
		      FLOAT ProjectedFieldRightEdge[],
		      int ProjectedFieldDims[], float *ProjectedField[], 
		      int ProjectionDimension, int ProjectionSmooth,
                      int NumberOfProjectedFields, int level,
		      int XrayUseLookupTable, float XrayLowerCutoffkeV,
		      float XrayUpperCutoffkeV, char *XrayFileName);

/* Set the fields to zero under the active region of the specified subgrid. */

   int ZeroSolutionUnderSubgrid(grid *Subgrid, int FieldsToZero, 
                                float Value = 1.0);

/* Convert the grid data to particle data for output. */

   int OutputAsParticleData(FLOAT RegionLeftEdge[], FLOAT RegionRightEdge[],
                           ListOfParticles *ParticleList[NUM_PARTICLE_TYPES],
                           float BaseRadius);

/* Output star particles to a binary file */

   int OutputStarParticleInformation(FILE *StarFile);

/* Return some information about the grid. */

   int CollectGridInformation(int &GridMemory, float &GridVolume, 
                              int &NumberOfCells, float &AxialRatio,
                              int &CellsTotal, int &Particles);

/* Output grid information (for movie generation). */

   int OutputGridMovieData(FILE *Gridfptr, FILE *DMfptr, FILE *Starfptr,
			   FLOAT RegionLeftEdge[], FLOAT RegionRightEdge[],
			   FLOAT WriteOutTime, int NumberOfPoints[3],
			   int NumberOfValuesPerPoint[3],
			   char *PointValueNames[3][20], float BaseRadius);

/* Output movie data (sequential format) */

   int WriteNewMovieData(FLOAT RegionLeftEdge[], FLOAT RegionRightEdge[], 
			 FLOAT StopTime, int lastMovieStep, int &cycle);

   int ReturnMovieTimestep() { return TimestepsSinceCreation; };

/* Output tracer particle information (interpolated from baryon grid). */

   int TracerParticleOutputData(FILE *ptr, FLOAT WriteOutTime);

// -------------------------------------------------------------------------
// Functions for radiative cooling and multi-species rate equations
//

/* Solve the radiative cooling/heating equations  */

   int SolveRadiativeCooling();

/* Solve the rate equations. */

   int SolveRateEquations();

/* Compute densities of various species for RadiationFieldUpdate. */

   int RadiationComputeDensities(int level);

// -------------------------------------------------------------------------
// Functions for Gadget cooling

/* O'SHEA: please add comment */

   int GadgetCalculateCooling(float *d, float *e, float *ge,
			      float *u, float *v, float *w,
			      int *in, int *jn, int *kn,
			      int *iexpand, hydro_method *imethod, 
			      int *idual, int *idim,
			      int *is, int *js, int *ks, int *ie, int *je,
			      int *ke, float *dt, float *aye,
			      float *fh, float *utem, float *uxyz,
			      float *uaye, float *urho, float *utim,
			      float *gamma);

/* O'SHEA: please add comment */

   int GadgetCoolingTime(float *d, float *e, float *ge,
			 float *u, float *v, float *w,
			 float *cooltime,
			 int *in, int *jn, int *kn,
			 int *iexpand, hydro_method *imethod, int *idual, int *idim,
			 int *is, int *js, int *ks, int *ie, int *je,
			 int *ke, float *dt, float *aye,
			 float *fh, float *utem, float *uxyz,
			 float *uaye, float *urho, float *utim,
			 float *gamma);


/* O'SHEA: please rename and add comment */

   void Gadgetfind_abundances_and_rates(float logT, float rho, float *ne_guess);

/* O'SHEA: please rename and add comment */

   float Gadgetconvert_u_to_temp(float u, float rho, float *ne_guess);

/* O'SHEA: please add comment */

   float GadgetCoolingRate(float logT, float rho, float *nelec, float redshift);

/* O'SHEA: please add comment */

   float Gadget_EquilibriumCooling(float u_old, float rho, float dt,
				   float *ne_guess, float *utem, float *uxyz, 
				   float *uaye, float *urho,
				   float *utim, float redshift);

/* O'SHEA: please add comment */

   float GadgetCoolingRateFromU(float u, float rho, float *ne_guess, 
				float redshift);

// -------------------------------------------------------------------------
// Functions for grid (re)generation.
//

/* Remove un-needed arrays before rebuilding. */

   void CleanUp();

/* Delete all the fields, but leave other grid data. */

   void DeleteAllFields();

/* Clear mass flagging field (gg #1) */

   void ClearMassFlaggingField();

/* Clear boolean flagging field (gg #0) */

   void ClearFlaggingField();

/* Set boolean flagging field */

   int SetFlaggingField(int &NumberOfFlaggedCells, int level);

/* Set flagging field from static regions */

   int SetFlaggingFieldStaticRegions(int level, int &NumberOfFlaggedCells);

/* Delete flagging field */

   void DeleteFlaggingField();

/* Particles: deposit particles living in this grid into the Mass Flagging
             field (gg #2) */

   void DepositParticlesToMassFlaggingField() {};

/* baryons: add baryon density to mass flaggin field (so the mass flagging
            field contains the mass in the cell (not the density) 
            (gg #3) */

   int AddFieldMassToMassFlaggingField();

/* Flag all points that require refining  (and delete Mass Flagging Field).
     Returns the number of flagged cells.  Returns the number of flagged cells
     (gg #4) */

   int FlagCellsToBeRefinedByMass(int level, int method);

/* Flag all points that require refining by their slope.
     Returns the number of flagged cells.  Returns the number of flagged cells
     (gg #4) */

   int FlagCellsToBeRefinedBySlope();

/* Flag all points that require refinging by the presence of shocks.
     Returns the number of flagged cells.  Returns the number of flagged cells
     (gg #4) */

   int FlagCellsToBeRefinedByShocks();

/* Flag all points that require refining by the Jean's length criterion. */

   int FlagCellsToBeRefinedByJeansLength();

/* Flag all points that require refining by Shear. */

   int FlagCellsToBeRefinedByShear();

/* Flag all cells for which tcool < dx/sound_speed. */

   int FlagCellsToBeRefinedByCoolingTime();

/* Flag all cells which are near a must-refine particle. */

   int FlagCellsToBeRefinedByMustRefineParticles();

/* Flagging all cell adjacent to a previous flagged cell.  Also, remove all
   Flagged cells in the boundary zones and within one zone of the boundary. */

   int FlagBufferZones();

/* Identify new subgrids for this grid (and prove Fermat's last theorem too)
   (gg #5) */

   void IdentifyNewSubgrids(GridList &list);

/* Identify new subgrids for this grid (1x1x1 subgrids).
   (gg #5) */

   void IdentifyNewSubgridsSmall(GridList &list);

/* Coalesce neighbouring subgrids */

   // void CoalesceSubgrids(GridList &list);

/* Inherit properties (rank, baryon field types, etc.) from ParentGrid
   (gg # 5,6) */

   void InheritProperties(grid *ParentGrid);

/* set the grid dimensions, left, right edges and cell quantities based
   on arguments (gg #5,6) */

   void PrepareGrid(int Rank, int Dimensions[], 
		    FLOAT LeftEdge[], FLOAT RightEdge[], int NumParticles);

/* Allocates space for grids (dims and NumberOfBaryonFields must be set). */

   void AllocateGrids();

/* set the grid derived quantites (CellLeftEdge, CellWidth & BoundaryFluxes) */

   void PrepareGridDerivedQuantities();

/* baryons: interpolate field values from the Parent Grid (gg #6).
            Returns SUCCESS or FAIL. */

   int InterpolateFieldValues(grid *ParentGrid);

/* baryons: check for coincident zones between grids & copy if found.
            (correctly includes periodic boundary conditions). */

   int CheckForOverlap(grid *OtherGrid,
		       boundary_type LeftFaceBoundaryCondition[],
		       boundary_type RightFaceBoundaryCondition[],
		       int (grid::*CopyFunction)(grid *OtherGrid,
						 FLOAT EdgeOffset[]));

   /* baryons: check for subgrids adjacent to external boundary with reflecting BCs. */

   int CheckForExternalReflections(
				   boundary_type LeftFaceBoundaryCondition[],
				   boundary_type RightFaceBoundaryCondition[]);

/* David Collins flux correction - July 2005 */
#ifdef FLUX_FIX
   int CheckForSharedFace(grid *OtherGrid,
			       boundary_type LeftFaceBoundaryCondition[],
			       boundary_type RightFaceBoundaryCondition[]);

   int CheckForSharedFaceHelper(grid *OtherGrid,
				     FLOAT EdgeOffset[MAX_DIMENSION]);
#endif

/* baryons: check for overlap between grids & return TRUE if it exists
            (correctly includes periodic boundary conditions). */

   int CheckForPossibleOverlap(grid *OtherGrid,
                       boundary_type LeftFaceBoundaryCondition[],
                       boundary_type RightFaceBoundaryCondition[]);
   int CheckForPossibleOverlapHelper(grid *OtherGrid,
                                        FLOAT EdgeOffset[MAX_DIMENSION]);

/* baryons: copy coincident zone from the (old) grid in the argument
            (gg #7).  Return SUCCESS or FAIL. */

   int CopyZonesFromGrid(grid *GridOnSameLevel, 
			 FLOAT EdgeOffset[MAX_DIMENSION]);

/* gravity: copy coincident potential field zones from grid in the argument
            (gg #7).  Return SUCCESS or FAIL. */

   int CopyPotentialField(grid *GridOnSameLevel, 
			  FLOAT EdgeOffset[MAX_DIMENSION]);

/* baryons: check for coincident zone from the (old) grid in the argument
            (gg #7).  Return SUCCESS or FAIL. */

   int CopyZonesFromGridCountOnly(grid *GridOnSameLevel, int &Overlap);

/* Returns whether or not the subgrids of this grid are static. */

   int AreSubgridsStatic() {return SubgridsAreStatic;};

/* Check the energy conservation. */

   int ComputeEnergy(float EnergySum[]);

/* These two routines add grids to the chaining mesh used in the
   FastSiblingLocator method and use the chaining mesh to find
   possible siblings. */

#ifdef DC_OPT_SIBSUB
   int FastSiblingLocatorAddGrid(ChainingMeshStructure *mesh,
			  HierarchyEntry *HE_for_this = NULL);
#else
   int FastSiblingLocatorAddGrid(ChainingMeshStructure *mesh);
#endif

   int FastSiblingLocatorFindSiblings(ChainingMeshStructure *mesh,
                          SiblingGridList *list,
                          boundary_type LeftBoundaryCondition[],
                          boundary_type RightBoundaryCondition[]);

   /* hack: add density squared field to grid (used in ExtractSection). */

   void CreateDensitySquaredField() {
     int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
     BaryonField[NumberOfBaryonFields] = new float[size];
     for (int i = 0; i < size; i++)
       BaryonField[NumberOfBaryonFields][i] = 
	 BaryonField[0][i]*BaryonField[0][i];
     FieldType[NumberOfBaryonFields++] = Density;
   };

// -------------------------------------------------------------------------
// Functions for use with gravity.
//

/* Set the gravity boundary type of a grid. */

   void SetGravityParameters(gravity_boundary_type Boundary) {
     GravityBoundaryType = Boundary;};
   gravity_boundary_type ReturnGravityBoundaryType() 
     {return GravityBoundaryType;};

/* Gravity: Initialize, the gravitating Mass Field
    (for steps #5, #6). */

   int InitializeGravitatingMassField(int RefinementFactor);

/* Gravity: Initialize, the particle component of the mass field. */

   int InitializeGravitatingMassFieldParticles(int RefinementFactor);

/* Gravity: allocate & clear the GravitatingMassField. */

   int ClearGravitatingMassField();

/* Gravity & baryons: Copy the parent density field to the extra boundary
      region of GravitatingMassField (if any). */

   int CopyParentToGravitatingFieldBoundary(grid *ParentGrid);

/* Gravity & Particles: allocate & clear the GravitatingMassFieldParticles. */

   int ClearGravitatingMassFieldParticles();

/* Baryons: add the baryon mass to the GravitatingMassField. */

   int AddBaryonsToGravitatingMassField();

/* Generic deposit particles/grids to grid (either GravitatingMassField or
   GravitatingMassFieldParticles depending on the value of DepositField). */

   int DepositPositions(FLOAT *Positions[], float *Mass, int Number, 
			int DepositField);

/* deposit particles/grids to grid (if they are on the grid). */

/* int DepositPositionsEdgeOff(float *Positions[], float *Mass, int Number);*/

/* Gravity: Difference potential to get acceleration field. */

   int ComputeAccelerationField(int DifferenceType, int level);

/* Gravity: Interpolate accelerations from other grid. */

   int InterpolateAccelerations(grid *FromGrid);

/* Gravity: Compute particle and grid accelerations. */

   int ComputeAccelerations(int level);

/* Particles: add overlapping ParticleMassField to Target's 
   GravitatingMassField. */

   int CopyOverlappingMassField(grid *TargetGrid, 
				FLOAT EdgeOffset[MAX_DIMENSION]);

/* Gravity: Allocate and make initial guess for PotentialField. */

   int PreparePotentialField(grid *ParentGrid);

/* Gravity: Allocate and make initial guess for PotentialField. */

   int SolveForPotential(int &Done, int level, FLOAT PotentialTime = -1);

/* Gravity: Prepare the Greens Function. */

   int PrepareGreensFunction();
   int PreparePeriodicGreensFunction(region *GreensRegion);

/* Gravity: Copy potential/density into/out of FFT regions. */

   int PrepareFFT(region *InitialRegion, int Field, int DomainDim[]);
   int FinishFFT(region *InitialRegion, int Field, int DomainDim[]);

/* Gravity: set the potential boundary for isolated BC's */

   // int SetIsolatedPotentialBoundary();

/* Gravity: Set the external acceleration fields. */

   int ComputeAccelerationFieldExternal();

/* Particles + Gravity: Clear ParticleAccleration. */

   int ClearParticleAccelerations();

/* Baryons + Gravity: Interpolate the AccelerationField in FromGrid to
             AccelerationFieldForCells at the GridPositions in this grid. */

   int InterpolateGridPositions(grid *FromGrid);

/* Particles + Gravity: Interpolate the AccelerationField in FromGrid to
             ParticleAcceleration at the ParticlePositions in this grid. */

   int InterpolateParticlePositions(grid *FromGrid, int DifferenceType);

/* Generic routine for interpolating particles/grid. */

   int InterpolatePositions(FLOAT *Positions[], int dim, float *Field, 
			    int Number);

/* Gravity: Delete GravitatingMassField. */

   void DeleteGravitatingMassField() {
     delete GravitatingMassField; 
     GravitatingMassField = NULL;
   };

/* Gravity: Delete AccelerationField. */

   void DeleteAccelerationField() {
     for (int dim = 0; dim < GridRank; dim++) {
       delete AccelerationField[dim];
       AccelerationField[dim] = NULL;
     }
   };

/* Gravity: Add fixed, external acceleration to baryons & particles. */

   int AddExternalAcceleration();

/* Gravity: deposit baryons into target GravitatingMassField. */

   int DepositBaryons(grid *TargetGrid, FLOAT DepositTime);

#ifdef ISO_GRAV
/* Gravity: return the PotentialField (for isolating solver) */

   float* GetPotentialField() {return PotentialField;}

/* Gravity: set the PotentialField with a given pointer (for isolating) */

   void SetPotentialField(float *PotentialPtr) {
     if (PotentialField == NULL)
       PotentialField = PotentialPtr;
     else
       fprintf(stderr,"SetPotentialField Error: Potential already allocated\n");
   }

/* Gravity: return the GravitatingMassField (for isolating solver) */

   float* GetGravitatingMassField() {return GravitatingMassField;}
#endif


// -------------------------------------------------------------------------
// Functions for accessing specific baryon fields (all sources combined 
// in Grid_AccessBaryonFields.C)
//
   float* AccessDensity();
   float* AccessTotalEnergy();
   float* AccessGasEnergy();
   float* AccessVelocity1();
   float* AccessVelocity2();
   float* AccessVelocity3();
   float* AccessElectronDensity();
   float* AccessHIDensity();
   float* AccessHeIDensity();
   float* AccessHeIIDensity();
   float* AccessRadiationFrequency0();
   float* AccessRadiationFrequency1();
   float* AccessRadiationFrequency2();
   float* AccessRadiationFrequency3();
   float* AccessRadiationFrequency4();
   float* AccessRadiationFrequency5();


// -------------------------------------------------------------------------
// Functions for accessing various grid-based information
//
   int GetGridRank() {return GridRank;}
   int GetGridDimension(int Dimension) {return GridDimension[Dimension];}
   int GetGridStartIndex(int Dimension) {return GridStartIndex[Dimension];}
   int GetGridEndIndex(int Dimension) {return GridEndIndex[Dimension];}
   FLOAT GetGridLeftEdge(int Dimension) {return GridLeftEdge[Dimension];}
   FLOAT GetGridRightEdge(int Dimension) {return GridRightEdge[Dimension];}
   int GetGravitatingMassFieldDimension(int Dimension) {
     return GravitatingMassFieldDimension[Dimension];}
   FLOAT GetGravitatingMassFieldLeftEdge(int Dimension) {
     return GravitatingMassFieldLeftEdge[Dimension];}
   FLOAT GetGravitatingMassFieldCellSize() {
     return GravitatingMassFieldCellSize;}
   int GetGravityGhostZones() {return GravityGhostZones;}


// -------------------------------------------------------------------------
// Functions for accessing top-grid parallelism information
// (note: information only available/valid for this level)
//

/* Processor layout: get and set the number of procs in each 
   dim within the cartesian processor grid
   (1-based, i.e. {1 1 1} defines a single-processor layout) */ 
   int GetProcessorLayout(int Dimension) {return ProcLayout[Dimension];}
   void SetProcessorLayout(int Dimension, int procs) {
     if (Dimension < 0 || Dimension > MAX_DIMENSION)
       fprintf(stderr,"SetProcessorLayout: invalid Dimension.\n");
     else
       if (procs > 0)  ProcLayout[Dimension] = procs; 
       else fprintf(stderr,"SetProcessorLayout: invalid procs value.\n");
   }

/* Processor location: get and set the location of this grid's proc
   within the cartesian processor grid defined in ProcLayout
   (0-based, i.e. {0 0 0} defines the 1st proc in each dimension) */
   int GetProcessorLocation(int Dimension) {return ProcLocation[Dimension];}
   void SetProcessorLocation(int Dimension, int location) {
     if (Dimension < 0 || Dimension > MAX_DIMENSION)
       fprintf(stderr,"SetProcessorLocation: invalid Dimension.\n");
     else
       if (location >= 0)  ProcLocation[Dimension] = location; 
       else fprintf(stderr,"SetProcessorLocation: invalid location.\n");     
   }

/* Processor neighbors: get and set the grid IDs (not MPI process IDs) of this
   grid's neighbors within the cartesian processor grid defined in ProcLayout. 
     Get... returns the {left=0,right=1} neighbor grid ID in a given dim
     Set... provides access to set neighbor information into the grid */
   int GetProcessorNeighbors(int Dimension, int LR) {
     return ProcNeighbors[Dimension][LR];}
   void SetProcessorNeighbors(int Dimension, int LR, int NBid) { 
     if (Dimension < 0 || Dimension > MAX_DIMENSION)
       fprintf(stderr,"SetProcessorNeighbors: invalid Dimension.\n");
     else
       if (LR < 0 || LR > 1) 
	 fprintf(stderr,"SetProcessorNeighbors: invalid neighbor.\n");    
       else
	 if (NBid >= 0)  ProcNeighbors[Dimension][LR] = NBid; 
	 else fprintf(stderr,"SetProcessorNeighbors: invalid grid ID.\n");    
   }

// -------------------------------------------------------------------------
// Functions for use with particles.
//

/* Particles: Deposit particles in the specified field (DepositField) of the
              TargetGrid at the given time. */

   int DepositParticlePositions(grid *TargetGrid, FLOAT DepositTime, 
				int DepositField);

/* Particles: add overlapping ParticleMassField to Target's 
   GravitatingMassField. */

   int AddOverlappingParticleMassField(grid *TargetGrid, 
				       FLOAT EdgeOffset[MAX_DIMENSION]);

/* Particles: Apply particle acceleration to velocity for particles in this 
              grid
    (for step #9) */

   int UpdateParticleVelocity(float TimeStep);

/* Particles: Update particle positions (push)
    (for step #13) */

   int UpdateParticlePosition(float TimeStep);

/* Particles: Move particles from TargetGrid to this grid. */

   int MoveAllParticles(int NumberOfGrids, grid* TargetGrids[]);

/* Particles: Move particles that lie within this grid from the TargetGrid
              to this grid. */

//   int MoveSubgridParticles(grid *TargetGrid);


   int MoveSubgridParticles(grid *TargetGrid,
                            int *Counter,
                            int *Number,
                            int *Type,
                            float *Mass,
                            FLOAT *Position[],
                            float *Velocity[],
                            float *Attribute[]);


/* Particles: same as above, but a version that is much more efficient. */

   int MoveSubgridParticlesFast(int NumberOfSubgrids, grid *ToGrids[],
				int AllLocal);

/* Particles: Clean up moved particles (from MoveSubgridParticles). */

   int CleanUpMovedParticles();

/* Particles: delete accleration fields. */

   void DeleteParticleAcceleration() {
     for (int dim = 0; dim < GridRank+ComputePotential; dim++) {
       delete ParticleAcceleration[dim];
       ParticleAcceleration[dim] = NULL;
     }
   };

/* Particles & Gravity: Delete GravitatingMassField. */

   void DeleteGravitatingMassFieldParticles() {
     delete GravitatingMassFieldParticles; 
     GravitatingMassFieldParticles = NULL;
     GravitatingMassFieldParticlesCellSize = FLOAT_UNDEFINED;
   };

/* Particles: return number of particles. */

   int ReturnNumberOfParticles() {return NumberOfParticles;};

/* Particles: set number of particles. */

   void SetNumberOfParticles(int num) {NumberOfParticles = num;};

/* Particles: delete particle fields and set null. */

   void DeleteParticles() {
     delete ParticleMass;
     delete ParticleNumber;
     delete ParticleType;
     ParticleMass = NULL;
     ParticleNumber = NULL;
     ParticleType = NULL;
     for (int dim = 0; dim < GridRank; dim++) {
       delete ParticlePosition[dim];
       delete ParticleVelocity[dim];
       ParticlePosition[dim] = NULL;
       ParticleVelocity[dim] = NULL;
     }
     for (int i = 0; i < NumberOfParticleAttributes; i++) {
       delete ParticleAttribute[i];
       ParticleAttribute[i] = NULL;
     }   
   };

/* Particles: allocate new particle fields. */

   void AllocateNewParticles(int NumberOfNewParticles) {
     ParticleMass = new float[NumberOfNewParticles];
     ParticleNumber = new int[NumberOfNewParticles];
     ParticleType = new int[NumberOfNewParticles];
     for (int dim = 0; dim < GridRank; dim++) {
       ParticlePosition[dim] = new FLOAT[NumberOfNewParticles];
       ParticleVelocity[dim] = new float[NumberOfNewParticles];
     }
     for (int i = 0; i < NumberOfParticleAttributes; i++)
       ParticleAttribute[i] = new float[NumberOfNewParticles];
   };

/* Particles: Copy pointers passed into into grid. */

   void SetParticlePointers(float *Mass, int *Number, int *Type,
                            FLOAT *Position[], 
			    float *Velocity[], float *Attribute[]) {
    ParticleMass   = Mass;
    ParticleNumber = Number;
    ParticleType   = Type;
    for (int dim = 0; dim < GridRank; dim++) {
      ParticlePosition[dim] = Position[dim];
      ParticleVelocity[dim] = Velocity[dim];
    }
    for (int i = 0; i < NumberOfParticleAttributes; i++)
      ParticleAttribute[i] = Attribute[i];
   };

/* Particles: Set new star particle index. */

   void SetNewParticleIndex(int &NumberCount, int BaseNumber) {
    for (int n = 0; n < NumberOfParticles; n++)
      if (ParticleNumber[n] == INT_UNDEFINED)
	ParticleNumber[n] = BaseNumber + NumberCount++;
   };

/* Particles: Add given number to particle index. */

   void AddToParticleNumber(int *Count) {
     if (MyProcessorNumber == ProcessorNumber)
       for (int n = 0; n < NumberOfParticles; n++)
	 ParticleNumber[n] += *Count;
     *Count += NumberOfParticles;
   }

/* Particles: return particle type (1 - dm, 2 - star).   Note that this
   has now been superceded by a real particle type field. */

   int ReturnParticleType(int index) {
     if (NumberOfParticleAttributes > 0 && StarParticleCreation > 0)
       if (ParticleAttribute[0][index] > 0)
         return PARTICLE_TYPE_STAR;
     return PARTICLE_TYPE_DARK_MATTER;
   }

/* Particles: sort particle data in ascending order by number (id). */

void SortParticlesByNumber();

// -------------------------------------------------------------------------
// Communication functions
//

/* Set grid's processor number. */

  void SetProcessorNumber(int Proc) {
    ProcessorNumber = Proc;
  };

/* Return grid's processor number. */

  int ReturnProcessorNumber() {
    return ProcessorNumber;
  }

/* Send a region from a real grid to a 'fake' grid on another processor. */

  int CommunicationSendRegion(grid *ToGrid, int ToProcessor, int SendField, 
			     int NewOrOld, int RegionStart[], int RegionDim[]);

/* Send a region from a 'fake' grid to a real grid on another processor. */

  int CommunicationReceiveRegion(grid *ToGrid, int ToProcessor, 
				 int SendField, int NewOrOld, 
				 int RegionStart[], int RegionDim[],
				 int IncludeBoundary);

/* Move a grid from one processor to another. */

  int CommunicationMoveGrid(int ToProcessor);

/* Send particles from one grid to another. */

  int CommunicationSendParticles(grid *ToGrid, int ToProcessor, 
				int FromStart, int FromNumber, int ToStart);

/* Transfer particle amount level 0 grids. */

  int CommunicationTransferParticles(grid* Grids[], int NumberOfSubgrids, 
		 int ToGrid[6], int NumberToMove[6], 
		 float_int *ParticleData[6], int CopyDirection);

// -------------------------------------------------------------------------
// Helper functions (should be made private)
//

/* Baryons: find certain commonly used variables from the list of fields. */

  int IdentifyPhysicalQuantities(int &DensNum, int &GENum,   int &Vel1Num, 
				 int &Vel2Num, int &Vel3Num, int &TENum);

/* Identify Multi-species fields. */

  int IdentifySpeciesFields(int &DeNum, int &HINum, int &HIINum, 
			    int &HeINum, int &HeIINum, int &HeIIINum,
			    int &HMNum, int &H2INum, int &H2IINum,
                            int &DINum, int &DIINum, int &HDINum);


// AccelerationHack

  int AccelerationHack;

  //These should be moved later.
  //Used for boundary condition set of AccelerationField.
  int ActualNumberOfBaryonFields;
  int AttachAcceleration();
  int DetachAcceleration();

  int    ActualFieldType[MAX_NUMBER_OF_BARYON_FIELDS];
  float *ActualBaryonField[MAX_NUMBER_OF_BARYON_FIELDS];
  float *ActualOldBaryonField[MAX_NUMBER_OF_BARYON_FIELDS];
  float *OldAccelerationField[3];

// -------------------------------------------------------------------------
// Functions for Specific problems (usually problem generator functions).
//

/* Radiation Hydrodynamics test problem: initialize grid (SUCCESS or FAIL) */
  int RadHydroConstTestInitializeGrid(int NumChemicals,
				      float DensityConstant, 
				      float V0Constant,
				      float V1Constant, 
				      float V2Constant,
				      float TemperatureConstant, 
				      float EgConstant,
				      float HydrogenMassFraction, 
				      float InitialFractionHII, 
				      float InitialFractionHeII, 
				      float InitialFractionHeIII, 
				      float OmegaBaryonNow);

/* Radiation Hydrodynamics test problem: stream test (SUCCESS or FAIL) */
  int RadHydroStreamTestInitializeGrid(
				      float DensityConstant, 
				      float EgConstant);

/* Protostellar Collapse problem: initialize grid (returns SUCCESS or FAIL) */
  int ProtostellarCollapseInitializeGrid(float CoreDensity,
					 float CoreEnergy,
					 float CoreRadius,
					 float AngularVelocity);

/* ShockTube problem: initialize grid (returns SUCCESS or FAIL) */

  int ShockTubeInitializeGrid(int Direction, float Boundary, float Density[],
			      float Pressure[], float Velocity[]);

/* Initialize for a uniform grid (returns SUCCESS or FAIL) */

  int InitializeUniformGrid(float UniformDensity, float UniformTotalEnergy,
			    float UniformGasEnergy, float UniformVelocity[]);


/* Initialize a grid for the Double Mach reflection problem. */

  int DoubleMachInitializeGrid(float d0, float e0, float u0,float v0,float w0);


/* Initialize a grid for Implosion test problem */

  int ImplosionInitializeGrid(float ImplosionDiamondDensity,
			      float ImplosionDiamondTotalEnergy);


/* Initialize a grid for Sedov Explosion */

  int SedovBlastInitializeGrid(float SedovBlastInitialRadius,
                               float SedovBlastInnerTotalEnergy);

  int SedovBlastInitializeGrid3D(char * fname);

  /* Initialize a grid for the KH instability problem. */

  int KHInitializeGrid(float KHInnerDensity,
                       float KHInnerInternalEnergy,
                       float KHOuterInternalEnergy,
                       float KHPerturbationAmplitude,
                       float KHInnerVx, float KHOuterVx);

  /* Initialize a grid and set boundary for the 2D/3D Noh problem. */

  int NohInitializeGrid(float d0, float p0, float u0);
  int ComputeExternalNohBoundary();


/* Zeldovich Pancake: initial grid (returns SUCCESS or FAIL). */

  int ZeldovichPancakeInitializeGrid(int   ZeldovichPancakeDirection,
				     float ZeldovichPancakeCentralOffset,
				     float ZeldovichPancakeOmegaBaryonNow,
				     float ZeldovichPancakeOmegaCDMNow,
				     float ZeldovichPancakeCollapseRedshift,
				     float ZeldovichPancakeInitialTemperature);

/* 1D Pressureless Collapse: initialize grid. */

  int PressurelessCollapseInitializeGrid(int PressurelessCollapseDirection,
				    float PressurelessCollapseInitialDensity,
				      int PressurelessCollapseNumberOfCells);

/* Gravity Test: initialize grid. */

  int TestGravityInitializeGrid(float CentralDensity, 
				int NumberOfNewParticles, int UseBaryons);

/* Gravity Test: check results. */

  int TestGravityCheckResults(FILE *fptr, grid *TopGrid);

/* Gravity Test Motion: initialize grid. */

  int TestGravityMotionInitializeGrid(float InitialVelocity);

/* Gravity Test (Sphere): initialize grid. */

  int TestGravitySphereInitializeGrid(float InteriorDensity, 
				      float ExteriorDensity,
				      float SphereRadius, 
				      int SphereType, int UseBaryons,
				      FLOAT SphereCenter[]);

/* Gravity Test (Sphere): check results. */

  int TestGravitySphereCheckResults(FILE *fptr);

/* Spherical Infall Test: initialize grid. */

  int SphericalInfallInitializeGrid(float InitialPerturbation, int UseBaryons,
				    float SphericalInfallOmegaBaryonNow,
				    float SphericalInfallOmegaCDMNow,
				    int SubgridIsStatic);


/* Spherical Infall Test: get the profile from the center. */

  int SphericalInfallGetProfile(int level, int ReportLevel);

/* GravityEquilibriumTest: initialize grid. */

  int GravityEquilibriumTestInitializeGrid(float ScaleHeight);

/* CollapseTest: initialize grid. */

#define MAX_SPHERES 10

  int CollapseTestInitializeGrid(int NumberOfSpheres,
			     FLOAT SphereRadius[MAX_SPHERES],
			     FLOAT SphereCoreRadius[MAX_SPHERES],
			     float SphereDensity[MAX_SPHERES],
			     float SphereTemperature[MAX_SPHERES],
			     FLOAT SpherePosition[MAX_SPHERES][MAX_DIMENSION],
			     float SphereVelocity[MAX_SPHERES][MAX_DIMENSION],
			     int   SphereType[MAX_SPHERES],
			     int   SphereUseParticles, 
			     float UniformVelocity[MAX_DIMENSION],
			     int   SphereUseColour,
			     float InitialTemperature, int level);

/* CosmologySimulation: initialize grid. */

  int CosmologySimulationInitializeGrid(
			  float CosmologySimulationOmegaBaryonNow,
			  float CosmologySimulationOmegaCDMNow,
			  float CosmologySimulationInitialTemperature,
			  char *CosmologySimulationDensityName,
			  char *CosmologySimulationTotalEnergyName,
			  char *CosmologySimulationGasEnergyName,
			  char *CosmologySimulationVelocityNames[],
			  char *CosmologySimulationParticlePositionName,
			  char *CosmologySimulationParticleVelocityName,
			  char *CosmologySimulationParticleMassName,
			  char *CosmologySimulationParticleTypeName,
			  int   CosmologySimulationSubgridsAreStatic,
			  int   TotalRefinement,
			  float CosmologySimulationInitialFrctionHII,
			  float CosmologySimulationInitialFrctionHeII,
			  float CosmologySimulationInitialFrctionHeIII,
			  float CosmologySimulationInitialFrctionHM,
			  float CosmologySimulationInitialFrctionH2I,
			  float CosmologySimulationInitialFrctionH2II,
			  int   CosmologySimulationUseMetallicityField,
			  int  &CurrentNumberOfParticles);

/* CosmologySimulation: initialize partitioned nested grids. */

  int NestedCosmologySimulationInitializeGrid(
			  float CosmologySimulationOmegaBaryonNow,
			  float CosmologySimulationOmegaCDMNow,
			  float CosmologySimulationInitialTemperature,
			  char *CosmologySimulationDensityName,
			  char *CosmologySimulationTotalEnergyName,
			  char *CosmologySimulationGasEnergyName,
			  char *CosmologySimulationVelocityNames[],
			  char *CosmologySimulationParticlePositionName,
			  char *CosmologySimulationParticleVelocityName,
			  char *CosmologySimulationParticleMassName,
			  char* CosmologySimulationParticleTypeName,
			  int   CosmologySimulationSubgridsAreStatic,
			  int   TotalRefinement,
			  float CosmologySimulationInitialFrctionHII,
			  float CosmologySimulationInitialFrctionHeII,
			  float CosmologySimulationInitialFrctionHeIII,
			  float CosmologySimulationInitialFrctionHM,
			  float CosmologySimulationInitialFrctionH2I,
			  float CosmologySimulationInitialFrctionH2II,
			  int   CosmologySimulationUseMetallicityField,
			  int  &CurrentNumberOfParticles);

/* Supernova restart initialize grid. */

  int SupernovaRestartInitialize(float EjectaDensity, float EjectaRadius,
				 float EjectaThermalEnergy, 
				 FLOAT EjectaCenter[3], int ColourField,
				 int *NumberOfCellsSet);

/* Tricks for Random Forcing. */

  int ReturnNumberOfBaryonFields(){return NumberOfBaryonFields;};
  int SetNumberOfBaryonFields(int &number)
    {NumberOfBaryonFields = number; return 0;};
  int AppendForcingToBaryonFields();
  int DetachForcingFromBaryonFields();
  int RemoveForcingFromBaryonFields();
  int AddRandomForcing(float * norm, float dtTopGrid);
  int PrepareRandomForcingNormalization(float * GlobVal, int GlobNum);
  inline int ReadRandomForcingFields(FILE *main_file_pointer);

  inline bool isLocal () {return MyProcessorNumber == ProcessorNumber; };
 private:
  int ReadRandomForcingFieldsHDF4(FILE *main_file_pointer);
  int ReadRandomForcingFieldsHDF5(FILE *main_file_pointer);
 public:

/* TurbulenceSimulation: initialize grid. */

#define TURBULENCE_INIT_PARAMETERS_DECL \
     float TurbulenceSimulationInitialDensity, \
     float TurbulenceSimulationInitialTemperature, \
     char *TurbulenceSimulationDensityName, \
     char *TurbulenceSimulationTotalEnergyName, \
     char *TurbulenceSimulationGasEnergyName, \
     char *TurbulenceSimulationVelocityNames[], \
     char *TurbulenceSimulationRandomForcingNames[], \
     int   TurbulenceSimulationSubgridsAreStatic, \
     int   TotalRefinement


#define TURBULENCE_INIT_PARAMETERS \
     TurbulenceSimulationInitialDensity, \
     TurbulenceSimulationInitialTemperature, \
     TurbulenceSimulationDensityName, \
     TurbulenceSimulationTotalEnergyName, \
     TurbulenceSimulationGasEnergyName, \
     TurbulenceSimulationVelocityNames, \
     TurbulenceSimulationRandomForcingNames, \
     TurbulenceSimulationSubgridsAreStatic, \
     TotalRefinement


inline int TurbulenceSimulationInitializeGrid(TURBULENCE_INIT_PARAMETERS_DECL);

  // The following are private since they should only be called by
  // TurbulenceSimulationInitializeGrid()

 private:
  int TurbulenceSimulationInitializeGridHDF4(TURBULENCE_INIT_PARAMETERS_DECL);
  int TurbulenceSimulationInitializeGridHDF5(TURBULENCE_INIT_PARAMETERS_DECL);
 public:


/* Comoving coordinate expansion terms. */

  int ComovingExpansionTerms();

/* Adjust the gravity source terms for comoving coordinates. */

  int ComovingGravitySourceTerm();

/* Star Particle handler routine. */

  int StarParticleHandler(int level);

/* Apply a time-action to a grid. */

  int ApplyTimeAction(int Type, float Parameter);

/* Routine to set the tracer particle velocities from the grid velocity. */

  int TracerParticleSetVelocity();

/* Creates tracer particles in this grid. */

  int TracerParticleCreateParticles(FLOAT LeftEdge[], FLOAT RightEdge[],
                                    FLOAT Spacing, int &TotalParticleCount);

  /* Includes for analysis tools */

#ifdef ANALYSIS_TOOLS
#   include "Grid_AnalyzeClusters.h"
#endif

};

#define USE_HDF5

/***********************************************************************/
inline int grid::ReadRandomForcingFields (FILE *main_file_pointer)
{
#if defined (USE_HDF4)
  return ReadRandomForcingFieldsHDF4 (main_file_pointer);
#elif defined (USE_HDF5)
  return ReadRandomForcingFieldsHDF5 (main_file_pointer);
#else
  WARNING_MESSAGE;
  return FAIL; // Fail if neither USE_HDF4 nor USE_HDF5 defined
#endif
}

/***********************************************************************/

inline int 
grid::TurbulenceSimulationInitializeGrid (TURBULENCE_INIT_PARAMETERS_DECL)
{
#if defined (USE_HDF4)
  return TurbulenceSimulationInitializeGridHDF4 (TURBULENCE_INIT_PARAMETERS);
#elif defined (USE_HDF5)
  return TurbulenceSimulationInitializeGridHDF5 (TURBULENCE_INIT_PARAMETERS);
#else
  WARNING_MESSAGE;
  return FAIL; // Fail if neither USE_HDF4 nor USE_HDF5 defined
#endif
}



#endif
