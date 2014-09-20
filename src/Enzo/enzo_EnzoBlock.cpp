// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar  3 23:02:02 PST 2011
/// @brief    Implementation of the EnzoBlock class

#include "cello.hpp"

#include "enzo.hpp"

//======================================================================

/// Boundary

int  EnzoBlock::BoundaryRank;
int  EnzoBlock::BoundaryDimension[MAX_DIMENSION];
int  EnzoBlock::BoundaryFieldType[MAX_NUMBER_OF_BARYON_FIELDS];
bc_enum *EnzoBlock::BoundaryType[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2];
enzo_float *EnzoBlock::BoundaryValue[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2];

int EnzoBlock::ComovingCoordinates;
int EnzoBlock::UseMinimumPressureSupport;
enzo_float EnzoBlock::MinimumPressureSupportParameter;
enzo_float EnzoBlock::ComovingBoxSize;
enzo_float EnzoBlock::HubbleConstantNow;
enzo_float EnzoBlock::OmegaMatterNow;
enzo_float EnzoBlock::OmegaLambdaNow;
enzo_float EnzoBlock::MaxExpansionRate;

// Chemistry

int EnzoBlock::MultiSpecies;

// Gravity

int EnzoBlock::GravityOn;

// Physics

int EnzoBlock::PressureFree;
enzo_float EnzoBlock::Gamma;
enzo_float EnzoBlock::GravitationalConstant;

// Problem-specific

int EnzoBlock::ProblemType;

// Method PPM

int EnzoBlock::PPMFlatteningParameter;
int EnzoBlock::PPMDiffusionParameter;
int EnzoBlock::PPMSteepeningParameter;

// Numerics

int EnzoBlock::DualEnergyFormalism;
enzo_float EnzoBlock::DualEnergyFormalismEta1;
enzo_float EnzoBlock::DualEnergyFormalismEta2;

enzo_float EnzoBlock::pressure_floor;
enzo_float EnzoBlock::density_floor;
enzo_float EnzoBlock::number_density_floor;
enzo_float EnzoBlock::temperature_floor;

enzo_float EnzoBlock::CourantSafetyNumber;
enzo_float EnzoBlock::InitialRedshift;
enzo_float EnzoBlock::InitialTimeInCodeUnits;

// Domain

enzo_float EnzoBlock::DomainLeftEdge [MAX_DIMENSION];
enzo_float EnzoBlock::DomainRightEdge[MAX_DIMENSION];

// PPM

int EnzoBlock::GridRank;

int EnzoBlock::ghost_depth[MAX_DIMENSION];

// Fields

int EnzoBlock::NumberOfBaryonFields;

//----------------------------------------------------------------------

// STATIC
void EnzoBlock::initialize(EnzoConfig * enzo_config,
			   FieldDescr * field_descr)
{

  GridRank = 0;
  NumberOfBaryonFields = 0;

  int i,j,k;

  for (i=0; i<MAX_DIMENSION; i++) {
    DomainLeftEdge [i] = 0;
    DomainRightEdge[i] = 0;
    ghost_depth[i] = 0;
    for (j=0; j<MAX_NUMBER_OF_BARYON_FIELDS; j++) {
      for (k=0; k<2; k++) {
	BoundaryDimension[i]=0;
	BoundaryFieldType[j] = 0;
	BoundaryType[j][i][k] = 0;
	BoundaryValue[j][i][k]= 0;
      }
    }
  }

  ComovingCoordinates = enzo_config->physics_cosmology;
  Gamma               = enzo_config->field_gamma;

  GridRank            = enzo_config->mesh_root_rank;

  BoundaryRank = GridRank;

  // Chemistry parameters

  MultiSpecies = 0;    // 0:0 1:6 2:9 3:12

  // Gravity parameters

  GravityOn                       = 0;    // Whether gravity is included
  GravitationalConstant           = 1.0;  // used only in SetMinimumSupport()

  //Problem specific parameter

  ProblemType = 0;

  // PPM parameters

  InitialRedshift   = enzo_config->physics_cosmology_initial_redshift;
  HubbleConstantNow = enzo_config->physics_cosmology_hubble_constant_now;
  OmegaLambdaNow    = enzo_config->physics_cosmology_omega_lamda_now;
  OmegaMatterNow    = enzo_config->physics_cosmology_omega_matter_now;
  MaxExpansionRate  = enzo_config->physics_cosmology_max_expansion_rate;
  ComovingBoxSize   = enzo_config->physics_cosmology_comoving_box_size;

  PressureFree              = enzo_config->ppm_pressure_free;
  UseMinimumPressureSupport = enzo_config->ppm_use_minimum_pressure_support;
  MinimumPressureSupportParameter = 
    enzo_config->ppm_minimum_pressure_support_parameter;
  PPMFlatteningParameter    = enzo_config->ppm_flattening;
  PPMDiffusionParameter     = enzo_config->ppm_diffusion;
  PPMSteepeningParameter    = enzo_config->ppm_steepening;
  pressure_floor            = enzo_config->ppm_pressure_floor;
  density_floor             = enzo_config->ppm_density_floor;
  temperature_floor         = enzo_config->ppm_temperature_floor;
  number_density_floor      = enzo_config->ppm_number_density_floor;
  DualEnergyFormalism       = enzo_config->ppm_dual_energy;
  DualEnergyFormalismEta1   = enzo_config->ppm_dual_energy_eta_1;
  DualEnergyFormalismEta2   = enzo_config->ppm_dual_energy_eta_2;

  int gx = enzo_config->field_ghosts[0];
  int gy = enzo_config->field_ghosts[1];
  int gz = enzo_config->field_ghosts[2];

  if (GridRank < 1) gx = 0;
  if (GridRank < 2) gy = 0;
  if (GridRank < 3) gz = 0;

  ghost_depth[0] = gx;
  ghost_depth[1] = gy;
  ghost_depth[2] = gz;

  NumberOfBaryonFields = enzo_config->field_list.size();

  // Check NumberOfBaryonFields

  if (NumberOfBaryonFields == 0) {
    ERROR ("EnzoBlock::initialize",
	   "List parameter 'Field fields' must have length greater than zero");
  } else if (NumberOfBaryonFields > MAX_NUMBER_OF_BARYON_FIELDS) {
    ERROR2 ("EnzoBlock::initialize",
	    "MAX_NUMBER_OF_BARYON_FIELDS = %d is too small for %d fields",
	    MAX_NUMBER_OF_BARYON_FIELDS,NumberOfBaryonFields );
  }

  DomainLeftEdge [0] = enzo_config->domain_lower[0];
  DomainLeftEdge [1] = enzo_config->domain_lower[1];
  DomainLeftEdge [2] = enzo_config->domain_lower[2];

  DomainRightEdge[0] = enzo_config->domain_upper[0];
  DomainRightEdge[1] = enzo_config->domain_upper[1];
  DomainRightEdge[2] = enzo_config->domain_upper[2];

  CourantSafetyNumber = enzo_config->field_courant;

  double time  = enzo_config->initial_time;

  InitialTimeInCodeUnits = time;

} // void initialize()

//----------------------------------------------------------------------

EnzoBlock::EnzoBlock
(
 Index index,
 int nx, int ny, int nz,
 int num_field_blocks,
 int count_adapt,
 int cycle, double time, double dt,
 int narray, char * array, int op_array,
 int num_face_level, int * face_level,
 bool testing
) throw()
  : CommBlock 
    (
     index,
     nx,ny,nz,
     num_field_blocks,
     count_adapt,
     cycle, time, dt,
     narray,  array, op_array,
     num_face_level, face_level,
     testing),
    Time_(time),
    CycleNumber(cycle),
    OldTime(0),
    dt(dt),
    SubgridFluxes(0)
{
  initialize_enzo_();
  initialize();
}

//----------------------------------------------------------------------

void EnzoBlock::initialize_enzo_()
{
  for (int i=0; i<MAX_DIMENSION; i++) {
    AccelerationField[i] = 0;
    GridLeftEdge[i] = 0;
    GridDimension[i] = 0;
    GridStartIndex[i] = 0;
    GridEndIndex[i] = 0;
    CellWidth[i] = 0;
  }

}

//----------------------------------------------------------------------

EnzoBlock::~EnzoBlock() throw ()
{
}

//----------------------------------------------------------------------

void EnzoBlock::pup(PUP::er &p)
{ 

  TRACEPUP;
  TRACE ("BEGIN EnzoBlock::pup()");

  CommBlock::pup(p);

  p | Time_;

  p | CycleNumber;
  p | OldTime;
  p | dt;
  static bool warn0 = true;
  if (warn0) {
    WARNING("EnzoBlock::pup()", "skipping AccelerationField_ (not used)");
    warn0=false;
  }
  p | Time_;

  static bool warn1 = true;
  if (warn1) {
    WARNING("EnzoBlock::pup()", "skipping AccelerationField_ (not used)");
    warn1=false;
  }
  static bool warn2 = true;
  if (warn2) {
    WARNING("EnzoBlock::pup()", "skipping SubgridFluxes (not used)");
    warn2=false;
  }

  PUParray(p,GridLeftEdge,MAX_DIMENSION); 
  PUParray(p,GridDimension,MAX_DIMENSION); 
  PUParray(p,GridStartIndex,MAX_DIMENSION); 
  PUParray(p,GridEndIndex,MAX_DIMENSION); 
  PUParray(p,CellWidth,MAX_DIMENSION);

  static bool warn3 = true;
  if (warn3) {
    warn3=false;
    WARNING("EnzoBlock::pup()", "skipping OldBaryonField[] [not used]");
  }

  PUParray(p,method_turbulence_data,9);

  TRACE ("END EnzoBlock::pup()");

}

//======================================================================

void EnzoBlock::write(FILE * fp) throw ()
{
  fprintf (fp,"EnzoBlock: ComovingCoordinates %d\n",
	   ComovingCoordinates);
  fprintf (fp,"EnzoBlock: UseMinimumPressureSupport %d\n",
	   UseMinimumPressureSupport);
  fprintf (fp,"EnzoBlock: MinimumPressureSupportParameter %g\n",
	   MinimumPressureSupportParameter);
  fprintf (fp,"EnzoBlock: ComovingBoxSize %g\n",
	   ComovingBoxSize);
  fprintf (fp,"EnzoBlock: HubbleConstantNow %g\n",
	   HubbleConstantNow);
  fprintf (fp,"EnzoBlock: OmegaLambdaNow %g\n",
	   OmegaLambdaNow);
  fprintf (fp,"EnzoBlock: OmegaMatterNow %g\n",
	   OmegaMatterNow);
  fprintf (fp,"EnzoBlock: MaxExpansionRate %g\n",
	   MaxExpansionRate);

  // Chemistry

  fprintf (fp,"EnzoBlock: MultiSpecies %d\n",
	   MultiSpecies);

  // Gravity

  fprintf (fp,"EnzoBlock: GravityOn %d\n",
	   GravityOn);
  //  fprintf (fp,"EnzoBlock: *AccelerationField %g\n",
  //           *AccelerationField[MAX_DIMENSION)];

  // Physics

  fprintf (fp,"EnzoBlock: PressureFree %d\n",
	   PressureFree);
  fprintf (fp,"EnzoBlock: Gamma %g\n",
	   Gamma);
  fprintf (fp,"EnzoBlock: GravitationalConstant %g\n",
	   GravitationalConstant);

  // Problem-specific

  fprintf (fp,"EnzoBlock: ProblemType %d\n",
	   ProblemType);

  // Method PPM

  fprintf (fp,"EnzoBlock: PPMFlatteningParameter %d\n",
	   PPMFlatteningParameter);
  fprintf (fp,"EnzoBlock: PPMDiffusionParameter %d\n",
	   PPMDiffusionParameter);
  fprintf (fp,"EnzoBlock: PPMSteepeningParameter %d\n",
	   PPMSteepeningParameter);

  // Numerics

  fprintf (fp,"EnzoBlock: DualEnergyFormalism %d\n",
	   DualEnergyFormalism);
  fprintf (fp,"EnzoBlock: DualEnergyFormalismEta1 %g\n",
	   DualEnergyFormalismEta1);
  fprintf (fp,"EnzoBlock: DualEnergyFormalismEta2 %g\n",
	   DualEnergyFormalismEta2);
  fprintf (fp,"EnzoBlock: pressure_floor %g\n",
	   pressure_floor);
  fprintf (fp,"EnzoBlock: density_density_floor %g\n",
	   density_floor);
  fprintf (fp,"EnzoBlock: number_density_floor %g\n",
	   number_density_floor);
  fprintf (fp,"EnzoBlock: temperature_floor %g\n",
	   temperature_floor);

  fprintf (fp,"EnzoBlock: CourantSafetyNumber %g\n",
	   CourantSafetyNumber);
  fprintf (fp,"EnzoBlock: InitialRedshift %g\n",
	   InitialRedshift);
  fprintf (fp,"EnzoBlock: InitialTimeInCodeUnits %g\n",
	   InitialTimeInCodeUnits);
  fprintf (fp,"EnzoBlock: Time %g\n",
	   Time());
  fprintf (fp,"EnzoBlock: OldTime %g\n",
	   OldTime);

  // Domain

  fprintf (fp,"EnzoBlock: DomainLeftEdge %g %g %g\n",
	   DomainLeftEdge [0],DomainLeftEdge [0],DomainLeftEdge [0]);
  fprintf (fp,"EnzoBlock: DomainRightEdge %g %g %g\n",
	   DomainRightEdge[0],DomainRightEdge[1],DomainRightEdge[2]);

  // Fields

  // Grid

  fprintf (fp,"EnzoBlock: GridRank %d\n",    GridRank);
  fprintf (fp,"EnzoBlock: GridDimension %d %d %d\n",
	   GridDimension[0],GridDimension[1],GridDimension[2]);
  fprintf (fp,"EnzoBlock: GridStartIndex %d %d %d\n",
	   GridStartIndex[0],GridStartIndex[1],GridStartIndex[2]);
  fprintf (fp,"EnzoBlock: GridEndIndex %d %d %d\n",
	   GridEndIndex[0],GridEndIndex[1],GridEndIndex[2]);
  fprintf (fp,"EnzoBlock: GridLeftEdge %g %g %g\n",
	   GridLeftEdge[0],GridLeftEdge[1],GridLeftEdge[2]);

  fprintf (fp,"EnzoBlock: CellWidth %g %g %g\n", 
	   CellWidth[0], CellWidth[1], CellWidth[2] );

  fprintf (fp,"EnzoBlock: ghost %d %d %d\n",
	   ghost_depth[0],ghost_depth[1],ghost_depth[2]);


  fprintf (fp,"EnzoBlock: NumberOfBaryonFields %d\n",
	   NumberOfBaryonFields);

  fprintf (fp,"EnzoBlock: BoundaryRank %d\n", BoundaryRank);
  fprintf (fp,"EnzoBlock: BoundaryDimension %d %d %d\n",
	   BoundaryDimension[0],BoundaryDimension[1],BoundaryDimension[2]);

  for (int i=0; i<NumberOfBaryonFields; i++) {

    fprintf (fp,"EnzoBlock: BoundaryFieldType[%d] %d\n", 
	     i, BoundaryFieldType[i]);

    fprintf (fp,"EnzoBlock: BoundaryType[%d] %p %p %p %p %p %p\n", i, 
	     BoundaryType[i][0][0],
	     BoundaryType[i][0][1],
	     BoundaryType[i][1][0],
	     BoundaryType[i][1][1],
	     BoundaryType[i][2][0],
	     BoundaryType[i][2][1]);

    fprintf (fp,"EnzoBlock: BoundaryValue[%d] %p %p %p %p %p %p\n",i,
	     BoundaryValue[i][0][0],
	     BoundaryValue[i][0][1],
	     BoundaryValue[i][1][0],
	     BoundaryValue[i][1][1],
	     BoundaryValue[i][2][0],
	     BoundaryValue[i][2][1]);  
  }

  // problem

  fprintf (fp,"EnzoBlock: CycleNumber %d\n",   CycleNumber);
  fprintf (fp,"EnzoBlock: dt %g\n", dt);

  // fluxes

  fprintf (fp,"EnzoBlock: SubgridFluxes %p\n", SubgridFluxes);
  
}

//----------------------------------------------------------------------

void EnzoBlock::set_cycle (int cycle_start) throw ()
{
  CommBlock::set_cycle (cycle_start);
  TRACE2("%p EnzoBlock::set_cycle(%d)",this,cycle_start);

  CycleNumber = cycle_start;
}

//----------------------------------------------------------------------

void EnzoBlock::set_time (double time) throw ()
{
  CommBlock::set_time (time);

  //  if (cycle_ > 460) {
  //    char buffer[255];
  //    sprintf (buffer,"set_time %15.12f %15.12f",Time_,time);
  //    index_.print(buffer,-1,2,false,simulation());
  //  }

  //  if (! (Time_ == 0 || Time_ < time)) {
  if (! (Time_ == 0 || Time_ < time)) {
    index_.print("ERROR",-1,2,false,simulation());
    ASSERT2("EnzoBlock::set_time()",
   	    "set_time() may be called more than once or dt = 0.0\n"
   	    "Time_ = %15.8f time = %15.8f",
   	    Time_,time,
   	    Time_ == 0 || Time_ < time);
  }

  OldTime   = Time_;
  Time_     = time;

}

//----------------------------------------------------------------------

void EnzoBlock::set_dt (double dt_param) throw ()
{
  CommBlock::set_dt (dt_param);

  dt = dt_param;
}

//----------------------------------------------------------------------

void EnzoBlock::set_stop (bool stop) throw ()
{
  CommBlock::set_stop (stop);
}

//----------------------------------------------------------------------

void EnzoBlock::initialize () throw()
{
  TRACE ("EnzoBlock::initialize()\n");

  double xm,ym,zm;

  block()->lower(&xm,&ym,&zm);

  GridLeftEdge[0]  = xm;
  GridLeftEdge[1]  = ym;
  GridLeftEdge[2]  = zm;

  // Grid dimensions

  Field field = block()->field();

  int nx,ny,nz;
  field.size (&nx,&ny,&nz);

  int gx,gy,gz;

  gx = EnzoBlock::ghost_depth[0];
  gy = EnzoBlock::ghost_depth[1];
  gz = EnzoBlock::ghost_depth[2];

  GridDimension[0]  = nx + 2*gx;
  GridDimension[1]  = ny + 2*gy;
  GridDimension[2]  = nz + 2*gz;

  TRACE("Initializing GridStartIndex");

  GridStartIndex[0] = gx;
  GridStartIndex[1] = gy;
  GridStartIndex[2] = gz;

  GridEndIndex[0] = gx + nx - 1;
  GridEndIndex[1] = gy + ny - 1;
  GridEndIndex[2] = gz + nz - 1;

  // Initialize CellWidth

  double xp,yp,zp;
  block()->upper(&xp,&yp,&zp);
  double hx,hy,hz;
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);

  CellWidth[0] = hx;
  CellWidth[1] = hy;
  CellWidth[2] = hz;

  TRACE ("Exit  EnzoBlock::initialize()\n");
}

//----------------------------------------------------------------------

int EnzoBlock::CosmologyComputeExpansionFactor
(enzo_float time, enzo_float *a, enzo_float *dadt)
{
 
  /* Error check. */

  if (InitialTimeInCodeUnits == 0) {
    
    char error_message[ERROR_LENGTH];
    sprintf(error_message, "The cosmology parameters seem to be improperly set");
    ERROR("CosmologyComputeExpansionFactor",error_message);
  }
 
  *a = ENZO_FLOAT_UNDEFINED;
 
  /* Find Omega due to curvature. */
 
  enzo_float OmegaCurvatureNow = 1 - OmegaMatterNow - OmegaLambdaNow;
 
  /* Convert the time from code units to Time * H0 (c.f. CosmologyGetUnits). */
 
  enzo_float TimeUnits = 2.52e17/sqrt(OmegaMatterNow)/HubbleConstantNow/
                    pow(1 + InitialRedshift,enzo_float(1.5));
 
  enzo_float TimeHubble0 = time * TimeUnits * (HubbleConstantNow*3.24e-18);
 
  /* 1) For a flat universe with OmegaMatterNow = 1, it's easy. */
 
  if (fabs(OmegaMatterNow-1) < OMEGA_TOLERANCE &&
      OmegaLambdaNow < OMEGA_TOLERANCE)
    *a    = pow(time/InitialTimeInCodeUnits, enzo_float(2.0/3.0));
 
#define INVERSE_HYPERBOLIC_EXISTS
 
#ifdef INVERSE_HYPERBOLIC_EXISTS
 
 
  /* 2) For OmegaMatterNow < 1 and OmegaLambdaNow == 0 see
        Peebles 1993, eq. 13-3, 13-10.
	Actually, this is a little tricky since we must solve an equation
	of the form eta - sinh(eta) + x = 0..*/
 
  if (OmegaMatterNow < 1 && OmegaLambdaNow < OMEGA_TOLERANCE) {
 
    enzo_float eta, eta_old, x;
    int i;

    x = 2*TimeHubble0*pow(1.0 - OmegaMatterNow, 1.5) / OmegaMatterNow;
 
    /* Compute eta in a three step process, first from a third-order
       Taylor expansion of the formula above, then use that in a fifth-order
       approximation.  Then finally, iterate on the formula itself, solving for
       eta.  This works well because parts 1 & 2 are an excellent approximation
       when x is small and part 3 converges quickly when x is large. */
 
    eta = pow(6*x, enzo_float(1.0/3.0));                     // part 1
    eta = pow(120*x/(20+eta*eta), enzo_float(1.0/3.0));      // part 2
    for (i = 0; i < 40; i++) {                          // part 3
      eta_old = eta;
      eta = asinh(eta + x);
      if (fabs(eta-eta_old) < ETA_TOLERANCE) break;
    }
    if (i == 40) {
      fprintf(stderr, "Case 2 -- no convergence after %"ISYM" iterations.\n", i);
      return ENZO_FAIL;
    }
 
    /* Now use eta to compute the expansion factor (eq. 13-10, part 2). */
 
    *a = OmegaMatterNow/(2*(1 - OmegaMatterNow))*(cosh(eta) - 1);
    *a *= (1 + InitialRedshift);    // to convert to code units, divide by [a]
  }
 
  /* 3) For OmegaMatterNow > 1 && OmegaLambdaNow == 0, use sin/cos.
        Easy, but skip it for now. */
 
  if (OmegaMatterNow > 1 && OmegaLambdaNow < OMEGA_TOLERANCE) {
  }
 
  /* 4) For flat universe, with non-zero OmegaLambdaNow, see eq. 13-20. */
 
  if (fabs(OmegaCurvatureNow) < OMEGA_TOLERANCE &&
      OmegaLambdaNow > OMEGA_TOLERANCE) {
    *a = pow(enzo_float(OmegaMatterNow/(1 - OmegaMatterNow)), enzo_float(1.0/3.0)) *
         pow(enzo_float(sinh(1.5 * sqrt(1.0 - OmegaMatterNow)*TimeHubble0)),
	     enzo_float(2.0/3.0));
    *a *= (1 + InitialRedshift);    // to convert to code units, divide by [a]
  }
 
#endif /* INVERSE_HYPERBOLIC_EXISTS */
 
  /* Compute the derivative of the expansion factor (Peebles93, eq. 13.3). */
 
  enzo_float TempVal = (*a)/(1 + InitialRedshift);
  *dadt = sqrt( 2.0/(3.0*OmegaMatterNow*(*a)) *
	       (OmegaMatterNow + OmegaCurvatureNow*TempVal +
		OmegaLambdaNow*TempVal*TempVal*TempVal));
 
  /* Someday, we'll implement the general case... */
 
  if ((*a) == ENZO_FLOAT_UNDEFINED) {
    fprintf(stderr, "Cosmology selected is not implemented.\n");
    return ENZO_FAIL;
  }
 
  return ENZO_SUCCESS;
}

//---------------------------------------------------------------------- 
 
int EnzoBlock::CosmologyComputeExpansionTimestep
(enzo_float time, enzo_float *dtExpansion)
{
 
  /* Error check. */
 
  if (InitialTimeInCodeUnits == 0) {
    fprintf(stderr, "The cosmology parameters seem to be improperly set.\n");
    return ENZO_FAIL;
  }
 
  /* Compute the expansion factors. */
 
  enzo_float a, dadt;
  if (CosmologyComputeExpansionFactor(time, &a, &dadt) == ENZO_FAIL) {
    fprintf(stderr, "Error in ComputeExpnasionFactors.\n");
    return ENZO_FAIL;
  }
 
  /* Compute the maximum allwed timestep given the maximum allowed
     expansion factor. */
 
  *dtExpansion = MaxExpansionRate*a/dadt;
 
  return ENZO_SUCCESS;
}
