// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar  3 23:02:02 PST 2011
/// @brief    Implementation of the EnzoBlock class

#include "cello.hpp"

#include "enzo.hpp"

// #define DEBUG_ENZO_BLOCK

//======================================================================

int EnzoBlock::UseMinimumPressureSupport[CONFIG_NODE_SIZE];
enzo_float EnzoBlock::MinimumPressureSupportParameter[CONFIG_NODE_SIZE];
enzo_float EnzoBlock::ComovingBoxSize[CONFIG_NODE_SIZE];
enzo_float EnzoBlock::HubbleConstantNow[CONFIG_NODE_SIZE];
enzo_float EnzoBlock::OmegaMatterNow[CONFIG_NODE_SIZE];
enzo_float EnzoBlock::OmegaLambdaNow[CONFIG_NODE_SIZE];
enzo_float EnzoBlock::MaxExpansionRate[CONFIG_NODE_SIZE];

// Chemistry

int EnzoBlock::MultiSpecies[CONFIG_NODE_SIZE];

// Physics

int EnzoBlock::PressureFree[CONFIG_NODE_SIZE];
enzo_float EnzoBlock::Gamma[CONFIG_NODE_SIZE];
enzo_float EnzoBlock::GravitationalConstant[CONFIG_NODE_SIZE];

// Problem-specific

int EnzoBlock::ProblemType[CONFIG_NODE_SIZE];

// Method PPM

int EnzoBlock::PPMFlatteningParameter[CONFIG_NODE_SIZE];
int EnzoBlock::PPMDiffusionParameter[CONFIG_NODE_SIZE];
int EnzoBlock::PPMSteepeningParameter[CONFIG_NODE_SIZE];

// Numerics

int EnzoBlock::DualEnergyFormalism[CONFIG_NODE_SIZE];
enzo_float EnzoBlock::DualEnergyFormalismEta1[CONFIG_NODE_SIZE];
enzo_float EnzoBlock::DualEnergyFormalismEta2[CONFIG_NODE_SIZE];

enzo_float EnzoBlock::pressure_floor[CONFIG_NODE_SIZE];
enzo_float EnzoBlock::density_floor[CONFIG_NODE_SIZE];
enzo_float EnzoBlock::number_density_floor[CONFIG_NODE_SIZE];
enzo_float EnzoBlock::temperature_floor[CONFIG_NODE_SIZE];

enzo_float EnzoBlock::InitialRedshift[CONFIG_NODE_SIZE];
enzo_float EnzoBlock::InitialTimeInCodeUnits[CONFIG_NODE_SIZE];

// Domain

enzo_float EnzoBlock::DomainLeftEdge [3*CONFIG_NODE_SIZE];
enzo_float EnzoBlock::DomainRightEdge[3*CONFIG_NODE_SIZE];

// PPM

int EnzoBlock::GridRank[CONFIG_NODE_SIZE];

int EnzoBlock::ghost_depth[3*CONFIG_NODE_SIZE];

// Fields

int EnzoBlock::NumberOfBaryonFields[CONFIG_NODE_SIZE];

//----------------------------------------------------------------------

// STATIC
void EnzoBlock::initialize(EnzoConfig * enzo_config,
			   FieldDescr * field_descr)
{
#ifdef DEBUG_ENZO_BLOCK
  CkPrintf ("%d DEBUG_ENZO_BLOCK EnzoBlock::initialize\n",CkMyPe());
#endif
  int gx = enzo_config->field_ghost_depth[0];
  int gy = enzo_config->field_ghost_depth[1];
  int gz = enzo_config->field_ghost_depth[2];

  const int rank = enzo_config->mesh_root_rank;

  if (rank < 1) gx = 0;
  if (rank < 2) gy = 0;
  if (rank < 3) gz = 0;

  double time  = enzo_config->initial_time;

  for (int in=0; in<CONFIG_NODE_SIZE; in++) {

    GridRank[in] = 0;
    NumberOfBaryonFields[in] = 0;

    int i;

    for (i=0; i<MAX_DIMENSION; i++) {
      DomainLeftEdge [in*3+i] = 0;
      DomainRightEdge[in*3+i] = 0;
      ghost_depth[in*3+i] = 0;
    }

    Gamma[in]               = enzo_config->field_gamma;

    GridRank[in]            = enzo_config->mesh_root_rank;

    // Chemistry parameters

    MultiSpecies[in] = 0;    // 0:0 1:6 2:9 3:12

    // Gravity parameters

    GravitationalConstant[in]           = 1.0;  // used only in SetMinimumSupport()

    //Problem specific parameter

    ProblemType[in] = 0;

    // PPM parameters

    InitialRedshift[in]   = enzo_config->physics_cosmology_initial_redshift;
    HubbleConstantNow[in] = enzo_config->physics_cosmology_hubble_constant_now;
    OmegaLambdaNow[in]    = enzo_config->physics_cosmology_omega_lamda_now;
    OmegaMatterNow[in]    = enzo_config->physics_cosmology_omega_matter_now;
    MaxExpansionRate[in]  = enzo_config->physics_cosmology_max_expansion_rate;
    ComovingBoxSize[in]   = enzo_config->physics_cosmology_comoving_box_size;

    PressureFree[in]              = enzo_config->ppm_pressure_free;
    UseMinimumPressureSupport[in] = enzo_config->ppm_use_minimum_pressure_support;
    MinimumPressureSupportParameter[in] = 
      enzo_config->ppm_minimum_pressure_support_parameter;
    PPMFlatteningParameter[in]    = enzo_config->ppm_flattening;
    PPMDiffusionParameter[in]     = enzo_config->ppm_diffusion;
    PPMSteepeningParameter[in]    = enzo_config->ppm_steepening;
    pressure_floor[in]            = enzo_config->ppm_pressure_floor;
    density_floor[in]             = enzo_config->ppm_density_floor;
    temperature_floor[in]         = enzo_config->ppm_temperature_floor;
    number_density_floor[in]      = enzo_config->ppm_number_density_floor;
    DualEnergyFormalism[in]       = enzo_config->ppm_dual_energy;
    DualEnergyFormalismEta1[in]   = enzo_config->ppm_dual_energy_eta_1;
    DualEnergyFormalismEta2[in]   = enzo_config->ppm_dual_energy_eta_2;

    ghost_depth[in*3+0] = gx;
    ghost_depth[in*3+1] = gy;
    ghost_depth[in*3+2] = gz;

    NumberOfBaryonFields[in] = enzo_config->field_list.size();

    // Check NumberOfBaryonFields

    if (NumberOfBaryonFields[in] > MAX_NUMBER_OF_BARYON_FIELDS) {
      ERROR2 ("EnzoBlock::initialize",
	      "MAX_NUMBER_OF_BARYON_FIELDS = %d is too small for %d fields",
	      MAX_NUMBER_OF_BARYON_FIELDS,NumberOfBaryonFields[in] );
    }

    DomainLeftEdge [in*3+0] = enzo_config->domain_lower[0];
    DomainLeftEdge [in*3+1] = enzo_config->domain_lower[1];
    DomainLeftEdge [in*3+2] = enzo_config->domain_lower[2];

    DomainRightEdge[in*3+0] = enzo_config->domain_upper[0];
    DomainRightEdge[in*3+1] = enzo_config->domain_upper[1];
    DomainRightEdge[in*3+2] = enzo_config->domain_upper[2];

    InitialTimeInCodeUnits[in] = time;

  }

} // void initialize()

//----------------------------------------------------------------------

EnzoBlock::EnzoBlock
( MsgRefine * msg )
  : BASE_ENZO_BLOCK ( msg ),
    mg_iter_(0),
    mg_sync_(),
    dt(dt_),
    SubgridFluxes(NULL)
{
  initialize_enzo_();
  initialize();
#ifdef DEBUG_ENZO_BLOCK
  CkPrintf ("%d %p TRACE_BLOCK EnzoBlock(msg)\n",CkMyPe(),this);
  print();
#endif
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

EnzoBlock::~EnzoBlock()
{
#ifdef DEBUG_ENZO_BLOCK
  CkPrintf ("%d %p TRACE_BLOCK ~EnzoBlock(...)\n",CkMyPe(),this);
  print();
#endif
}

//----------------------------------------------------------------------

void EnzoBlock::pup(PUP::er &p)
{ 

  TRACEPUP;
  TRACE ("BEGIN EnzoBlock::pup()");

  BASE_ENZO_BLOCK::pup(p);

  p | dt;

  const int in = cello::index_static();

  static bool warn0[CONFIG_NODE_SIZE] = {true};
  if (warn0[in]) {
    warn0[in] = false;
    WARNING("EnzoBlock::pup()", "skipping AccelerationField_ (not used)");
  }

  static bool warn1[CONFIG_NODE_SIZE] = {true};
  if (warn1[in]) {
    warn1[in] = false;
    WARNING("EnzoBlock::pup()", "skipping SubgridFluxes (not used)");
  }

  PUParray(p,GridLeftEdge,MAX_DIMENSION); 
  PUParray(p,GridDimension,MAX_DIMENSION); 
  PUParray(p,GridStartIndex,MAX_DIMENSION); 
  PUParray(p,GridEndIndex,MAX_DIMENSION); 
  PUParray(p,CellWidth,MAX_DIMENSION);

  PUParray(p,method_turbulence_data,max_turbulence_array);

  p | mg_iter_;
  p | mg_sync_;

  TRACE ("END EnzoBlock::pup()");

}

//======================================================================

void EnzoBlock::write(FILE * fp) throw ()
{
  const int in = cello::index_static();

  fprintf (fp,"EnzoBlock: UseMinimumPressureSupport %d\n",
	   UseMinimumPressureSupport[in]);
  fprintf (fp,"EnzoBlock: MinimumPressureSupportParameter %g\n",
	   MinimumPressureSupportParameter[in]);
  fprintf (fp,"EnzoBlock: ComovingBoxSize %g\n",
	   ComovingBoxSize[in]);
  fprintf (fp,"EnzoBlock: HubbleConstantNow %g\n",
	   HubbleConstantNow[in]);
  fprintf (fp,"EnzoBlock: OmegaLambdaNow %g\n",
	   OmegaLambdaNow[in]);
  fprintf (fp,"EnzoBlock: OmegaMatterNow %g\n",
	   OmegaMatterNow[in]);
  fprintf (fp,"EnzoBlock: MaxExpansionRate %g\n",
	   MaxExpansionRate[in]);

  // Chemistry

  fprintf (fp,"EnzoBlock: MultiSpecies %d\n",
	   MultiSpecies[in]);

  // Physics

  fprintf (fp,"EnzoBlock: PressureFree %d\n",
	   PressureFree[in]);
  fprintf (fp,"EnzoBlock: Gamma %g\n",
	   Gamma[in]);
  fprintf (fp,"EnzoBlock: GravitationalConstant %g\n",
	   GravitationalConstant[in]);

  // Problem-specific

  fprintf (fp,"EnzoBlock: ProblemType %d\n",
	   ProblemType[in]);

  // Method PPM

  fprintf (fp,"EnzoBlock: PPMFlatteningParameter %d\n",
	   PPMFlatteningParameter[in]);
  fprintf (fp,"EnzoBlock: PPMDiffusionParameter %d\n",
	   PPMDiffusionParameter[in]);
  fprintf (fp,"EnzoBlock: PPMSteepeningParameter %d\n",
	   PPMSteepeningParameter[in]);

  // Numerics

  fprintf (fp,"EnzoBlock: DualEnergyFormalism %d\n",
	   DualEnergyFormalism[in]);
  fprintf (fp,"EnzoBlock: DualEnergyFormalismEta1 %g\n",
	   DualEnergyFormalismEta1[in]);
  fprintf (fp,"EnzoBlock: DualEnergyFormalismEta2 %g\n",
	   DualEnergyFormalismEta2[in]);
  fprintf (fp,"EnzoBlock: pressure_floor %g\n",
	   pressure_floor[in]);
  fprintf (fp,"EnzoBlock: density_density_floor %g\n",
	   density_floor[in]);
  fprintf (fp,"EnzoBlock: number_density_floor %g\n",
	   number_density_floor[in]);
  fprintf (fp,"EnzoBlock: temperature_floor %g\n",
	   temperature_floor[in]);

  fprintf (fp,"EnzoBlock: InitialRedshift %g\n",
	   InitialRedshift[in]);
  fprintf (fp,"EnzoBlock: InitialTimeInCodeUnits %g\n",
	   InitialTimeInCodeUnits[in]);

  // Domain

  fprintf (fp,"EnzoBlock: DomainLeftEdge %g %g %g\n",
	   DomainLeftEdge [in*3+0],
	   DomainLeftEdge [in*3+1],
	   DomainLeftEdge [in*3+2]);
  fprintf (fp,"EnzoBlock: DomainRightEdge %g %g %g\n",
	   DomainRightEdge[in*3+0],
	   DomainRightEdge[in*3+1],
	   DomainRightEdge[in*3+2]);

  // Fields

  // Grid

  fprintf (fp,"EnzoBlock: GridRank %d\n",    GridRank[in]);
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
	   ghost_depth[in*3+0],
	   ghost_depth[in*3+1],
	   ghost_depth[in*3+2]);


  fprintf (fp,"EnzoBlock: NumberOfBaryonFields %d\n",
	   NumberOfBaryonFields[in]);

  // problem

  fprintf (fp,"EnzoBlock: dt %g\n", dt);

  // fluxes

  fprintf (fp,"EnzoBlock: SubgridFluxes %p\n", SubgridFluxes);
  
}

//----------------------------------------------------------------------

void EnzoBlock::set_dt (double dt_param) throw ()
{
  Block::set_dt (dt_param);

  dt = dt_param;
}

//----------------------------------------------------------------------

void EnzoBlock::set_stop (bool stop) throw ()
{
  Block::set_stop (stop);
}

//----------------------------------------------------------------------

void EnzoBlock::initialize () throw()
{
  TRACE ("EnzoBlock::initialize()\n");

  double xm,ym,zm;

  data()->lower(&xm,&ym,&zm);

  GridLeftEdge[0]  = xm;
  GridLeftEdge[1]  = ym;
  GridLeftEdge[2]  = zm;

  // Grid dimensions

  Field field = data()->field();

  int nx,ny,nz;
  field.size (&nx,&ny,&nz);

  int gx,gy,gz;

  const int in = cello::index_static();

  gx = EnzoBlock::ghost_depth[in*3+0];
  gy = EnzoBlock::ghost_depth[in*3+1];
  gz = EnzoBlock::ghost_depth[in*3+2];

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
  data()->upper(&xp,&yp,&zp);
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
 
  const int in = cello::index_static();

  /* Error check. */

  if (InitialTimeInCodeUnits[in] == 0) {
    
    char error_message[ERROR_LENGTH];
    sprintf(error_message, "The cosmology parameters seem to be improperly set");
    ERROR("CosmologyComputeExpansionFactor",error_message);
  }
 
  *a = ENZO_FLOAT_UNDEFINED;
 
  /* Find Omega due to curvature. */
 
  enzo_float OmegaCurvatureNow = 1 - OmegaMatterNow[in] - OmegaLambdaNow[in];
 
  /* Convert the time from code units to Time * H0 (c.f. CosmologyGetUnits). */
 
  enzo_float TimeUnits = 2.52e17/sqrt(OmegaMatterNow[in])/HubbleConstantNow[in]/
                    pow(1 + InitialRedshift[in],enzo_float(1.5));
 
  enzo_float TimeHubble0 = time * TimeUnits * (HubbleConstantNow[in]*3.24e-18);
 
  /* 1) For a flat universe with OmegaMatterNow = 1, it's easy. */
 
  if (fabs(OmegaMatterNow[in]-1) < OMEGA_TOLERANCE &&
      OmegaLambdaNow[in] < OMEGA_TOLERANCE)
    *a    = pow(time/InitialTimeInCodeUnits[in], enzo_float(2.0/3.0));
 
#define INVERSE_HYPERBOLIC_EXISTS
 
#ifdef INVERSE_HYPERBOLIC_EXISTS
 
 
  /* 2) For OmegaMatterNow < 1 and OmegaLambdaNow == 0 see
        Peebles 1993, eq. 13-3, 13-10.
	Actually, this is a little tricky since we must solve an equation
	of the form eta - sinh(eta) + x = 0..*/
 
  if (OmegaMatterNow[in] < 1 && OmegaLambdaNow[in] < OMEGA_TOLERANCE) {
 
    enzo_float eta, eta_old, x;
    int i;

    x = 2*TimeHubble0*pow(1.0 - OmegaMatterNow[in], 1.5) / OmegaMatterNow[in];
 
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
      fprintf(stderr, "Case 2 -- no convergence after %" ISYM " iterations.\n", i);
      return ENZO_FAIL;
    }
 
    /* Now use eta to compute the expansion factor (eq. 13-10, part 2). */
 
    *a = OmegaMatterNow[in]/(2*(1 - OmegaMatterNow[in]))*(cosh(eta) - 1);
    *a *= (1 + InitialRedshift[in]);    // to convert to code units, divide by [a]
  }
 
  /* 3) For OmegaMatterNow > 1 && OmegaLambdaNow == 0, use sin/cos.
        Easy, but skip it for now. */
 
  if (OmegaMatterNow[in] > 1 && OmegaLambdaNow[in] < OMEGA_TOLERANCE) {
  }
 
  /* 4) For flat universe, with non-zero OmegaLambdaNow, see eq. 13-20. */
 
  if (fabs(OmegaCurvatureNow) < OMEGA_TOLERANCE &&
      OmegaLambdaNow[in] > OMEGA_TOLERANCE) {
    *a = pow(enzo_float(OmegaMatterNow[in]/(1 - OmegaMatterNow[in])), enzo_float(1.0/3.0)) *
         pow(enzo_float(sinh(1.5 * sqrt(1.0 - OmegaMatterNow[in])*TimeHubble0)),
	     enzo_float(2.0/3.0));
    *a *= (1 + InitialRedshift[in]);    // to convert to code units, divide by [a]
  }
 
#endif /* INVERSE_HYPERBOLIC_EXISTS */
 
  /* Compute the derivative of the expansion factor (Peebles93, eq. 13.3). */
 
  enzo_float TempVal = (*a)/(1 + InitialRedshift[in]);
  *dadt = sqrt( 2.0/(3.0*OmegaMatterNow[in]*(*a)) *
	       (OmegaMatterNow[in] + OmegaCurvatureNow*TempVal +
		OmegaLambdaNow[in]*TempVal*TempVal*TempVal));
 
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
 
  const int in = cello::index_static();

  if (InitialTimeInCodeUnits[in] == 0) {
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
 
  *dtExpansion = MaxExpansionRate[in]*a/dadt;
 
  return ENZO_SUCCESS;
}

