// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoDescr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Aug 31 15:38:36 PDT 2010
/// @brief    Implementation of EnzoDescr methods

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoDescr::EnzoDescr() throw ()
  : ComovingCoordinates(0),
    UseMinimumPressureSupport(0),
    MinimumPressureSupportParameter(0),
    ComovingBoxSize(0),
    HubbleConstantNow(0),
    OmegaMatterNow(0),
    OmegaLambdaNow(0),
    MaxExpansionRate(0),
    MultiSpecies(0),
    GravityOn(0),
    PressureFree(0),
    Gamma(0),
    GravitationalConstant(0),
    ProblemType(0),
    PPMFlatteningParameter(0),
    PPMDiffusionParameter(0),
    PPMSteepeningParameter(0),
    ProcessorNumber(0),
    DualEnergyFormalism(0),
    DualEnergyFormalismEta1(0),
    DualEnergyFormalismEta2(0),
    pressure_floor(0),
    density_floor(0),
    number_density_floor(0),
    temperature_floor(0),
    CourantSafetyNumber(0),
    InitialRedshift(0),
    InitialTimeInCodeUnits(0),
    Time(0),
    OldTime(0),
    field_density(-1),
    field_total_energy(-1),
    field_internal_energy(-1),
    field_velocity_x(-1),
    field_velocity_y(-1),
    field_velocity_z(-1),
    field_color(-1),
    field_magnetic_x(-1),
    field_magnetic_y(-1),
    field_magnetic_z(-1),
    field_density_xp(-1),
    field_velocity_x_xp(-1),
    field_velocity_y_xp(-1),
    field_velocity_z_xp(-1),
    field_magnetic_x_xp(-1),
    field_magnetic_y_xp(-1),
    field_magnetic_z_xp(-1),
    field_density_yp(-1),
    field_velocity_x_yp(-1),
    field_velocity_y_yp(-1),
    field_velocity_z_yp(-1),
    field_magnetic_x_yp(-1),
    field_magnetic_y_yp(-1),
    field_magnetic_z_yp(-1),
    field_density_zp(-1),
    field_velocity_x_zp(-1),
    field_velocity_y_zp(-1),
    field_velocity_z_zp(-1),
    field_magnetic_x_zp(-1),
    field_magnetic_y_zp(-1),
    field_magnetic_z_zp(-1),
    GridRank(0),
    NumberOfBaryonFields(0),
    BoundaryRank(0),
    CycleNumber(0),
    dt(0),
    SubgridFluxes(0)

{

  int i,j,k;

  for (i=0; i<MAX_DIMENSION; i++) {
    AccelerationField[i] = 0;
    DomainLeftEdge [i] = 0;
    DomainRightEdge[i] = 0;
    GridDimension[i] = 0;
    GridStartIndex[i] = 0;
    GridEndIndex[i] = 0;
    GridLeftEdge[i] = 0;
    CellWidth[i] = 0;
    ghost_depth[i] = 0;
  }

  for (j=0; j<MAX_NUMBER_OF_BARYON_FIELDS; j++) {
    BaryonField[j] = 0;
    OldBaryonField[j] = 0;
    FieldType[j] = 0;
  }

  for (i=0; i<MAX_DIMENSION; i++) {
    for (j=0; j<MAX_NUMBER_OF_BARYON_FIELDS; j++) {
      for (k=0; k<2; k++) {
	BoundaryDimension[i]=0;
	BoundaryFieldType[j] = 0;
	BoundaryType[j][i][k] = 0;
	BoundaryValue[j][i][k]= 0; 
      }
    }
  }
}

//----------------------------------------------------------------------

void
EnzoDescr::initialize(Parameters * parameters) throw ()
{

  //--------------------------------------------------
  parameters->set_current_group ("Physics");
  //--------------------------------------------------

  ComovingCoordinates  = parameters->value_logical ("cosmology",false);
  Gamma                = parameters->value_scalar  ("gamma",5.0/3.0);
  CycleNumber = 0;

  // PPM parameters

  //--------------------------------------------------
  parameters->set_current_subgroup ("cosmology");
  //--------------------------------------------------

  InitialRedshift   = parameters->value_scalar ("initial_redshift",  20.0);
  HubbleConstantNow = parameters->value_scalar ("hubble_constant_now",0.701);
  OmegaLambdaNow    = parameters->value_scalar ("omega_lambda_now",   0.721);
  OmegaMatterNow    = parameters->value_scalar ("omega_matter_now",   0.279);
  MaxExpansionRate  = parameters->value_scalar ("max_expansion_rate", 0.01);
  ComovingBoxSize   = parameters->value_scalar ("comoving_box_size", 64.0);

  //--------------------------------------------------
  parameters->set_current_group ("Method","ppm");
  //--------------------------------------------------

  PressureFree = parameters->value_scalar("pressure_free",false);
  UseMinimumPressureSupport 
    =              parameters->value_logical("use_minimum_pressure_support",false);
  MinimumPressureSupportParameter 
    =              parameters->value_integer("minimum_pressure_support_parameter",100);
  PPMFlatteningParameter = parameters->value_logical ("flattening", false);
  PPMDiffusionParameter  = parameters->value_logical ("diffusion",  false);
  PPMSteepeningParameter = parameters->value_logical ("steepening", false);

  double floor_default = 1e-6;
  pressure_floor       = parameters->value_scalar("pressure_floor",      floor_default);
  density_floor        = parameters->value_scalar("density_floor",       floor_default);
  temperature_floor    = parameters->value_scalar("temperature_floor",   floor_default);
  number_density_floor = parameters->value_scalar("number_density_floor",floor_default);

  DualEnergyFormalism     = parameters->value_logical ("dual_energy",false);
  DualEnergyFormalismEta1 = parameters->value_scalar  ("dual_energy_eta_1",0.001);
  DualEnergyFormalismEta2 = parameters->value_scalar  ("dual_energy_eta_1",0.1);

  //--------------------------------------------------
  parameters->set_current_group ("Physics");
  //--------------------------------------------------

  BoundaryRank = parameters->value_integer("dimensions",0);

  //--------------------------------------------------
  parameters->set_current_group ("Mesh");
  //--------------------------------------------------

  int nx = parameters->list_value_integer(0,"block_size",1);
  int ny = parameters->list_value_integer(1,"block_size",1);
  int nz = parameters->list_value_integer(2,"block_size",1);

  //--------------------------------------------------
  parameters->set_current_group ("Field");
  //--------------------------------------------------

  int gx = parameters->list_value_integer(0,"ghosts",0);
  int gy = parameters->list_value_integer(1,"ghosts",0);
  int gz = parameters->list_value_integer(2,"ghosts",0);

  BoundaryDimension[0] = nx + 2*gx;
  BoundaryDimension[1] = ny + 2*gy;
  BoundaryDimension[2] = nz + 2*gz;


  // Initialize Enzo field-related attributes

  //--------------------------------------------------
  parameters->set_current_group("Field");
  //--------------------------------------------------

  NumberOfBaryonFields = parameters->list_length("fields");

  if (NumberOfBaryonFields == 0) {
    ERROR_MESSAGE ("EnzoDescr::EnzoDescr",
		   "List parameter 'Field fields' must have length greater than zero");
  }

  for (int field_index=0; field_index<NumberOfBaryonFields; field_index++) {

    std::string name = 
      parameters->list_value_string(field_index,"fields");

    if        (name == "density") {
      field_density          = field_index;
      FieldType[field_index] = Density;
    } else if (name == "velocity_x") {
      field_velocity_x       = field_index;
      FieldType[field_index] = Velocity1;
    } else if (name == "velocity_y") {
      field_velocity_y       = field_index;
      FieldType[field_index] = Velocity2;
    } else if (name == "velocity_z") {
      field_velocity_z       = field_index;
      FieldType[field_index] = Velocity3;
    } else if (name == "total_energy") {
      field_total_energy     = field_index;
      FieldType[field_index] = TotalEnergy;
    } else if (name == "internal_energy") {
      field_internal_energy  = field_index;
      FieldType[field_index] = InternalEnergy;
    } else if (name == "electron_density") {
      field_color            = field_index;
      FieldType[field_index] = ElectronDensity;
    } else {
      char error_message[ERROR_MESSAGE_LENGTH];
      sprintf (error_message,"Unknown field %s",name.c_str());
      ERROR_MESSAGE ("EnzoDescr::EnzoDescr", error_message);
    }
  }

  //--------------------------------------------------
  parameters->set_current_group("Physics");
  //--------------------------------------------------

  GridRank = parameters->value_integer ("dimensions",0);

  // Chemistry parameters

  MultiSpecies = 0;    // 0:0 1:6 2:9 3:12

  // Gravity parameters

  GravityOn                       = 0;    // Whether gravity is included
  GravitationalConstant           = 1.0;  // used only in SetMinimumSupport()
  AccelerationField[0]            = NULL;
  AccelerationField[1]            = NULL;
  AccelerationField[2]            = NULL;

  //Problem specific parameter

  ProblemType = 0;

  // Field parameters

  //--------------------------------------------------
  parameters->set_current_group ("Field");
  //--------------------------------------------------

  ghost_depth[0] = (GridRank >= 1) ? 
    parameters->list_value_integer(0,"ghosts",3) : 0;
  ghost_depth[1] = (GridRank >= 2) ? 
    parameters->list_value_integer(1,"ghosts",3) : 0;
  ghost_depth[2] = (GridRank >= 3) ? 
    parameters->list_value_integer(2,"ghosts",3) : 0;

  printf ("NumberOfBaryonFields = %d\n",NumberOfBaryonFields );
  ASSERT ("initialize",
	  "MAX_NUMBER_OF_BARYON_FIELDS is too small",
	  NumberOfBaryonFields <= MAX_NUMBER_OF_BARYON_FIELDS);

  Time                   = 0;

  // Domain parameters

  //--------------------------------------------------
  parameters->set_current_group ("Domain");
  //--------------------------------------------------
  
  DomainLeftEdge [0] = parameters->list_value_scalar(0,"extent",0.0);
  DomainRightEdge[0] = parameters->list_value_scalar(1,"extent",1.0);
  DomainLeftEdge [1] = parameters->list_value_scalar(2,"extent",0.0);
  DomainRightEdge[1] = parameters->list_value_scalar(3,"extent",1.0);
  DomainLeftEdge [2] = parameters->list_value_scalar(4,"extent",0.0);
  DomainRightEdge[2] = parameters->list_value_scalar(5,"extent",1.0);

  // Initial conditions

  //--------------------------------------------------
  parameters->set_current_group ("Initial");
  //--------------------------------------------------

  InitialTimeInCodeUnits = parameters->value_scalar ("time",0.0);
  Time = InitialTimeInCodeUnits;
  OldTime = Time;

  // Parallel parameters

  //--------------------------------------------------
  parameters->set_current_group ("Parallel");
  //--------------------------------------------------

  std::string parallel_method = 
    parameters->list_value_string(0,"method","serial");

  GroupProcess * parallel = 0;

  parallel = GroupProcess::create();

  ProcessorNumber = parallel->rank();

  delete parallel;

  //--------------------------------------------------
  parameters->set_current_group ("Field");
  //--------------------------------------------------
  
  CourantSafetyNumber = parameters->value_scalar ("courant",0.6);

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // END: Moved from Enzo MethodControl::initialize
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

}

EnzoDescr::~EnzoDescr() throw ()
{
}

//----------------------------------------------------------------------

void EnzoDescr::write(FILE * fp) throw ()
{
  fprintf (fp,"EnzoDescr: ComovingCoordinates %d\n",
	   ComovingCoordinates);
  fprintf (fp,"EnzoDescr: UseMinimumPressureSupport %d\n",
	   UseMinimumPressureSupport);
  fprintf (fp,"EnzoDescr: MinimumPressureSupportParameter %g\n",
	   MinimumPressureSupportParameter);
  fprintf (fp,"EnzoDescr: ComovingBoxSize %g\n",
	   ComovingBoxSize);
  fprintf (fp,"EnzoDescr: HubbleConstantNow %g\n",
	   HubbleConstantNow);
  fprintf (fp,"EnzoDescr: OmegaLambdaNow %g\n",
	   OmegaLambdaNow);
  fprintf (fp,"EnzoDescr: OmegaMatterNow %g\n",
	   OmegaMatterNow);
  fprintf (fp,"EnzoDescr: MaxExpansionRate %g\n",
	   MaxExpansionRate);

  // Chemistry

  fprintf (fp,"EnzoDescr: MultiSpecies %d\n",
	   MultiSpecies);

  // Gravity

  fprintf (fp,"EnzoDescr: GravityOn %d\n",
	   GravityOn);
  //  fprintf (fp,"EnzoDescr: *AccelerationField %g\n",
  //           *AccelerationField[MAX_DIMENSION)];

  // Physics

  fprintf (fp,"EnzoDescr: PressureFree %d\n",
	   PressureFree);
  fprintf (fp,"EnzoDescr: Gamma %g\n",
	   Gamma);
  fprintf (fp,"EnzoDescr: GravitationalConstant %g\n",
	   GravitationalConstant);

  // Problem-specific

  fprintf (fp,"EnzoDescr: ProblemType %d\n",
	   ProblemType);

  // Method PPM

  fprintf (fp,"EnzoDescr: PPMFlatteningParameter %d\n",
	   PPMFlatteningParameter);
  fprintf (fp,"EnzoDescr: PPMDiffusionParameter %d\n",
	   PPMDiffusionParameter);
  fprintf (fp,"EnzoDescr: PPMSteepeningParameter %d\n",
	   PPMSteepeningParameter);

  // Parallel

  fprintf (fp,"EnzoDescr: ProcessorNumber %d\n",
	   ProcessorNumber);

  // Numerics

  fprintf (fp,"EnzoDescr: DualEnergyFormalism %d\n",
	   DualEnergyFormalism);
  fprintf (fp,"EnzoDescr: DualEnergyFormalismEta1 %g\n",
	   DualEnergyFormalismEta1);
  fprintf (fp,"EnzoDescr: DualEnergyFormalismEta2 %g\n",
	   DualEnergyFormalismEta2);
  fprintf (fp,"EnzoDescr: pressure_floor %g\n",
	   pressure_floor);
  fprintf (fp,"EnzoDescr: density_density_floor %g\n",
	   density_floor);
  fprintf (fp,"EnzoDescr: number_density_floor %g\n",
	   number_density_floor);
  fprintf (fp,"EnzoDescr: temperature_floor %g\n",
	   temperature_floor);

  fprintf (fp,"EnzoDescr: CourantSafetyNumber %g\n",
	   CourantSafetyNumber);
  fprintf (fp,"EnzoDescr: InitialRedshift %g\n",
	   InitialRedshift);
  fprintf (fp,"EnzoDescr: InitialTimeInCodeUnits %g\n",
	   InitialTimeInCodeUnits);
  fprintf (fp,"EnzoDescr: Time %g\n",
	   Time);
  fprintf (fp,"EnzoDescr: OldTime %g\n",
	   OldTime);

  // Domain

  fprintf (fp,"EnzoDescr: DomainLeftEdge %g %g %g\n",
	   DomainLeftEdge [0],DomainLeftEdge [0],DomainLeftEdge [0]);
  fprintf (fp,"EnzoDescr: DomainRightEdge %g %g %g\n",
	   DomainRightEdge[0],DomainRightEdge[1],DomainRightEdge[2]);

  // Fields

  if (field_density != -1) 
    fprintf (fp,"EnzoDescr: field_density %d\n", field_density);
  if (field_total_energy != -1) 
    fprintf (fp,"EnzoDescr: field_total_energy %d\n", field_total_energy);
  if (field_internal_energy != -1) 
    fprintf (fp,"EnzoDescr: field_internal_energy %d\n", field_internal_energy);
  if (field_velocity_x != -1) 
    fprintf (fp,"EnzoDescr: field_velocity_x %d\n", field_velocity_x);
  if (field_velocity_y != -1) 
    fprintf (fp,"EnzoDescr: field_velocity_y %d\n", field_velocity_y);
  if (field_velocity_z != -1) 
    fprintf (fp,"EnzoDescr: field_velocity_z %d\n", field_velocity_z);
  if (field_color != -1) 
    fprintf (fp,"EnzoDescr: field_color %d\n", field_color);

  if (field_magnetic_x != -1) 
    fprintf (fp,"EnzoDescr: field_magnetic_x %d\n", field_magnetic_x);
  if (field_magnetic_y != -1) 
    fprintf (fp,"EnzoDescr: field_magnetic_y %d\n", field_magnetic_y);
  if (field_magnetic_z != -1) 
    fprintf (fp,"EnzoDescr: field_magnetic_z %d\n", field_magnetic_z);

  if (field_density_xp != -1) 
    fprintf (fp,"EnzoDescr: field_density_xp %d\n", field_density_xp);
  if (field_velocity_x_xp != -1) 
    fprintf (fp,"EnzoDescr: field_velocity_x_xp %d\n", field_velocity_x_xp);
  if (field_velocity_y_xp != -1) 
    fprintf (fp,"EnzoDescr: field_velocity_y_xp %d\n", field_velocity_y_xp);
  if (field_velocity_z_xp != -1) 
    fprintf (fp,"EnzoDescr: field_velocity_z_xp %d\n", field_velocity_z_xp);
  if (field_magnetic_x_xp != -1) 
    fprintf (fp,"EnzoDescr: field_magnetic_x_xp %d\n", field_magnetic_x_xp);
  if (field_magnetic_y_xp != -1) 
    fprintf (fp,"EnzoDescr: field_magnetic_y_xp %d\n", field_magnetic_y_xp);
  if (field_magnetic_z_xp != -1) 
    fprintf (fp,"EnzoDescr: field_magnetic_z_xp %d\n", field_magnetic_z_xp);

  if (field_density_yp != -1) 
    fprintf (fp,"EnzoDescr: field_density_yp %d\n", field_density_yp);
  if (field_velocity_x_yp != -1) 
    fprintf (fp,"EnzoDescr: field_velocity_x_yp %d\n", field_velocity_x_yp);
  if (field_velocity_y_yp != -1) 
    fprintf (fp,"EnzoDescr: field_velocity_y_yp %d\n", field_velocity_y_yp);
  if (field_velocity_z_yp != -1) 
    fprintf (fp,"EnzoDescr: field_velocity_z_yp %d\n", field_velocity_z_yp);
  if (field_magnetic_x_yp != -1) 
    fprintf (fp,"EnzoDescr: field_magnetic_x_yp %d\n", field_magnetic_x_yp);
  if (field_magnetic_y_yp != -1) 
    fprintf (fp,"EnzoDescr: field_magnetic_y_yp %d\n", field_magnetic_y_yp);
  if (field_magnetic_z_yp != -1) 
    fprintf (fp,"EnzoDescr: field_magnetic_z_yp %d\n", field_magnetic_z_yp);

  if (field_density_zp != -1) 
    fprintf (fp,"EnzoDescr: field_density_zp %d\n", field_density_zp);
  if (field_velocity_x_zp != -1) 
    fprintf (fp,"EnzoDescr: field_velocity_x_zp %d\n", field_velocity_x_zp);
  if (field_velocity_y_zp != -1) 
    fprintf (fp,"EnzoDescr: field_velocity_y_zp %d\n", field_velocity_y_zp);
  if (field_velocity_z_zp != -1) 
    fprintf (fp,"EnzoDescr: field_velocity_z_zp %d\n", field_velocity_z_zp);
  if (field_magnetic_x_zp != -1) 
    fprintf (fp,"EnzoDescr: field_magnetic_x_zp %d\n", field_magnetic_x_zp);
  if (field_magnetic_y_zp != -1) 
    fprintf (fp,"EnzoDescr: field_magnetic_y_zp %d\n", field_magnetic_y_zp);
  if (field_magnetic_z_zp != -1) 
    fprintf (fp,"EnzoDescr: field_magnetic_z_zp %d\n", field_magnetic_z_zp);

  // Grid

  fprintf (fp,"EnzoDescr: GridRank %d\n",    GridRank);
  fprintf (fp,"EnzoDescr: GridDimension %d %d %d\n",
	   GridDimension[0],GridDimension[1],GridDimension[2]);
  fprintf (fp,"EnzoDescr: GridStartIndex %d %d %d\n",
	   GridStartIndex[0],GridStartIndex[1],GridStartIndex[2]);
  fprintf (fp,"EnzoDescr: GridEndIndex %d %d %d\n",
	   GridEndIndex[0],GridEndIndex[1],GridEndIndex[2]);
  fprintf (fp,"EnzoDescr: GridLeftEdge %g %g %g\n",
	   GridLeftEdge[0],GridLeftEdge[1],GridLeftEdge[2]);

  fprintf (fp,"EnzoDescr: CellWidth %p %p %p\n", 
	   CellWidth[0], CellWidth[1], CellWidth[2] );
  if (CellWidth[0] != NULL) 
    fprintf (fp,"EnzoDescr: CellWidth[0][0] %g\n", CellWidth[0][0]);
  if (CellWidth[1] != NULL) 
    fprintf (fp,"EnzoDescr: CellWidth[1][0] %g\n", CellWidth[1][0]);
  if (CellWidth[2] != NULL) 
    fprintf (fp,"EnzoDescr: CellWidth[2][0] %g\n", CellWidth[2][0]);

  fprintf (fp,"EnzoDescr: ghost %d %d %d\n",
	   ghost_depth[0],ghost_depth[1],ghost_depth[2]);


  fprintf (fp,"EnzoDescr: NumberOfBaryonFields %d\n",
	   NumberOfBaryonFields);
  int i;
  for (i=0; i<NumberOfBaryonFields; i++) {
    fprintf (fp,"EnzoDescr: BaryonField[%d] %p\n",
	     i, BaryonField[i]);
    fprintf (fp,"EnzoDescr: OldBaryonField[%d] %p\n",
	     i, OldBaryonField[i]);
    fprintf (fp,"EnzoDescr: FieldType[%d] %d\n",
	     i, FieldType[i]);
  }

  fprintf (fp,"EnzoDescr: BoundaryRank %d\n", BoundaryRank);
  fprintf (fp,"EnzoDescr: BoundaryDimension %d %d %d\n",
	   BoundaryDimension[0],BoundaryDimension[1],BoundaryDimension[2]);

  // unknown, reflecting, outflow, inflow, periodic
  //  const char * bc_string[] = 
  //    {"unknown", "reflecting", "outflow", "inflow", "periodic"};

  for (i=0; i<NumberOfBaryonFields; i++) {

    fprintf (fp,"EnzoDescr: BoundaryFieldType[%d] %d\n", 
	     i, BoundaryFieldType[i]);

    fprintf (fp,"EnzoDescr: BoundaryType[%d] %p %p %p %p %p %p\n", i, 
	     BoundaryType[i][0][0],
	     BoundaryType[i][0][1],
	     BoundaryType[i][1][0],
	     BoundaryType[i][1][1],
	     BoundaryType[i][2][0],
	     BoundaryType[i][2][1]);

    fprintf (fp,"EnzoDescr: BoundaryValue[%d] %p %p %p %p %p %p\n",i,
	     BoundaryValue[i][0][0],
	     BoundaryValue[i][0][1],
	     BoundaryValue[i][1][0],
	     BoundaryValue[i][1][1],
	     BoundaryValue[i][2][0],
	     BoundaryValue[i][2][1]);  
  }

  // problem

  fprintf (fp,"EnzoDescr: CycleNumber %d\n",   CycleNumber);
  fprintf (fp,"EnzoDescr: dt %g\n", dt);

  // fluxes

  fprintf (fp,"EnzoDescr: SubgridFluxes %p\n", SubgridFluxes);
  

}

//----------------------------------------------------------------------
