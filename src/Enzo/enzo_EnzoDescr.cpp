// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoDescr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Aug 31 15:38:36 PDT 2010
/// @brief    Implementation of EnzoDescr methods

#include "cello.hpp"

#include "field.hpp"
#include "parameters.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoDescr::EnzoDescr(Global * global) throw ()
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
    field_density(0),
    field_total_energy(0),
    field_internal_energy(0),
    field_velocity_x(0),
    field_velocity_y(0),
    field_velocity_z(0),
    field_color(0),
    field_magnetic_x(0),
    field_magnetic_y(0),
    field_magnetic_z(0),
    field_density_xp(0),
    field_velocity_x_xp(0),
    field_velocity_y_xp(0),
    field_velocity_z_xp(0),
    field_magnetic_x_xp(0),
    field_magnetic_y_xp(0),
    field_magnetic_z_xp(0),
    field_density_yp(0),
    field_velocity_x_yp(0),
    field_velocity_y_yp(0),
    field_velocity_z_yp(0),
    field_magnetic_x_yp(0),
    field_magnetic_y_yp(0),
    field_magnetic_z_yp(0),
    field_density_zp(0),
    field_velocity_x_zp(0),
    field_velocity_y_zp(0),
    field_velocity_z_zp(0),
    field_magnetic_x_zp(0),
    field_magnetic_y_zp(0),
    field_magnetic_z_zp(0),
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


  // --------------------------------------------------
  // Initiazize user descriptors
  // --------------------------------------------------

  // --------------------------------------------------
  // Initiazize user method descriptors
  // --------------------------------------------------

  Parameters * parameters = global->parameters();

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

    std::string method_name = 
      parameters->list_value_string(field_index,"fields");

    if        (method_name == "density") {
      field_density          = field_index;
      FieldType[field_index] = Density;
    } else if (method_name == "velocity_x") {
      field_velocity_x       = field_index;
      FieldType[field_index] = Velocity1;
    } else if (method_name == "velocity_y") {
      field_velocity_y       = field_index;
      FieldType[field_index] = Velocity2;
    } else if (method_name == "velocity_z") {
      field_velocity_z       = field_index;
      FieldType[field_index] = Velocity3;
    } else if (method_name == "total_energy") {
      field_total_energy     = field_index;
      FieldType[field_index] = TotalEnergy;
    } else if (method_name == "internal_energy") {
      field_internal_energy  = field_index;
      FieldType[field_index] = InternalEnergy;
    } else if (method_name == "electron_density") {
      field_color            = field_index;
      FieldType[field_index] = ElectronDensity;
    } else {
      char error_message[ERROR_MESSAGE_LENGTH];
      sprintf (error_message,"Unknown field %s",method_name.c_str());
      ERROR_MESSAGE ("EnzoDescr::EnzoDescr", error_message);
    }
  }

  //--------------------------------------------------
  parameters->set_current_group("Physics");
  //--------------------------------------------------

  GridRank = parameters->value_integer ("dimensions",0);

  if (GridRank >= 1) {
  }

  if (GridRank >= 2) {
  }

  if (GridRank >= 3) {
  }


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

//----------------------------------------------------------------------

void
EnzoDescr::read_parameters_(Parameters * parameters) throw ()
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

}

//----------------------------------------------------------------------

EnzoDescr::~EnzoDescr() throw ()
{
}

//----------------------------------------------------------------------
