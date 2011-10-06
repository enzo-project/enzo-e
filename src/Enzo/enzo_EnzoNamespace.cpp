// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoNamespace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Aug 31 15:38:36 PDT 2010
/// @brief    "Global" Enzo data

#include "cello.hpp"

#include "enzo.hpp"

namespace enzo {

  /// @namespace  enzo
  /// @brief      Namespace for Enzo functions and global variables

  /// Boundary

  int  BoundaryRank;
  int  BoundaryDimension[MAX_DIMENSION];
  int  BoundaryFieldType[MAX_NUMBER_OF_BARYON_FIELDS];
  bc_enum *BoundaryType[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2];
  enzo_float *BoundaryValue[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2];

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

  //  int ProcessorNumber;

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

  int ghost_depth[MAX_DIMENSION];

  // Fields

  int NumberOfBaryonFields;      // active baryon fields

  int FieldType[MAX_NUMBER_OF_BARYON_FIELDS];

  void initialize(Parameters * parameters)

  {

    BoundaryRank = 0;
    ComovingCoordinates = 0;
    UseMinimumPressureSupport = 0;
    MinimumPressureSupportParameter = 0;
    ComovingBoxSize = 0;
    HubbleConstantNow = 0;
    OmegaMatterNow = 0;
    OmegaLambdaNow = 0;
    MaxExpansionRate = 0;
    MultiSpecies = 0;
    GravityOn = 0;
    PressureFree = 0;
    Gamma = 0;
    GravitationalConstant = 0;
    ProblemType = 0;
    PPMFlatteningParameter = 0;
    PPMDiffusionParameter = 0;
    PPMSteepeningParameter = 0;
    //    ProcessorNumber = 0;
    DualEnergyFormalism = 0;
    DualEnergyFormalismEta1 = 0;
    DualEnergyFormalismEta2 = 0;
    pressure_floor = 0;
    density_floor = 0;
    number_density_floor = 0;
    temperature_floor = 0;
    CourantSafetyNumber = 0;
    InitialRedshift = 0;
    InitialTimeInCodeUnits = 0;
    //    Time = 0;
    //    OldTime = 0;
    field_density = -1;
    field_total_energy = -1;
    field_internal_energy = -1;
    field_velocity_x = -1;
    field_velocity_y = -1;
    field_velocity_z = -1;
    field_color = -1;
    field_magnetic_x = -1;
    field_magnetic_y = -1;
    field_magnetic_z = -1;
    field_density_xp = -1;
    field_velocity_x_xp = -1;
    field_velocity_y_xp = -1;
    field_velocity_z_xp = -1;
    field_magnetic_x_xp = -1;
    field_magnetic_y_xp = -1;
    field_magnetic_z_xp = -1;
    field_density_yp = -1;
    field_velocity_x_yp = -1;
    field_velocity_y_yp = -1;
    field_velocity_z_yp = -1;
    field_magnetic_x_yp = -1;
    field_magnetic_y_yp = -1;
    field_magnetic_z_yp = -1;
    field_density_zp = -1;
    field_velocity_x_zp = -1;
    field_velocity_y_zp = -1;
    field_velocity_z_zp = -1;
    field_magnetic_x_zp = -1;
    field_magnetic_y_zp = -1;
    field_magnetic_z_zp = -1;
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

    for (j=0; j<MAX_NUMBER_OF_BARYON_FIELDS; j++) {
      FieldType[j] = 0;
    }

    //--------------------------------------------------
    // parameter: Physics : cosmology
    // parameter: Physics : gamma
    // parameter: Physics : dimensions
    //--------------------------------------------------

    ComovingCoordinates = parameters->value_logical ("Physics:cosmology",false);
    Gamma               = parameters->value_float   ("Physics:gamma",5.0/3.0);
    GridRank            = parameters->value_integer ("Physics:dimensions",0);
    BoundaryRank = GridRank;

    // Chemistry parameters

    MultiSpecies = 0;    // 0:0 1:6 2:9 3:12

    // Gravity parameters

    GravityOn                       = 0;    // Whether gravity is included
    GravitationalConstant           = 1.0;  // used only in SetMinimumSupport()

    //Problem specific parameter

    ProblemType = 0;

    // PPM parameters

    //--------------------------------------------------
    // parameter: Physics:cosmology:initial_redshift
    // parameter: Physics:cosmology:hubble_constant_now
    // parameter: Physics:cosmology:omega_lambda_now
    // parameter: Physics:cosmology:omega_matter_now
    // parameter: Physics:cosmology:max_expansion_rate
    // parameter: Physics:cosmology:comoving_box_size
    //--------------------------------------------------

    InitialRedshift   = parameters->value_float
      ("Physics:cosmology:initial_redshift",  20.0);

    HubbleConstantNow = parameters->value_float
      ("Physics:cosmology:hubble_constant_now",0.701);

    OmegaLambdaNow    = parameters->value_float
      ("Physics:cosmology:omega_lambda_now",   0.721);

    OmegaMatterNow    = parameters->value_float
      ("Physics:cosmology:omega_matter_now",   0.279);

    MaxExpansionRate  = parameters->value_float
      ("Physics:cosmology:max_expansion_rate", 0.01);

    ComovingBoxSize   = parameters->value_float
      ("Physics:cosmology:comoving_box_size", 64.0);

    //--------------------------------------------------
    // parameter: Method:ppm:pressure_free
    // parameter: Method:ppm:use_minimum_pressure_support
    // parameter: Method:ppm:minimum_pressure_support_parameter
    // parameter: Method:ppm:flattening
    // parameter: Method:ppm:diffusion
    // parameter: Method:ppm:steepening
    // parameter: Method:ppm:pressure_floor
    // parameter: Method:ppm:density_floor
    // parameter: Method:ppm:temperature_floor
    // parameter: Method:ppm:number_density_floor
    // parameter: Method:ppm:dual_energy
    // parameter: Method:ppm:dual_energy_eta_1
    // parameter: Method:ppm:dual_energy_eta_2
    //--------------------------------------------------

    PressureFree = parameters->value_logical
      ("Method:ppm:pressure_free",false);

    UseMinimumPressureSupport = parameters->value_logical
      ("Method:ppm:use_minimum_pressure_support",false);

    MinimumPressureSupportParameter = parameters->value_integer
      ("Method:ppm:minimum_pressure_support_parameter",100);

    PPMFlatteningParameter = parameters->value_integer
      ("Method:ppm:flattening", 3);

    PPMDiffusionParameter  = parameters->value_logical 
      ("Method:ppm:diffusion",  false);

    PPMSteepeningParameter = parameters->value_logical 
      ("Method:ppm:steepening", false);

    double floor_default = 1e-6;

    pressure_floor       = parameters->value_float
      ("Method:ppm:pressure_floor", floor_default);

    density_floor        = parameters->value_float
      ("Method:ppm:density_floor",  floor_default);

    temperature_floor    = parameters->value_float
      ("Method:ppm:temperature_floor", floor_default);

    number_density_floor = parameters->value_float
      ("Method:ppm:number_density_floor", floor_default);

    DualEnergyFormalism     = parameters->value_logical 
      ("Method:ppm:dual_energy",false);

    DualEnergyFormalismEta1 = parameters->value_float
      ("Method:ppm:dual_energy_eta_1", 0.001);

    DualEnergyFormalismEta2 = parameters->value_float
      ("Method:ppm:dual_energy_eta_2", 0.1);

    //--------------------------------------------------
    // parameter: Field::ghosts
    // parameter: Field::fields
    //--------------------------------------------------

    int gx = 1;
    int gy = 1;
    int gz = 1;

    if (parameters->type("Field:ghosts") == parameter_integer) {
      gx = gy = gz = parameters->value_integer("Field:ghosts",1);
    } else if (parameters->type("Field:ghosts") == parameter_list) {
      gx = parameters->list_value_integer(0,"Field:ghosts",1);
      gy = parameters->list_value_integer(1,"Field:ghosts",1);
      gz = parameters->list_value_integer(2,"Field:ghosts",1);
    }

    if (GridRank < 1) gx = 0;
    if (GridRank < 2) gy = 0;
    if (GridRank < 3) gz = 0;

    ghost_depth[0] = gx;
    ghost_depth[1] = gy;
    ghost_depth[2] = gz;

    NumberOfBaryonFields = parameters->list_length("Field:fields");

    // Check NumberOfBaryonFields

    if (NumberOfBaryonFields == 0) {
      ERROR ("EnzoBlock::initialize",
	     "List parameter 'Field fields' must have length greater than zero");
    } else if (NumberOfBaryonFields > MAX_NUMBER_OF_BARYON_FIELDS) {
      char buffer[ERROR_LENGTH];
      sprintf (buffer,
	       "MAX_NUMBER_OF_BARYON_FIELDS = %d is too small for %d fields",
	       NumberOfBaryonFields,NumberOfBaryonFields );
      ERROR ("EnzoBlock::initialize",  buffer);
    }

    for (int field_index=0; field_index<NumberOfBaryonFields; field_index++) {

      std::string name = 
	parameters->list_value_string(field_index,"Field:fields");

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
	char error_message[ERROR_LENGTH];
	sprintf (error_message,"Unknown field '%s'",name.c_str());
	ERROR ("EnzoBlock::EnzoBlock", error_message);
      }
    }

    // BoundaryDimension[0] = nx + 2*ghost_depth[0];
    // BoundaryDimension[1] = ny + 2*ghost_depth[1];
    // BoundaryDimension[2] = nz + 2*ghost_depth[2];

    //--------------------------------------------------
    // parameter: Domain::lower
    // parameter: Domain::upper
    //--------------------------------------------------
  
    // (ALREADY READ BY Simulation::initialize_simulation_())

    DomainLeftEdge [0] = parameters->list_value_float (0,"Domain:lower",0.0);
    DomainLeftEdge [1] = parameters->list_value_float (1,"Domain:lower",0.0);
    DomainLeftEdge [2] = parameters->list_value_float (2,"Domain:lower",0.0);

    DomainRightEdge[0] = parameters->list_value_float (0,"Domain:upper",1.0);
    DomainRightEdge[1] = parameters->list_value_float (1,"Domain:upper",1.0);
    DomainRightEdge[2] = parameters->list_value_float (2,"Domain:upper",1.0);

    //--------------------------------------------------
    // parameter: Field::courant
    //--------------------------------------------------

    CourantSafetyNumber = parameters->value_float  ("Field:courant",0.6);

    //--------------------------------------------------
    // parameter: Initial::time
    //--------------------------------------------------

    double time  = parameters->value_float  ("Initial:time",0.0);

    InitialTimeInCodeUnits = time;

  } // void initialize()
} // namespace enzo

