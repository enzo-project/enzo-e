
/// @namespace  enzo
/// @brief      Namespace for Enzo functions and global variables

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

// Parallel

//  int ProcessorNumber;

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

int EnzoBlock::field_density;
int EnzoBlock::field_total_energy;
int EnzoBlock::field_internal_energy;
int EnzoBlock::field_velocity_x;
int EnzoBlock::field_velocity_y;
int EnzoBlock::field_velocity_z;

int EnzoBlock::field_color;

// PPM

int EnzoBlock::field_velox;
int EnzoBlock::field_veloy;
int EnzoBlock::field_veloz;
int EnzoBlock::field_bfieldx;
int EnzoBlock::field_bfieldy;
int EnzoBlock::field_bfieldz;

int EnzoBlock::field_dens_rx;
int EnzoBlock::field_velox_rx;
int EnzoBlock::field_veloy_rx;
int EnzoBlock::field_veloz_rx;
int EnzoBlock::field_bfieldx_rx;
int EnzoBlock::field_bfieldy_rx;
int EnzoBlock::field_bfieldz_rx;

int EnzoBlock::field_dens_ry;
int EnzoBlock::field_velox_ry;
int EnzoBlock::field_veloy_ry;
int EnzoBlock::field_veloz_ry;
int EnzoBlock::field_bfieldx_ry;
int EnzoBlock::field_bfieldy_ry;
int EnzoBlock::field_bfieldz_ry;

int EnzoBlock::field_dens_rz;
int EnzoBlock::field_velox_rz;
int EnzoBlock::field_veloy_rz;
int EnzoBlock::field_veloz_rz;
int EnzoBlock::field_bfieldx_rz;
int EnzoBlock::field_bfieldy_rz;
int EnzoBlock::field_bfieldz_rz;

int EnzoBlock::GridRank;

int EnzoBlock::ghost_depth[MAX_DIMENSION];

// Fields

int EnzoBlock::NumberOfBaryonFields;      // active baryon fields

int EnzoBlock::FieldType[MAX_NUMBER_OF_BARYON_FIELDS];

void EnzoBlock::initialize(Parameters * parameters,
			   FieldDescr * field_descr)

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

  // PPM

  field_density = field_undefined;
  field_total_energy = field_undefined;
  field_internal_energy = field_undefined;
  field_velocity_x = field_undefined;
  field_velocity_y = field_undefined;
  field_velocity_z = field_undefined;

  field_color = field_undefined;

  // PPML

  field_velox = field_undefined;
  field_veloy = field_undefined;
  field_veloz = field_undefined;
  field_bfieldx = field_undefined;
  field_bfieldy = field_undefined;
  field_bfieldz = field_undefined;

  field_dens_rx = field_undefined;
  field_velox_rx = field_undefined;
  field_veloy_rx = field_undefined;
  field_veloz_rx = field_undefined;
  field_bfieldx_rx = field_undefined;
  field_bfieldy_rx = field_undefined;
  field_bfieldz_rx = field_undefined;
  field_dens_ry = field_undefined;

  field_velox_ry = field_undefined;
  field_veloy_ry = field_undefined;
  field_veloz_ry = field_undefined;
  field_bfieldx_ry = field_undefined;
  field_bfieldy_ry = field_undefined;
  field_bfieldz_ry = field_undefined;

  field_dens_rz = field_undefined;
  field_velox_rz = field_undefined;
  field_veloy_rz = field_undefined;
  field_veloz_rz = field_undefined;
  field_bfieldx_rz = field_undefined;
  field_bfieldy_rz = field_undefined;
  field_bfieldz_rz = field_undefined;

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
  // parameter: Mesh : root_rank
  //--------------------------------------------------

  ComovingCoordinates = parameters->value_logical ("Physics:cosmology",false);
  Gamma               = parameters->value_float   ("Physics:gamma",5.0/3.0);
  GridRank            = parameters->value_integer ("Mesh:root_rank",0);

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
  // parameter: Physics : cosmology : initial_redshift
  // parameter: Physics : cosmology : hubble_constant_now
  // parameter: Physics : cosmology : omega_lambda_now
  // parameter: Physics : cosmology : omega_matter_now
  // parameter: Physics : cosmology : max_expansion_rate
  // parameter: Physics : cosmology : comoving_box_size
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
  // parameter: Enzo : ppm : pressure_free
  // parameter: Enzo : ppm : use_minimum_pressure_support
  // parameter: Enzo : ppm : minimum_pressure_support_parameter
  // parameter: Enzo : ppm : flattening
  // parameter: Enzo : ppm : diffusion
  // parameter: Enzo : ppm : steepening
  // parameter: Enzo : ppm : pressure_floor
  // parameter: Enzo : ppm : density_floor
  // parameter: Enzo : ppm : temperature_floor
  // parameter: Enzo : ppm : number_density_floor
  // parameter: Enzo : ppm : dual_energy
  // parameter: Enzo : ppm : dual_energy_eta_1
  // parameter: Enzo : ppm : dual_energy_eta_2
  //--------------------------------------------------

  PressureFree = parameters->value_logical
    ("Enzo:ppm:pressure_free",false);

  UseMinimumPressureSupport = parameters->value_logical
    ("Enzo:ppm:use_minimum_pressure_support",false);

  MinimumPressureSupportParameter = parameters->value_integer
    ("Enzo:ppm:minimum_pressure_support_parameter",100);

  PPMFlatteningParameter = parameters->value_integer
    ("Enzo:ppm:flattening", 3);

  PPMDiffusionParameter  = parameters->value_logical 
    ("Enzo:ppm:diffusion",  false);

  PPMSteepeningParameter = parameters->value_logical 
    ("Enzo:ppm:steepening", false);

  double floor_default = 1e-6;

  pressure_floor       = parameters->value_float
    ("Enzo:ppm:pressure_floor", floor_default);

  density_floor        = parameters->value_float
    ("Enzo:ppm:density_floor",  floor_default);

  temperature_floor    = parameters->value_float
    ("Enzo:ppm:temperature_floor", floor_default);

  number_density_floor = parameters->value_float
    ("Enzo:ppm:number_density_floor", floor_default);

  DualEnergyFormalism     = parameters->value_logical 
    ("Enzo:ppm:dual_energy",false);

  DualEnergyFormalismEta1 = parameters->value_float
    ("Enzo:ppm:dual_energy_eta_1", 0.001);

  DualEnergyFormalismEta2 = parameters->value_float
    ("Enzo:ppm:dual_energy_eta_2", 0.1);

  //--------------------------------------------------
  // parameter: Field : ghosts
  // parameter: Field : fields
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
    ERROR2 ("EnzoBlock::initialize",
	    "MAX_NUMBER_OF_BARYON_FIELDS = %d is too small for %d fields",
	    MAX_NUMBER_OF_BARYON_FIELDS,NumberOfBaryonFields );
  }

  // [ Ideally this should not be a loop ]

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
    } else if (name == "velox") {
      field_velox            = field_index;
      FieldType[field_index] = Velocity1;
    } else if (name == "veloy") {
      field_veloy           = field_index;
      FieldType[field_index] = Velocity2;
    } else if (name == "veloz") {
      field_veloz           = field_index;
      FieldType[field_index] = Velocity3;
    } else if (name == "bfieldx") {
      field_bfieldx           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "bfieldy") {
      field_bfieldy           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "bfieldz") {
      field_bfieldz           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "dens_rx") {
      field_dens_rx           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "velox_rx") {
      field_velox_rx           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "veloy_rx") {
      field_veloy_rx           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "veloz_rx") {
      field_veloz_rx           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "bfieldx_rx") {
      field_bfieldx_rx           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "bfieldy_rx") {
      field_bfieldy_rx           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "bfieldz_rx") {
      field_bfieldz_rx           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "dens_ry") {
      field_dens_ry           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "velox_ry") {
      field_velox_ry           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "veloy_ry") {
      field_veloy_ry           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "veloz_ry") {
      field_veloz_ry           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "bfieldx_ry") {
      field_bfieldx_ry           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "bfieldy_ry") {
      field_bfieldy_ry           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "bfieldz_ry") {
      field_bfieldz_ry           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "dens_rz") {
      field_dens_rz           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "velox_rz") {
      field_velox_rz           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "veloy_rz") {
      field_veloy_rz           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "veloz_rz") {
      field_veloz_rz           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "bfieldx_rz") {
      field_bfieldx_rz           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "bfieldy_rz") {
      field_bfieldy_rz           = field_index;
      FieldType[field_index] = 0;
    } else if (name == "bfieldz_rz") {
      field_bfieldz_rz           = field_index;
      FieldType[field_index] = 0;
    } else {
      FieldType[field_index] = 0;
      WARNING1 ("EnzoBlock::EnzoBlock", 
		"Unknown field type for field %s",
		name.c_str());
    }
  }

  // BoundaryDimension[0] = nx + 2*ghost_depth[0];
  // BoundaryDimension[1] = ny + 2*ghost_depth[1];
  // BoundaryDimension[2] = nz + 2*ghost_depth[2];

  //--------------------------------------------------
  // parameter: Domain : lower
  // parameter: Domain : upper
  //--------------------------------------------------
  
  // (ALREADY READ BY Simulation::initialize_simulation_())

  DomainLeftEdge [0] = parameters->list_value_float (0,"Domain:lower",0.0);
  DomainLeftEdge [1] = parameters->list_value_float (1,"Domain:lower",0.0);
  DomainLeftEdge [2] = parameters->list_value_float (2,"Domain:lower",0.0);

  DomainRightEdge[0] = parameters->list_value_float (0,"Domain:upper",1.0);
  DomainRightEdge[1] = parameters->list_value_float (1,"Domain:upper",1.0);
  DomainRightEdge[2] = parameters->list_value_float (2,"Domain:upper",1.0);

  //--------------------------------------------------
  // parameter: Field : courant
  //--------------------------------------------------

  CourantSafetyNumber = parameters->value_float  ("Field:courant",0.6);

  //--------------------------------------------------
  // parameter: Initial : time
  //--------------------------------------------------

  double time  = parameters->value_float  ("Initial:time",0.0);

  InitialTimeInCodeUnits = time;

} // void initialize()

