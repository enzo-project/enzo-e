
  //----------------------------------------------------------------------
  // functions

  static void initialize (Parameters * parameters, FieldDescr *);

  //----------------------------------------------------------------------
  // variables

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

  static int field_density;
  static int field_total_energy;
  static int field_internal_energy;
  static int field_velocity_x;
  static int field_velocity_y;
  static int field_velocity_z;

  static int field_color;

  // PPM

  static int field_velox;
  static int field_veloy;
  static int field_veloz;
  static int field_bfieldx;
  static int field_bfieldy;
  static int field_bfieldz;

  static int field_dens_rx;
  static int field_velox_rx;
  static int field_veloy_rx;
  static int field_veloz_rx;
  static int field_bfieldx_rx;
  static int field_bfieldy_rx;
  static int field_bfieldz_rx;

  static int field_dens_ry;
  static int field_velox_ry;
  static int field_veloy_ry;
  static int field_veloz_ry;
  static int field_bfieldx_ry;
  static int field_bfieldy_ry;
  static int field_bfieldz_ry;

  static int field_dens_rz;
  static int field_velox_rz;
  static int field_veloy_rz;
  static int field_veloz_rz;
  static int field_bfieldx_rz;
  static int field_bfieldy_rz;
  static int field_bfieldz_rz;


  static int GridRank;

  static int ghost_depth[MAX_DIMENSION];

  // Fields

  static int NumberOfBaryonFields;      // active baryon fields

  static int FieldType[MAX_NUMBER_OF_BARYON_FIELDS];

