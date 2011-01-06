#include "enzo.hpp"

void EnzoDescr::write_hydro()
{
  printf ("write_hydro: ComovingCoordinates %d\n",    ComovingCoordinates);
  printf ("write_hydro: UseMinimumPressureSupport %d\n",    UseMinimumPressureSupport);
  printf ("write_hydro: MinimumPressureSupportParameter %g\n",  MinimumPressureSupportParameter);
  printf ("write_hydro: ComovingBoxSize %g\n",  ComovingBoxSize);
  printf ("write_hydro: HubbleConstantNow %g\n",  HubbleConstantNow);
  printf ("write_hydro: OmegaLambdaNow %g\n",  OmegaLambdaNow);
  printf ("write_hydro: OmegaMatterNow %g\n",  OmegaMatterNow);
  printf ("write_hydro: MaxExpansionRate %g\n",  MaxExpansionRate);

  // Chemistry

  printf ("write_hydro: MultiSpecies %d\n",    MultiSpecies);

  // Gravity

  printf ("write_hydro: GravityOn %d\n",    GravityOn);
  //  printf ("write_hydro: *AccelerationField %g\n", *AccelerationField[MAX_DIMENSION)];

  // Physics

  printf ("write_hydro: PressureFree %d\n",    PressureFree);
  printf ("write_hydro: Gamma %g\n",  Gamma);
  printf ("write_hydro: GravitationalConstant %g\n",  GravitationalConstant);

  // Problem-specific

  printf ("write_hydro: ProblemType %d\n",    ProblemType);

  // Method PPM

  printf ("write_hydro: PPMFlatteningParameter %d\n",    PPMFlatteningParameter);
  printf ("write_hydro: PPMDiffusionParameter %d\n",    PPMDiffusionParameter);
  printf ("write_hydro: PPMSteepeningParameter %d\n",    PPMSteepeningParameter);

  // Parallel

  printf ("write_hydro: ProcessorNumber %d\n",    ProcessorNumber);

  // Numerics

  printf ("write_hydro: DualEnergyFormalism %d\n",    DualEnergyFormalism);
  printf ("write_hydro: DualEnergyFormalismEta1 %g\n",  DualEnergyFormalismEta1);
  printf ("write_hydro: DualEnergyFormalismEta2 %g\n",  DualEnergyFormalismEta2);
  printf ("write_hydro: pressure %g\n",  pressure_floor);
  printf ("write_hydro: density %g\n",  density_floor);
  printf ("write_hydro: number %g\n",  number_density_floor);
  printf ("write_hydro: temperature %g\n",  temperature_floor);

  printf ("write_hydro: CourantSafetyNumber %g\n",  CourantSafetyNumber);
  printf ("write_hydro: InitialRedshift %g\n",  InitialRedshift);
  printf ("write_hydro: InitialTimeInCodeUnits %g\n",  InitialTimeInCodeUnits);
  printf ("write_hydro: Time %g\n",  Time);
  printf ("write_hydro: OldTime %g\n",  OldTime);

  // Domain

  printf ("write_hydro: DomainLeftEdge %g %g %g\n",  DomainLeftEdge [0],DomainLeftEdge [0],DomainLeftEdge [0]);
  printf ("write_hydro: %g %g %g\n",  DomainRightEdge[0],DomainRightEdge[1],DomainRightEdge[2]);

  // Grid

  printf ("write_hydro: GridRank %d\n",    GridRank);
  printf ("write_hydro: GridDimension %d %d %d\n",    GridDimension[0],GridDimension[1],GridDimension[2]);
  printf ("write_hydro: GridStartIndex %d %d %d\n",    GridStartIndex[0],GridStartIndex[1],GridStartIndex[2]);
  printf ("write_hydro: GridEndIndex %d %d %d\n",    GridEndIndex[0],GridEndIndex[1],GridEndIndex[2]);
  printf ("write_hydro: GridLeftEdge %g %g %g\n",  GridLeftEdge[0],GridLeftEdge[1],GridLeftEdge[2]);
  //  printf ("write_hydro: *CellWidth %g\n", *CellWidth[MAX_DIMENSION)];
  printf ("write_hydro: ghost %d %d %d\n",    ghost_depth[0],ghost_depth[1],ghost_depth[2]);

  // Fields

  printf ("write_hydro: field %d\n", field_density);
  printf ("write_hydro: field %d\n", field_total_energy);
  printf ("write_hydro: field %d\n", field_internal_energy);
  printf ("write_hydro: field %d\n", field_velocity_x);
  printf ("write_hydro: field %d\n", field_velocity_y);
  printf ("write_hydro: field %d\n", field_velocity_z);
  printf ("write_hydro: field %d\n", field_color);

  printf ("write_hydro: field %d\n", field_magnetic_x);
  printf ("write_hydro: field %d\n", field_magnetic_y);
  printf ("write_hydro: field %d\n", field_magnetic_z);

  printf ("write_hydro: field %d\n", field_density_xp);
  printf ("write_hydro: field %d\n", field_velocity_x_xp);
  printf ("write_hydro: field %d\n", field_velocity_y_xp);
  printf ("write_hydro: field %d\n", field_velocity_z_xp);
  printf ("write_hydro: field %d\n", field_magnetic_x_xp);
  printf ("write_hydro: field %d\n", field_magnetic_y_xp);
  printf ("write_hydro: field %d\n", field_magnetic_z_xp);

  printf ("write_hydro: field %d\n", field_density_yp);
  printf ("write_hydro: field %d\n", field_velocity_x_yp);
  printf ("write_hydro: field %d\n", field_velocity_y_yp);
  printf ("write_hydro: field %d\n", field_velocity_z_yp);
  printf ("write_hydro: field %d\n", field_magnetic_x_yp);
  printf ("write_hydro: field %d\n", field_magnetic_y_yp);
  printf ("write_hydro: field %d\n", field_magnetic_z_yp);

  printf ("write_hydro: field %d\n", field_density_zp);
  printf ("write_hydro: field %d\n", field_velocity_x_zp);
  printf ("write_hydro: field %d\n", field_velocity_y_zp);
  printf ("write_hydro: field %d\n", field_velocity_z_zp);
  printf ("write_hydro: field %d\n", field_magnetic_x_zp);
  printf ("write_hydro: field %d\n", field_magnetic_y_zp);
  printf ("write_hydro: field %d\n", field_magnetic_z_zp);


  printf ("write_hydro: NumberOfBaryonFields %d\n",    NumberOfBaryonFields);
  //  printf ("write_hydro: *BaryonField %g\n", *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS)];
  //  printf ("write_hydro: *OldBaryonField %g\n", *OldBaryonField[MAX_NUMBER_OF_BARYON_FIELDS)];
  //  printf ("write_hydro: FieldType %d\n",    FieldType[MAX_NUMBER_OF_BARYON_FIELDS)];

  printf ("write_hydro: BoundaryRank %d\n",      BoundaryRank);
  printf ("write_hydro: BoundaryDimension %d %d %d\n",      BoundaryDimension[0],BoundaryDimension[1],BoundaryDimension[2]);
  //  printf ("write_hydro: BoundaryFieldType %d\n",      BoundaryFieldType[MAX_NUMBER_OF_BARYON_FIELDS)];
  //  printf ("write_hydro: *BoundaryType bc\n", *BoundaryType[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2)];
  //  printf ("write_hydro: *BoundaryValue %g\n",   *BoundaryValue[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2)];  

  // problem

  printf ("write_hydro: CycleNumber %d\n",   CycleNumber);
  printf ("write_hydro: dt %g\n", dt);

}

//----------------------------------------------------------------------
