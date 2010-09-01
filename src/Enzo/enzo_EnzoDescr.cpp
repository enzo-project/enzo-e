// $Id: enzo_EnzoDescr.cpp 1688 2010-08-03 22:34:22Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoDescr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Aug 31 15:38:36 PDT 2010
/// @brief    Implementation of EnzoDescr methods


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
}

//----------------------------------------------------------------------

EnzoDescr::~EnzoDescr() throw ()
{
}

//----------------------------------------------------------------------
