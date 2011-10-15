// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoNamespace.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-07-31
/// @todo     Comment Enzo functions / attributes
/// @todo     Move boundary functionality to Block
/// @brief    [\ref Enzo] "Global" data for Enzo

#ifndef ENZO_ENZO_NAMESPACE_HPP
#define ENZO_ENZO_NAMESPACE_HPP

enum bc_enum 
  { // explicitly enumerated to match what Enzo expects
    bc_unknown    = 0, 
    bc_reflecting = 1, 
    bc_outflow    = 2, 
    bc_inflow     = 3, 
    bc_periodic   = 4 
  };

namespace enzo {

  //----------------------------------------------------------------------
  // functions

  extern void initialize (Parameters * parameters);

  //----------------------------------------------------------------------
  // variables

  /// Boundary

  extern int  BoundaryRank;
  extern int  BoundaryDimension[MAX_DIMENSION];
  extern int  BoundaryFieldType[MAX_NUMBER_OF_BARYON_FIELDS];
  extern bc_enum *BoundaryType[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2];
  extern enzo_float *BoundaryValue[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2]; 

  extern int ComovingCoordinates;
  extern int UseMinimumPressureSupport;
  extern enzo_float MinimumPressureSupportParameter;
  extern enzo_float ComovingBoxSize;
  extern enzo_float HubbleConstantNow;
  extern enzo_float OmegaMatterNow;
  extern enzo_float OmegaLambdaNow;
  extern enzo_float MaxExpansionRate;

  // Chemistry

  extern int MultiSpecies;

  // Gravity

  extern int GravityOn;

  // Physics

  extern int PressureFree;
  extern enzo_float Gamma;
  extern enzo_float GravitationalConstant;

  // Problem-specific

  extern int ProblemType;

  // Method PPM

  extern int PPMFlatteningParameter;
  extern int PPMDiffusionParameter;
  extern int PPMSteepeningParameter;

  // Parallel

  //  extern int ProcessorNumber;

  // Numerics

  extern int DualEnergyFormalism;
  extern enzo_float DualEnergyFormalismEta1;
  extern enzo_float DualEnergyFormalismEta2;

  extern enzo_float pressure_floor;
  extern enzo_float density_floor;
  extern enzo_float number_density_floor;
  extern enzo_float temperature_floor;

  extern enzo_float CourantSafetyNumber;
  extern enzo_float InitialRedshift;
  extern enzo_float InitialTimeInCodeUnits;

  // Domain

  extern enzo_float DomainLeftEdge [MAX_DIMENSION];
  extern enzo_float DomainRightEdge[MAX_DIMENSION];

  // Grid

  extern int field_density;
  extern int field_total_energy;
  extern int field_internal_energy;
  extern int field_velocity_x;
  extern int field_velocity_y;
  extern int field_velocity_z;
  extern int field_color;
  extern int field_magnetic_x;
  extern int field_magnetic_y;
  extern int field_magnetic_z;

  extern int field_density_xp;
  extern int field_velocity_x_xp;
  extern int field_velocity_y_xp;
  extern int field_velocity_z_xp;
  extern int field_magnetic_x_xp;
  extern int field_magnetic_y_xp;
  extern int field_magnetic_z_xp;

  extern int field_density_yp;
  extern int field_velocity_x_yp;
  extern int field_velocity_y_yp;
  extern int field_velocity_z_yp;
  extern int field_magnetic_x_yp;
  extern int field_magnetic_y_yp;
  extern int field_magnetic_z_yp;

  extern int field_density_zp;
  extern int field_velocity_x_zp;
  extern int field_velocity_y_zp;
  extern int field_velocity_z_zp;
  extern int field_magnetic_x_zp;
  extern int field_magnetic_y_zp;
  extern int field_magnetic_z_zp;

  extern int GridRank;

  extern int ghost_depth[MAX_DIMENSION];

  // Fields

  extern int NumberOfBaryonFields;      // active baryon fields

  extern int FieldType[MAX_NUMBER_OF_BARYON_FIELDS];

};

using namespace enzo;

#endif /* ENZO_ENZO_NAMESPACE_HPP */

