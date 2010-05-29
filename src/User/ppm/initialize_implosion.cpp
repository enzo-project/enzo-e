// $Id$
// See LICENSE_ENZO file for license and copyright information

/// @file      initialize_implosion.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Sat Aug 29 14:20:09 PDT 2009
/// @brief     Initialize variables in cello_hydro.h

#include "cello_hydro.h"

const bool debug = false;
 
void initialize_implosion (int size_param, int cycles_param)

{

  int grid_size [] = { size_param, size_param };
  float density_out = 1.0;
  float density_in  = 0.125;
  float pressure_out = 1.0;
  float pressure_in  = 0.14;
  float velocity_x = 0.0;
  float velocity_y = 0.0;
  float cutoff = 0.1517;

  // Physics

  Gamma                           = 1.4;

  // Method PPM

  PPMFlatteningParameter = 3;
  PPMDiffusionParameter  = 1;
  PPMSteepeningParameter = 1;

  // Control

  time_stop              = 2.5;
  cycle_stop             = cycles_param;

  CourantSafetyNumber    = 0.8;
  InitialRedshift        = 20;
  InitialTimeInCodeUnits = 0;
  Time                   = 0;
  OldTime                = 0;

  // Domain

  DomainLeftEdge[0]  = 0.0;
  DomainLeftEdge[1]  = 0.0;

  DomainRightEdge[0] = 0.3;
  DomainRightEdge[1] = 0.3;

  // Grid

  GridRank           =  2;
  GridDimension[0]   = grid_size[0] + 6;
  GridDimension[1]   = grid_size[1] + 6;
  GridDimension[2]   = 1;
  GridStartIndex[0]  = 3;
  GridStartIndex[1]  = 3;
  GridStartIndex[2]  = 0;
  GridEndIndex[0]    = grid_size[0] + 3 - 1;
  GridEndIndex[1]    = grid_size[1] + 3 - 1;
  GridEndIndex[2]    = 0;
  GridLeftEdge[0]    = 0.0;
  GridLeftEdge[1]    = 0.0;
  GridLeftEdge[2]    = 0.0;

  for (int dim=0; dim<GridRank; dim++) {
    CellWidth[dim] = new ENZO_FLOAT[GridDimension[dim]];
    float h = (DomainRightEdge[dim] - DomainLeftEdge[dim]) / 
      (GridEndIndex[dim] - GridStartIndex[dim] + 1);
    for (int i=0; i<GridDimension[dim]; i++) {
      CellWidth[dim][i] = h;
    }
  }

  // Grid variables
  

  int k = 0;

  FieldType[field_density      = k++] = Density;
  FieldType[field_total_energy = k++] = TotalEnergy;
  FieldType[field_velocity_x   = k++] = Velocity1;
  FieldType[field_velocity_y   = k++] = Velocity2;
  FieldType[field_color        = k++] = ElectronDensity;

  NumberOfBaryonFields = k;

  ASSERT ("initialize_implosion",
	  "MAX_NUMBER_OF_BARYON_FIELDS is too small",
	  NumberOfBaryonFields <= MAX_NUMBER_OF_BARYON_FIELDS);
  
  int nd = GridDimension[0] * GridDimension[1] * GridDimension[2];

  float * baryon_fields = new float [NumberOfBaryonFields * nd];
  for (int field = 0; field < NumberOfBaryonFields; field++) {
    BaryonField[field] = baryon_fields + field*nd;
  }

//   float * old_baryon_fields = new float [NumberOfBaryonFields * nd];
//   for (int field = 0; field < NumberOfBaryonFields; field++) {
//     OldBaryonField[field] = baryon_fields + field*nd;
//   }

  int ndx = GridDimension[0];
  int ndy = GridDimension[1];

  float hx = CellWidth[0][0];
  float hy = CellWidth[1][0];

  for (int iy = GridStartIndex[1]; iy<=GridEndIndex[1]; iy++) {

    float y = (iy - GridStartIndex[1] + 0.5)*hy;

    for (int ix = GridStartIndex[0]; ix<=GridEndIndex[0]; ix++) {

      float x = (ix - GridStartIndex[0] + 0.5)*hx;

      int i = ix + ndx * iy;

      // Initialize density

      if (x + y < cutoff) {
	BaryonField[ field_density ] [ i ] = density_in;
      } else {
	BaryonField[ field_density ] [ i ] = density_out;
      }

      // Initialize total energy

      if (x + y < cutoff) {
	BaryonField[ field_total_energy ][ i ] = 
	  pressure_in / ((Gamma - 1.0)*density_in);
      } else {
	BaryonField[ field_total_energy ][ i ] = 
	  pressure_out / ((Gamma - 1.0)*density_out);
      }

      // Initialize internal energy

      //      BaryonField[ field_internal_energy ][ i ] = 

      // Initialize velocity

      BaryonField[ field_velocity_x ][ i ] = velocity_x;
      BaryonField[ field_velocity_y ][ i ] = velocity_y;

      // Initialize color

      BaryonField[ field_color ][ i ] = 0.0;

    }
  }


  AccelerationField[0] = NULL;
  AccelerationField[1] = NULL;
  AccelerationField[2] = NULL;

  // Numerics

  DualEnergyFormalism             = 0;
  DualEnergyFormalismEta1         = 0.001;
  DualEnergyFormalismEta2         = 0.1;
  pressure_floor                  = 1e-6;
  number_density_floor            = 1e-6;
  density_floor                   = 1e-6;
  temperature_floor               = 1e-6;

  // boundary
 
  BoundaryRank = 2;
  BoundaryDimension[0] = GridDimension[0];
  BoundaryDimension[1] = GridDimension[1];

  for (int field=0; field<NumberOfBaryonFields; field++) {
    BoundaryFieldType[field] = FieldType[field];
    for (int dim = 0; dim < 3; dim++) {
      for (int face = 0; face < 2; face++) {
	int n1 = GridDimension[(dim+1)%3];
	int n2 = GridDimension[(dim+2)%3];
	int size = n1*n2;
	BoundaryType [field][dim][face] = new bc_type [size];
	BoundaryValue[field][dim][face] = NULL;
	for (int i2 = 0; i2<n2; i2++) {
	  for (int i1 = 0; i1<n1; i1++) {
	    int i = i1 + n1*i2;
	    BoundaryType[field][dim][face][i] = bc_reflecting;
	  }
	}
      }
    }
  }

}
    
