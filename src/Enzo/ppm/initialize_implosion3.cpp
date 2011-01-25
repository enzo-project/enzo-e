// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      initialize_implosion3.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Sat Aug 29 14:20:09 PDT 2009
/// @brief     Initialize variables in cello_hydro.h

#include "cello.hpp"

#include "enzo.hpp"

const bool debug = false;
 
void EnzoDescr::initialize_implosion3 (int size_param)

{

  int grid_size [] = { size_param, size_param, size_param };

  float density_out = 1.0;
  float density_in  = 0.125;
  float pressure_out = 1.0;
  float pressure_in  = 0.14;

  // Physics

  Gamma                           = 1.4;

  // Method PPM

  PPMFlatteningParameter = 3;
  PPMDiffusionParameter  = 1;
  PPMSteepeningParameter = 1;

  CourantSafetyNumber    = 0.8;
  InitialRedshift        = 20;
  InitialTimeInCodeUnits = 0;
  Time                   = 0;
  OldTime                = 0;

  // Domain

  DomainLeftEdge[0]  = 0.0;
  DomainLeftEdge[1]  = 0.0;
  DomainLeftEdge[2]  = 0.0;

  DomainRightEdge[0] = 0.3;
  DomainRightEdge[1] = 0.3;
  DomainRightEdge[2] = 0.3;

  // Grid

  GridRank           =  3;
  GridDimension[0]   = grid_size[0] + 6;
  GridDimension[1]   = grid_size[1] + 6;
  GridDimension[2]   = grid_size[2] + 6;
  GridStartIndex[0]  = 3;
  GridStartIndex[1]  = 3;
  GridStartIndex[2]  = 3;
  GridEndIndex[0]    = grid_size[0] + 3 - 1;
  GridEndIndex[1]    = grid_size[1] + 3 - 1;
  GridEndIndex[2]    = grid_size[2] + 3 - 1;
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

  FieldType[field_density         = k++] = Density;
  FieldType[field_total_energy    = k++] = TotalEnergy;
  FieldType[field_velocity_x      = k++] = Velocity1;
  FieldType[field_velocity_y      = k++] = Velocity2;
  FieldType[field_velocity_z      = k++] = Velocity3;
  FieldType[field_color        = k++] = ElectronDensity;

  //  FieldType[field_internal_energy = k++] = InternalEnergy;

  NumberOfBaryonFields = k;
  
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
  int ndz = GridDimension[2];

  float xd = (DomainRightEdge[0] - DomainLeftEdge[0]) ;
  float yd = (DomainRightEdge[1] - DomainLeftEdge[1]) ;
  float zd = (DomainRightEdge[2] - DomainLeftEdge[2]) ;
  int  ixg = (GridEndIndex[0] - GridStartIndex[0] + 1);
  int  iyg = (GridEndIndex[1] - GridStartIndex[1] + 1);
  int  izg = (GridEndIndex[2] - GridStartIndex[2] + 1);
  float hx = CellWidth[0][0];
  float hy = CellWidth[1][0];
  float hz = CellWidth[2][0];

  if (debug) printf ("Size = %d %d %d\n",ndx,ndy,ndz);
  if (debug) printf ("%g  %g %g  %g %g\n",
	  Gamma, 
	  pressure_out,density_out,
	  pressure_in,density_in);
  if (debug) printf ("total energy: %g %g\n",
	  pressure_out / ((Gamma - 1.0)*density_out),
	  pressure_in / ((Gamma - 1.0)*density_in));

  for (int iz = GridStartIndex[2]; iz<=GridEndIndex[2]; iz++) {

    float z = 0.5*hz + (iz - GridStartIndex[2]) * zd / izg;

    for (int iy = GridStartIndex[1]; iy<=GridEndIndex[1]; iy++) {

      float y = 0.5*hy + (iy - GridStartIndex[1]) * yd / iyg;

      for (int ix = GridStartIndex[0]; ix<=GridEndIndex[0]; ix++) {

	float x = 0.5*hx + (ix - GridStartIndex[0]) * xd / ixg;

	int i = ix + ndx * (iy + ndy * iz);

	// Initialize density

	if (x + y + z < 0.1517) {
	  BaryonField[ field_density ] [ i ] = density_in;
	} else {
	  BaryonField[ field_density ] [ i ] = density_out;
	}

	// Initialize total energy

	if (x + y + z < 0.1517) {
	  BaryonField[ field_total_energy ][ i ] = 
	    pressure_in / ((Gamma - 1.0)*density_in);
	} else {
	  BaryonField[ field_total_energy ][ i ] = 
	    pressure_out / ((Gamma - 1.0)*density_out);
	}

	// Initialize internal energy

	//      BaryonField[ field_internal_energy ][ i ] = 

	// Initialize velocity

	BaryonField[ field_velocity_x ][ i ] = 0;
	BaryonField[ field_velocity_y ][ i ] = 0;
	BaryonField[ field_velocity_z ][ i ] = 0;

	// Initialize color

	BaryonField[ field_color ][ i ] = 0.0;

      }
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

  BoundaryRank = 3;
  BoundaryDimension[0] = GridDimension[0];
  BoundaryDimension[1] = GridDimension[1];
  BoundaryDimension[2] = GridDimension[2];

  for (int field=0; field<NumberOfBaryonFields; field++) {
    BoundaryFieldType[field] = FieldType[field];
    for (int dim = 0; dim < 3; dim++) {
      for (int face = 0; face < 2; face++) {
	int n1 = GridDimension[(dim+1)%3];
	int n2 = GridDimension[(dim+2)%3];
	int size = n1*n2;
	BoundaryType [field][dim][face] = new bc_enum [size];
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
    
