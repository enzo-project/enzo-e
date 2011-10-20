// See LICENSE_CELLO file for license and copyright information

/// @file      initialize_image.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Sat Aug 29 14:20:09 PDT 2009
/// @brief     Initialize variables in cello_hydro.h
///
///    Initialize variables in cello_hydro.h.  Initial density and
///    pressure are given by an image saved using "gimp" with the ".h"
///    format.  This file is sym-linked or copied to image.h before
///    compiling.

#include "cello.hpp"

#include "enzo.hpp"


inline enzo_float color_value_png
(size_t nx, size_t ny,
 enzo_float x, enzo_float y, enzo_float lower[2], enzo_float upper[2],
 pngwriter * png)
// Return boolean flag whether point is inside the text "Enzo"
{
  if (x < lower[0] || x > upper[0]) return false;
  if (y < lower[1] || y > upper[1]) return false;

  int width  = png->getwidth();
  int height = png->getheight();

  int ix = width*(x - lower[0]) / (upper[0] - lower[0]);
  int iy = height*(y - lower[1]) / (upper[1] - lower[1]);
  iy = height - iy;
  if (ix == width) ix--;
  if (iy == height) iy--;
  ASSERT ("color_value","ix or iy out of range",ix < width && iy < height);

  return (1.0*png->read(ix+1,iy+1)/256);
} 

void EnzoBlock::initialize_image ()

{

  pngwriter png;

  png.readfromfile("input/ppm-image.png");
  int width  = png.getwidth();
  int height = png.getheight();

  
  int grid_size [] = { width, height };

  enzo_float density_out = 1.0;
  enzo_float density_in  = 0.125;
  enzo_float pressure_out = 1.0;
  enzo_float pressure_in  = 0.14;
  enzo_float velocity_x = 0.0;
  enzo_float velocity_y = 0.0;

  Gamma                           = 1.4;

  // Method PPM

  PPMFlatteningParameter = 3;
  PPMDiffusionParameter  = 1;
  PPMSteepeningParameter = 1;

  CourantSafetyNumber    = 0.8;
  InitialRedshift        = 20;
  InitialTimeInCodeUnits = 0;
  Time_                  = 0;
  OldTime                = 0;

  // Domain

  DomainLeftEdge[0]  = 0.0;
  DomainLeftEdge[1]  = 0.0;

  DomainRightEdge[0] = 1.0*width  / 2000.0;
  DomainRightEdge[1] = 1.0*height / 2000.0;

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
    enzo_float h = (DomainRightEdge[dim] - DomainLeftEdge[dim]) / 
      (GridEndIndex[dim] - GridStartIndex[dim] + 1);
    CellWidth[dim] = h;
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

  enzo_float * baryon_fields = new enzo_float [NumberOfBaryonFields * nd];
  for (int field = 0; field < NumberOfBaryonFields; field++) {
    BaryonField[field] = baryon_fields + field*nd;
  }

//   enzo_float * old_baryon_fields = new enzo_float [NumberOfBaryonFields * nd];
//   for (int field = 0; field < NumberOfBaryonFields; field++) {
//     OldBaryonField[field] = baryon_fields + field*nd;
//   }

  int ndx = GridDimension[0];
  //  int ndy = GridDimension[1];

  enzo_float hx = CellWidth[0];
  enzo_float hy = CellWidth[1];

  for (int iy = GridStartIndex[1]; iy<=GridEndIndex[1]; iy++) {

    enzo_float y = (iy - GridStartIndex[1] + 0.5)*hy;

    for (int ix = GridStartIndex[0]; ix<=GridEndIndex[0]; ix++) {

      enzo_float x = (ix - GridStartIndex[0] + 0.5)*hx;

      int i = ix + ndx*(GridEndIndex[1] - (iy-GridStartIndex[1]));

      // Initialize density and total energy

      enzo_float a;

      a = color_value_png(width,height,x,y,DomainLeftEdge,DomainRightEdge,
			  &png);

      enzo_float density  = a*density_in  + (1-a)*density_out;
      enzo_float pressure = a*pressure_in + (1-a)*pressure_out;

      BaryonField[ field_density ] [ i ] = density;
      BaryonField[ field_total_energy ][ i ] = 
	pressure / ((Gamma - 1.0)*density);

      // Initialize internal energy

      //      BaryonField[ field_internal_energy ][ i ] = 

      // Initialize velocity

      BaryonField[ field_velocity_x ][ i ] = velocity_x;
      BaryonField[ field_velocity_y ][ i ] = velocity_y;

      // Initialize color

      //      BaryonField[ field_color ][ i ] = BaryonField[ field_density ] [ i ];
      BaryonField[ field_color ][ i ] = 0;
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
	int n1 = GridDimension[(dim+2)%3];
	int n2 = GridDimension[(dim+1)%3];
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
    
