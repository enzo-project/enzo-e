//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2008 James Bordner
 * Copyright (C) 2008 Laboratory for Computational Astrophysics
 * Copyright (C) 2008 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */


/** 
*********************************************************************
*
* @file      initialize_color.cpp
* @brief     Initialize variables in cello_hydro.h
* @author    James Bordner
* @date      Sat Aug 29 14:20:09 PDT 2009

*
* DESCRIPTION 
* 
*    Initialize variables in cello_hydro.h.  Initial density and pressure
*    are given by an image saved using "gimp" with the ".h" format.  This
*    file is sym-linked or copied to image.h before compiling.
*
* PACKAGES
*
*    NONE
* 
* INCLUDES
*  
*    cello_hydro.h
*
* PUBLIC FUNCTIONS
*  
*    initialize_color ();
*
* PRIVATE FUCTIONS
*  
*    NONE
*
* $Id$
*
*********************************************************************
*/

#include "cello_hydro.h"
#include "assert.h"
#include "image.h"

inline float color_value 
(float * image, size_t nx, size_t ny,
 float x, float y, double enzo_lower[2], double enzo_upper[2])
// Return boolean flag whether point is inside the text "Enzo"
{
  if (x < enzo_lower[0] || x > enzo_upper[0]) return false;
  if (y < enzo_lower[1] || y > enzo_upper[1]) return false;

  size_t ix = width*(x - enzo_lower[0]) / (enzo_upper[0] - enzo_lower[0]);
  size_t iy = height*(y - enzo_lower[1]) / (enzo_upper[1] - enzo_lower[1]);
  if (ix == width) ix--;
  if (iy == height) iy--;
  assert (ix < width);
  assert (iy < height);
  return (image[ix + width*iy]);
} 

void initialize_color ()

{

  int grid_size [] = { width, height };

  float color_density_out = 1.0;
  float color_density_in  = 0.125;
  float color_pressure_out = 1.0;
  float color_pressure_in  = 0.14;
  float color_velocity_x = 10.0;
  float color_velocity_y = 0.0;
  float color_in = 1.0;
  float color_out= 0.0;

  int pixel[3];
  const char * data = header_data;

  float * image = new float [width*height];
  for (size_t iy=0; iy<height; iy++) {
    for (size_t ix=0; ix<width; ix++) {
      HEADER_PIXEL(data,pixel);
      int i=ix + width*iy;
      image [i] = 1.0*(pixel[0] + pixel[1] + pixel[2])/(255*3);
    }
  }

  // Physics

  Gamma                           = 1.4;

  // Method PPM

  PPMFlatteningParameter = 3;
  PPMDiffusionParameter  = 1;
  PPMSteepeningParameter = 1;

  // Control

  time_stop              = 2.5;
  cycle_stop             = 5000;

  CourantSafetyNumber    = 0.8;
  InitialRedshift        = 20;
  InitialTimeInCodeUnits = 0;
  Time                   = 0;
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
    CellWidth[dim] = new FLOAT[GridDimension[dim]];
    float h = double (DomainRightEdge[dim] - DomainLeftEdge[dim]) / 
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

  assert (NumberOfBaryonFields <= MAX_NUMBER_OF_BARYON_FIELDS);
  
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
//   int ndy = GridDimension[1];

//   float xd = (DomainRightEdge[0] - DomainLeftEdge[0]) ;
//   float yd = (DomainRightEdge[1] - DomainLeftEdge[1]) ;
//   int  ixg = (GridEndIndex[0] - GridStartIndex[0] + 1);
//   int  iyg = (GridEndIndex[1] - GridStartIndex[1] + 1);
  float hx = CellWidth[0][0];
  float hy = CellWidth[1][0];

  for (int iy = GridStartIndex[1]; iy<=GridEndIndex[1]; iy++) {

    float y = (iy - GridStartIndex[1] + 0.5)*hy;

    for (int ix = GridStartIndex[0]; ix<=GridEndIndex[0]; ix++) {

      float x = (ix - GridStartIndex[0] + 0.5)*hx;

      int i = ix + ndx*iy;

      // Initialize density and total energy

      float front = 0.2*(DomainRightEdge[0] - DomainLeftEdge[0]);

      if (x < front) {
	BaryonField[ field_density ] [ i ] = color_density_in;
      } else {
	BaryonField[ field_density ] [ i ] = color_density_out;
      }

      // Initialize total energy

      if (x < front) {
	BaryonField[ field_total_energy ][ i ] = 
	  color_pressure_in / ((Gamma - 1.0)*color_density_in);
      } else {
	BaryonField[ field_total_energy ][ i ] = 
	  color_pressure_out / ((Gamma - 1.0)*color_density_out);
      }

      // Initialize internal energy

      //      BaryonField[ field_internal_energy ][ i ] = 

      // Initialize velocity

      float vel_low = 0.4*(DomainRightEdge[1] - DomainLeftEdge[1]);
      float vel_high = 0.6*(DomainRightEdge[1] - DomainLeftEdge[1]);
      BaryonField[ field_velocity_x ][ i ] = 
	(y < vel_low || y > vel_high) ? 0.0 : color_velocity_x;
      BaryonField[ field_velocity_y ][ i ] = color_velocity_y;

      // Initialize color

      float a = color_value(image, width,height,x,y,
			    DomainLeftEdge,DomainRightEdge);

      float color  = a * color_in  + (1-a) * color_out;

      BaryonField[ field_color ][ i ] = color;

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
    
