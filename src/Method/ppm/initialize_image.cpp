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
* @file      initialize_image.cpp
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
*    initialize_image ();
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
(float * image, int nx, int ny,
 float x, float y, float enzo_lower[2], float enzo_upper[2])
// Return boolean flag whether point is inside the text "Enzo"
{
  if (x < enzo_lower[0] || x > enzo_upper[0]) return false;
  if (y < enzo_lower[1] || y > enzo_upper[1]) return false;

  int ix = width*(x - enzo_lower[0]) / (enzo_upper[0] - enzo_lower[0]);
  int iy = height*(y - enzo_lower[1]) / (enzo_upper[1] - enzo_lower[1]);
  if (ix == width) ix--;
  if (iy == height) iy--;
  assert (ix >= 0);
  assert (iy >= 0);
  assert (ix < width);
  assert (iy < height);
  return (image[ix + width*iy]);
} 

void initialize_image ()

{

  int grid_size [] = { width, height };

  float enzo_density_out = 1.0;
  float enzo_density_in  = 0.125;
  float enzo_pressure_out = 1.0;
  float enzo_pressure_in  = 0.14;
  float enzo_velocity_x = 0.0;
  float enzo_velocity_y = 0.0;
  int pixel[3];
  char * data = header_data;

  float * image = new float [width*height];
  for (int iy=0; iy<height; iy++) {
    for (int ix=0; ix<width; ix++) {
      HEADER_PIXEL(data,pixel);
      int i=ix + width*iy;
      image [i] = 1.0*(pixel[0] + pixel[1] + pixel[2])/(255*3);
    }
  }

  // Set extents of the E

  // 5/p   @ @ @
  // 4/4   @
  // 3/4   @ @
  // 2/4   @
  // 1/4   @ @ @
  //  0 

  float enzo_lower[2] = {0.0, 0.0};
  float enzo_upper[2] = {2.0, 2.0};

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

  DomainRightEdge[0] = 4;
  DomainRightEdge[1] = 2;

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

  FieldType[field_density         = k++] = Density;
  FieldType[field_total_energy    = k++] = TotalEnergy;
  FieldType[field_velocity_x      = k++] = Velocity1;
  FieldType[field_velocity_y      = k++] = Velocity2;
  //  FieldType[field_internal_energy = k++] = InternalEnergy;

  NumberOfBaryonFields = k;
  
  int nd = GridDimension[0] * GridDimension[1] * GridDimension[2];

  float * baryon_fields = new float [NumberOfBaryonFields * nd];
  for (int field = 0; field < NumberOfBaryonFields; field++) {
    BaryonField[field] = baryon_fields + field*nd;
  }

  float * old_baryon_fields = new float [NumberOfBaryonFields * nd];
  for (int field = 0; field < NumberOfBaryonFields; field++) {
    OldBaryonField[field] = baryon_fields + field*nd;
  }

  int ndx = GridDimension[0];
  int ndy = GridDimension[1];

  float xd = (DomainRightEdge[0] - DomainLeftEdge[0]) ;
  float yd = (DomainRightEdge[1] - DomainLeftEdge[1]) ;
  int  ixg = (GridEndIndex[0] - GridStartIndex[0] + 1);
  int  iyg = (GridEndIndex[1] - GridStartIndex[1] + 1);
  float hx = CellWidth[0][0];
  float hy = CellWidth[1][0];

  for (int iy = GridStartIndex[1]; iy<=GridEndIndex[1]; iy++) {

    float y = (iy - GridStartIndex[1] + 0.5)*hy;

    for (int ix = GridStartIndex[0]; ix<=GridEndIndex[0]; ix++) {

      float x = (ix - GridStartIndex[0] + 0.5)*hx;

      int i = ix + ndx*iy;

      // Initialize density and total energy

      float a = color_value(image, width,height,x,y,enzo_lower,enzo_upper);
      float density  = a*enzo_density_in  + (1-a)*enzo_density_out;
      float pressure = a*enzo_pressure_in + (1-a)*enzo_pressure_out;

      BaryonField[ field_density ] [ i ] = density;
      BaryonField[ field_total_energy ][ i ] = 
	pressure / ((Gamma - 1.0)*density);

      // Initialize internal energy

      //      BaryonField[ field_internal_energy ][ i ] = 

      // Initialize velocity

      BaryonField[ field_velocity_x ][ i ] = enzo_velocity_x;
      BaryonField[ field_velocity_y ][ i ] = enzo_velocity_y;

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
    
