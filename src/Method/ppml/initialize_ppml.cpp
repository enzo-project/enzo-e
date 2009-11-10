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
 * @file      initialize_ppml.cpp
 * @brief     Initialize variables in cello_hydro.h
 * @author    James Bordner
 * @date      Sat Aug 29 14:20:09 PDT 2009

 *
 * DESCRIPTION 
 * 
 *    Initialize variables in cello_hydro.h
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
 *    initialize_ppml ();
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

const bool debug = false;
 
void initialize_ppml (int size_param, int cycles_param)

{

  int grid_size [] = { size_param, size_param, size_param };

  float MHDBlastDensity = 100.0;
  float MHDBlastField = 10.0;

  float radiusBlast = 0.125;

  const float R2 (radiusBlast*radiusBlast);

  float density_out  = 1.0;
  float density_in   = MHDBlastDensity;
  float magnetic_out = MHDBlastField;
  float magnetic_in  = MHDBlastField;

  bool is_full = true;

  // Physics

  Gamma                           = 1.4;

  // Method PPM

  PPMFlatteningParameter = 3;
  PPMDiffusionParameter  = 1;
  PPMSteepeningParameter = 1;

  // Control

  time_stop              = 10000;
  cycle_stop             = cycles_param;

  CourantSafetyNumber    = 0.8;
  InitialRedshift        = 20;
  InitialTimeInCodeUnits = 0;
  Time                   = 0;
  OldTime                = 0;

  // Domain

  if (is_full) {
    DomainLeftEdge[0]  = -0.5;
    DomainLeftEdge[1]  = -0.5;
    DomainLeftEdge[2]  = -0.5;

    DomainRightEdge[0] = 0.5;
    DomainRightEdge[1] = 0.5;
    DomainRightEdge[2] = 0.5;
  } else {
    DomainLeftEdge[0]  = 0.0;
    DomainLeftEdge[1]  = 0.0;
    DomainLeftEdge[2]  = 0.0;

    DomainRightEdge[0] = 1;
    DomainRightEdge[1] = 1;
    DomainRightEdge[2] = 1;
  }

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
    CellWidth[dim] = new FLOAT[GridDimension[dim]];
    float h = (DomainRightEdge[dim] - DomainLeftEdge[dim]) / 
      (GridEndIndex[dim] - GridStartIndex[dim] + 1);
    for (int i=0; i<GridDimension[dim]; i++) {
      CellWidth[dim][i] = h;
    }
  }

  // Grid variables
  

  int k = 0;

  FieldType[field_density      = k] = k; k++;
  FieldType[field_velocity_x   = k] = k; k++;
  FieldType[field_velocity_y   = k] = k; k++;
  FieldType[field_velocity_z   = k] = k; k++;
  FieldType[field_magnetic_x   = k] = k; k++;
  FieldType[field_magnetic_y   = k] = k; k++;
  FieldType[field_magnetic_z   = k] = k; k++;

  FieldType[field_density_xp    = k] = k; k++;
  FieldType[field_velocity_x_xp = k] = k; k++;
  FieldType[field_velocity_y_xp = k] = k; k++;
  FieldType[field_velocity_z_xp = k] = k; k++;
  FieldType[field_magnetic_x_xp = k] = k; k++;
  FieldType[field_magnetic_y_xp = k] = k; k++;
  FieldType[field_magnetic_z_xp = k] = k; k++;

  FieldType[field_density_yp    = k] = k; k++;
  FieldType[field_velocity_x_yp = k] = k; k++;
  FieldType[field_velocity_y_yp = k] = k; k++;
  FieldType[field_velocity_z_yp = k] = k; k++;
  FieldType[field_magnetic_x_yp = k] = k; k++;
  FieldType[field_magnetic_y_yp = k] = k; k++;
  FieldType[field_magnetic_z_yp = k] = k; k++;

  FieldType[field_density_zp    = k] = k; k++;
  FieldType[field_velocity_x_zp = k] = k; k++;
  FieldType[field_velocity_y_zp = k] = k; k++;
  FieldType[field_velocity_z_zp = k] = k; k++;
  FieldType[field_magnetic_x_zp = k] = k; k++;
  FieldType[field_magnetic_y_zp = k] = k; k++;
  FieldType[field_magnetic_z_zp = k] = k; k++;



  //  FieldType[field_internal_energy = k++] = InternalEnergy;

  NumberOfBaryonFields = k;
  
  int nd = GridDimension[0] * GridDimension[1] * GridDimension[2];

  float * baryon_fields = new float [NumberOfBaryonFields * nd];
  for (int field = 0; field < NumberOfBaryonFields; field++) {
    BaryonField[field] = baryon_fields + field*nd;
  }

  int ndx = GridDimension[0];
  int ndy = GridDimension[1];
  int ndz = GridDimension[2];

  float hx = CellWidth[0][0];
  float hy = CellWidth[1][0];
  float hz = CellWidth[2][0];

  if (debug) printf ("Size = %d %d %d\n",ndx,ndy,ndz);

  // Clear all fields

  for (int field=0; field<NumberOfBaryonFields; field++) {
    for (int i=0; i<nd; i++) {
      BaryonField[field][i] = 0.0;
    }
  }

  // Initializize background density and magnetic_x
  for (int i=0; i<nd; i++) {
    BaryonField[field_density   ][i] = density_out;
    BaryonField[field_density_xp][i] = density_out;
    BaryonField[field_density_yp][i] = density_out;
    BaryonField[field_density_zp][i] = density_out;

    BaryonField[field_magnetic_x   ][i] = magnetic_out;
    BaryonField[field_magnetic_x_xp][i] = magnetic_out;
    BaryonField[field_magnetic_x_yp][i] = magnetic_out;
    BaryonField[field_magnetic_x_zp][i] = magnetic_out;
  }



  for (int iz = GridStartIndex[2]; iz<=GridEndIndex[2]; iz++) {

    float z0 = 0.5*hz + (iz - GridStartIndex[2]) * hz;
    float zp = 1.0*hz + (iz - GridStartIndex[2]) * hz;

    for (int iy = GridStartIndex[1]; iy<=GridEndIndex[1]; iy++) {

      float y0 = 0.5*hy + (iy - GridStartIndex[1]) * hy;
      float yp = 1.0*hy + (iy - GridStartIndex[1]) * hy;

      for (int ix = GridStartIndex[0]; ix<=GridEndIndex[0]; ix++) {

	float x0 = 0.5*hx + (ix - GridStartIndex[0]) * hx;
	float xp = 1.0*hx + (ix - GridStartIndex[0]) * hx;

	int i = ix + ndx * (iy + ndy * iz);


	float r2   = x0*x0 + y0*y0 + z0*z0;
	float rxp2 = xp*xp + y0*y0 + z0*z0;
	float ryp2 = x0*x0 + yp*yp + z0*z0;
	float rzp2 = x0*x0 + y0*y0 + zp*zp;

	// Initialize fields

	if (r2 < R2) {
	  BaryonField[field_density]    [i] = density_in;
	  BaryonField[field_magnetic_x] [i] = magnetic_in;
	} else {
	  BaryonField[field_density]    [i] = density_out;
	  BaryonField[field_magnetic_x] [i] = magnetic_in;
	}

	if (rxp2 < R2) {
	  BaryonField[field_density_xp]    [i] = density_in;
	  BaryonField[field_magnetic_x_xp] [i] = magnetic_in;
	} else {
	  BaryonField[field_density_xp]    [i] = density_out;
	  BaryonField[field_magnetic_x_xp] [i] = magnetic_in;
	}

	if (ryp2 < R2) {
	  BaryonField[field_density_yp]    [i] = density_in;
	  BaryonField[field_magnetic_x_yp] [i] = magnetic_in;
	} else {
	  BaryonField[field_density_yp]    [i] = density_out;
	  BaryonField[field_magnetic_x_yp] [i] = magnetic_in;
	}

	if (rzp2 < R2) {
	  BaryonField[field_density_zp]    [i] = density_in;
	  BaryonField[field_magnetic_x_zp] [i] = magnetic_in;
	} else {
	  BaryonField[field_density_zp]    [i] = density_out;
	  BaryonField[field_magnetic_x_zp] [i] = magnetic_in;
	}
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
	BoundaryType [field][dim][face] = new bc_type [size];
	BoundaryValue[field][dim][face] = NULL;
	for (int i2 = 0; i2<n2; i2++) {
	  for (int i1 = 0; i1<n1; i1++) {
	    int i = i1 + n1*i2;
	    BoundaryType[field][dim][face][i] = bc_periodic;
	  }
	}
      }
    }
  }

}
    
