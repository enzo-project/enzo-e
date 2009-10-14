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
 * @file      initialize_hydro.cpp
 * @brief     Initialize variables in cello_hydro.h
 * @author    James Bordner
 * @date      Sun Aug 23 12:22:10 PDT 2009
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
 *    initialize_hydro ();
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
 
void initialize_hydro ()

{

  // Cosmology

  ComovingCoordinates             = 0;    // Physics: Cosmology
  UseMinimumPressureSupport       = 0;    // call UseMinimumPressureSupport() ?
  MinimumPressureSupportParameter = 100;  // SetMinimumSupport() Enzo parameter
  ComovingBoxSize                 = 64;   // Physics cosmology: Mpc/h at z=0
  HubbleConstantNow               = 0.5;  // Physics: cosmology parameter
  OmegaLambdaNow                  = 0.0;  // Physics: cosmology parameter
  OmegaMatterNow                  = 1.0;  // Physics: cosmology parameter
  MaxExpansionRate                = 0.01; // Cosmology timestep constraint

  // Chemistry

  MultiSpecies                    = 0;    // 0:0 1:6 2:9 3:12

  // Gravity

  GravityOn                       = 0;    // Whether gravity is included
  AccelerationField[0]            = NULL;
  AccelerationField[1]            = NULL;
  AccelerationField[2]            = NULL;

  // Physics

  PressureFree                    = 0;    // Physics: passed into ppm_de
  Gamma                           = 5.0/3.0; // Physics: ideal gas law constant
  GravitationalConstant           = 1.0;  // used only in SetMinimumSupport()

  // Problem-specific

  ProblemType                     = 0;    //  7 (sedov) or 60 61 (turbulence)

  // Method PPM

  PPMFlatteningParameter          = 0;
  PPMDiffusionParameter           = 0;
  PPMSteepeningParameter          = 0;

  // Parallel

  ProcessorNumber = 0;

  // Numerics

  DualEnergyFormalism             = 0;    // Method: PPM parameter
  DualEnergyFormalismEta1         = 0.001;// Method: PPM parameter
  DualEnergyFormalismEta2         = 0.1;  // Method: PPM parameter
  pressure_floor                  = 1e-20; // Was "tiny_number"
  number_density_floor            = 1e-20; // Was "tiny_number"
  density_floor                   = 1e-20; // Was "tiny_number"
  temperature_floor               = 1e-20; // Was "tiny_number"

  // Boundary

  BoundaryRank         = 0;
  BoundaryDimension[0] = 1;
  BoundaryDimension[1] = 1;
  BoundaryDimension[2] = 1;
 
  /* Clear BoundaryType and BoundaryValue pointers. */
 
  for (int field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++) {
    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
      for (int i = 0; i < 2; i++) {
	BoundaryType[field][dim][i] = NULL;
	BoundaryValue[field][dim][i] = NULL;
      }
    }
  }
}
    
