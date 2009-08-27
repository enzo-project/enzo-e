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

  // Global parameters

  ComovingCoordinates             = 0;    // Physics: Cosmology
  DualEnergyFormalism             = 0;    // Method: PPM parameter
  MultiSpecies                    = 0;    // 0:0 1:6 2:9 3:12
  GravityOn                       = 0;    // Whether gravity is included
  PressureFree                    = 0;    // Physics: passed into ppm_de
  ProblemType                     = 0;    //  7 (sedov) or 60 61 (turbulence)
  UseMinimumPressureSupport       = 0;    // call UseMinimumPressureSupport() ?
  ComovingBoxSize                 = 64;   // Physics cosmology: Mpc/h at z=0
  DualEnergyFormalismEta1         = 0.001;// Method: PPM parameter
  DualEnergyFormalismEta2         = 0.1;  // Method: PPM parameter
  Gamma                           = 5.0/3.0; // Physics: ideal gas law constant
  GravitationalConstant           = 1.0;  // used only in SetMinimumSupport()
  HubbleConstantNow               = 0.5;  // Physics: cosmology parameter
  MinimumPressureSupportParameter = 100; // SetMinimumSupport() Enzo parameter
  OmegaLambdaNow                  = 0.0;  // Physics: cosmology parameter
  OmegaMatterNow                  = 1.0;  // Physics: cosmology parameter
  pressure_floor                  = 1e-20; // Was "tiny_number"
  number_density_floor            = 1e-20; // Was "tiny_number"
  density_floor                   = 1e-20; // Was "tiny_number"
  temperature_floor               = 1e-20; // Was "tiny_number"
  DomainLeftEdge[0]               = 0; // for computing index: SolveHydroEquations()
  DomainLeftEdge[1]               = 0; // for computing index: SolveHydroEquations()
  DomainLeftEdge[2]               = 0; // for computing index: SolveHydroEquations()
  InitialRedshift                 = 20;  // Physics: Cosmology parameter
  WARNING("InitialTimeInCodeUnits not initialized properly");
  InitialTimeInCodeUnits          = -1;        // CosmologyComputeExpansionFactor()

  // Grid parameters
  
  GridRank           =  2;  // number of dimensions
  GridDimension[0]   = 38;  // total dimensions
  GridDimension[1]   = 38;  // total dimensions
  GridDimension[2]   = 38;  // total dimensions
  GridStartIndex[0]  =  3;  // starting index of the active region
  GridStartIndex[1]  =  3;  // starting index of the active region
  GridStartIndex[2]  =  3;  // starting index of the active region
  GridEndIndex[0]    = 35;  // stoping index of the active region
  GridEndIndex[1]    = 35;  // stoping index of the active region
  GridEndIndex[2]    = 35;  // stoping index of the active region
  GridLeftEdge[0]    = 0.0; // starting pos (active problem space)
  GridLeftEdge[1]    = 0.0; // starting pos (active problem space)
  GridLeftEdge[2]    = 0.0; // starting pos (active problem space)
  CourantSafetyNumber = 0.8; // Courant safety factor
  MaxExpansionRate   = 0.01; // Cosmology timestep constraint
  // Grid variables

//     float  dtFixed;                        // current (fixed) timestep
//     FLOAT  Time;                           // current problem time
//     FLOAT  OldTime;                        // time corresponding to OldBaryonField

//     int    NumberOfBaryonFields;                        // active baryon fields
//     float *BaryonField[MAX_NUMBER_OF_BARYON_FIELDS];    // pointers to arrays
//     float *OldBaryonField[MAX_NUMBER_OF_BARYON_FIELDS]; // pointers to old arrays
//     int    FieldType[MAX_NUMBER_OF_BARYON_FIELDS];
//     FLOAT *CellWidth[MAX_DIMENSION];

//     int    PPMFlatteningParameter;
//     int    PPMDiffusionParameter;
//     int    PPMSteepeningParameter;
//     float *AccelerationField[MAX_DIMENSION]; // cell cntr acceleration at n+1/2
//     int    ProcessorNumber;

//   //----------------------------------------------------------------------
//   // STRUCTS
//   //----------------------------------------------------------------------

//   struct fluxes
//   {
//     global_index LeftFluxStartGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
//     global_index LeftFluxEndGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
//     global_index RightFluxStartGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
//     global_index RightFluxEndGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
//     float       *LeftFluxes[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION];
//     float       *RightFluxes[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION];
//   };

//   //----------------------------------------------------------------------
//   // CLASSES
//   //----------------------------------------------------------------------

//   class grid
//   {


//   public:

//     int ComputeGammaField(float *GammaField);
//     int ComputePressure(FLOAT time, float *pressure);
//     int ComputePressureDualEnergyFormalism(FLOAT time, float *pressure);
//     int ComputeTemperatureField(float *temperature);
//     int IdentifyPhysicalQuantities(int &DensNum, int &GENum,   int &Vel1Num, 
// 				   int &Vel2Num, int &Vel3Num, int &TENum);
//     int IdentifySpeciesFields(int &DeNum, int &HINum, int &HIINum, 
// 			      int &HeINum, int &HeIINum, int &HeIIINum,
// 			      int &HMNum, int &H2INum, int &H2IINum,
// 			      int &DINum, int &DIINum, int &HDINum);
//     int SetMinimumSupport(float &MinimumSupportEnergyCoefficient);
//     int SolveHydroEquations(int CycleNumber, int NumberOfSubgrids, 
// 			    fluxes *SubgridFluxes[], int level);
//   };


}
    
