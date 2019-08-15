// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodHydro.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Oct 18 12:35:40 PDT 2017
/// @brief    Implements the EnzoMethodHydro class

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodHydro::EnzoMethodHydro 
( std::string method,
  enzo_float gamma,
  bool gravity,
  bool comoving_coordinates,
  bool dual_energy,
  enzo_float dual_energy_eta1,
  enzo_float dual_energy_eta2,
  std::string reconstruct_method,
  bool reconstruct_conservative,
  bool reconstruct_positive,
  enzo_float ppm_density_floor,
  enzo_float ppm_pressure_floor,
  int ppm_pressure_free,
  int ppm_diffusion,
  int ppm_flattening,
  int ppm_steepening,
  std::string riemann_solver
  )
  : Method(),
    method_(method),
    gamma_(gamma),
    gravity_(gravity?1:0),
    comoving_coordinates_(comoving_coordinates),
    dual_energy_(dual_energy?1:0),
    dual_energy_eta1_(dual_energy_eta1),
    dual_energy_eta2_(dual_energy_eta2),
    reconstruct_method_(reconstruct_method),
    reconstruct_conservative_(reconstruct_conservative ? 1:0),
    reconstruct_positive_(reconstruct_positive ? 1:0),
    ppm_density_floor_(ppm_density_floor),
    ppm_pressure_floor_(ppm_pressure_floor),
    ppm_pressure_free_(ppm_pressure_free),
    ppm_diffusion_(ppm_diffusion),
    ppm_flattening_(ppm_flattening),
    ppm_steepening_(ppm_steepening),
    riemann_solver_(riemann_solver)
    
{
  // Initialize default Refresh object

#ifdef NEW_REFRESH
#else  
  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier,
			     enzo_sync_id_method_ppm);

  refresh(ir)->add_field("density");
  refresh(ir)->add_field("velocity_x");
  refresh(ir)->add_field("velocity_y");
  refresh(ir)->add_field("velocity_z");
  refresh(ir)->add_field("acceleration_x");
  refresh(ir)->add_field("acceleration_y");
  refresh(ir)->add_field("acceleration_z");
  refresh(ir)->add_field("internal_energy");
  refresh(ir)->add_field("total_energy");
  refresh(ir)->add_field("pressure");
#endif
}

//----------------------------------------------------------------------

void EnzoMethodHydro::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | method_;
  p | gamma_;
  p | gravity_;
  p | comoving_coordinates_;
  p | dual_energy_;
  p | dual_energy_eta1_;
  p | dual_energy_eta2_;
  p | reconstruct_method_;
  p | reconstruct_conservative_;
  p | reconstruct_positive_;
  p | ppm_density_floor_;
  p | ppm_pressure_floor_;
  p | ppm_pressure_free_;
  p | ppm_diffusion_;
  p | ppm_flattening_;
  p | ppm_steepening_;
  p | riemann_solver_;
}

//----------------------------------------------------------------------

void EnzoMethodHydro::compute ( Block * block) throw()
{

  Field field = block->data()->field();

  if ( ! block->is_leaf() || field.field_count() == 0 ) {
    block->compute_done();
    return;
  }

  /* initialize */

  //     // MAX_COLOR is defined in fortran.def
  //     int dim, i, j, field, size, subgrid, n, colnum[MAX_COLOR];
  //     Elong_int GridGlobalStart[MAX_DIMENSION];
  //     FLOAT a = 1, dadt;

  //     /* Compute size (in floats) of the current grid. */

  //     size = 1;
  //     for (dim = 0; dim < GridRank; dim++)
  //       size *= GridDimension[dim];

  //     /* If multi-species being used, then treat them as colour variables
  //        (note: the solver has been modified to treat these as density vars). */

  //     int NumberOfColours = 0, ColourNum;

  //     // use different color fields for RadiativeTransferFLD problems
  //     //   first, the standard Enzo color field advection
  //     if (MultiSpecies > 0 && RadiativeTransferFLD != 2) {
  //       NumberOfColours = 6 + 3*(MultiSpecies-1);

  //       if ((ColourNum =
  //            FindField(ElectronDensity, FieldType, NumberOfBaryonFields)) < 0) {
  //         ENZO_FAIL("Could not find ElectronDensity.");
  //       }

  //       /* Generate an array of field numbers corresponding to the colour fields
  // 	 (here assumed to start with ElectronDensity and continue in order). */

  //       for (i = 0; i < NumberOfColours; i++)
  //         colnum[i] = ColourNum+i;

  //     }
  //     // second, the color field advection if using RadiativeTransferFLD for 
  //     // a radiation propagation problem (i.e. not using ray-tracing)
  //     if (RadiativeTransferFLD == 2) {
  //       if (ImplicitProblem < 4)  {  // grey radiation problem
	
  // 	// set the grey radiation field (required)
  // 	if ((ColourNum =
  // 	     FindField(RadiationFreq0, FieldType, NumberOfBaryonFields)) < 0) 
  // 	  ENZO_FAIL("Could not find RadiationFreq0.");
  // 	colnum[0] = ColourNum;

  // 	// check for other chemistry fields; add if they're present
  // 	//   ElectronDensity
  // 	if ((ColourNum =
  // 	     FindField(ElectronDensity, FieldType, NumberOfBaryonFields)) >= 0) 
  // 	  colnum[++NumberOfColours] = ColourNum;
  // 	//   HIDensity
  // 	if ((ColourNum =
  // 	     FindField(HIDensity, FieldType, NumberOfBaryonFields)) >= 0) 
  // 	  colnum[++NumberOfColours] = ColourNum;
  // 	//   HIIDensity
  // 	if ((ColourNum =
  // 	     FindField(HIIDensity, FieldType, NumberOfBaryonFields)) >= 0) 
  // 	  colnum[++NumberOfColours] = ColourNum;
  // 	//   HeIDensity
  // 	if ((ColourNum =
  // 	     FindField(HeIDensity, FieldType, NumberOfBaryonFields)) >= 0) 
  // 	  colnum[++NumberOfColours] = ColourNum;
  // 	//   HeIIDensity
  // 	if ((ColourNum =
  // 	     FindField(HeIIDensity, FieldType, NumberOfBaryonFields)) >= 0) 
  // 	  colnum[++NumberOfColours] = ColourNum;
  // 	//   HeIIIDensity
  // 	if ((ColourNum =
  // 	     FindField(HeIIIDensity, FieldType, NumberOfBaryonFields)) >= 0) 
  // 	  colnum[++NumberOfColours] = ColourNum;
  // 	//   HMDensity
  // 	if ((ColourNum =
  // 	     FindField(HMDensity, FieldType, NumberOfBaryonFields)) >= 0) 
  // 	  colnum[++NumberOfColours] = ColourNum;
  // 	//   H2IDensity
  // 	if ((ColourNum =
  // 	     FindField(H2IDensity, FieldType, NumberOfBaryonFields)) >= 0) 
  // 	  colnum[++NumberOfColours] = ColourNum;
  // 	//   H2IIDensity
  // 	if ((ColourNum =
  // 	     FindField(H2IIDensity, FieldType, NumberOfBaryonFields)) >= 0) 
  // 	  colnum[++NumberOfColours] = ColourNum;
  // 	//   DIDensity
  // 	if ((ColourNum =
  // 	     FindField(DIDensity, FieldType, NumberOfBaryonFields)) >= 0) 
  // 	  colnum[++NumberOfColours] = ColourNum;
  // 	//   DIIDensity
  // 	if ((ColourNum =
  // 	     FindField(DIIDensity, FieldType, NumberOfBaryonFields)) >= 0) 
  // 	  colnum[++NumberOfColours] = ColourNum;
  // 	//   HDIDensity
  // 	if ((ColourNum =
  // 	     FindField(HDIDensity, FieldType, NumberOfBaryonFields)) >= 0) 
  // 	  colnum[++NumberOfColours] = ColourNum;
  //       }
  //     }

  //     /* Add "real" colour fields (metallicity, etc.) as colour variables. */

  //     int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum,
  //       MetalIaNum, MetalIINum; 

  //     if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum,
  //                 MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum) == FAIL)
  //       ENZO_FAIL("Error in grid->IdentifyColourFields.\n");

  //     if (MetalNum != -1) {
  //       colnum[NumberOfColours++] = MetalNum;
  //       if (MultiMetals || TestProblemData.MultiMetals) {
  // 	colnum[NumberOfColours++] = MetalNum+1; //ExtraType0
  // 	colnum[NumberOfColours++] = MetalNum+2; //ExtraType1
  //       }
  //     }

  //     if (MetalIaNum       != -1) colnum[NumberOfColours++] = MetalIaNum;
  //     if (MetalIINum       != -1) colnum[NumberOfColours++] = MetalIINum;
  //     if (SNColourNum      != -1) colnum[NumberOfColours++] = SNColourNum;
  //     if (MBHColourNum     != -1) colnum[NumberOfColours++] = MBHColourNum;
  //     if (Galaxy1ColourNum != -1) colnum[NumberOfColours++] = Galaxy1ColourNum;
  //     if (Galaxy2ColourNum != -1) colnum[NumberOfColours++] = Galaxy2ColourNum;


  //     /* Add Simon Glover's chemistry species as color fields */

  //     if(TestProblemData.GloverChemistryModel){

  //       // Declarations for Simon Glover's cooling.
  //       int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
  // 	DINum, DIINum, HDINum;

  //       int CINum,CIINum,OINum,OIINum,SiINum,SiIINum,SiIIINum,CHINum,CH2INum,
  // 	CH3IINum,C2INum,COINum,HCOIINum,OHINum,H2OINum,O2INum;

  //       int GCM = TestProblemData.GloverChemistryModel;  // purely for convenience

  //       if (IdentifyGloverSpeciesFields(HIINum,HINum,H2INum,DINum,DIINum,HDINum,
  // 				      HeINum,HeIINum,HeIIINum,CINum,CIINum,OINum,
  // 				      OIINum,SiINum,SiIINum,SiIIINum,CHINum,CH2INum,
  // 				      CH3IINum,C2INum,COINum,HCOIINum,OHINum,H2OINum,
  // 				      O2INum) == FAIL) {
  // 	ENZO_FAIL("Error in IdentifyGloverSpeciesFields.");
  //       }

  //       colnum[NumberOfColours++] = HIINum;
  //       colnum[NumberOfColours++] = HINum;
  //       colnum[NumberOfColours++] = H2INum;

  //       if( (GCM==1) || (GCM==2) || (GCM==3) || (GCM==7) ){
  // 	colnum[NumberOfColours++] = DINum;
  // 	colnum[NumberOfColours++] = DIINum;
  // 	colnum[NumberOfColours++] = HDINum;
  // 	colnum[NumberOfColours++] = HeINum;
  // 	colnum[NumberOfColours++] = HeIINum;
  // 	colnum[NumberOfColours++] = HeIIINum;
  //       }

  //       if( (GCM==3) || (GCM==5) || (GCM==7) ){
  // 	colnum[NumberOfColours++] = COINum;
  //       }

  //       if( (GCM==2) || (GCM==3) || (GCM==7) ){
  // 	colnum[NumberOfColours++] = CINum;
  // 	colnum[NumberOfColours++] = CIINum;
  // 	colnum[NumberOfColours++] = OINum;
  // 	colnum[NumberOfColours++] = OIINum;
  //       }

  //       if( (GCM==2) || (GCM==3) ){
  // 	colnum[NumberOfColours++] = SiINum;
  // 	colnum[NumberOfColours++] = SiIINum;
  // 	colnum[NumberOfColours++] = SiIIINum;
  //       }

  //       if( (GCM==3) || (GCM==7) ){
  // 	colnum[NumberOfColours++] = CHINum;
  // 	colnum[NumberOfColours++] = CH2INum;
  // 	colnum[NumberOfColours++] = CH3IINum;
  // 	colnum[NumberOfColours++] = C2INum;
  // 	colnum[NumberOfColours++] = HCOIINum;
  // 	colnum[NumberOfColours++] = OHINum;
  // 	colnum[NumberOfColours++] = H2OINum;
  // 	colnum[NumberOfColours++] = O2INum;
  //       }
      
  //     } // if(TestProblemData.GloverChemistryModel)


  //     /* Add Cosmic Ray Energy Density as a colour variable. */
  //     if(CRModel){
  //       int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, CRNum;
  //       if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
  //                        Vel3Num, TENum, CRNum ) == FAIL )
  //         ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  //       colnum[NumberOfColours++] = CRNum;
  //     } // end CR if

  //     /* Add shock variables as a colour variable. */

  //     if(ShockMethod){
  //       int MachNum, PSTempNum,PSDenNum;
      
  //       if (IdentifyShockSpeciesFields(MachNum,PSTempNum,PSDenNum) == FAIL) {
  // 	ENZO_FAIL("Error in IdentifyShockSpeciesFields.")
  //       }
      
  //       colnum[NumberOfColours++] = MachNum;
  //       if(StorePreShockFields){
  // 	colnum[NumberOfColours++] = PSTempNum;
  // 	colnum[NumberOfColours++] = PSDenNum;
  //       }
  //     }
  //     /* Determine if Gamma should be a scalar or a field. */
    
  //     int UseGammaField = FALSE;
  //     float *GammaField = NULL;
  //     if (HydroMethod == Zeus_Hydro && MultiSpecies > 1) {
  //       UseGammaField = TRUE;
  //       GammaField = new float[size];
  //       if (this->ComputeGammaField(GammaField) == FAIL) {
  // 	ENZO_FAIL("Error in grid->ComputeGammaField.");
  //       }
  //     } else {
  //       GammaField = new float[1];
  //       GammaField[0] = gamma_;

  //     }
    
  //     /* Set lowest level flag (used on Zeus hydro). */

  //     int LowestLevel = (level > MaximumRefinementLevel-1) ? TRUE : FALSE;

  //     /* Set minimum support (used natively in zeus hydro). */

  //     float MinimumSupportEnergyCoefficient = 0;
  //     if (UseMinimumPressureSupport == TRUE && level > MaximumRefinementLevel-1)
  //       if (this->SetMinimumSupport(MinimumSupportEnergyCoefficient) == FAIL) {
  // 	ENZO_FAIL("Error in grid->SetMinimumSupport,");
  //       }

  //     /* allocate space for fluxes */

  //     /* Set up our restart dump fluxes container */
  //     this->SubgridFluxStorage = SubgridFluxes;
  //     this->NumberOfSubgrids = NumberOfSubgrids;

  //     for (i = 0; i < NumberOfSubgrids; i++) {
  //       for (dim = 0; dim < GridRank; dim++)  {

  // 	/* compute size (in floats) of flux storage */

  //         size = 1;
  //         for (j = 0; j < GridRank; j++)
  //           size *= SubgridFluxes[i]->LeftFluxEndGlobalIndex[dim][j] -
  //                   SubgridFluxes[i]->LeftFluxStartGlobalIndex[dim][j] + 1;

  // 	/* set unused dims (for the solver, which is hardwired for 3d). */

  //         for (j = GridRank; j < 3; j++) {
  //           SubgridFluxes[i]->LeftFluxStartGlobalIndex[dim][j] = 0;
  //           SubgridFluxes[i]->LeftFluxEndGlobalIndex[dim][j] = 0;
  //           SubgridFluxes[i]->RightFluxStartGlobalIndex[dim][j] = 0;
  //           SubgridFluxes[i]->RightFluxEndGlobalIndex[dim][j] = 0;
  //         }

  // 	/* Allocate space (if necessary). */

  //         for (field = 0; field < NumberOfBaryonFields; field++) {
  // 	  //
  // 	  if (SubgridFluxes[i]->LeftFluxes[field][dim] == NULL)
  // 	    SubgridFluxes[i]->LeftFluxes[field][dim]  = new float[size];
  // 	  if (SubgridFluxes[i]->RightFluxes[field][dim] == NULL)
  // 	    SubgridFluxes[i]->RightFluxes[field][dim] = new float[size];
  // 	  for (n = 0; n < size; n++) {
  // 	    SubgridFluxes[i]->LeftFluxes[field][dim][n] = 0;
  // 	    SubgridFluxes[i]->RightFluxes[field][dim][n] = 0;
  // 	  }
  //         }

  // 	for (field = NumberOfBaryonFields; field < MAX_NUMBER_OF_BARYON_FIELDS;
  // 	     field++) {
  //           SubgridFluxes[i]->LeftFluxes[field][dim] = NULL;
  //           SubgridFluxes[i]->RightFluxes[field][dim] = NULL;
  // 	}

  //       }  // next dimension

  //       /* make things pretty */

  //       for (dim = GridRank; dim < 3; dim++)
  //         for (field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++) {
  //           SubgridFluxes[i]->LeftFluxes[field][dim] = NULL;
  //           SubgridFluxes[i]->RightFluxes[field][dim] = NULL;
  // 	}

  //     } // end of loop over subgrids

  //     /* compute global start index for left edge of entire grid 
  //        (including boundary zones) */

  //     for (dim = 0; dim < GridRank; dim++)
  //       GridGlobalStart[dim] = nlongint((GridLeftEdge[dim]-DomainLeftEdge[dim])/(*(CellWidth[dim]))) -
  // 	GridStartIndex[dim];

  //     /* fix grid quantities so they are defined to at least 3 dims */

  //     for (i = GridRank; i < 3; i++) {
  //       GridDimension[i]   = 1;
  //       GridStartIndex[i]  = 0;
  //       GridEndIndex[i]    = 0;
  //       GridVelocity[i]    = 0.0;
  //       GridGlobalStart[i] = 0;
  //     }

  //     /* If using comoving coordinates, multiply dx by a(n+1/2).
  //        In one fell swoop, this recasts the equations solved by solver
  //        in comoving form (except for the expansion terms which are taken
  //        care of elsewhere). */


  //     /* Prepare Gravity. */

  //     int GravityOn = 0, FloatSize = sizeof(float);
  //     if (SelfGravity || UniformGravity || PointSourceGravity || DiskGravity || ExternalGravity )
  //       GravityOn = 1;
  // #ifdef TRANSFER
  //     if (RadiationPressure)
  //       GravityOn = 1;
  // #endif    

  //     //Some setup for MHDCT

  //     float *MagneticFlux[3][2];
  //     if ( UseMHDCT ) {
  //         MHDCT_ConvertEnergyToConservedS();  //Energy toggle.  Probably will be removed soon.
  //         for(field=0;field<3;field++){
  //             if(ElectricField[field] == NULL ) 
  //                 ElectricField[field] = new float[ElectricSize[field]];
  //             for(i=0;i<ElectricSize[field]; i++) ElectricField[field][i] = 0.0;
  //         }
  //         for(field=0;field<3;field++){
  //           MagneticFlux[field][0] = new float[2*MagneticSize[field]];
  //           MagneticFlux[field][1] =  MagneticFlux[field][0] +MagneticSize[field];
  //           for (i=0; i< 2*MagneticSize[field]; i++) MagneticFlux[field][0][i] = 0.0;
  //         }
  //         CenterMagneticField();
  // #ifdef BIERMANN
        
  //         /* Compute Units. */
        
  //         float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1, TimeUnits = 1,
  //           VelocityUnits = 1, BFieldUnits = 1;
        
  //         if(ComovingCoordinates){
  //           if (MHDCosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
  //           		     &TimeUnits, &VelocityUnits, Time,&BFieldUnits) == FAIL) {
  //             fprintf(stderr, "Error in MHD CosmologyGetUnits.\n");
  //             return FAIL;
  //           }
  //           /* Transform speed of light, hydrogen mass and electron change in units of ENZO */
  //           // Biermann battery constants in cgs units, convert to enzo units later
  //           float speedoflight = clight/VelocityUnits;
  //           float hydrogenmass = mh/DensityUnits*POW(LengthUnits,3);
  //           float electroncharge = 4.803e-10*
  //                 TimeUnits*BFieldUnits/(speedoflight*DensityUnits*POW(LengthUnits,3));
  //           float chi=1.0;
  //         }
  // #endif //BIERMANN
  //     }//UseMHDCT



  //     float* Fluxes[3] = {MagneticFlux[0][0],MagneticFlux[1][0],MagneticFlux[2][0]};
  //     int CurlStart[3] = {0,0,0}, 
  //     CurlEnd[3] = {mx-1,my-1,mz-1};
  //     if ( UseMHDCT ){
  //         if (HydroMethod == MHD_Li){
  //           this->SolveMHD_Li(CycleNumber, NumberOfSubgrids, SubgridFluxes, 
  //                 CellWidthTemp, GridGlobalStart, gravity_, 
  //                 NumberOfColours, colnum, Fluxes);
  //         }

  //         if( HydroMethod != NoHydro )
  //             switch( MHD_CT_Method ){
  //                 case CT_BalsaraSpicer: //1
  //                 case CT_Athena_LF:     //2
  //                 case CT_Athena_Switch: //3
  //                     ComputeElectricField(dtFixed, Fluxes);
  //                     break;
  //                 case CT_Biermann:      //4
  // #ifdef BIERMANN
  //                     FORTRAN_NAME(woc_create_e_biermann)(MagneticFlux[0][0], MagneticFlux[1][0], MagneticFlux[2][0],
  //                                  MagneticFlux[0][1], MagneticFlux[1][1], MagneticFlux[2][1],
  //                                  ElectricField[0], ElectricField[1], ElectricField[2],
  //                                  hxa, CellWidthTemp[1], CellWidthTemp[2],
  //                                  de,ei,
  //                                  GridDimension, GridDimension +1, GridDimension +2,
  //                                  GridStartIndex, GridEndIndex,
  //                                  GridStartIndex+1, GridEndIndex+1,
  //                                  GridStartIndex+2, GridEndIndex+2, &dtFixed, &UseDT,
  //                                  &gamma_, &speedoflight,&hydrogenmass, &electroncharge, &chi, a);
  //                     break;
  // #endif //BIERMANN
  //                 case CT_None:
  //                 default:
  //                     if(MyProcessorNumber == ROOT_PROCESSOR )
  //                         fprintf(stderr, "Warning: No CT method used with MHD_Li.\n");
  //                 break;
  //             }
  //         MHD_UpdateMagneticField(level, NULL, TRUE);
  //         CenterMagneticField();

  //         MHDCT_ConvertEnergyToSpecificS();
  //         for(field=0;field<3;field++){
  //           delete [] MagneticFlux[field][0];
  //         }
  //     }
  //     /* Call Solver on this grid.
  //        Notice that it is hard-wired for three dimensions, but it does
  //        the right thing for < 3 dimensions. */
  //     /* note: Start/EndIndex are zero based */
        
  if (method_ == "ppm") {

    ppm_method_ (block);
    
  }

  //     /* PPM LR has been withdrawn. */

  //     if (HydroMethod == PPM_LagrangeRemap) {
  // #ifdef PPM_LR
  //       FORTRAN_NAME(woc_ppm_lr)(
  // 			density, totalenergy, velocity1, velocity2, velocity3,
  //                           gasenergy,
  // 			&gravity_, AccelerationField[0],
  //                            AccelerationField[1],
  //                            AccelerationField[2],
  // 			&gamma_, &dtFixed, &CycleNumber,
  //                           hxa, CellWidthTemp[1], CellWidthTemp[2],
  // 			&GridRank, &mx, &my,
  //                            &mz, GridStartIndex, GridEndIndex,
  // 			GridVelocity, &PPMFlatteningParameter,
  //                            &ppm_pressure_free_,
  // 			&PPMDiffusionParameter, &PPMSteepeningParameter,
  //                            &dual_energy_, &dual_energy_eta1_,
  //                            &dual_energy_eta2_,
  // 			&NumberOfSubgrids, leftface, rightface,
  // 			istart, iend, jstart, jend,
  // 			standard, dindex, Eindex, uindex, vindex, windex,
  // 			  geindex, temp,
  //                         &NumberOfColours, colourpt, coloff, colindex);
  // #else /* PPM_LR */
  //       ENZO_FAIL("PPM LR is not supported.");
  // #endif /* PPM_LR */
  //     }

  //     if (HydroMethod == Zeus_Hydro)
  //       if (this->ZeusSolver(gamma_Field, UseGammaField, CycleNumber, 
  //                hxa, CellWidthTemp[1], CellWidthTemp[2],
  //                gravity_, NumberOfSubgrids, GridGlobalStart,
  //                SubgridFluxes,
  //                NumberOfColours, colnum, LowestLevel,
  //                MinimumSupportEnergyCoefficient) == FAIL)
  // 	ENZO_FAIL("ZeusSolver() failed!\n");
	


  //     /* Clean up allocated fields. */

  //     delete [] GammaField;   

  //     for (dim = 0; dim < MAX_DIMENSION; dim++)
  //       delete [] CellWidthTemp[dim];

  //   /* If we're supposed to be outputting on Density, we need to update
  //      the current maximum value of that Density. */

  //     if (OutputOnDensity == 1 || 
  //         StopFirstTimeAtDensity > 0. || 
  //         StopFirstTimeAtMetalEnrichedDensity > 0.) {
  //       int DensNum = FindField(Density, FieldType, NumberOfBaryonFields);

  //       int MetalNum = 0, SNColourNum = 0;
  //       int MetalFieldPresent = FALSE;
  //       float *MetalPointer;
  //       float *TotalMetals = NULL;

  //       if (StopFirstTimeAtMetalEnrichedDensity > 0.) {

  //         // First see if there's a metal field (so we can conserve species in
  //         // the solver)
  //         MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);
  //         SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields);
  //         MetalFieldPresent = (MetalNum != -1 || SNColourNum != -1);

  //         // Double check if there's a metal field when we have metal cooling
  //         if (MetalFieldPresent == FALSE) {
  //           ENZO_FAIL("StopFirstTimeAtMetalEnrichedDensity is set, but no metal field is present.\n");
  //         }

  //         /* If both metal fields (Pop I/II and III) exist, create a field
  //            that contains their sum */

  //         if (MetalNum != -1 && SNColourNum != -1) {
  //           TotalMetals = new float[size];
  //           for (int i = 0; i < size; i++)
  //             TotalMetals[i] = BaryonField[MetalNum][i] + BaryonField[SNColourNum][i];
  //           MetalPointer = TotalMetals;
  //         } // ENDIF both metal types
  //         else {
  //           if (MetalNum != -1)
  //             MetalPointer = BaryonField[MetalNum];
  //           else if (SNColourNum != -1)
  //             MetalPointer = BaryonField[SNColourNum];
  //         } // ENDELSE both metal types
 
  //       }

  //       for (i = 0; i < size; i++) {

  //         CurrentMaximumDensity =
  //             max(de[i], CurrentMaximumDensity);

  //         if (StopFirstTimeAtMetalEnrichedDensity > 0. &&
  //             (MetalPointer[i] / de[i]) > EnrichedMetalFraction) {
  //           CurrentMaximumMetalEnrichedDensity = 
  //             max(de[i], CurrentMaximumMetalEnrichedDensity);
  //         }

  //       }

  //       delete [] TotalMetals;

  //     } // end: if (OutputOnDensity == 1 || ...

  block->compute_done();


 
}

//----------------------------------------------------------------------

void EnzoMethodHydro::ppm_method_ ( Block * block )
{

  Field field = block->data()->field();
  
  // Assume gravity if acceleration fields exist
  
  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);

  // compute pressure

  EnzoComputePressure compute_pressure (gamma_,	comoving_coordinates_);
  compute_pressure.compute(block);

  const int cycle = block->cycle();
  const int rank  = cello::rank();

  for (int i0=0; i0<3; i0++) {
    int i = (i0 + cycle) % rank;

    // update in x-direction
    if ((mx > 1) && (i % rank == 0)) {
      for (int iz=0; iz<mz; iz++) {
	ppm_euler_x_(block,iz);    }
    }
    // update in y-direction
    if ((my > 1) && (i % rank == 1)) {
      for (int ix=0; ix<mx; ix++) {
	ppm_euler_y_(block,ix);
      }
    }
    // update in z-direction
    if ((mz > 1) && (i % rank == 2 )) {
      for (int iy=0; iy<my; iy++) {
	ppm_euler_z_(block,iy);
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoMethodHydro::ppm_euler_x_(Block * block, int iz)
{
  // int dim = 0, idim = 1, jdim = 2;
  // int dim_p1 = dim+1;   // To match definition in calcdiss
  // int ierr = 0;

  // /* Find fields: density, total energy, velocity1-3. */

  // int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  
  // this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
  // 				   Vel3Num, TENum);

  // int nxz, nyz, nzz, ixyz;

  // nxz = GridEndIndex[0] - GridStartIndex[0] + 1;
  // nyz = GridEndIndex[1] - GridStartIndex[1] + 1;
  // nzz = GridEndIndex[2] - GridStartIndex[2] + 1;

  // float MinimumPressure = tiny_number;
  
  // // Copy from field to slice

  // float *dslice, *eslice, *uslice, *vslice, *wslice, *grslice, *geslice, 
  //   *colslice, *pslice;


  Field field = block->data()->field();
  
  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);

  int nx,ny,nz;
  field.size (&nx,&ny,&nz);

  int gx,gy,gz;
  field.ghost_depth (0,&gx,&gy,&gz);

  // determine temporary array size

  //   ... compute slice size
  const int ns = mx*my;

  int na = 0;

  // ...density, total energy, velocities, pressure
  na += 6*ns;

  // ...add slice for gravity if needed
  if (gravity_) na += 1*ns;

  // ... add slice for gas energy if needed
  if (dual_energy_) na += 1*ns;

  // ... add slices for colour fields
  Grouping * field_groups = block->data()->field().groups();
  int nc = field_groups->size("colour");
  na += nc*ns;

  // allocate array
  enzo_float * slice_array = new enzo_float [na];

  // initialize array of slices
  
  enzo_float * pa = slice_array;

  // ...density, pressure, and velocities
  enzo_float * dslice = pa; pa += ns;
  enzo_float * eslice = pa; pa += ns;
  enzo_float * pslice = pa; pa += ns;
  enzo_float * uslice = pa; pa += ns;
  enzo_float * vslice = pa; pa += ns;
  enzo_float * wslice = pa; pa += ns;

  // ...add slice for gravity if needed
  enzo_float * grslice = NULL;
  if (gravity_) {
    grslice = pa; pa +=ns;
  }
  
  // ... add slice for gas energy if needed
  enzo_float * geslice = NULL;
  if (dual_energy_) {
    geslice = pa; pa +=ns;
  }

  // add slices for colour fields
  enzo_float * colslice = NULL;
  if (nc > 0) {
    colslice = pa; pa += nc*ns;
  }

  enzo_float * de = (enzo_float *) field.values("density");
  enzo_float * te = (enzo_float *) field.values("total_energy");
  enzo_float * vx = (enzo_float *) field.values("velocity_x");
  enzo_float * vy = (enzo_float *) field.values("velocity_y");
  enzo_float * vz = (enzo_float *) field.values("velocity_z");
  enzo_float * pr = (enzo_float *) field.values("pressure");

  const int rank = cello::rank();
  
  for (int iy=0; iy<my; iy++) {
    for (int ix=0; ix<mx; ix++) {
      int i  = ix + mx*(iy + my*iz);
      int is = ix + mx*iy;
      dslice[is] = de[i];
      eslice[is] = te[i];
      pslice[is] = pr[i];
      uslice[is] = vx[i];
      vslice[is] = (rank >= 2) ? vy[i] : 0.0;
      wslice[is] = (rank >= 3) ? vz[i] : 0.0;
    } // ENDFOR i
  }

  if (gravity_) {
    enzo_float * ax = (enzo_float *) field.values("acceleration_x");
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {
	int i  = ix + mx*(iy + my*iz);
	int is = ix + mx*iy;
	grslice[is] = ax[i];
      }
    }
  }

  if (dual_energy_) {
    enzo_float * ei = (enzo_float *) field.values("internal_energy");
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {
	int i  = ix + mx*(iy + my*iz);
	int is = ix + mx*iy;
	geslice[is] = ei[i];
      }
    }
  }

  for (int ic=0; ic<nc; ic++) {
    enzo_float * c = (enzo_float *)
      field.values(field_groups->item("colour",ic));
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {
	int i = ix + mx*(iy + my*iz);
	int k = ix + mx*(iy + my*ic);
	colslice[k] = c[i];
      }
    }
  }
    
  ASSERT2("EnzoMethodHydro::ppm_xeuler_x",
	  "temporary slice array actual size %d differs from expected size %d",
	  (pa-slice_array), na,
	  ((pa -slice_array) == na));

  // Allocate memory for fluxes

  enzo_float *dls, *drs, *flatten, *pbar, *pls, *prs, *ubar, *uls, *urs, *vls, 
    *vrs, *gels, *gers, *wls, *wrs, *diffcoef, *df, *ef, *uf, *vf, *wf, *gef,
    *ges, *colf, *colls, *colrs;

  int nf = (23 + 3*nc)*ns;
  
  enzo_float * fluxes_array = new enzo_float [nf];

  enzo_float * pf = fluxes_array;
  
  dls      = pf; pf += ns;
  drs      = pf; pf += ns;	
  flatten  = pf; pf += ns;	
  pbar     = pf; pf += ns;	
  pls      = pf; pf += ns;	
  prs      = pf; pf += ns;	
  ubar     = pf; pf += ns;	
  uls      = pf; pf += ns;	
  urs      = pf; pf += ns;	
  vls      = pf; pf += ns;	
  vrs      = pf; pf += ns;	
  gels     = pf; pf += ns;	
  gers     = pf; pf += ns;	
  wls      = pf; pf += ns;	
  wrs      = pf; pf += ns;	
  diffcoef = pf; pf += ns;	
  df       = pf; pf += ns;		
  ef       = pf; pf += ns;		
  uf       = pf; pf += ns;		
  vf       = pf; pf += ns;		
  wf       = pf; pf += ns;		
  gef      = pf; pf += ns;	
  ges      = pf; pf += ns;	
  colf     = pf; pf += nc*ns;
  colls    = pf; pf += nc*ns;
  colrs    = pf; pf += nc*ns;

  ASSERT2("EnzoMethodHydro::ppm_xeuler_x",
	  "temporary fluxes array actual size %d differs from expected size %d",
	  (pf-fluxes_array), nf,
	  ((pf -fluxes_array) == nf));

  // Convert start and end indexes into 1-based for FORTRAN

  // int is, ie, js, je, is_m3, ie_p3, ie_p1, k_p1;

  // is = GridStartIndex[0] + 1;
  // ie = GridEndIndex[0] + 1;
  // js = 1;
  // je = my;
  // is_m3 = is - 3;
  // ie_p1 = ie + 1;
  // ie_p3 = ie + 3;
  // k_p1 = k + 1;

  // Compute the pressure on a slice

  int is = 1 + gx;
  int ie = 1 + gx+nx;
  int js = 1;
  int je = my;
  int ie_p1 = ie + 1;
  int k_p1  = iz + 1;
  

  if (dual_energy_) {

    FORTRAN_NAME(woc_pgas2d_dual)
      (dslice, eslice, geslice, pslice, uslice, vslice, 
       wslice, &dual_energy_eta1_,
       &dual_energy_eta2_,
       &mx,&my, &is, &ie, &js, &je, 
       &gamma_, &ppm_pressure_floor_);

  } else {

    FORTRAN_NAME(woc_pgas2d)
      (dslice, eslice, pslice, uslice, vslice, 
       wslice, &mx, &my,
       &is, &ie, &js, &je, &gamma_, &ppm_pressure_floor_);
  }

  
  // If requested, compute diffusion and slope flattening coefficients

  // Adjust cell widths for cosmological expansion if needed
  
  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  enzo_float cosmo_a    = 1.0;
  enzo_float cosmo_dadt = 0.0;
  
  if (cosmology) {
    cosmology->compute_expansion_factor
      (&cosmo_a, &cosmo_dadt, (enzo_float)block->time());
  }

  double h = 0.0;
  block->cell_width(&h);

  enzo_float hxa = cosmo_a*h;
  enzo_float hya = cosmo_a*h;
  enzo_float hza = cosmo_a*h;
  
  enzo_float dt = block->dt();
  
  enzo_float * flatten_array = new enzo_float[ns];

  int riemann_solver_fallback = 1;
  
  if (ppm_diffusion_ || ppm_flattening_) {

    int dim_p1 = 1;

    FORTRAN_NAME(woc_calcdiss)(dslice, eslice, uslice,
			   vy, vz,
			   pslice,
			   &hxa, &hya, &hza,
			   &mx,  &my, &mz,
  			   &is, &ie, &js, &je, &k_p1,
  			   &nz, &dim_p1, &mx, &my, &mz,
  			   &dt, &gamma_, &ppm_diffusion_,
  			   &ppm_flattening_, diffcoef, flatten_array);

  }

  // Compute Eulerian left and right states at zone edges via interpolation

  if (reconstruct_method_ == "ppm") {

    FORTRAN_NAME(woc_inteuler)
      (dslice, pslice, &gravity_, grslice, geslice, uslice,
       vslice, wslice, &hxa, flatten,
       &mx, &my,
       &is, &ie, &js, &je, &dual_energy_,
       &dual_energy_eta1_, &dual_energy_eta2_,
       &ppm_steepening_, &ppm_flattening_,
       &reconstruct_conservative_, &reconstruct_positive_,
       &dt, &gamma_, &ppm_pressure_free_, 
       dls, drs, pls, prs, gels, gers, uls, urs, vls, vrs,
       wls, wrs, &nc, colslice, colls, colrs);
  }

  // Compute (Lagrangian part of the) Riemann problem at each zone boundary

  if (riemann_solver_ == "two_shock") {

    FORTRAN_NAME(woc_twoshock)
      (dls, drs, pls, prs, uls, urs,
       &mx, &my,
       &is, &ie_p1, &js, &je,
       &dt, &gamma_, &ppm_pressure_floor_, &ppm_pressure_free_,
       pbar, ubar, &gravity_, grslice,
       &dual_energy_, &dual_energy_eta1_);
    
    FORTRAN_NAME(woc_flux_twoshock)
      (dslice, eslice, geslice, uslice, vslice, wslice,
       &hxa, diffcoef, 
       &mx, &my,
       &is, &ie, &js, &je, &dt, &gamma_,
       &ppm_diffusion_, &dual_energy_,
       &dual_energy_eta1_,
       &riemann_solver_fallback,
       dls, drs, pls, prs, gels, gers, uls, urs,
       vls, vrs, wls, wrs, pbar, ubar,
       df, ef, uf, vf, wf, gef, ges,
       &nc, colslice, colls, colrs, colf);
    
  } else if (riemann_solver_ == "hll") {

    FORTRAN_NAME(woc_flux_hll)
      (dslice, eslice, geslice, uslice, vslice, wslice,
       &hxa, diffcoef, 
       &mx, &my,
       &is, &ie, &js, &je, &dt, &gamma_,
       &ppm_diffusion_, &dual_energy_,
       &dual_energy_eta1_,
       &riemann_solver_fallback,
       dls, drs, pls, prs, uls, urs,
       vls, vrs, wls, wrs, gels, gers,
       df, uf, vf, wf, ef, gef, ges,
       &nc, colslice, colls, colrs, colf);

  } else if (riemann_solver_ == "hllc") {

    FORTRAN_NAME(woc_flux_hllc)
      (dslice, eslice, geslice, uslice, vslice, wslice,
       &hxa, diffcoef, 
       &mx, &my,
       &is, &ie, &js, &je, &dt, &gamma_,
       &ppm_diffusion_, &dual_energy_,
       &dual_energy_eta1_,
       &riemann_solver_fallback,
       dls, drs, pls, prs, uls, urs,
       vls, vrs, wls, wrs, gels, gers,
       df, uf, vf, wf, ef, gef, ges,
       &nc, colslice, colls, colrs, colf);

  } else {

    for (int i = 0; i < ns; i++) {
      df[i] = 0;
      ef[i] = 0;
      uf[i] = 0;
      vf[i] = 0;
      wf[i] = 0;
      gef[i] = 0;
      ges[i] = 0;
    }

  }

  // Compute Eulerian fluxes and update zone-centered quantities

  FORTRAN_NAME(woc_euler)
    (dslice, eslice, grslice, geslice, uslice, vslice, wslice,
     &hxa, diffcoef, 
     &mx, &my, 
     &is, &ie, &js, &je, &dt, &gamma_, 
     &ppm_diffusion_, &gravity_, &dual_energy_, 
     &dual_energy_eta1_, &dual_energy_eta2_,
     df, ef, uf, vf, wf, gef, ges,
     &nc, colslice, colf, &ppm_density_floor_);

  // /* If necessary, recompute the pressure to correctly set ge and e */

  // if (dual_energy_)
  //   FORTRAN_NAME(woc_pgas2d_dual)(dslice, eslice, geslice, pslice, uslice, vslice, 
  // 			      wslice, &dual_energy_eta1_, 
  // 			      &dual_energy_eta2_, &mx, 
  // 			      &my, &is_m3, &ie_p3, &js, &je, 
  // 			      &gamma_, &ppm_pressure_floor_);

  // /* Check this slice against the list of subgrids (all subgrid
  //    quantities are zero-based) */

  // int jstart, jend, offset, nfi, lface, rface, lindex, rindex, 
  //   fistart, fiend, fjstart, fjend, clindex, crindex;
  
  // for (n = 0; n < NumberOfSubgrids; n++) {

  //   fistart = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][idim] - 
  //     GridGlobalStart[idim];
  //   fiend = SubgridFluxes[n]->RightFluxEndGlobalIndex[dim][idim] -
  //     GridGlobalStart[idim];
  //   fjstart = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][jdim] - 
  //     GridGlobalStart[jdim];
  //   fjend = SubgridFluxes[n]->RightFluxEndGlobalIndex[dim][jdim] -
  //     GridGlobalStart[jdim];

  //   if (k >= fjstart && k <= fjend) {

  //     nfi = fiend - fistart + 1;
  //     for (j = fistart; j <= fiend; j++) {

  // 	offset = (j-fistart) + (k-fjstart)*nfi;

  // 	lface = SubgridFluxes[n]->LeftFluxStartGlobalIndex[dim][dim] -
  // 	  GridGlobalStart[dim];
  // 	lindex = j * GridDimension[dim] + lface;

  // 	rface = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][dim] -
  // 	  GridGlobalStart[dim] + 1;
  // 	rindex = j * GridDimension[dim] + rface;	

  // 	SubgridFluxes[n]->LeftFluxes [DensNum][dim][offset] = df[lindex];
  // 	SubgridFluxes[n]->RightFluxes[DensNum][dim][offset] = df[rindex];
  // 	SubgridFluxes[n]->LeftFluxes [TENum][dim][offset]   = ef[lindex];
  // 	SubgridFluxes[n]->RightFluxes[TENum][dim][offset]   = ef[rindex];
  // 	SubgridFluxes[n]->LeftFluxes [Vel1Num][dim][offset] = uf[lindex];
  // 	SubgridFluxes[n]->RightFluxes[Vel1Num][dim][offset] = uf[rindex];

  // 	if (nyz > 1) {
  // 	  SubgridFluxes[n]->LeftFluxes [Vel2Num][dim][offset] = vf[lindex];
  // 	  SubgridFluxes[n]->RightFluxes[Vel2Num][dim][offset] = vf[rindex];
  // 	} // ENDIF y-data

  // 	if (nzz > 1) {
  // 	  SubgridFluxes[n]->LeftFluxes [Vel3Num][dim][offset] = wf[lindex];
  // 	  SubgridFluxes[n]->RightFluxes[Vel3Num][dim][offset] = wf[rindex];
  // 	} // ENDIF z-data

  // 	if (dual_energy_) {
  // 	  SubgridFluxes[n]->LeftFluxes [GENum][dim][offset] = gef[lindex];
  // 	  SubgridFluxes[n]->RightFluxes[GENum][dim][offset] = gef[rindex];
  // 	} // ENDIF dual_energy_

  // 	for (ncolour = 0; ncolour < nc; ncolour++) {
  // 	  clindex = (j + ncolour * my) * GridDimension[dim] +
  // 	    lface;
  // 	  crindex = (j + ncolour * my) * GridDimension[dim] +
  // 	    rface;

  // 	  SubgridFluxes[n]->LeftFluxes [colnum[ncolour]][dim][offset] = 
  // 	    colf[clindex];
  // 	  SubgridFluxes[n]->RightFluxes[colnum[ncolour]][dim][offset] = 
  // 	    colf[crindex];
  // 	} // ENDFOR ncolour

  //     } // ENDFOR J

  //   } // ENDIF k inside

  // } // ENDFOR n

  // /* Copy from slice to field */

  // for (j = 0; j < my; j++) {

  //   index2 = j * mx;

  //   for (i = 0; i < mx; i++) {
  //     index3 = (k*my + j)*mx + i;
  //     de[index3] = dslice[index2+i];
  //     et[index3] = eslice[index2+i];
  //     vx[index3] = uslice[index2+i];
  //   } // ENDFOR i

  //   if (GridRank > 1)
  //     for (i = 0; i < mx; i++) {
  // 	index3 = (k*my + j)*mx + i;
  // 	vy[index3] = vslice[index2+i];
  //     }

  //   if (GridRank > 2)
  //     for (i = 0; i < mx; i++) {
  // 	index3 = (k*my + j)*mx + i;
  // 	vz[index3] = wslice[index2+i];
  //     }

  //   if (dual_energy_)
  //     for (i = 0; i < mx; i++) {
  // 	index3 = (k*my + j)*mx + i;
  // 	ei[index3] = geslice[index2+i];
  //     }

  //   for (n = 0; n < nc; n++) {
  //     index2 = (n*my + j) * mx;
  //     for (i = 0; i < mx; i++) {
  // 	index3 = (k*my + j) * mx + i;
  // 	BaryonField[colnum[n]][index3] = colslice[index2+i];
  //     }
  //   } // ENDFOR colours
  // } // ENDFOR j

  // deallocate array
  delete [] flatten_array;
  delete [] fluxes_array;
  delete [] slice_array;
}

//----------------------------------------------------------------------

void EnzoMethodHydro::ppm_euler_y_ (Block * block, int ix)
{

  // int dim = 1, idim = 0, jdim = 2;
  // int dim_p1 = dim+1;   // To match definition in calcdiss

  // /* Find fields: density, total energy, velocity1-3. */
  
  // int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;

  // this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
  // 				   Vel3Num, TENum);

  // int nxz, nyz, nzz, ixyz;

  // nxz = GridEndIndex[0] - GridStartIndex[0] + 1;
  // nyz = GridEndIndex[1] - GridStartIndex[1] + 1;
  // nzz = GridEndIndex[2] - GridStartIndex[2] + 1;

  // float MinimumPressure = tiny_number;
  
  // // Copy from field to slice

  // float *dslice, *eslice, *uslice, *vslice, *wslice, *grslice, *geslice, 
  //   *colslice, *pslice;

  // int size = my * mz;
  // dslice = new float[size];
  // eslice = new float[size];
  // uslice = new float[size];
  // vslice = new float[size];
  // wslice = new float[size];
  // pslice = new float[size];
  // if (gravity_) {
  //   grslice = new float[size];  
  // }
  // if (dual_energy_) {
  //   geslice = new float[size];  
  // }
  // if (nc > 0) {
  //   colslice = new float[nc * size];  
  // }

  // int j, k, n, ncolour, index2, index3;
  // for (k = 0; k < mz; k++) {

  //   index2 = k * my;

  //   for (j = 0; j < my; j++) {
  //     index3 = (k*my + j) * mx + i;
  //     dslice[index2+j] = de[index3];
  //     eslice[index2+j] = et[index3];
  //     pslice[index2+j] = pressure[index3];
  //     wslice[index2+j] = vx[index3];
  //   } // ENDFOR i

  //   // Set velocities to zero if rank < 3 since hydro routines are
  //   // hard-coded for 3-d

  //   if (GridRank > 1) 
  //     for (j = 0; j < my; j++) {
  // 	index3 = (k*my + j) * mx + i;
  // 	uslice[index2+j] = vy[index3];
  //     }
  //   else
  //     for (j = 0; j < my; j++)
  // 	uslice[index2+j] = 0;
  
  //   if (GridRank > 2)
  //     for (j = 0; j < my; j++) {
  // 	index3 = (k*my + j) * mx + i;
  // 	vslice[index2+j] = vz[index3];
  //     }
  //   else
  //     for (j = 0; j < my; j++)
  // 	vslice[index2+j] = 0;
    
  //   if (gravity_)
  //     for (j = 0; j < my; j++) {
  // 	index3 = (k*my + j) * mx + i;
  // 	grslice[index2+j] = AccelerationField[dim][index3];
  //     }

  //   if (dual_energy_)
  //     for (j = 0; j < my; j++) {
  // 	index3 = (k*my + j) * mx + i;
  // 	geslice[index2+j] = ei[index3];
  //     }

  //   for (n = 0; n < nc; n++) {
  //     index2 = (n*mz + k) * my;
  //     for (j = 0; j < my; j++) {
  // 	index3 = (k*my + j) * mx + i;
  // 	colslice[index2+j] = BaryonField[colnum[n]][index3];
  //     }
  //   } // ENDFOR colours

  // } // ENDFOR j

  // /* Allocate memory for fluxes */

  // float *dls, *drs, *flatten, *pbar, *pls, *prs, *ubar, *uls, *urs, *vls, 
  //   *vrs, *gels, *gers, *wls, *wrs, *diffcoef, *df, *ef, *uf, *vf, *wf, *gef,
  //   *ges, *colf, *colls, *colrs;

  // dls = new float[size];	
  // drs = new float[size];	
  // flatten = new float[size];	
  // pbar = new float[size];	
  // pls = new float[size];	
  // prs = new float[size];	
  // ubar = new float[size];	
  // uls = new float[size];	
  // urs = new float[size];	
  // vls = new float[size];	
  // vrs = new float[size];	
  // gels = new float[size];	
  // gers = new float[size];	
  // wls = new float[size];	
  // wrs = new float[size];	
  // diffcoef = new float[size];	
  // df = new float[size];		
  // ef = new float[size];		
  // uf = new float[size];		
  // vf = new float[size];		
  // wf = new float[size];		
  // gef = new float[size];	
  // ges = new float[size];
  // colf = new float[nc*size];  
  // colls = new float[nc*size];  
  // colrs = new float[nc*size];  

  // /* Convert start and end indexes into 1-based for FORTRAN */

  // int is, ie, js, je, is_m3, ie_p3, ie_p1, k_p1;

  // is = GridStartIndex[1] + 1;
  // ie = GridEndIndex[1] + 1;
  // js = 1;
  // je = mz;
  // is_m3 = is - 3;
  // ie_p1 = ie + 1;
  // ie_p3 = ie + 3;
  // k_p1 = i + 1;

  // /* Compute the pressure on a slice */

  // /*
  // if (dual_energy_)
  //   FORTRAN_NAME(woc_pgas2d_dual)(dslice, eslice, geslice, pslice, uslice, vslice, 
  // 			      wslice, &dual_energy_eta1_, 
  // 			      &dual_energy_eta2_, &my, 
  // 			      &mz, &is_m3, &ie_p3, &js, &je, 
  // 			      &gamma_, &ppm_pressure_floor_);
  // else
  //   FORTRAN_NAME(woc_pgas2d)(dslice, eslice, pslice, uslice, vslice,
  // 			 wslice, &my, &mz, 
  // 			 &is_m3, &ie_p3, &js, &je, &gamma_, &ppm_pressure_floor_);
  // */
  // /* If requested, compute diffusion and slope flattening coefficients */

  // if (ppm_diffusion_ != 0 || PPMFlatteningParameter != 0)
  //   FORTRAN_NAME(woc_calcdiss)(dslice, eslice, uslice, vz,
  // 			   vx, pslice, CellWidthTemp[1],
  // 			   CellWidthTemp[2], hxa, 
  // 			   &my, &mz,
  // 			   &mx, &is, &ie, &js, &je, &k_p1,
  // 			   &nxz, &dim_p1, &mx,
  // 			   &my, &mz,
  // 			   &dt, &gamma_, &ppm_diffusion_,
  // 			   &PPMFlatteningParameter, diffcoef, flatten);

  // /* Compute Eulerian left and right states at zone edges via interpolation */

  // if (ReconstructionMethod == PPM)
  //   FORTRAN_NAME(woc_inteuler)(dslice, pslice, &gravity_, grslice, geslice, uslice,
  // 			   vslice, wslice, CellWidthTemp[1], flatten,
  // 			   &my, &mz,
  // 			   &is, &ie, &js, &je, &dual_energy_, 
  // 			   &dual_energy_eta1_, &dual_energy_eta2_,
  // 			   &PPMSteepeningParameter, &PPMFlatteningParameter,
  // 			   &ConservativeReconstruction, &PositiveReconstruction,
  // 			   &dt, &gamma_, &ppm_pressure_free_, 
  // 			   dls, drs, pls, prs, gels, gers, uls, urs, vls, vrs,
  // 			   wls, wrs, &nc, colslice, colls, colrs);

  // /* Compute (Lagrangian part of the) Riemann problem at each zone boundary */

  // switch (RiemannSolver) {
  // case TwoShock:
  //   FORTRAN_NAME(woc_twoshock)(dls, drs, pls, prs, uls, urs,
  // 			   &my, &mz,
  // 			   &is, &ie_p1, &js, &je,
  // 			   &dt, &gamma_, &ppm_pressure_floor_, &ppm_pressure_free_,
  // 			   pbar, ubar, &gravity_, grslice,
  // 			   &dual_energy_, &dual_energy_eta1_);
    
  //   FORTRAN_NAME(woc_flux_twoshock)(dslice, eslice, geslice, uslice, vslice, wslice,
  // 				CellWidthTemp[1], diffcoef, 
  // 				&my, &mz,
  // 				&is, &ie, &js, &je, &dt, &gamma_,
  // 				&ppm_diffusion_, &dual_energy_,
  // 				&dual_energy_eta1_,
  // 				&riemann_solver_fallback,
  // 				dls, drs, pls, prs, gels, gers, uls, urs,
  // 				vls, vrs, wls, wrs, pbar, ubar,
  // 				df, ef, uf, vf, wf, gef, ges,
  // 				&nc, colslice, colls, colrs, colf);
  //   break;

  // case HLL:
  //   FORTRAN_NAME(woc_flux_hll)(dslice, eslice, geslice, uslice, vslice, wslice,
  // 			   CellWidthTemp[1], diffcoef, 
  // 			   &my, &mz,
  // 			   &is, &ie, &js, &je, &dt, &gamma_,
  // 			   &ppm_diffusion_, &dual_energy_,
  // 			   &dual_energy_eta1_,
  // 			   &riemann_solver_fallback,
  // 			   dls, drs, pls, prs, uls, urs,
  // 			   vls, vrs, wls, wrs, gels, gers,
  // 			   df, uf, vf, wf, ef, gef, ges,
  // 			   &nc, colslice, colls, colrs, colf);
  //   break;

  // case HLLC:
  //   FORTRAN_NAME(woc_flux_hllc)(dslice, eslice, geslice, uslice, vslice, wslice,
  // 			    CellWidthTemp[1], diffcoef, 
  // 			    &my, &mz,
  // 			    &is, &ie, &js, &je, &dt, &gamma_,
  // 			    &ppm_diffusion_, &dual_energy_,
  // 			    &dual_energy_eta1_,
  // 			    &riemann_solver_fallback,
  // 			    dls, drs, pls, prs, uls, urs,
  // 			    vls, vrs, wls, wrs, gels, gers,
  // 			    df, uf, vf, wf, ef, gef, ges,
  // 			    &nc, colslice, colls, colrs, colf);
  //   break;

  // default:
  //   for (int index = 0; index < size; index++) {
  //     df[index] = 0;
  //     ef[index] = 0;
  //     uf[index] = 0;
  //     vf[index] = 0;
  //     wf[index] = 0;
  //     gef[index] = 0;
  //     ges[index] = 0;
  //   }
  //   break;

  // } // ENDCASE


  // /* Compute Eulerian fluxes and update zone-centered quantities */

  // FORTRAN_NAME(woc_euler)(dslice, eslice, grslice, geslice, uslice, vslice, wslice,
  // 		      CellWidthTemp[1], diffcoef, 
  // 		      &my, &mz, 
  // 		      &is, &ie, &js, &je, &dt, &gamma_, 
  // 		      &ppm_diffusion_, &gravity_, &dual_energy_, 
  // 		      &dual_energy_eta1_, &dual_energy_eta2_,
  // 		      df, ef, uf, vf, wf, gef, ges,
  // 		      &nc, colslice, colf, &density_floor_);

  // /* If necessary, recompute the pressure to correctly set ge and e */

  // if (dual_energy_)
  //   FORTRAN_NAME(woc_pgas2d_dual)(dslice, eslice, geslice, pslice, uslice, vslice, 
  // 			      wslice, &dual_energy_eta1_, 
  // 			      &dual_energy_eta2_, &my, 
  // 			      &mz, &is_m3, &ie_p3, &js, &je, 
  // 			      &gamma_, &ppm_pressure_floor_);

  // /* Check this slice against the list of subgrids (all subgrid
  //    quantities are zero-based) */

  // int jstart, jend, offset, nfi, lface, rface, lindex, rindex, 
  //   fistart, fiend, fjstart, fjend, clindex, crindex;
  
  // for (n = 0; n < NumberOfSubgrids; n++) {

  //   fistart = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][idim] - 
  //     GridGlobalStart[idim];
  //   fiend = SubgridFluxes[n]->RightFluxEndGlobalIndex[dim][idim] -
  //     GridGlobalStart[idim];
  //   fjstart = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][jdim] - 
  //     GridGlobalStart[jdim];
  //   fjend = SubgridFluxes[n]->RightFluxEndGlobalIndex[dim][jdim] -
  //     GridGlobalStart[jdim];

  //   if (i >= fistart && i <= fiend) {

  //     nfi = fiend - fistart + 1;
  //     for (k = fjstart; k <= fjend; k++) {

  // 	offset = (i-fistart) + (k-fjstart)*nfi;

  // 	lface = SubgridFluxes[n]->LeftFluxStartGlobalIndex[dim][dim] -
  // 	  GridGlobalStart[dim];
  // 	lindex = k * GridDimension[dim] + lface;

  // 	rface = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][dim] -
  // 	  GridGlobalStart[dim] + 1;
  // 	rindex = k * GridDimension[dim] + rface;	

  // 	SubgridFluxes[n]->LeftFluxes [DensNum][dim][offset] = df[lindex];
  // 	SubgridFluxes[n]->RightFluxes[DensNum][dim][offset] = df[rindex];
  // 	SubgridFluxes[n]->LeftFluxes [TENum][dim][offset]   = ef[lindex];
  // 	SubgridFluxes[n]->RightFluxes[TENum][dim][offset]   = ef[rindex];

  // 	if (nxz > 1) {
  // 	  SubgridFluxes[n]->LeftFluxes [Vel1Num][dim][offset] = wf[lindex];
  // 	  SubgridFluxes[n]->RightFluxes[Vel1Num][dim][offset] = wf[rindex];
  // 	} // ENDIF x-data

  // 	SubgridFluxes[n]->LeftFluxes [Vel2Num][dim][offset] = uf[lindex];
  // 	SubgridFluxes[n]->RightFluxes[Vel2Num][dim][offset] = uf[rindex];

  // 	if (nzz > 1) {
  // 	  SubgridFluxes[n]->LeftFluxes [Vel3Num][dim][offset] = vf[lindex];
  // 	  SubgridFluxes[n]->RightFluxes[Vel3Num][dim][offset] = vf[rindex];
  // 	} // ENDIF z-data

  // 	if (dual_energy_) {
  // 	  SubgridFluxes[n]->LeftFluxes [GENum][dim][offset] = gef[lindex];
  // 	  SubgridFluxes[n]->RightFluxes[GENum][dim][offset] = gef[rindex];
  // 	} // ENDIF dual_energy_

  // 	for (ncolour = 0; ncolour < nc; ncolour++) {
  // 	  clindex = (k + ncolour * mz) * GridDimension[dim] +
  // 	    lface;
  // 	  crindex = (k + ncolour * mz) * GridDimension[dim] +
  // 	    rface;

  // 	  SubgridFluxes[n]->LeftFluxes [colnum[ncolour]][dim][offset] = 
  // 	    colf[clindex];
  // 	  SubgridFluxes[n]->RightFluxes[colnum[ncolour]][dim][offset] = 
  // 	    colf[crindex];
  // 	} // ENDFOR ncolour

  //     } // ENDFOR J

  //   } // ENDIF k inside

  // } // ENDFOR n

  // /* Copy from slice to field */

  // for (k = 0; k < mz; k++) {
  //   index2 = k * my;
  //   for (j = 0; j < my; j++) {
  //     index3 = (k*my + j)*mx + i;
  //     de[index3] = dslice[index2+j];
  //     et[index3] = eslice[index2+j];
  //     vx[index3] = wslice[index2+j];
  //   } // ENDFOR i

  //   if (GridRank > 1)
  //     for (j = 0; j < my; j++) {
  // 	index3 = (k*my + j)*mx + i;
  // 	vy[index3] = uslice[index2+j];
  //     }

  //   if (GridRank > 2)
  //     for (j = 0; j < my; j++) {
  // 	index3 = (k*my + j)*mx + i;
  // 	vz[index3] = vslice[index2+j];
  //     }

  //   if (dual_energy_)
  //     for (j = 0; j < my; j++) {
  // 	index3 = (k*my + j)*mx + i;
  // 	ei[index3] = geslice[index2+j];
  //     }

  //   for (n = 0; n < nc; n++) {
  //     index2 = (n*mz + k) * my;
  //     for (j = 0; j < my; j++) {
  // 	index3 = (k*my + j) * mx + i;
  // 	BaryonField[colnum[n]][index3] = colslice[index2+j];
  //     }
  //   } // ENDFOR colours    

  // } // ENDFOR j

  // /* Delete all temporary slices */

  // delete [] dslice;
  // delete [] eslice;
  // delete [] uslice;
  // delete [] vslice;
  // delete [] wslice;
  // delete [] pslice;
  // if (gravity_)
  //   delete [] grslice;
  // if (dual_energy_)
  //   delete [] geslice;
  // if (nc > 0)
  //   delete [] colslice;

  // delete [] dls;
  // delete [] drs;
  // delete [] flatten;
  // delete [] pbar;
  // delete [] pls;
  // delete [] prs;
  // delete [] ubar;
  // delete [] uls;
  // delete [] urs;
  // delete [] vls;
  // delete [] vrs;
  // delete [] gels;
  // delete [] gers;
  // delete [] wls;
  // delete [] wrs;
  // delete [] diffcoef;
  // delete [] df;
  // delete [] ef;
  // delete [] uf;
  // delete [] vf;
  // delete [] wf;
  // delete [] gef;
  // delete [] ges;
  // delete [] colf;
  // delete [] colls;
  // delete [] colrs;

  // return SUCCESS;

}

//----------------------------------------------------------------------

void EnzoMethodHydro::ppm_euler_z_ (Block * block, int iy)
{

  // int dim = 2, idim = 0, jdim = 1;
  // int dim_p1 = dim+1;   // To match definition in calcdiss

  // /* Find fields: density, total energy, velocity1-3. */
  
  // int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;

  // this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
  // 				   Vel3Num, TENum);

  // int nxz, nyz, nzz, ixyz;

  // nxz = GridEndIndex[0] - GridStartIndex[0] + 1;
  // nyz = GridEndIndex[1] - GridStartIndex[1] + 1;
  // nzz = GridEndIndex[2] - GridStartIndex[2] + 1;

  // float MinimumPressure = tiny_number;
  
  // // Copy from field to slice

  // float *dslice, *eslice, *uslice, *vslice, *wslice, *grslice, *geslice, 
  //   *colslice, *pslice;

  // int size = mz * mx;
  // dslice = new float[size];  
  // eslice = new float[size];  
  // uslice = new float[size];  
  // vslice = new float[size];  
  // wslice = new float[size];  
  // pslice = new float[size];  
  // if (gravity_) {
  //   grslice = new float[size];  
  // }
  // if (dual_energy_) {
  //   geslice = new float[size];  
  // }
  // if (NumberOfColours > 0) {
  //   colslice = new float[NumberOfColours * size];  
  // }

  // int i, k, n, ncolour, index2, index3;

  // for (i = 0; i < mx; i++) {
  //   index2 = i * mz;
  //   for (k = 0; k < mz; k++) {
  //     index3 = (k*my + j) * mx + i;
  //     dslice[index2+k] = de[index3];
  //     eslice[index2+k] = et[index3];
  //     pslice[index2+k] = pressure[index3];
  //     vslice[index2+k] = vx[index3];
  //   } // ENDFOR i

  //   // Set velocities to zero if rank < 3 since hydro routines are
  //   // hard-coded for 3-d

  //   if (GridRank > 1)
  //     for (k = 0; k < mz; k++) {
  // 	index3 = (k*my + j) * mx + i;
  // 	wslice[index2+k] = vy[index3];
  //     }
  //   else
  //     for (k = 0; k < mz; k++)
  // 	wslice[index2+k] = 0;
  
  //   if (GridRank > 2)
  //     for (k = 0; k < mz; k++) {
  // 	index3 = (k*my + j) * mx + i;
  // 	uslice[index2+k] = vz[index3];
  //     }
  //   else
  //     for (k = 0; k < mz; k++)
  // 	uslice[index2+k] = 0;

  //   if (gravity_)
  //     for (k = 0; k < mz; k++) {
  // 	index3 = (k*my + j) * mx + i;
  // 	grslice[index2+k] = AccelerationField[dim][index3];
  //     }

  //   if (dual_energy_)
  //     for (k = 0; k < mz; k++) {
  // 	index3 = (k*my + j) * mx + i;
  // 	geslice[index2+k] = ei[index3];
  //     }

  //   for (n = 0; n < NumberOfColours; n++) {
  //     index2 = (n*mx + i) * mz;
  //     for (k = 0; k < mz; k++) {
  // 	index3 = (k*my + j) * mx + i;
  // 	colslice[index2+k] = BaryonField[colnum[n]][index3];
  //     }
  //   } // ENDFOR colours
  // } // ENDFOR j

  // /* Allocate memory for temporaries used in solver */

  // float *dls, *drs, *flatten, *pbar, *pls, *prs, *ubar, *uls, *urs, *vls, 
  //   *vrs, *gels, *gers, *wls, *wrs, *diffcoef, *df, *ef, *uf, *vf, *wf, *gef,
  //   *ges, *colf, *colls, *colrs;

  // dls = new float[size];	
  // drs = new float[size];	
  // flatten = new float[size];	
  // pbar = new float[size];	
  // pls = new float[size];	
  // prs = new float[size];	
  // ubar = new float[size];	
  // uls = new float[size];	
  // urs = new float[size];	
  // vls = new float[size];	
  // vrs = new float[size];	
  // gels = new float[size];	
  // gers = new float[size];	
  // wls = new float[size];	
  // wrs = new float[size];	
  // diffcoef = new float[size];	
  // df = new float[size];		
  // ef = new float[size];		
  // uf = new float[size];		
  // vf = new float[size];		
  // wf = new float[size];		
  // gef = new float[size];	
  // ges = new float[size];
  // colf = new float[NumberOfColours*size];  
  // colls = new float[NumberOfColours*size];  
  // colrs = new float[NumberOfColours*size];  

  // /* Convert start and end indexes into 1-based for FORTRAN */

  // int is, ie, js, je, is_m3, ie_p3, ie_p1, k_p1;

  // is = GridStartIndex[2] + 1;
  // ie = GridEndIndex[2] + 1;
  // js = 1;
  // je = mx;
  // is_m3 = is - 3;
  // ie_p1 = ie + 1;
  // ie_p3 = ie + 3;
  // k_p1 = j + 1;

  // /* Compute the pressure on a slice */
  // /*
  // if (dual_energy_)
  //   FORTRAN_NAME(woc_pgas2d_dual)(dslice, eslice, geslice, pslice, uslice, vslice, 
  // 			      wslice, &dual_energy_eta1_, 
  // 			      &dual_energy_eta2_, &mz, 
  // 			      &mx, &is_m3, &ie_p3, &js, &je, 
  // 			      &gamma_, &ppm_pressure_floor_);
  // else
  //   FORTRAN_NAME(woc_pgas2d)(dslice, eslice, pslice, uslice, vslice, 
  // 			 wslice, &mz, &mx, 
  // 			 &is_m3, &ie_p3, &js, &je, &gamma_, &ppm_pressure_floor_);
  // */
  // /* If requested, compute diffusion and slope flattening coefficients */

  // if (ppm_diffusion_ != 0 || PPMFlatteningParameter != 0)
  //   FORTRAN_NAME(woc_calcdiss)(dslice, eslice, uslice, vx,
  // 			   vy, pslice, CellWidthTemp[2],
  // 			   hxa, CellWidthTemp[1], 
  // 			   &mz, &mx, 
  // 			   &my, &is, &ie, &js, &je, &k_p1,
  // 			   &nyz, &dim_p1, &mx,
  // 			   &my, &mz,
  // 			   &dt, &gamma_, &ppm_diffusion_,
  // 			   &PPMFlatteningParameter, diffcoef, flatten);

  // /* Compute Eulerian left and right states at zone edges via interpolation */

  // if (ReconstructionMethod == PPM)
  //   FORTRAN_NAME(woc_inteuler)(dslice, pslice, &gravity_, grslice, geslice, uslice,
  // 			   vslice, wslice, CellWidthTemp[2], flatten,
  // 			   &mz, &mx,
  // 			   &is, &ie, &js, &je, &dual_energy_, 
  // 			   &dual_energy_eta1_, &dual_energy_eta2_,
  // 			   &PPMSteepeningParameter, &PPMFlatteningParameter,
  // 			   &ConservativeReconstruction, &PositiveReconstruction,
  // 			   &dt, &gamma_, &ppm_pressure_free_, 
  // 			   dls, drs, pls, prs, gels, gers, uls, urs, vls, vrs,
  // 			   wls, wrs, &NumberOfColours, colslice, colls, colrs);

  // /* Compute (Lagrangian part of the) Riemann problem at each zone boundary */

  // switch (RiemannSolver) {
  // case TwoShock:
  //   FORTRAN_NAME(woc_twoshock)(dls, drs, pls, prs, uls, urs,
  // 			   &mz, &mx,
  // 			   &is, &ie_p1, &js, &je,
  // 			   &dt, &gamma_, &ppm_pressure_floor_, &ppm_pressure_free_,
  // 			   pbar, ubar, &gravity_, grslice,
  // 			   &dual_energy_, &dual_energy_eta1_);
    
  //   FORTRAN_NAME(woc_flux_twoshock)(dslice, eslice, geslice, uslice, vslice, wslice,
  // 				CellWidthTemp[2], diffcoef, 
  // 				&mz, &mx,
  // 				&is, &ie, &js, &je, &dt, &gamma_,
  // 				&ppm_diffusion_, &dual_energy_,
  // 				&dual_energy_eta1_,
  // 				&riemann_solver_fallback,
  // 				dls, drs, pls, prs, gels, gers, uls, urs,
  // 				vls, vrs, wls, wrs, pbar, ubar,
  // 				df, ef, uf, vf, wf, gef, ges,
  // 				&NumberOfColours, colslice, colls, colrs, colf);
  //   break;

  // case HLL:
  //   FORTRAN_NAME(woc_flux_hll)(dslice, eslice, geslice, uslice, vslice, wslice,
  // 			   CellWidthTemp[2], diffcoef, 
  // 			   &mz, &mx,
  // 			   &is, &ie, &js, &je, &dt, &gamma_,
  // 			   &ppm_diffusion_, &dual_energy_,
  // 			   &dual_energy_eta1_,
  // 			   &riemann_solver_fallback,
  // 			   dls, drs, pls, prs, uls, urs,
  // 			   vls, vrs, wls, wrs, gels, gers,
  // 			   df, uf, vf, wf, ef, gef, ges,
  // 			   &NumberOfColours, colslice, colls, colrs, colf);
  //   break;

  // case HLLC:
  //   FORTRAN_NAME(woc_flux_hllc)(dslice, eslice, geslice, uslice, vslice, wslice,
  // 			    CellWidthTemp[2], diffcoef, 
  // 			    &mz, &mx,
  // 			    &is, &ie, &js, &je, &dt, &gamma_,
  // 			    &ppm_diffusion_, &dual_energy_,
  // 			    &dual_energy_eta1_,
  // 			    &riemann_solver_fallback,
  // 			    dls, drs, pls, prs, uls, urs,
  // 			    vls, vrs, wls, wrs, gels, gers,
  // 			    df, uf, vf, wf, ef, gef, ges,
  // 			    &NumberOfColours, colslice, colls, colrs, colf);
  //   break;

  // default:
  //   for (int index = 0; index < size; index++) {
  //     df[index] = 0;
  //     ef[index] = 0;
  //     uf[index] = 0;
  //     vf[index] = 0;
  //     wf[index] = 0;
  //     gef[index] = 0;
  //     ges[index] = 0;
  //   }
  //   break;

  // } // ENDCASE

  // /* Compute Eulerian fluxes and update zone-centered quantities */

  // FORTRAN_NAME(woc_euler)(dslice, eslice, grslice, geslice, uslice, vslice, wslice,
  // 		      CellWidthTemp[2], diffcoef, 
  // 		      &mz, &mx, 
  // 		      &is, &ie, &js, &je, &dt, &gamma_, 
  // 		      &ppm_diffusion_, &gravity_, &dual_energy_, 
  // 		      &dual_energy_eta1_, &dual_energy_eta2_,
  // 		      df, ef, uf, vf, wf, gef, ges,
  // 		      &NumberOfColours, colslice, colf, &density_floor_);

  // /* If necessary, recompute the pressure to correctly set ge and e */

  // if (dual_energy_)
  //   FORTRAN_NAME(woc_pgas2d_dual)(dslice, eslice, geslice, pslice, uslice, vslice, 
  // 			      wslice, &dual_energy_eta1_, 
  // 			      &dual_energy_eta2_, &mz, 
  // 			      &mx, &is_m3, &ie_p3, &js, &je, 
  // 			      &gamma_, &ppm_pressure_floor_);

  // /* Check this slice against the list of subgrids (all subgrid
  //    quantities are zero-based) */

  // int jstart, jend, offset, nfi, lface, rface, lindex, rindex, 
  //   fistart, fiend, fjstart, fjend, clindex, crindex;
  
  // for (n = 0; n < NumberOfSubgrids; n++) {

  //   fistart = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][idim] - 
  //     GridGlobalStart[idim];
  //   fiend = SubgridFluxes[n]->RightFluxEndGlobalIndex[dim][idim] -
  //     GridGlobalStart[idim];
  //   fjstart = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][jdim] - 
  //     GridGlobalStart[jdim];
  //   fjend = SubgridFluxes[n]->RightFluxEndGlobalIndex[dim][jdim] -
  //     GridGlobalStart[jdim];

  //   if (j >= fjstart && j <= fjend) {

  //     nfi = fiend - fistart + 1;
  //     for (i = fistart; i <= fiend; i++) {

  // 	offset = (i-fistart) + (j-fjstart)*nfi;

  // 	lface = SubgridFluxes[n]->LeftFluxStartGlobalIndex[dim][dim] -
  // 	  GridGlobalStart[dim];
  // 	lindex = i * GridDimension[dim] + lface;

  // 	rface = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][dim] -
  // 	  GridGlobalStart[dim] + 1;
  // 	rindex = i * GridDimension[dim] + rface;	

  // 	SubgridFluxes[n]->LeftFluxes [DensNum][dim][offset] = df[lindex];
  // 	SubgridFluxes[n]->RightFluxes[DensNum][dim][offset] = df[rindex];
  // 	SubgridFluxes[n]->LeftFluxes [TENum][dim][offset]   = ef[lindex];
  // 	SubgridFluxes[n]->RightFluxes[TENum][dim][offset]   = ef[rindex];

  // 	if (nxz > 1) {
  // 	  SubgridFluxes[n]->LeftFluxes [Vel1Num][dim][offset] = vf[lindex];
  // 	  SubgridFluxes[n]->RightFluxes[Vel1Num][dim][offset] = vf[rindex];
  // 	} // ENDIF x-data

  // 	if (nyz > 1) {
  // 	  SubgridFluxes[n]->LeftFluxes [Vel2Num][dim][offset] = wf[lindex];
  // 	  SubgridFluxes[n]->RightFluxes[Vel2Num][dim][offset] = wf[rindex];
  // 	} // ENDIF y-data

  // 	SubgridFluxes[n]->LeftFluxes [Vel3Num][dim][offset] = uf[lindex];
  // 	SubgridFluxes[n]->RightFluxes[Vel3Num][dim][offset] = uf[rindex];

  // 	if (dual_energy_) {
  // 	  SubgridFluxes[n]->LeftFluxes [GENum][dim][offset] = gef[lindex];
  // 	  SubgridFluxes[n]->RightFluxes[GENum][dim][offset] = gef[rindex];
  // 	} // ENDIF dual_energy_

  // 	for (ncolour = 0; ncolour < NumberOfColours; ncolour++) {
  // 	  clindex = (i + ncolour * mx) * GridDimension[dim] +
  // 	    lface;
  // 	  crindex = (i + ncolour * mx) * GridDimension[dim] +
  // 	    rface;

  // 	  SubgridFluxes[n]->LeftFluxes [colnum[ncolour]][dim][offset] = 
  // 	    colf[clindex];
  // 	  SubgridFluxes[n]->RightFluxes[colnum[ncolour]][dim][offset] = 
  // 	    colf[crindex];
  // 	} // ENDFOR ncolour

  //     } // ENDFOR J

  //   } // ENDIF k inside

  // } // ENDFOR n

  // /* Copy from slice to field */

  // for (i = 0; i < mx; i++) {
  //   index2 = i * mz;
  //   for (k = 0; k < mz; k++) {
  //     index3 = (k*my + j)*mx + i;
  //     de[index3] = dslice[index2+k];
  //     et[index3] = eslice[index2+k];
  //     vx[index3] = vslice[index2+k];
  //   } // ENDFOR i

  //   if (GridRank > 1)
  //     for (k = 0; k < mz; k++) {
  // 	index3 = (k*my + j)*mx + i;
  // 	vy[index3] = wslice[index2+k];
  //     }

  //   if (GridRank > 2)
  //     for (k = 0; k < mz; k++) {
  // 	index3 = (k*my + j)*mx + i;
  // 	vz[index3] = uslice[index2+k];
  //     }

  //   if (dual_energy_)
  //     for (k = 0; k < mz; k++) {
  // 	index3 = (k*my + j)*mx + i;
  // 	ei[index3] = geslice[index2+k];
  //     }

  //   for (n = 0; n < NumberOfColours; n++) {
  //     index2 = (n*mx + i) * mz;
  //     for (k = 0; k < mz; k++) {
  // 	index3 = (k*my + j) * mx + i;
  // 	BaryonField[colnum[n]][index3] = colslice[index2+k];
  //     }
  //   } // ENDFOR colours

  // } // ENDFOR j

  // /* Delete all temporary slices */

  // delete [] dslice;
  // delete [] eslice;
  // delete [] uslice;
  // delete [] vslice;
  // delete [] wslice;
  // delete [] pslice;
  // if (gravity_)
  //   delete [] grslice;
  // if (dual_energy_)
  //   delete [] geslice;
  // if (NumberOfColours > 0)
  //   delete [] colslice;

  // delete [] dls;
  // delete [] drs;
  // delete [] flatten;
  // delete [] pbar;
  // delete [] pls;
  // delete [] prs;
  // delete [] ubar;
  // delete [] uls;
  // delete [] urs;
  // delete [] vls;
  // delete [] vrs;
  // delete [] gels;
  // delete [] gers;
  // delete [] wls;
  // delete [] wrs;
  // delete [] diffcoef;
  // delete [] df;
  // delete [] ef;
  // delete [] uf;
  // delete [] vf;
  // delete [] wf;
  // delete [] gef;
  // delete [] ges;
  // delete [] colf;
  // delete [] colls;
  // delete [] colrs;

  // return SUCCESS;

  
}

//----------------------------------------------------------------------

double EnzoMethodHydro::timestep ( Block * block ) const throw()
{

  double dt = std::numeric_limits<double>::max();
  
  return dt;
}
