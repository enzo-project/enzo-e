// $Id$
// See LICENSE_ENZO file for license and copyright information

/// @file      Grid_SolveHydroEquations.cpp
/// @author    Greg Bryan
/// @date      November, 1994
/// @brief     Solve the hydro equations, saving subgrid fluxes

#include "cello.hpp"

#include "enzo.hpp"

static bool first_time = true;

int EnzoDescr::SolveHydroEquations (DataBlock * data_block,
				    int CycleNumber, enzo_float dt)
{
  if (first_time && data_block) {
    first_time = false;
    WARNING("EnzoDescr::SolveHydroEquations",
	    "Ignoring data_block input parameter");
  }

  // @@@@ ASSUME UNIGIRD PROBLEM @@@@

  int NumberOfSubgrids = 0;

  /* initialize */

  int dim, i,j,  size;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, coloff[MAX_COLOR];
  long long GridGlobalStart[MAX_DIMENSION];
  enzo_float a = 1, dadt;

  enzo_float *colourpt = NULL;

  /* Compute size (in enzo_floats) of the current grid. */

  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Find fields: density, total energy, velocity1-3. */

  if (IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				 Vel3Num, TENum) == ENZO_FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return ENZO_FAIL;
  }

  /* If multi-species being used, then treat them as colour variables
     (note: the solver has been modified to treat these as density vars). */

  int NumberOfColours = 0;
  int ColourNum = -1;

  if (MultiSpecies > 0) {

    NumberOfColours = 6 + 3*(MultiSpecies-1);

    if ((ColourNum =
	 FindField(ElectronDensity, FieldType, NumberOfBaryonFields)) < 0) {
      fprintf(stderr, "Could not find ElectronDensity.\n");
      return ENZO_FAIL;
    }

    /* Set Offsets from the first field (here assumed to be the electron
       density) to the other multi-species fields. */

    colourpt = BaryonField[ColourNum];

    for (i = 0; i < NumberOfColours; i++)
      coloff[i] = BaryonField[ColourNum+i] - colourpt;

  }

  // Initialize color field

  //  NumberOfColours = 1;
  //  ColourNum       = field_color;
  //  colourpt        = BaryonField[field_color];
  //  coloff[0]       = 0;

  /* Add metallicity as a colour variable. */

  int MetalNum;

  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) != -1) {
    if (colourpt == NULL) {
      colourpt = BaryonField[MetalNum];
      ColourNum = MetalNum;
    }

    coloff[NumberOfColours++] = BaryonField[MetalNum  ] - colourpt;
    coloff[NumberOfColours++] = BaryonField[MetalNum+1] - colourpt;
    coloff[NumberOfColours++] = BaryonField[MetalNum+2] - colourpt;

  }

  /* Get easy to handle pointers for each variable. */

  enzo_float *density     = BaryonField[DensNum];
  enzo_float *totalenergy = BaryonField[TENum];
  enzo_float *gasenergy   = BaryonField[GENum];

  /* Velocity1 must exist, but if 2 & 3 aren't present, then create blank
     buffers for them (since the solver needs to advect something). */

  enzo_float *velocity1, *velocity2, *velocity3;
  velocity1 = BaryonField[Vel1Num];

  if (GridRank > 1)
    velocity2 = BaryonField[Vel2Num];
  else {
    velocity2 = new enzo_float[size];
    for (i = 0; i < size; i++)
      velocity2[i] = 0.0;
  }

  if (GridRank > 2)
    velocity3 = BaryonField[Vel3Num];
  else {
    velocity3 = new enzo_float[size];
    for (i = 0; i < size; i++)
      velocity3[i] = 0.0;
  }

  /* Determine if Gamma should be a scalar or a field. */

  enzo_float *GammaField;
  GammaField = new enzo_float[1];
  GammaField[0] = Gamma;

  /* Set minimum support. */

  enzo_float MinimumSupportEnergyCoefficient = 0;
  if (UseMinimumPressureSupport == TRUE)
    if (SetMinimumSupport(MinimumSupportEnergyCoefficient) == ENZO_FAIL) {
      fprintf(stderr, "Error in grid->SetMinimumSupport,\n");
      return ENZO_FAIL;
    }

  /* allocate space for fluxes */

  SubgridFluxes = new fluxes *[NumberOfSubgrids];

  for (i = 0; i < NumberOfSubgrids; i++) {
    SubgridFluxes[i] = new fluxes;
      
    for (dim = 0; dim < GridRank; dim++)  {

      /* compute size (in enzo_floats) of flux storage */

      size = 1;
      for (j = 0; j < GridRank; j++)
	size *= SubgridFluxes[i]->LeftFluxEndGlobalIndex[dim][j] -
	  SubgridFluxes[i]->LeftFluxStartGlobalIndex[dim][j] + 1;

      /* set unused dims (for the solver, which is hardwired for 3d). */

      for (j = GridRank; j < 3; j++) {
	SubgridFluxes[i]->LeftFluxStartGlobalIndex[dim][j] = 0;
	SubgridFluxes[i]->LeftFluxEndGlobalIndex[dim][j] = 0;
	SubgridFluxes[i]->RightFluxStartGlobalIndex[dim][j] = 0;
	SubgridFluxes[i]->RightFluxEndGlobalIndex[dim][j] = 0;
      }

      /* Allocate space (if necessary). */

      for (int field = 0; field < NumberOfBaryonFields; field++) {
	if (SubgridFluxes[i]->LeftFluxes[field][dim] == NULL)
	  SubgridFluxes[i]->LeftFluxes[field][dim]  = new enzo_float[size];
	if (SubgridFluxes[i]->RightFluxes[field][dim] == NULL)
	  SubgridFluxes[i]->RightFluxes[field][dim] = new enzo_float[size];
	for (int n = 0; n < size; n++) {
	  SubgridFluxes[i]->LeftFluxes[field][dim][n] = 0;
	  SubgridFluxes[i]->RightFluxes[field][dim][n] = 0;
	}
      }

      for (int field = NumberOfBaryonFields; field < MAX_NUMBER_OF_BARYON_FIELDS;
	   field++) {
	SubgridFluxes[i]->LeftFluxes[field][dim] = NULL;
	SubgridFluxes[i]->RightFluxes[field][dim] = NULL;
      }

    }  // next dimension

    /* make things pretty */

    for (dim = GridRank; dim < 3; dim++)
      for (int field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++) {
	SubgridFluxes[i]->LeftFluxes[field][dim] = NULL;
	SubgridFluxes[i]->RightFluxes[field][dim] = NULL;
      }

  } // end of loop over subgrids

    /* compute global start index for left edge of entire grid
       (including boundary zones) */

  for (dim = 0; dim < GridRank; dim++)
    GridGlobalStart[dim] =
      NINT((GridLeftEdge[dim] - DomainLeftEdge[dim])/CellWidth[dim]) -
      GridStartIndex[dim];

  /* fix grid quantities so they are defined to at least 3 dims */

  for (i = GridRank; i < 3; i++) {
    GridDimension[i]   = 1;
    GridStartIndex[i]  = 0;
    GridEndIndex[i]    = 0;
    //      GridGlobalStart[i] = 0;
  }

  /* allocate temporary space for solver (enough to fit 31 of the largest
     possible 2d slices plus 4*NumberOfColours). */

  int tempsize = MAX(MAX(GridDimension[0]*GridDimension[1],
			 GridDimension[1]*GridDimension[2]),
		     GridDimension[2]*GridDimension[0]  );
  enzo_float *temp = new enzo_float[tempsize*(31+NumberOfColours*4)];

  /* create and fill in arrays which are easier for the solver to
     understand. */

  size = NumberOfSubgrids*3*(18+2*NumberOfColours) + 1;
  int * array = new int[size];
  for (int i=0; i<size; i++) array[i] = 0;

  int *leftface  = array + NumberOfSubgrids*3*0;
  int *rightface = array + NumberOfSubgrids*3*1;
  int *istart    = array + NumberOfSubgrids*3*2;
  int *jstart    = array + NumberOfSubgrids*3*3;
  int *iend      = array + NumberOfSubgrids*3*4;
  int *jend      = array + NumberOfSubgrids*3*5;
  int *dindex    = array + NumberOfSubgrids*3*6;
  int *Eindex    = array + NumberOfSubgrids*3*8;
  int *uindex    = array + NumberOfSubgrids*3*10;
  int *vindex    = array + NumberOfSubgrids*3*12;
  int *windex    = array + NumberOfSubgrids*3*14;
  int *geindex   = array + NumberOfSubgrids*3*16;
  int *colindex  = array + NumberOfSubgrids*3*18;

  enzo_float standard[1];
  //    enzo_float *standard = SubgridFluxes[0]->LeftFluxes[0][0];

  for (int subgrid = 0; subgrid < NumberOfSubgrids; subgrid++)
    for (dim = 0; dim < GridRank; dim++) {

      /* Set i,j dimensions of 2d flux slice (this works even if we
	 are in 1 or 2d) the correspond to the dimensions of the global
	 indicies.  I.e. for dim = 0, the plane is dims 1,2
	 for dim = 1, the plane is dims 0,2
	 for dim = 2, the plane is dims 0,1 */

      int idim = (dim == 0) ? 1 : 0;
      int jdim = (dim == 2) ? 1 : 2;

      /* Set the index (along the dimension perpindicular to the flux
	 plane) of the left and right flux planes.  The index is zero
	 based from the left side of the entire grid. */

      leftface[subgrid*3+dim] =
	SubgridFluxes[subgrid]->LeftFluxStartGlobalIndex[dim][dim] -
	GridGlobalStart[dim];
      rightface[subgrid*3+dim] =
	SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[dim][dim] -
	GridGlobalStart[dim];   // (+1 done by fortran code)

      /* set the start and end indicies (zero based on entire grid)
	 of the 2d flux plane. */

      istart[subgrid*3+dim] =
	SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[dim][idim] -
	GridGlobalStart[idim];
      jstart[subgrid*3+dim] =
	SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[dim][jdim] -
	GridGlobalStart[jdim];
      iend[subgrid*3+dim] =
	SubgridFluxes[subgrid]->RightFluxEndGlobalIndex[dim][idim] -
	GridGlobalStart[idim];
      jend[subgrid*3+dim] =
	SubgridFluxes[subgrid]->RightFluxEndGlobalIndex[dim][jdim] -
	GridGlobalStart[jdim];

      /* Compute offset from the standard pointer to the start of
	 each set of flux data.  This is done to compensate for
	 fortran's inability to handle arrays of pointers or structs.
	 NOTE: This pointer arithmatic is illegal; some other way should
	 be found to do it (like write higher level ppm stuff in c++). */

      dindex[subgrid*6+dim*2] =
	SubgridFluxes[subgrid]->LeftFluxes[DensNum][dim] - standard;
      dindex[subgrid*6+dim*2+1] =
	SubgridFluxes[subgrid]->RightFluxes[DensNum][dim] - standard;
      Eindex[subgrid*6+dim*2] =
	SubgridFluxes[subgrid]->LeftFluxes[TENum][dim] - standard;
      Eindex[subgrid*6+dim*2+1] =
	SubgridFluxes[subgrid]->RightFluxes[TENum][dim] - standard;
      uindex[subgrid*6+dim*2] =
	SubgridFluxes[subgrid]->LeftFluxes[Vel1Num][dim] - standard;
      uindex[subgrid*6+dim*2+1] =
	SubgridFluxes[subgrid]->RightFluxes[Vel1Num][dim] - standard;

      if (GridRank > 1) {
	vindex[subgrid*6+dim*2] =
	  SubgridFluxes[subgrid]->LeftFluxes[Vel2Num][dim] - standard;
	vindex[subgrid*6+dim*2+1] =
	  SubgridFluxes[subgrid]->RightFluxes[Vel2Num][dim] - standard;
      }
      if (GridRank > 2) {
	windex[subgrid*6+dim*2] =
	  SubgridFluxes[subgrid]->LeftFluxes[Vel3Num][dim] - standard;
	windex[subgrid*6+dim*2+1] =
	  SubgridFluxes[subgrid]->RightFluxes[Vel3Num][dim] - standard;
      }

      if (DualEnergyFormalism) {
	geindex[subgrid*6+dim*2] =
	  SubgridFluxes[subgrid]->LeftFluxes[GENum][dim] - standard;
	geindex[subgrid*6+dim*2+1] =
	  SubgridFluxes[subgrid]->RightFluxes[GENum][dim] - standard;
      }

      for (i = 0; i < NumberOfColours; i++) {
	colindex[i*NumberOfSubgrids*6 + subgrid*6 + dim*2] =
	  SubgridFluxes[subgrid]->LeftFluxes[ColourNum+i][dim] - standard;
	colindex[i*NumberOfSubgrids*6 + subgrid*6 + dim*2 + 1] =
	  SubgridFluxes[subgrid]->RightFluxes[ColourNum+i][dim] - standard;
      }

    }

  /* If using comoving coordinates, multiply dx by a(n+1/2).
     In one fell swoop, this recasts the equations solved by solver
     in comoving form (except for the expansion terms which are taken
     care of elsewhere). */

  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(Time+0.5*dt, &a, &dadt)
	== ENZO_FAIL) {
      fprintf(stderr, "Error in CsomologyComputeExpansionFactors.\n");
      return ENZO_FAIL;
    }

  /* Create a cell width array to pass (and convert to absolute coords). */

  enzo_float CellWidthTemp[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    if (dim < GridRank)
      CellWidthTemp[dim] = enzo_float(a*CellWidth[dim]);
    else
      CellWidthTemp[dim] = 1.0;
  }


  /* call a Fortran routine to actually compute the hydro equations
     on this grid.
     Notice that it is hard-wired for three dimensions, but it does
     the right thing for < 3 dimensions. */
  /* note: Start/EndIndex are zero based */

  //     printf (" density = %g\n", density[0]);
  //     printf (" totalenergy = %g\n", totalenergy[0]);
  //     printf (" velocity1 = %g\n", velocity1[0]);
  //     printf (" velocity2 = %g\n", velocity2[0]);
  //     printf (" velocity3 = %g\n", velocity3[0]);
  //     printf (" gasenergy = %g\n", gasenergy[0]);
  //     printf (" GravityOn = %d\n", GravityOn);
  //     AccelerationField[0] = density;
  //     AccelerationField[1] = density;
  //     AccelerationField[2] = density;
  //     printf (" AccelerationField = %g\n", AccelerationField[0][0]);
  //     printf (" Gamma = %g\n", Gamma);
  //     printf (" dt = %g\n", dt);
  //     printf (" CycleNumber = %d\n", CycleNumber);
  //     printf (" CellWidthTemp = %g\n", CellWidthTemp[0][0]);
  //     printf (" GridRank = %d\n", GridRank);
  //     printf (" GridDimension = %d\n", GridDimension[0]);
  //     printf (" GridStartIndex = %d\n", GridStartIndex[0]);
  //     printf (" GridEndIndex = %d\n", GridEndIndex[0]);
  //     printf (" PPMFlatteningParameter = %d\n", PPMFlatteningParameter);
  //     printf (" PressureFree = %d\n", PressureFree);
  //     printf (" PPMDiffusionParameter = %d\n", PPMDiffusionParameter);
  //     printf (" PPMSteepeningParameter = %d\n", PPMSteepeningParameter);
  //     printf (" DualEnergyFormalism = %d\n", DualEnergyFormalism);
  //     printf (" DualEnergyFormalismEta1 = %g\n", DualEnergyFormalismEta1);
  //     printf (" DualEnergyFormalismEta2 = %g\n", DualEnergyFormalismEta2);
  //     printf (" NumberOfSubgrids = %d\n", NumberOfSubgrids);
  //     printf (" leftface = %d\n", leftface[0]);
  //     printf (" rightface = %d\n", rightface[0]);
  //     printf (" istart = %d\n", istart[0]);
  //     printf (" iend = %d\n", iend[0]);
  //     printf (" jstart = %d\n", jstart[0]);
  //     printf (" jend = %d\n", jend[0]);
  //     printf (" standard = %g\n", standard[0]);
  //     printf (" dindex = %d\n", dindex[0]);
  //     printf (" Eindex = %d\n", Eindex[0]);
  //     printf (" uindex = %d\n", uindex[0]);
  //     printf (" vindex = %d\n", vindex[0]);
  //     printf (" windex = %d\n", windex[0]);
  //     printf (" geindex = %d\n", geindex[0]);
  //     printf (" temp = %g\n", temp[0]);
  //     printf (" NumberOfColours = %d\n", NumberOfColours);
  //     printf (" colourpt = %g\n", colourpt[0]);
  //     printf (" coloff = %d\n", coloff[0]);
  //     printf (" colindex = %d\n", colindex[0]);
  FORTRAN_NAME(ppm_de)
    (
     density, totalenergy, velocity1, velocity2, velocity3,
     gasenergy,
     &GravityOn, AccelerationField[0],
     AccelerationField[1],
     AccelerationField[2],
     &Gamma, &dt, &CycleNumber,
     &CellWidthTemp[0], &CellWidthTemp[1], &CellWidthTemp[2],
     &GridRank, &GridDimension[0], &GridDimension[1],
     &GridDimension[2], GridStartIndex, GridEndIndex,
     &PPMFlatteningParameter,
     &PressureFree,
     &PPMDiffusionParameter, &PPMSteepeningParameter,
     &DualEnergyFormalism, &DualEnergyFormalismEta1,
     &DualEnergyFormalismEta2,
     &NumberOfSubgrids, leftface, rightface,
     istart, iend, jstart, jend,
     standard, dindex, Eindex, uindex, vindex, windex,
     geindex, temp,
     &NumberOfColours, colourpt, coloff, colindex
     );

  /* deallocate temporary space for solver */

  delete [] temp;
  if (GridRank < 3) delete [] velocity3;
  if (GridRank < 2) delete [] velocity2;

  delete [] leftface;
  delete [] GammaField;

  for (int i=0; i<NumberOfSubgrids; i++) {
    delete SubgridFluxes[i];
  }
  delete [] SubgridFluxes;

  return ENZO_SUCCESS;

}
