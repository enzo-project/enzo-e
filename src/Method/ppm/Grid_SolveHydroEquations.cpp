/***********************************************************************
/
/  GRID CLASS (SOLVE THE HYDRO EQUATIONS, SAVING SUBGRID FLUXES)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
// Solve the hydro equations with the solver, saving the subgrid fluxes
//
 
#include "cello_hydro.h"

/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int FindField(int f, int farray[], int n);
extern "C" void FORTRAN_NAME(ppm_de)(
			  float *d, float *E, float *u, float *v, float *w,
			    float *ge,
                          int *grav, float *gr_ax, float *gr_ay, float *gr_az,
			  float *gamma, float *dt, int *cycle_number,
                            float dx[], float dy[], float dz[],
			  int *rank, int *in, int *jn, int *kn,
                            int is[], int ie[],
			  int *flatten, int *ipresfree,
			  int *diff, int *steepen, int *idual,
                            float *eta1, float *eta2,
			  int *num_subgrids, int leftface[], int rightface[],
			  int istart[], int iend[], int jstart[], int jend[],
			  float *standard, int dindex[], int Eindex[],
			  int uindex[], int vindex[], int windex[],
			    int geindex[], float *temp,
                          int *ncolour, float *colourpt, int *coloff,
                            int colindex[]);
 
int SolveHydroEquations
(
 int   CycleNumber, 
 float dt)
{

  int NumberOfSubgrids = 0;

  if (NumberOfBaryonFields > 0) {
 
    /* initialize */
 
    int dim, i, idim, j, jdim, field, size, subgrid, n;
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, coloff[MAX_COLOR];
    long long GridGlobalStart[MAX_DIMENSION];
    FLOAT a = 1, dadt;
 
    float *colourpt = NULL;
 
    /* Compute size (in floats) of the current grid. */
 
    size = 1;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];
 
    /* Find fields: density, total energy, velocity1-3. */
 
    if (IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == FAIL) {
      fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
      return FAIL;
    }
 
    /* If multi-species being used, then treat them as colour variables
       (note: the solver has been modified to treat these as density vars). */
 
    int NumberOfColours = 0, ColourNum;
 
    if (MultiSpecies > 0) {
 
      NumberOfColours = 6 + 3*(MultiSpecies-1);
 
      if ((ColourNum =
	   FindField(ElectronDensity, FieldType, NumberOfBaryonFields)) < 0) {
	fprintf(stderr, "Could not find ElectronDensity.\n");
	return FAIL;
      }
 
      /* Set Offsets from the first field (here assumed to be the electron
	 density) to the other multi-species fields. */
 
      colourpt = BaryonField[ColourNum];
 
      for (i = 0; i < NumberOfColours; i++)
	coloff[i] = BaryonField[ColourNum+i] - colourpt;
 
    }

    // Initialize color field
 
    NumberOfColours = 1;
    ColourNum       = field_color;
    colourpt        = BaryonField[field_color];
    coloff[0]       = 0;

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
 
    float *density     = BaryonField[DensNum];
    float *totalenergy = BaryonField[TENum];
    float *gasenergy   = BaryonField[GENum];
 
    /* Velocity1 must exist, but if 2 & 3 aren't present, then create blank
       buffers for them (since the solver needs to advect something). */
 
    float *velocity1, *velocity2, *velocity3;
    velocity1 = BaryonField[Vel1Num];
 
    if (GridRank > 1)
      velocity2 = BaryonField[Vel2Num];
    else {
      velocity2 = new float[size];
      for (i = 0; i < size; i++)
        velocity2[i] = 0.0;
    }
 
    if (GridRank > 2)
      velocity3 = BaryonField[Vel3Num];
    else {
      velocity3 = new float[size];
      for (i = 0; i < size; i++)
        velocity3[i] = 0.0;
    }
 
    /* Determine if Gamma should be a scalar or a field. */
 
    int UseGammaField = FALSE;
    float *GammaField;
    GammaField = new float[1];
    GammaField[0] = Gamma;
 
    /* Set minimum support. */
 
    float MinimumSupportEnergyCoefficient = 0;
    if (UseMinimumPressureSupport == TRUE)
      if (SetMinimumSupport(MinimumSupportEnergyCoefficient) == FAIL) {
	fprintf(stderr, "Error in grid->SetMinimumSupport,\n");
	return FAIL;
      }
 
//     /* allocate space for fluxes */
 
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
 
    /* compute global start index for left edge of entire grid
       (including boundary zones) */
 
    for (dim = 0; dim < GridRank; dim++)
      GridGlobalStart[dim] =
	nint((GridLeftEdge[dim] - DomainLeftEdge[dim])/(*(CellWidth[dim]))) -
	GridStartIndex[dim];
 
    /* fix grid quantities so they are defined to at least 3 dims */
 
    for (i = GridRank; i < 3; i++) {
      GridDimension[i]   = 1;
      GridStartIndex[i]  = 0;
      GridEndIndex[i]    = 0;
      GridGlobalStart[i] = 0;
    }
 
    /* allocate temporary space for solver (enough to fit 31 of the largest
       possible 2d slices plus 4*NumberOfColours). */
 
    int tempsize = max(max(GridDimension[0]*GridDimension[1],
                           GridDimension[1]*GridDimension[2]),
		           GridDimension[2]*GridDimension[0]  );
    float *temp = new float[tempsize*(31+NumberOfColours*4)];
 
    /* create and fill in arrays which are easiler for the solver to
       understand. */
 
    int *leftface  = new int[NumberOfSubgrids*3*(18+2*NumberOfColours)];
    int *rightface = leftface + NumberOfSubgrids*3*1;
    int *istart    = leftface + NumberOfSubgrids*3*2;
    int *jstart    = leftface + NumberOfSubgrids*3*3;
    int *iend      = leftface + NumberOfSubgrids*3*4;
    int *jend      = leftface + NumberOfSubgrids*3*5;
    int *dindex    = leftface + NumberOfSubgrids*3*6;
    int *Eindex    = leftface + NumberOfSubgrids*3*8;
    int *uindex    = leftface + NumberOfSubgrids*3*10;
    int *vindex    = leftface + NumberOfSubgrids*3*12;
    int *windex    = leftface + NumberOfSubgrids*3*14;
    int *geindex   = leftface + NumberOfSubgrids*3*16;
    int *colindex  = leftface + NumberOfSubgrids*3*18;
    float *standard = NULL;
//     if (NumberOfSubgrids > 0) standard = SubgridFluxes[0]->LeftFluxes[0][0];
 
//     for (subgrid = 0; subgrid < NumberOfSubgrids; subgrid++)
//       for (dim = 0; dim < GridRank; dim++) {
 
//         /* Set i,j dimensions of 2d flux slice (this works even if we
//            are in 1 or 2d) the correspond to the dimensions of the global
//            indicies.  I.e. for dim = 0, the plane is dims 1,2
//                            for dim = 1, the plane is dims 0,2
//                            for dim = 2, the plane is dims 0,1 */
 
// 	idim = (dim == 0) ? 1 : 0;
// 	jdim = (dim == 2) ? 1 : 2;
 
//         /* Set the index (along the dimension perpindicular to the flux
//            plane) of the left and right flux planes.  The index is zero
//            based from the left side of the entire grid. */
 
// 	leftface[subgrid*3+dim] =
// 	  SubgridFluxes[subgrid]->LeftFluxStartGlobalIndex[dim][dim] -
// 	    GridGlobalStart[dim];
// 	rightface[subgrid*3+dim] =
// 	  SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[dim][dim] -
// 	    GridGlobalStart[dim];   // (+1 done by fortran code)
 
//         /* set the start and end indicies (zero based on entire grid)
//            of the 2d flux plane. */
 
// 	istart[subgrid*3+dim] =
// 	  SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[dim][idim] -
// 	    GridGlobalStart[idim];
// 	jstart[subgrid*3+dim] =
// 	  SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[dim][jdim] -
// 	    GridGlobalStart[jdim];
// 	iend[subgrid*3+dim] =
// 	  SubgridFluxes[subgrid]->RightFluxEndGlobalIndex[dim][idim] -
// 	    GridGlobalStart[idim];
// 	jend[subgrid*3+dim] =
// 	  SubgridFluxes[subgrid]->RightFluxEndGlobalIndex[dim][jdim] -
// 	    GridGlobalStart[jdim];
 
//         /* Compute offset from the standard pointer to the start of
//            each set of flux data.  This is done to compensate for
//            fortran's inability to handle arrays of pointers or structs.
// 	   NOTE: This pointer arithmatic is illegal; some other way should
// 	   be found to do it (like write higher level ppm stuff in c++). */
 
// 	dindex[subgrid*6+dim*2] =
// 	  SubgridFluxes[subgrid]->LeftFluxes[DensNum][dim] - standard;
// 	dindex[subgrid*6+dim*2+1] =
// 	  SubgridFluxes[subgrid]->RightFluxes[DensNum][dim] - standard;
// 	Eindex[subgrid*6+dim*2] =
// 	  SubgridFluxes[subgrid]->LeftFluxes[TENum][dim] - standard;
// 	Eindex[subgrid*6+dim*2+1] =
// 	  SubgridFluxes[subgrid]->RightFluxes[TENum][dim] - standard;
// 	uindex[subgrid*6+dim*2] =
// 	  SubgridFluxes[subgrid]->LeftFluxes[Vel1Num][dim] - standard;
// 	uindex[subgrid*6+dim*2+1] =
// 	  SubgridFluxes[subgrid]->RightFluxes[Vel1Num][dim] - standard;
 
// 	if (GridRank > 1) {
//           vindex[subgrid*6+dim*2] =
// 	    SubgridFluxes[subgrid]->LeftFluxes[Vel2Num][dim] - standard;
//           vindex[subgrid*6+dim*2+1] =
// 	    SubgridFluxes[subgrid]->RightFluxes[Vel2Num][dim] - standard;
//         }
// 	if (GridRank > 2) {
//           windex[subgrid*6+dim*2] =
// 	    SubgridFluxes[subgrid]->LeftFluxes[Vel3Num][dim] - standard;
//           windex[subgrid*6+dim*2+1] =
// 	    SubgridFluxes[subgrid]->RightFluxes[Vel3Num][dim] - standard;
//         }
 
//         if (DualEnergyFormalism) {
//           geindex[subgrid*6+dim*2] =
//             SubgridFluxes[subgrid]->LeftFluxes[GENum][dim] - standard;
//           geindex[subgrid*6+dim*2+1] =
//             SubgridFluxes[subgrid]->RightFluxes[GENum][dim] - standard;
//         }
 
// 	for (i = 0; i < NumberOfColours; i++) {
// 	  colindex[i*NumberOfSubgrids*6 + subgrid*6 + dim*2] =
// 	    SubgridFluxes[subgrid]->LeftFluxes[ColourNum+i][dim] - standard;
// 	  colindex[i*NumberOfSubgrids*6 + subgrid*6 + dim*2 + 1] =
// 	    SubgridFluxes[subgrid]->RightFluxes[ColourNum+i][dim] - standard;
// 	}
 
//       }
 
    /* If using comoving coordinates, multiply dx by a(n+1/2).
       In one fell swoop, this recasts the equations solved by solver
       in comoving form (except for the expansion terms which are taken
       care of elsewhere). */
 
    if (ComovingCoordinates)
      if (CosmologyComputeExpansionFactor(Time+0.5*dt, &a, &dadt)
	  == FAIL) {
	fprintf(stderr, "Error in CsomologyComputeExpansionFactors.\n");
	return FAIL;
      }
 
    /* Create a cell width array to pass (and convert to absolute coords). */
 
    float *CellWidthTemp[MAX_DIMENSION];
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      CellWidthTemp[dim] = new float[GridDimension[dim]];
      for (i = 0; i < GridDimension[dim]; i++)
	if (dim < GridRank)
	  CellWidthTemp[dim][i] = float(a*CellWidth[dim][i]);
	else
	  CellWidthTemp[dim][i] = 1.0;
    }
 
 
    /* call a Fortran routine to actually compute the hydro equations
       on this grid.
       Notice that it is hard-wired for three dimensions, but it does
       the right thing for < 3 dimensions. */
    /* note: Start/EndIndex are zero based */
 
    FORTRAN_NAME(ppm_de)
      (
       density, totalenergy, velocity1, velocity2, velocity3,
       gasenergy,
       &GravityOn, AccelerationField[0],
       AccelerationField[1],
       AccelerationField[2],
       &Gamma, &dt, &CycleNumber,
       CellWidthTemp[0], CellWidthTemp[1], CellWidthTemp[2],
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
 
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      delete [] CellWidthTemp[dim];
 
  }  // end: if (NumberOfBaryonFields > 0)
 
  return SUCCESS;
 
}
