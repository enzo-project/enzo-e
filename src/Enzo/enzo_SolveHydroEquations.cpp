// See LICENSE_ENZO file for license and copyright information

/// @file      enzo_SolveHydroEquations.cpp
/// @author    Greg Bryan
/// @date      November, 1994
/// @brief     Solve the hydro equations, saving subgrid fluxes

#include "cello.hpp"

#include "enzo.hpp"

int EnzoBlock::SolveHydroEquations 
(
 enzo_float time,
 enzo_float dt,
 int comoving_coordinates
 )
{
  int NumberOfSubgrids = 0;

  /* initialize */

  int dim, i,j,  size;
  enzo_float a = 1, dadt;

  Field field = data()->field();

  //------------------------------
  // Prepare colour field parameters
  //------------------------------

  // ncolour: number of colour fields

  int    ncolour  = field.groups()->size("colour");

  // colourpt: the color 'array' (contains all color fields)
  enzo_float * colourpt = (enzo_float *) field.permanent();

  // coloff: offsets into the color array (for each color field)
  int * coloff   = new int [ncolour];
  int index_colour = 0;
  for (int index_field = 0;
       index_field < field.field_count();
       index_field++) {
    std::string name = field.field_name(index_field);
    if (field.groups()->is_in(name,"colour")) {
      coloff[index_colour++] 
	= (enzo_float *)(field.values(index_field)) - colourpt;
    }

  }

  // No subgrids, so colindex is NULL
  int *colindex         = NULL;

  /* Compute size (in enzo_floats) of the current grid. */

  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Find fields: density, total energy, velocity1-3. */

  // enzo_float * species_De    = (enzo_float *) field.values("species_De");
  // enzo_float * species_HI    = (enzo_float *) field.values("species_HI");
  // enzo_float * species_HII   = (enzo_float *) field.values("species_HII");
  // enzo_float * species_HeI   = (enzo_float *) field.values("species_HeI");
  // enzo_float * species_HeII  = (enzo_float *) field.values("species_HeII");
  // enzo_float * species_HeIII = (enzo_float *) field.values("species_HeIII");
  // enzo_float * species_HM    = (enzo_float *) field.values("species_HM");
  // enzo_float * species_H2I   = (enzo_float *) field.values("species_H2I");
  // enzo_float * species_H2II  = (enzo_float *) field.values("species_H2II");
  // enzo_float * species_DI    = (enzo_float *) field.values("species_DI");
  // enzo_float * species_DII   = (enzo_float *) field.values("species_DII");
  // enzo_float * species_HDI   = (enzo_float *) field.values("species_HDI");


  enzo_float * density         = (enzo_float*) field.values("density");
  enzo_float * total_energy    = (enzo_float *)field.values("total_energy");
  enzo_float * internal_energy = (enzo_float *)field.values("internal_energy");
  enzo_float * velocity_x      = (enzo_float*) field.values("velocity_x");
  enzo_float * velocity_y      = (enzo_float*) field.values("velocity_y");
  enzo_float * velocity_z      = (enzo_float*) field.values("velocity_z");
  enzo_float * acceleration_x  = (enzo_float*) field.values("acceleration_x");
  enzo_float * acceleration_y  = (enzo_float*) field.values("acceleration_y");
  enzo_float * acceleration_z  = (enzo_float*) field.values("acceleration_z");

  /* velocity_x must exist, but if y & z aren't present, then create blank
     buffers for them (since the solver needs to advect something). */

  if (GridRank < 2) {
    velocity_y = new enzo_float[size];
    for (int i=0; i<size; i++) velocity_y[i] = 0.0;
  }
  if (GridRank < 3) {
    velocity_z = new enzo_float[size];
    for (int i=0; i<size; i++) velocity_z[i] = 0.0;
  }

  /* Determine if Gamma should be a scalar or a field. */

  enzo_float *GammaField;
  GammaField = new enzo_float[1];
  GammaField[0] = Gamma;

  /* Set minimum support. */

  enzo_float MinimumSupportEnergyCoefficient = 0;
  if (UseMinimumPressureSupport == TRUE)
    if (SetMinimumSupport(MinimumSupportEnergyCoefficient,
			  comoving_coordinates) == ENZO_FAIL) {
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

  /* fix grid quantities so they are defined to at least 3 dims */

  for (i = GridRank; i < 3; i++) {
    GridDimension[i]   = 1;
    GridStartIndex[i]  = 0;
    GridEndIndex[i]    = 0;
    //      GridGlobalStart[i] = 0;
  }

  /* allocate temporary space for solver (enough to fit 31 of the largest
     possible 2d slices plus 4*ncolour). */

  int tempsize = MAX(MAX(GridDimension[0]*GridDimension[1],
			 GridDimension[1]*GridDimension[2]),
		     GridDimension[2]*GridDimension[0]);
  enzo_float *temp = new enzo_float[tempsize*(31+ncolour*4)];

  /* create and fill in arrays which are easier for the solver to
     understand. */

  size = NumberOfSubgrids*3*(18+2*ncolour) + 1;
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

  enzo_float standard[1];
  //    enzo_float *standard = SubgridFluxes[0]->LeftFluxes[0][0];

  /* If using comoving coordinates, multiply dx by a(n+1/2).
     In one fell swoop, this recasts the equations solved by solver
     in comoving form (except for the expansion terms which are taken
     care of elsewhere). */

  if (comoving_coordinates)
    if (CosmologyComputeExpansionFactor(time+0.5*dt, &a, &dadt)
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

  //     AccelerationField[0] = density;
  //     AccelerationField[1] = density;
  //     AccelerationField[2] = density;
  int gravity_on = (acceleration_x != NULL) ? 1 : 0;
  
  FORTRAN_NAME(ppm_de)
    (
     density, total_energy, velocity_x, velocity_y, velocity_z,
     internal_energy,
     &gravity_on, 
     acceleration_x,
     acceleration_y,
     acceleration_z,
     &Gamma, &dt, &cycle_,
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
     &ncolour, colourpt, coloff, colindex
     );

  /* deallocate temporary space for solver */

  delete [] temp;
  if (GridRank < 2) delete [] velocity_y;
  if (GridRank < 3) delete [] velocity_z;

  delete [] leftface;
  delete [] GammaField;

  for (int i=0; i<NumberOfSubgrids; i++) {
    delete SubgridFluxes[i];
  }
  delete [] SubgridFluxes;
  delete [] coloff;

  return ENZO_SUCCESS;

}
