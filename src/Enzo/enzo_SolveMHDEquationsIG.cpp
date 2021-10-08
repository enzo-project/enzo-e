// See LICENSE_ENZO file for license and copyright information

/***********************************************************************
/
/  GRID CLASS (SOLVE THE MHD EQUATIONS (with PPML for Ideal Gas), SAVING SUBGRID FLUXES)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Alexei Kritsuk, April 2009 (PPML), Sept 2018 (PPML-IG)
/
/  PURPOSE:
/
/  RETURNS:
/    ENZO_SUCCESS or ENZO_FAIL
/
************************************************************************/
 
// Solve the MHD equations with the solver, saving the subgrid fluxes

#include "cello.hpp"
#include "enzo.hpp" 

//----------------------------------------------------------------------

#define DEBUG_FIELDS

#ifdef DEBUG_FIELDS
#   define CHECK_FIELD(VALUES,NAME)             \
  ASSERT1("CHECK_FIELD",                        \
          "Field %s must be defined",           \
          NAME,                                 \
          (VALUES != nullptr));

#   define FIELD_STATS(NAME,VALUES,mx,my,mz,gx,gy,gz)           \
  {                                                             \
   double avg=0.0, max=-1.0, min=1e9;                           \
   int count=0;                                                 \
   for (int iz=gz; iz<mz-gz; iz++) {                            \
     for (int iy=gy; iy<my-gy; iy++) {                          \
       for (int ix=gx; ix<mx-gx; ix++) {                        \
         const int i=ix+mx*(iy+my*iz);                          \
         avg += VALUES[i];                                      \
         max = std::max(max,VALUES[i]);                         \
         min = std::min(min,VALUES[i]);                         \
         count++;                                               \
       }                                                        \
     }                                                          \
   }                                                            \
   avg /= count;                                                \
   CkPrintf ("FIELD_STATS %s  %g %g %g\n",NAME,min,avg,max);    \
  }
#else
#   define CHECK_FIELD(VALUES,NAME) /* ... */
#   define FIELD_STATS(NAME,VALUES,mx,my,mz,gx,gy,gz) /* ... */
#endif

//----------------------------------------------------------------------

int EnzoBlock::SolveMHDEquationsIG ( enzo_float dt, enzo_float gamma )
{
 
  /* exit if not 3D */

  // @@ assert GridRank == 3
  //  if (GridRank != 3) 
  //    my_exit(EXIT_ENZO_FAILURE);

  if (NumberOfBaryonFields > 0) {
 
    /* initialize */
 
    int dim, i,  size;
    Elong_int GridGlobalStart[MAX_DIMENSION];
 
    /* Compute size (in floats) of the current grid. */
 
    const int in = cello::index_static();

    size = 1;
    for (dim = 0; dim < GridRank[in]; dim++)
      size *= GridDimension[dim];
 
    /* Get easy to handle pointers for each variable. */


    Field field = data()->field();

    enzo_float *density    = (enzo_float *) field.values ("density");
    enzo_float *velox      = (enzo_float *) field.values ("velox");
    enzo_float *veloy      = (enzo_float *) field.values ("veloy");
    enzo_float *veloz      = (enzo_float *) field.values ("veloz");
    enzo_float *bfieldx    = (enzo_float *) field.values ("bfieldx");
    enzo_float *bfieldy    = (enzo_float *) field.values ("bfieldy");
    enzo_float *bfieldz    = (enzo_float *) field.values ("bfieldz");
    enzo_float *pressure   = (enzo_float *) field.values ("pressure");

    enzo_float *dens_rx    = (enzo_float *) field.values ("dens_rx");
    enzo_float *velox_rx   = (enzo_float *) field.values ("velox_rx");
    enzo_float *veloy_rx   = (enzo_float *) field.values ("veloy_rx");
    enzo_float *veloz_rx   = (enzo_float *) field.values ("veloz_rx");
    enzo_float *bfieldx_rx = (enzo_float *) field.values ("bfieldx_rx");
    enzo_float *bfieldy_rx = (enzo_float *) field.values ("bfieldy_rx");
    enzo_float *bfieldz_rx = (enzo_float *) field.values ("bfieldz_rx");
    enzo_float *press_rx   = (enzo_float *) field.values ("press_rx");

    enzo_float *dens_ry    = (enzo_float *) field.values ("dens_ry");
    enzo_float *velox_ry   = (enzo_float *) field.values ("velox_ry");
    enzo_float *veloy_ry   = (enzo_float *) field.values ("veloy_ry");
    enzo_float *veloz_ry   = (enzo_float *) field.values ("veloz_ry");
    enzo_float *bfieldx_ry = (enzo_float *) field.values ("bfieldx_ry");
    enzo_float *bfieldy_ry = (enzo_float *) field.values ("bfieldy_ry");
    enzo_float *bfieldz_ry = (enzo_float *) field.values ("bfieldz_ry");
    enzo_float *press_ry   = (enzo_float *) field.values ("press_ry");

    enzo_float *dens_rz    = (enzo_float *) field.values ("dens_rz");
    enzo_float *velox_rz   = (enzo_float *) field.values ("velox_rz");
    enzo_float *veloy_rz   = (enzo_float *) field.values ("veloy_rz");
    enzo_float *veloz_rz   = (enzo_float *) field.values ("veloz_rz");
    enzo_float *bfieldx_rz = (enzo_float *) field.values ("bfieldx_rz");
    enzo_float *bfieldy_rz = (enzo_float *) field.values ("bfieldy_rz");
    enzo_float *bfieldz_rz = (enzo_float *) field.values ("bfieldz_rz");
    enzo_float *press_rz   = (enzo_float *) field.values ("press_rz");

    CHECK_FIELD(density,"density");
    CHECK_FIELD(velox,"velox");
    CHECK_FIELD(veloy,"veloy");
    CHECK_FIELD(veloz,"veloz");
    CHECK_FIELD(bfieldx,"bfieldx");
    CHECK_FIELD(bfieldy,"bfieldy");
    CHECK_FIELD(bfieldz,"bfieldz");
    CHECK_FIELD(pressure,"pressure");

    CHECK_FIELD(dens_rx,"dens_rx");
    CHECK_FIELD(velox_rx,"velox_rx");
    CHECK_FIELD(veloy_rx,"veloy_rx");
    CHECK_FIELD(veloz_rx,"veloz_rx");
    CHECK_FIELD(bfieldx_rx,"bfieldx_rx");
    CHECK_FIELD(bfieldy_rx,"bfieldy_rx");
    CHECK_FIELD(bfieldz_rx,"bfieldz_rx");
    CHECK_FIELD(press_rx,"press_rx");

    CHECK_FIELD(dens_ry,"dens_ry");
    CHECK_FIELD(velox_ry,"velox_ry");
    CHECK_FIELD(veloy_ry,"veloy_ry");
    CHECK_FIELD(veloz_ry,"veloz_ry");
    CHECK_FIELD(bfieldx_ry,"bfieldx_ry");
    CHECK_FIELD(bfieldy_ry,"bfieldy_ry");
    CHECK_FIELD(bfieldz_ry,"bfieldz_ry");
    CHECK_FIELD(press_ry,"press_ry");

    CHECK_FIELD(dens_rz,"dens_rz");
    CHECK_FIELD(velox_rz,"velox_rz");
    CHECK_FIELD(veloy_rz,"veloy_rz");
    CHECK_FIELD(veloz_rz,"veloz_rz");
    CHECK_FIELD(bfieldx_rz,"bfieldx_rz");
    CHECK_FIELD(bfieldy_rz,"bfieldy_rz");
    CHECK_FIELD(bfieldz_rz,"bfieldz_rz");
    CHECK_FIELD(press_rz,"press_rz");

    /* allocate space for fluxes */
 
//     for (i = 0; i < NumberOfSubgrids; i++) {
//       for (dim = 0; dim < GridRank; dim++)  {
 
// 	/* compute size (in floats) of flux storage */
 
//         size = 1;
//         for (j = 0; j < GridRank; j++)
//           size *= SubgridFluxes[i]->LeftFluxEndGlobalIndex[dim][j] -
//                   SubgridFluxes[i]->LeftFluxStartGlobalIndex[dim][j] + 1;
 
// 	/* set unused dims (for the solver, which is hardwired for 3d). */
 
//         for (j = GridRank; j < 3; j++) {
//           SubgridFluxes[i]->LeftFluxStartGlobalIndex[dim][j]  = 0;
//           SubgridFluxes[i]->LeftFluxEndGlobalIndex[dim][j]    = 0;
//           SubgridFluxes[i]->RightFluxStartGlobalIndex[dim][j] = 0;
//           SubgridFluxes[i]->RightFluxEndGlobalIndex[dim][j]   = 0;
//         }
 
// 	/* Allocate space (if necessary). */
 
//         for (field = 0; field < NumberOfBaryonFields; field++) {
// 	  if (SubgridFluxes[i]->LeftFluxes[field][dim] == NULL)
// 	    SubgridFluxes[i]->LeftFluxes[field][dim]  = new enzo_float[size];
// 	  if (SubgridFluxes[i]->RightFluxes[field][dim] == NULL)
// 	    SubgridFluxes[i]->RightFluxes[field][dim] = new enzo_float[size];
// 	  for (n = 0; n < size; n++) {
// 	    SubgridFluxes[i]->LeftFluxes[field][dim][n]  = 0;
// 	    SubgridFluxes[i]->RightFluxes[field][dim][n] = 0;
// 	  }
//         }
 
// 	for (field = NumberOfBaryonFields; field < MAX_NUMBER_OF_BARYON_FIELDS;
// 	     field++) {
//           SubgridFluxes[i]->LeftFluxes[field][dim]  = NULL;
//           SubgridFluxes[i]->RightFluxes[field][dim] = NULL;
// 	}
 
//       }  // next dimension
 
//       /* make things pretty */
 
//       for (dim = GridRank; dim < 3; dim++)
//         for (field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++) {
//           SubgridFluxes[i]->LeftFluxes[field][dim]  = NULL;
//           SubgridFluxes[i]->RightFluxes[field][dim] = NULL;
// 	}
 
//     } // end of loop over subgrids
 
    /* compute global start index for left edge of entire grid
       (including boundary zones) */
 
     for (dim = 0; dim < GridRank[in]; dim++)
       GridGlobalStart[dim] =
     	NINT((GridLeftEdge[dim] - DomainLeftEdge[in*3+dim])/CellWidth[dim]) -
     	GridStartIndex[dim];
 
    /* fix grid quantities so they are defined to at least 3 dims */
 
    for (i = GridRank[in]; i < 3; i++) {
      GridDimension[i]   = 1;
      GridStartIndex[i]  = 0;
      GridEndIndex[i]    = 0;
      GridGlobalStart[i] = 0;
    }
 
    /* allocate temporary space for solver */

    int k = 0;
    enzo_float *temp = new enzo_float[size*(35)];
    enzo_float *f1 = &temp[k*size];  k++;
    enzo_float *f2 = &temp[k*size];  k++;
    enzo_float *f3 = &temp[k*size];  k++;
    enzo_float *f4 = &temp[k*size];  k++;
    enzo_float *f5 = &temp[k*size];  k++;
    enzo_float *f6 = &temp[k*size];  k++;
    enzo_float *f7 = &temp[k*size];  k++;
    enzo_float *f8 = &temp[k*size];  k++;
    enzo_float *g1 = &temp[k*size];  k++;
    enzo_float *g2 = &temp[k*size];  k++;
    enzo_float *g3 = &temp[k*size];  k++;
    enzo_float *g4 = &temp[k*size];  k++;
    enzo_float *g5 = &temp[k*size];  k++;
    enzo_float *g6 = &temp[k*size];  k++;
    enzo_float *g7 = &temp[k*size];  k++;
    enzo_float *g8 = &temp[k*size];  k++;
    enzo_float *h1 = &temp[k*size];  k++;
    enzo_float *h2 = &temp[k*size];  k++;
    enzo_float *h3 = &temp[k*size];  k++;
    enzo_float *h4 = &temp[k*size];  k++;
    enzo_float *h5 = &temp[k*size];  k++;
    enzo_float *h6 = &temp[k*size];  k++;
    enzo_float *h7 = &temp[k*size];  k++;
    enzo_float *h8 = &temp[k*size];  k++;
    enzo_float *ex = &temp[k*size];  k++;
    enzo_float *ey = &temp[k*size];  k++;
    enzo_float *ez = &temp[k*size];  k++;
    enzo_float *qu1 = &temp[k*size];  k++;
    enzo_float *qu2 = &temp[k*size];  k++;
    enzo_float *qu3 = &temp[k*size];  k++;
    enzo_float *qu4 = &temp[k*size];  k++;
    enzo_float *qu5 = &temp[k*size];  k++;
    enzo_float *qu6 = &temp[k*size];  k++;
    enzo_float *qu7 = &temp[k*size];  k++;
    enzo_float *qu8 = &temp[k*size];  k++;

    ASSERT ("EnzoBlock::SolveMHDEquationsIG",
	    "Insufficient temporary storage",
	    k <= 35);
    /* create and fill in arrays which are easiler for the solver to
       understand. */

    int NumberOfSubgrids = 0; // JB

    // int *leftface  = new int[NumberOfSubgrids*3*22];
    // int *rightface = leftface + NumberOfSubgrids*3*1;
    // int *istart    = leftface + NumberOfSubgrids*3*2;
    // int *jstart    = leftface + NumberOfSubgrids*3*3;
    // int *iend      = leftface + NumberOfSubgrids*3*4;
    // int *jend      = leftface + NumberOfSubgrids*3*5;
    // int *dnindex   = leftface + NumberOfSubgrids*3*6;
    // int *vxindex   = leftface + NumberOfSubgrids*3*8;
    // int *vyindex   = leftface + NumberOfSubgrids*3*10;
    // int *vzindex   = leftface + NumberOfSubgrids*3*12;
    // int *bxindex   = leftface + NumberOfSubgrids*3*14;
    // int *byindex   = leftface + NumberOfSubgrids*3*16;
    // int *bzindex   = leftface + NumberOfSubgrids*3*18;
    // int *pindex    = leftface + NumberOfSubgrids*3*20;

    enzo_float *standard = NULL;

    int *leftface  =  NULL;
    int *rightface =  NULL;
    int *istart    =  NULL;
    int *jstart    =  NULL;
    int *iend      =  NULL;
    int *jend      =  NULL;
    int *dnindex   =  NULL;
    int *vxindex   =  NULL;
    int *vyindex   =  NULL;
    int *vzindex   =  NULL;
    int *bxindex   =  NULL;
    int *byindex   =  NULL;
    int *bzindex   =  NULL;
    int *prindex   =  NULL;

    // enzo_float *standard = NULL;
    //    if (NumberOfSubgrids > 0) standard = SubgridFluxes[0]->LeftFluxes[0][0];
 
    // for (subgrid = 0; subgrid < NumberOfSubgrids; subgrid++)
    //   for (dim = 0; dim < GridRank[in]; dim++) {
 
        /* Set i,j dimensions of 2d flux slice (this works even if we
           are in 1 or 2d) the correspond to the dimensions of the global
           indicies.  I.e. for dim = 0, the plane is dims 1,2
                           for dim = 1, the plane is dims 0,2
                           for dim = 2, the plane is dims 0,1 */
 
// 	idim = (dim == 0) ? 1 : 0;
// 	jdim = (dim == 2) ? 1 : 2;
 
        /* Set the index (along the dimension perpindicular to the flux
           plane) of the left and right flux planes.  The index is zero
           based from the left side of the entire grid. */
 
 	// leftface[subgrid*3+dim] =
 	//   SubgridFluxes[subgrid]->LeftFluxStartGlobalIndex[dim][dim] -
 	//     GridGlobalStart[dim];
 	// rightface[subgrid*3+dim] =
 	//   SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[dim][dim] -
 	//     GridGlobalStart[dim];   // (+1 done by fortran code)
 
        /* set the start and end indicies (zero based on entire grid)
           of the 2d flux plane. */
 
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
// 	   NOTE: This pointer arithmetic is illegal; some other way should
// 	   be found to do it (like write higher level ppm stuff in c++). */
 
// 	dnindex[subgrid*6+dim*2] =
// 	  SubgridFluxes[subgrid]->LeftFluxes[0][dim]  - standard;
// 	dnindex[subgrid*6+dim*2+1] =
// 	  SubgridFluxes[subgrid]->RightFluxes[0][dim] - standard;
// 	vxindex[subgrid*6+dim*2] =
// 	  SubgridFluxes[subgrid]->LeftFluxes[1][dim]  - standard;
// 	vxindex[subgrid*6+dim*2+1] =
// 	  SubgridFluxes[subgrid]->RightFluxes[1][dim] - standard;
// 	vyindex[subgrid*6+dim*2] =
// 	  SubgridFluxes[subgrid]->LeftFluxes[2][dim]  - standard;
// 	vyindex[subgrid*6+dim*2+1] =
// 	  SubgridFluxes[subgrid]->RightFluxes[2][dim] - standard;
// 	vzindex[subgrid*6+dim*2] =
// 	  SubgridFluxes[subgrid]->LeftFluxes[3][dim]  - standard;
// 	vzindex[subgrid*6+dim*2+1] =
// 	  SubgridFluxes[subgrid]->RightFluxes[3][dim] - standard;
// 	bxindex[subgrid*6+dim*2] =
// 	  SubgridFluxes[subgrid]->LeftFluxes[4][dim]  - standard;
// 	bxindex[subgrid*6+dim*2+1] =
// 	  SubgridFluxes[subgrid]->RightFluxes[4][dim] - standard;
// 	byindex[subgrid*6+dim*2] =
// 	  SubgridFluxes[subgrid]->LeftFluxes[5][dim]  - standard;
// 	byindex[subgrid*6+dim*2+1] =
// 	  SubgridFluxes[subgrid]->RightFluxes[5][dim] - standard;
// 	bzindex[subgrid*6+dim*2] =
// 	  SubgridFluxes[subgrid]->LeftFluxes[6][dim]  - standard;
// 	bzindex[subgrid*6+dim*2+1] =
// 	  SubgridFluxes[subgrid]->RightFluxes[6][dim] - standard;
 
      // }
 
    /* If using comoving coordinates, multiply dx by a(n+1/2).
       In one fell swoop, this recasts the equations solved by solver
       in comoving form (except for the expansion terms which are taken
       care of elsewhere). */
 
    /* Create a cell width array to pass (and convert to absolute coords). */
    // this is not going to work for cosmology right away !AK

    enzo_float a = 1.0;
    enzo_float CellWidthTemp[MAX_DIMENSION];
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      if (dim < GridRank[in])
	CellWidthTemp[dim] = enzo_float(a*CellWidth[dim]);
      else
	CellWidthTemp[dim] = 1.0;
    }
 
    /* Prepare Gravity. */
 
    /* call a Fortran routine to actually compute the hydro equations
       on this grid.
       Notice that it is hard-wired for three dimensions, but it does
       the right thing for < 3 dimensions. */
    /* note: Start/EndIndex are zero based */


    /* current PPML-IG implementation only supports 3D and does not
       support color fields */	

    int mx,my,mz;
    field.dimensions(0,&mx,&my,&mz);
    const int m = mx*my*mz;

    enzo_float *velocity_x      = (enzo_float *) field.values ("velocity_x");
    enzo_float *velocity_y      = (enzo_float *) field.values ("velocity_y");
    enzo_float *velocity_z      = (enzo_float *) field.values ("velocity_z");
    bool have_velocity = (velocity_x != nullptr);

    if (have_velocity) {
      std::copy_n(velocity_x,m,velox);
      std::copy_n(velocity_y,m,veloy);
      std::copy_n(velocity_z,m,veloz);
    }

    std::copy_n(pressure,m,press_rx);
    std::copy_n(pressure,m,press_ry);
    std::copy_n(pressure,m,press_rz);
    
    enzo_float b0[3] = { 0.0, 0.0, 0.0 };
      FORTRAN_NAME(ppml_ig)
	(density,velox,   veloy,   veloz,   bfieldx,   bfieldy,   bfieldz,   pressure,
	 dens_rx,velox_rx,veloy_rx,veloz_rx,bfieldx_rx,bfieldy_rx,bfieldz_rx,press_rx,
	 dens_ry,velox_ry,veloy_ry,veloz_ry,bfieldx_ry,bfieldy_ry,bfieldz_ry,press_ry,
	 dens_rz,velox_rz,veloy_rz,veloz_rz,bfieldx_rz,bfieldy_rz,bfieldz_rz,press_rz,
	 //	 gravity,gx,gy,gz,
	 b0, &gamma,
	 &dt, &CellWidthTemp[0], &CellWidthTemp[1], &CellWidthTemp[2],
	 &GridDimension[0], &GridDimension[1], &GridDimension[2], 	 
	 GridStartIndex, GridEndIndex,
	 &NumberOfSubgrids, leftface, rightface,
	 istart, iend, jstart, jend,
	 standard, dnindex,
	 vxindex, vyindex, vzindex,
	 bxindex, byindex, bzindex,
	 prindex,
	 f1,f2,f3,f4,f5,f6,f7,f8,
	 g1,g2,g3,g4,g5,g6,g7,g8,
	 h1,h2,h3,h4,h5,h6,h7,h8,
	 ex,ey,ez,
	 qu1,qu2,qu3,qu4,qu5,qu6,qu7,qu8);

    if (have_velocity) {
      std::copy_n(velox,m,velocity_x);
      std::copy_n(veloy,m,velocity_y);
      std::copy_n(veloz,m,velocity_z);
    }
    std::copy_n(press_rx,m,pressure);
    std::copy_n(press_rx,m,pressure);
    std::copy_n(press_rx,m,pressure);

    /* deallocate temporary space for solver */
 
    delete [] temp;
 
    delete [] leftface;
 
  }  // end: if (NumberOfBaryonFields > 0)
 
  return ENZO_SUCCESS;
 
}
