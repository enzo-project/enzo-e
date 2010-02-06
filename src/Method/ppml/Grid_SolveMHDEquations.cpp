/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @bug       
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */
/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @bug       
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */
/***********************************************************************
/
/  GRID CLASS (SOLVE THE MHD EQUATIONS (with PPML), SAVING SUBGRID FLUXES)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Alexei Kritsuk, April 2009 (PPML)
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
// Solve the MHD equations with the solver, saving the subgrid fluxes

#include "cello_hydro.h"
 
  
/* function prototypes */
 
void my_exit(int status);

/*
!        i1,j1,k1 - initial point of grid
!        i2,j2,k2 - last point of grid
!    ihost,jhost,khost - number of host zones    
!    dn,vx,vy,vz,bx,by,bz - central primitive variables
!    dnrx,vxrx,vyrx,vzrx,bxrx,byrx,bzrx - x-right primitive variables
!    dnry,vxry,vyry,vzry,bxry,byry,bzry - y-right primitive variables
!    dnrz,vxrz,vyrz,vzrz,bxrz,byrz,bzrz - z-right primitive variables 
!    f1,f2,f3,f4,f5,f6,f7 - x-fluxes
!        g1,g2,g3,g4,g5,g6,g7 - y-fluxes
!    h1,h2,h3,h4,h5,h6,h7 - z-fluxes 
*/

extern "C" void FORTRAN_NAME(ppml)(float *dn,   float *vx,   float *vy,   float *vz,
			           float *bx,   float *by,   float *bz,
			           float *dnrx, float *vxrx, float *vyrx, float *vzrx,
			           float *bxrx, float *byrx, float *bzrx,
			           float *dnry, float *vxry, float *vyry, float *vzry,
			           float *bxry, float *byry, float *bzry,
			           float *dnrz, float *vxrz, float *vyrz, float *vzrz,
			           float *bxrz, float *byrz, float *bzrz,
                                   float *dt, float dx[], float dy[], float dz[],
                                   int *in, int *jn, int *kn,
				   int is[], int ie[],
                                   int *num_subgrids, int leftface[], int rightface[],
                                   int istart[], int iend[], int jstart[], int jend[],
                                   float *standard, int dnindex[], 
                                   int vxindex[], int vyindex[], int vzindex[],
                                   int bxindex[], int byindex[], int bzindex[]);

/*extern "C" void FORTRAN_NAME(ppm_de)(
			  float *d, float *E, float *u, float *v, float *w,
			    float *ge,
                          int *grav, float *gr_ax, float *gr_ay, float *gr_az,
			  float *gamma, float *dt,
                            float dx[], float dy[], float dz[],
			  int *rank, int *in, int *jn, int *kn,
                            int is[], int ie[],
			  float gridvel[], int *flatten, int *ipresfree,
			  int *diff, int *steepen, int *idual,
                            float *eta1, float *eta2,
			  int *num_subgrids, int leftface[], int rightface[],
			  int istart[], int iend[], int jstart[], int jend[],
			  float *standard, int dindex[], int Eindex[],
			  int uindex[], int vindex[], int windex[],
			    int geindex[], float *temp,
                          int *ncolour, float *colourpt, int *coloff,
                            int colindex[]);
*/
int SolveMHDEquations(int cycle, float dt)
{
 
  /* exit if not 3D */

  // @@ assert GridRank == 3
  //  if (GridRank != 3) 
  //    my_exit(EXIT_FAILURE);

  if (NumberOfBaryonFields > 0) {
 
    /* initialize */
 
    int dim, i, idim,  jdim,  size, subgrid;
    Elong_int GridGlobalStart[MAX_DIMENSION];
 
    /* Compute size (in floats) of the current grid. */
 
    size = 1;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];
 
    /* Get easy to handle pointers for each variable. */
 

    float *density    = BaryonField[ 0];
    float *velox      = BaryonField[ 1];
    float *veloy      = BaryonField[ 2];
    float *veloz      = BaryonField[ 3];
    float *bfieldx    = BaryonField[ 4];
    float *bfieldy    = BaryonField[ 5];
    float *bfieldz    = BaryonField[ 6];

    float *dens_rx    = BaryonField[ 7];
    float *velox_rx   = BaryonField[ 8];
    float *veloy_rx   = BaryonField[ 9];
    float *veloz_rx   = BaryonField[10];
    float *bfieldx_rx = BaryonField[11];
    float *bfieldy_rx = BaryonField[12];
    float *bfieldz_rx = BaryonField[13];

    float *dens_ry    = BaryonField[14];
    float *velox_ry   = BaryonField[15];
    float *veloy_ry   = BaryonField[16];
    float *veloz_ry   = BaryonField[17];
    float *bfieldx_ry = BaryonField[18];
    float *bfieldy_ry = BaryonField[19];
    float *bfieldz_ry = BaryonField[20];

    float *dens_rz    = BaryonField[21];
    float *velox_rz   = BaryonField[22];
    float *veloy_rz   = BaryonField[23];
    float *veloz_rz   = BaryonField[24];
    float *bfieldx_rz = BaryonField[25];
    float *bfieldy_rz = BaryonField[26];
    float *bfieldz_rz = BaryonField[27];

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
// 	    SubgridFluxes[i]->LeftFluxes[field][dim]  = new float[size];
// 	  if (SubgridFluxes[i]->RightFluxes[field][dim] == NULL)
// 	    SubgridFluxes[i]->RightFluxes[field][dim] = new float[size];
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
    // no colours here for now !AK 

    int tempsize = max(max(GridDimension[0]*GridDimension[1],
                           GridDimension[1]*GridDimension[2]),
		           GridDimension[2]*GridDimension[0]  );
    float *temp = new float[tempsize*(31)];
 
    /* create and fill in arrays which are easiler for the solver to
       understand. */

    int NumberOfSubgrids = 0; // JB

    int *leftface  = new int[NumberOfSubgrids*3*20];
    int *rightface = leftface + NumberOfSubgrids*3*1;
    int *istart    = leftface + NumberOfSubgrids*3*2;
    int *jstart    = leftface + NumberOfSubgrids*3*3;
    int *iend      = leftface + NumberOfSubgrids*3*4;
    int *jend      = leftface + NumberOfSubgrids*3*5;
    int *dnindex   = leftface + NumberOfSubgrids*3*6;
    int *vxindex   = leftface + NumberOfSubgrids*3*8;
    int *vyindex   = leftface + NumberOfSubgrids*3*10;
    int *vzindex   = leftface + NumberOfSubgrids*3*12;
    int *bxindex   = leftface + NumberOfSubgrids*3*14;
    int *byindex   = leftface + NumberOfSubgrids*3*16;
    int *bzindex   = leftface + NumberOfSubgrids*3*18;

    float *standard = NULL;
    //    if (NumberOfSubgrids > 0) standard = SubgridFluxes[0]->LeftFluxes[0][0];
 
    for (subgrid = 0; subgrid < NumberOfSubgrids; subgrid++)
      for (dim = 0; dim < GridRank; dim++) {
 
        /* Set i,j dimensions of 2d flux slice (this works even if we
           are in 1 or 2d) the correspond to the dimensions of the global
           indicies.  I.e. for dim = 0, the plane is dims 1,2
                           for dim = 1, the plane is dims 0,2
                           for dim = 2, the plane is dims 0,1 */
 
	idim = (dim == 0) ? 1 : 0;
	jdim = (dim == 2) ? 1 : 2;
 
        /* Set the index (along the dimension perpindicular to the flux
           plane) of the left and right flux planes.  The index is zero
           based from the left side of the entire grid. */
 
// 	leftface[subgrid*3+dim] =
// 	  SubgridFluxes[subgrid]->LeftFluxStartGlobalIndex[dim][dim] -
// 	    GridGlobalStart[dim];
// 	rightface[subgrid*3+dim] =
// 	  SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[dim][dim] -
// 	    GridGlobalStart[dim];   // (+1 done by fortran code)
 
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
 
      }
 
    /* If using comoving coordinates, multiply dx by a(n+1/2).
       In one fell swoop, this recasts the equations solved by solver
       in comoving form (except for the expansion terms which are taken
       care of elsewhere). */
 
    /* Create a cell width array to pass (and convert to absolute coords). */
    // this is not going to work for cosmology right away !AK

    float a = 1.0;
    float *CellWidthTemp[MAX_DIMENSION];
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      CellWidthTemp[dim] = new float[GridDimension[dim]];
      for (i = 0; i < GridDimension[dim]; i++)
	if (dim < GridRank)
	  CellWidthTemp[dim][i] = float(a*CellWidth[dim][i]);
	else
	  CellWidthTemp[dim][i] = 1.0;
    }
 
    /* Prepare Gravity. */
 
//     int GravityOn = 0;
//     if (SelfGravity || UniformGravity || PointSourceGravity)
//       GravityOn = 1;
 
    /* call a Fortran routine to actually compute the hydro equations
       on this grid.
       Notice that it is hard-wired for three dimensions, but it does
       the right thing for < 3 dimensions. */
    /* note: Start/EndIndex are zero based */


    /* current PPML implementation only supports 3D and does not support color fields */	

      FORTRAN_NAME(ppml)
	(density,velox,   veloy,   veloz,   bfieldx,   bfieldy,   bfieldz,
	 dens_rx,velox_rx,veloy_rx,veloz_rx,bfieldx_rx,bfieldy_rx,bfieldz_rx,
	 dens_ry,velox_ry,veloy_ry,veloz_ry,bfieldx_ry,bfieldy_ry,bfieldz_ry,
	 dens_rz,velox_rz,veloy_rz,veloz_rz,bfieldx_rz,bfieldy_rz,bfieldz_rz,
	 &dt, CellWidthTemp[0], CellWidthTemp[1], CellWidthTemp[2],
	 &GridDimension[0], &GridDimension[1], &GridDimension[2], 
	 GridStartIndex, GridEndIndex,
	 &NumberOfSubgrids, leftface, rightface,
	 istart, iend, jstart, jend,
	 standard, dnindex,
	 vxindex, vyindex, vzindex,
	 bxindex, byindex, bzindex);

    /* deallocate temporary space for solver */
 
    delete [] temp;
 
    delete [] leftface;
 
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      delete [] CellWidthTemp[dim];
 
  }  // end: if (NumberOfBaryonFields > 0)
 
  return SUCCESS;
 
}
