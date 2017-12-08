// See LICENSE_ENZO file for license and copyright information

/// @file      enzo_SolveHydroEquations.cpp
/// @author    Greg Bryan
/// @date      November, 1994
/// @brief     Solve the hydro equations, saving subgrid fluxes

#include "cello.hpp"
#include "enzo.hpp"
#include <stdio.h>
// #define DEBUG_TRACE_PPM
// #define DEBUG_READ_FIELDS
// #define DEBUG_WRITE_FIELDS

#ifdef DEBUG_TRACE_PPM
#  define TRACE_PPM(MESSAGE)						\
  CkPrintf ("%s:%d %s %s\n",						\
	    __FILE__,__LINE__,block->name().c_str(),MESSAGE);				
#else
#  define TRACE_PPM(MESSAGE) /* ... */
#endif

#ifdef DEBUG_READ_FIELDS
#   define READ_FIELD(FIELD,NAME,CYCLE,field,gx,gy,gz,mx,my,mz)	\
  {								\
    if (CYCLE==0) {						\
      enzo_float * array = (enzo_float*) field.values(FIELD);	\
      char buffer[80];						\
      sprintf (buffer,NAME,CYCLE);				\
      printf ("READ_FIELD %s cycle=%d\n",FIELD,CYCLE);		\
      FILE * fp = fopen(buffer,"r");				\
      for (int iz=gz; iz<mz-gz; iz++) {				\
	for (int iy=gy; iy<my-gy; iy++) {			\
	  for (int ix=gx; ix<mx-gx; ix++) {			\
	    int i=ix+mx*(iy+my*iz);				\
	    int jx,jy,jz;					\
	    double value;					\
	    fscanf (fp,"%d %d %d %lf\n",&jx,&jy,&jz,&value);	\
	    array[i] = value;					\
	  }							\
	}							\
      }								\
      fclose(fp);						\
      fp=NULL;							\
    }								\
  }
#  define READ_STATE(NAME,CYCLE,TIME,DT,DX)		\
  {							\
    char buffer[80];					\
    sprintf (buffer,NAME,CYCLE);			\
    FILE * fp = fopen(buffer,"r");			\
    printf ("fscanf returned %d\n",fscanf (fp,"%f %f %f",&TIME,&DT,&DX)); \
    printf ("DEBUG_TIME %20.15f %20.15f %20.15f\n",TIME,DT,DX);		\
    fclose(fp);								\
}
#else
#   define READ_FIELD(FIELD,NAME,CYCLE,field,gx,gy,gz,mx,my,mz)	\
  /* EMPTY */
#  define READ_STATE(NAME,CYCLE,TIME,DT,DX)	\
  /* EMPTY */
#endif

#ifdef DEBUG_WRITE_FIELDS
#   define WRITE_FIELD(FIELD,NAME,CYCLE,field,gx,gy,gz,mx,my,mz)	\
  {								\
    enzo_float * array = (enzo_float*) field.values(FIELD);	\
    char buffer[80];						\
    sprintf (buffer,NAME,CYCLE);				\
    printf ("WRITE_FIELD %s cycle=%d\n",FIELD,CYCLE);		\
    FILE * fp = fopen(buffer,"w");				\
    for (int iz=gz; iz<mz-gz; iz++) {				\
      for (int iy=gy; iy<my-gy; iy++) {				\
  	for (int ix=gx; ix<mx-gx; ix++) {			\
  	  int i=ix+mx*(iy+my*iz);				\
  	  fprintf (fp,"%d %d %d %20.16g\n",ix-gx,iy-gy,iz-gz,array[i]);	\
  	}							\
      }								\
    }								\
    fclose(fp);							\
    fp=NULL;							\
  }
#else
#   define WRITE_FIELD(FIELD,FILE,CYCLE,field,gx,gy,gz,mx,my,mz)	\
  /* EMPTY */
#endif

//----------------------------------------------------------------------

int EnzoBlock::SolveHydroEquations 
(
 enzo_float time,
 enzo_float dt,
 bool comoving_coordinates
 )
{
  /* initialize */

  int dim, size;

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

  int rank = this->rank();

  size = 1;
  for (dim = 0; dim < rank; dim++)
    size *= GridDimension[dim];

  enzo_float * density         = (enzo_float *) field.values("density");
  enzo_float * total_energy    = (enzo_float *) field.values("total_energy");
  enzo_float * internal_energy = (enzo_float *) field.values("internal_energy");

  /* velocity_x must exist, but if y & z aren't present, then create blank
     buffers for them (since the solver needs to advect something). */

  enzo_float * velocity_x = NULL;
  enzo_float * velocity_y = NULL;
  enzo_float * velocity_z = NULL;

  velocity_x = (enzo_float *) field.values("velocity_x");

  if (rank >= 2) {
    velocity_y = (enzo_float *) field.values("velocity_y");
  } else {
    velocity_y = new enzo_float[size];
    for (int i=0; i<size; i++) velocity_y[i] = 0.0;
  }

    if (rank >= 3) {
    velocity_z = (enzo_float *) field.values("velocity_z");
  } else {
    velocity_z = new enzo_float[size];
    for (int i=0; i<size; i++) velocity_z[i] = 0.0;
  }

  enzo_float * acceleration_x  = field.is_field("acceleration_x") ? 
    (enzo_float *) field.values("acceleration_x") : NULL;
  enzo_float * acceleration_y  = field.is_field("acceleration_y") ? 
    (enzo_float *)field.values("acceleration_y") : NULL;
  enzo_float * acceleration_z  = field.is_field("acceleration_z") ? 
    (enzo_float *)field.values("acceleration_z") : NULL;


  /* Determine if Gamma should be a scalar or a field. */

  const int in = cello::index_static();

  /* Set minimum support. */

  enzo_float MinimumSupportEnergyCoefficient = 0;
  if (UseMinimumPressureSupport[in] == TRUE) {
    if (SetMinimumSupport(MinimumSupportEnergyCoefficient,
			  comoving_coordinates) == ENZO_FAIL) {
      ERROR("EnzoBlock::SolveHydroEquations()",
	    "Grid::SetMinimumSupport() returned ENZO_FAIL");
    }
  }
  /* allocate space for fluxes */

  int NumberOfSubgrids = 0;

  //  SubgridFluxes = new fluxes *[NumberOfSubgrids];
  SubgridFluxes = NULL;

  // for (i = 0; i < NumberOfSubgrids; i++) {
  //   SubgridFluxes[i] = new fluxes;
      
  //   for (dim = 0; dim < rank; dim++)  {

  //     /* compute size (in enzo_floats) of flux storage */

  //     size = 1;
  //     for (j = 0; j < rank; j++)
  // 	size *= SubgridFluxes[i]->LeftFluxEndGlobalIndex[dim][j] -
  // 	  SubgridFluxes[i]->LeftFluxStartGlobalIndex[dim][j] + 1;

  //     /* set unused dims (for the solver, which is hardwired for 3d). */

  //     for (j = rank; j < 3; j++) {
  // 	SubgridFluxes[i]->LeftFluxStartGlobalIndex[dim][j] = 0;
  // 	SubgridFluxes[i]->LeftFluxEndGlobalIndex[dim][j] = 0;
  // 	SubgridFluxes[i]->RightFluxStartGlobalIndex[dim][j] = 0;
  // 	SubgridFluxes[i]->RightFluxEndGlobalIndex[dim][j] = 0;
  //     }

  //     /* Allocate space (if necessary). */

  //     for (int field = 0; field < NumberOfBaryonFields[in]; field++) {
  // 	if (SubgridFluxes[i]->LeftFluxes[field][dim] == NULL)
  // 	  SubgridFluxes[i]->LeftFluxes[field][dim]  = new enzo_float[size];
  // 	if (SubgridFluxes[i]->RightFluxes[field][dim] == NULL)
  // 	  SubgridFluxes[i]->RightFluxes[field][dim] = new enzo_float[size];
  // 	for (int n = 0; n < size; n++) {
  // 	  SubgridFluxes[i]->LeftFluxes[field][dim][n] = 0;
  // 	  SubgridFluxes[i]->RightFluxes[field][dim][n] = 0;
  // 	}
  //     }

  //     for (int field = NumberOfBaryonFields[in]; 
  // 	   field < MAX_NUMBER_OF_BARYON_FIELDS;
  // 	   field++) {
  // 	SubgridFluxes[i]->LeftFluxes[field][dim] = NULL;
  // 	SubgridFluxes[i]->RightFluxes[field][dim] = NULL;
  //     }

  //   }  // next dimension

  //   /* make things pretty */

  //   for (dim = rank; dim < 3; dim++)
  //     for (int field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++) {
  // 	SubgridFluxes[i]->LeftFluxes[field][dim] = NULL;
  // 	SubgridFluxes[i]->RightFluxes[field][dim] = NULL;
  //     }

  // } // end of loop over subgrids

  /* fix grid quantities so they are defined to at least 3 dims */

  for (int i = rank; i < 3; i++) {
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

  enzo_float *temp = new enzo_float[tempsize*(32+ncolour*4)];

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

  enzo_float a_cosmo = 1, dadt_cosmo;

  EnzoPhysicsCosmology * cosmology = (EnzoPhysicsCosmology * )
    simulation()->problem()->physics("cosmology");

  if (comoving_coordinates) {
    cosmology->compute_expansion_factor(&a_cosmo, &dadt_cosmo, time+0.5*dt);
  }

  ASSERT ("EnzoBlock::SolveHydroEquations()",
	  "comoving_coordinates enabled but missing EnzoPhysicsCosmology",
	  ! (comoving_coordinates && (cosmology == NULL)) );

  /* Create a cell width array to pass (and convert to absolute coords). */

  enzo_float * CellWidthTemp[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    CellWidthTemp[dim] = new enzo_float [GridDimension[dim]];
    if (dim < rank) {
      for (int i=0; i<GridDimension[dim]; i++) 
	CellWidthTemp[dim][i] = enzo_float(a_cosmo*CellWidth[dim]);
    } else {
      for (int i=0; i<GridDimension[dim]; i++) 
	CellWidthTemp[dim][i] = 1.0;
    }
  }

#ifdef DEBUG_READ_FIELDS
  if (cycle_ == 0) {
    enzo_float dx;
    READ_STATE("state-enzo-0-%03d.data",cycle_,time,dt,dx);
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      if (dim < rank) {
	for (int i=0; i<GridDimension[dim]; i++) 
	  CellWidthTemp[dim][i] = dx;
      }
    }
  }
#endif  

  /* call a Fortran routine to actually compute the hydro equations
     on this grid.
     Notice that it is hard-wired for three dimensions, but it does
     the right thing for < 3 dimensions. */
  /* note: Start/EndIndex are zero based */

  int gravity_on = (acceleration_x != NULL) ? 1 : 0;

#ifdef DEBUG_PPM  
  int mx,my,mz;
  int gx=0,gy=0,gz=0;
  field.dimensions(0,&mx,&my,&mz);

  READ_FIELD("density","de-enzo-0-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  READ_FIELD("velocity_x","vx-enzo-0-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  READ_FIELD("velocity_y","vy-enzo-0-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  READ_FIELD("velocity_z","vz-enzo-0-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  READ_FIELD("acceleration_x","ax-enzo-0-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  READ_FIELD("acceleration_y","ay-enzo-0-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  READ_FIELD("acceleration_z","az-enzo-0-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  READ_FIELD("total_energy","te-enzo-0-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  READ_FIELD("internal_energy","ie-enzo-0-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  
  WRITE_FIELD("density","de-enzop-0-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  WRITE_FIELD("velocity_x","vx-enzop-0-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  WRITE_FIELD("velocity_y","vy-enzop-0-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  WRITE_FIELD("velocity_z","vz-enzop-0-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  WRITE_FIELD("acceleration_x","ax-enzop-0-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  WRITE_FIELD("acceleration_y","ay-enzop-0-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  WRITE_FIELD("acceleration_z","az-enzop-0-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  WRITE_FIELD("total_energy","te-enzop-0-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  WRITE_FIELD("internal_energy","ie-enzop-0-%03d.data",cycle_,field,0,0,0,mx,my,mz);

  enzo_float * potential = (enzo_float *) field.values("potential_temp");
  enzo_float * density_total = (enzo_float *) field.values("B_temp");
  gx=gy=gz=3;
  TRACE_FIELD("ppm-0-density",density,1.0);
  TRACE_FIELD("ppm-0-potential",potential,1.0);
  TRACE_FIELD("ppm-0-density-total",density_total,1.0);
  TRACE_FIELD("ppm-0-total_energy",total_energy,1.0);
  TRACE_FIELD("ppm-0-velocity_x",velocity_x,1.0);
  TRACE_FIELD("ppm-0-velocity_y",velocity_y,1.0);
  TRACE_FIELD("ppm-0-velocity_z",velocity_z,1.0);
  PRINT_FIELD("ppm-0-acceleration_x",acceleration_x,1.0);
  TRACE_FIELD("ppm-0-acceleration_y",acceleration_y,1.0);
  TRACE_FIELD("ppm-0-acceleration_z",acceleration_z,1.0);
  TRACE_FIELD("ppm-0-internal_energy",internal_energy,1.0);


  Particle particle = data()->particle();
  TRACE_PARTICLE("ppm-0-velocity_x",particle,"dark","vx");
  TRACE_PARTICLE("ppm-0-velocity_y",particle,"dark","vy");
  TRACE_PARTICLE("ppm-0-velocity_z",particle,"dark","vz");
  TRACE_PARTICLE("ppm-0-acceleration_x",particle,"dark","ax");
  TRACE_PARTICLE("ppm-0-acceleration_y",particle,"dark","ay");
  TRACE_PARTICLE("ppm-0-acceleration_z",particle,"dark","az");
#endif
  
  int iconsrec = 0;
  int iposrec = 0;

  FORTRAN_NAME(ppm_de)
    (
     density, total_energy, velocity_x, velocity_y, velocity_z,
     internal_energy,
     &gravity_on, 
     acceleration_x,
     acceleration_y,
     acceleration_z,
     &Gamma[in], &dt, &cycle_,
     CellWidthTemp[0], CellWidthTemp[1], CellWidthTemp[2],
     &rank, &GridDimension[0], &GridDimension[1],
     &GridDimension[2], GridStartIndex, GridEndIndex,
     &PPMFlatteningParameter[in],
     &PressureFree[in],
     &iconsrec, &iposrec,
     &PPMDiffusionParameter[in], &PPMSteepeningParameter[in],
     &DualEnergyFormalism[in], &DualEnergyFormalismEta1[in],
     &DualEnergyFormalismEta2[in],
     &NumberOfSubgrids, leftface, rightface,
     istart, iend, jstart, jend,
     standard, dindex, Eindex, uindex, vindex, windex,
     geindex, temp,
     &ncolour, colourpt, coloff, colindex
     );

#ifdef DEBUG_READ_FIELDS
  READ_FIELD("density_diff","de-enzo-1-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  READ_FIELD("velocity_x_diff","vx-enzo-1-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  READ_FIELD("velocity_y_diff","vy-enzo-1-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  READ_FIELD("velocity_z_diff","vz-enzo-1-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  READ_FIELD("acceleration_x_diff","ax-enzo-1-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  READ_FIELD("acceleration_y_diff","ay-enzo-1-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  READ_FIELD("acceleration_z_diff","az-enzo-1-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  READ_FIELD("total_energy_diff","te-enzo-1-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  READ_FIELD("internal_energy_diff","ie-enzo-1-%03d.data",cycle_,field,0,0,0,mx,my,mz);

  enzo_float * density_diff    = (enzo_float *) field.values("density_diff");
  enzo_float * total_energy_diff    = (enzo_float *) field.values("total_energy_diff");
  enzo_float * internal_energy_diff = (enzo_float *) field.values("internal_energy_diff");
  enzo_float * velocity_x_diff = (enzo_float *) field.values("velocity_x_diff");
  enzo_float * velocity_y_diff = (enzo_float *) field.values("velocity_y_diff");
  enzo_float * velocity_z_diff = (enzo_float *) field.values("velocity_z_diff");
  enzo_float * acceleration_x_diff = (enzo_float *) field.values("acceleration_x_diff");
  enzo_float * acceleration_y_diff = (enzo_float *) field.values("acceleration_y_diff");
  enzo_float * acceleration_z_diff = (enzo_float *) field.values("acceleration_z_diff");

  for (int iz=0; iz<mz; iz++) {
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {
	const int i=ix+mx*(iy+my*iz);
	density_diff[i] -= density[i];
	total_energy_diff[i] -= total_energy[i];
	internal_energy_diff[i] -= internal_energy[i];
	velocity_x_diff[i] -= velocity_x[i];
	velocity_y_diff[i] -= velocity_y[i];
	velocity_z_diff[i] -= velocity_z[i];
	acceleration_x_diff[i] -= acceleration_x[i];
	acceleration_y_diff[i] -= acceleration_y[i];
	acceleration_z_diff[i] -= acceleration_z[i];
      }
    }
  }
#endif
  
  WRITE_FIELD("density","de-enzop-1-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  WRITE_FIELD("velocity_x","vx-enzop-1-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  WRITE_FIELD("velocity_y","vy-enzop-1-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  WRITE_FIELD("velocity_z","vz-enzop-1-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  WRITE_FIELD("acceleration_x","ax-enzop-1-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  WRITE_FIELD("acceleration_y","ay-enzop-1-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  WRITE_FIELD("acceleration_z","az-enzop-1-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  WRITE_FIELD("total_energy","te-enzop-1-%03d.data",cycle_,field,0,0,0,mx,my,mz);
  WRITE_FIELD("internal_energy","ie-enzop-1-%03d.data",cycle_,field,0,0,0,mx,my,mz);

  TRACE_FIELD("ppm-1-density",density,1.0);
  TRACE_FIELD("ppm-1-total_energy",total_energy,1.0);
  TRACE_FIELD("ppm-1-velocity_x",velocity_x,1.0);
  TRACE_FIELD("ppm-1-velocity_y",velocity_y,1.0);
  TRACE_FIELD("ppm-1-velocity_z",velocity_z,1.0);
  TRACE_FIELD("ppm-1-acceleration_x",acceleration_x,1.0);
  TRACE_FIELD("ppm-1-acceleration_y",acceleration_y,1.0);
  TRACE_FIELD("ppm-1-acceleration_z",acceleration_z,1.0);
  TRACE_FIELD("ppm-1-internal_energy",internal_energy,1.0);
  
  /* deallocate temporary space for solver */

  delete [] temp;
  if (rank < 2) delete [] velocity_y;
  if (rank < 3) delete [] velocity_z;

  delete [] array;
  
  for (int i=0; i<NumberOfSubgrids; i++) {
    delete SubgridFluxes[i];
  }
  delete [] SubgridFluxes;
  if (ncolour > 0) delete [] coloff;

  return ENZO_SUCCESS;

}
