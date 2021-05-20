// See LICENSE_ENZO file for license and copyright information

/// @file      enzo_SolveHydroEquations.cpp
/// @author    Greg Bryan
/// @date      November, 1994
/// @brief     Solve the hydro equations, saving subgrid fluxes

#include "cello.hpp"
#include "enzo.hpp"
#include <stdio.h>

// #define IE_ERROR_FIELD
// #define DEBUG_PPM
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
  int gx,gy,gz;
  int mx,my,mz;
  field.ghost_depth(0,&gx,&gy,&gz);
  field.dimensions(0,&mx,&my,&mz);

#ifdef IE_ERROR_FIELD
  int num_ie_error = 0;
  int *ie_error_x = new int [mx*my*mz];
  int *ie_error_y = new int [mx*my*mz];
  int *ie_error_z = new int [mx*my*mz];
#else
  int num_ie_error = -1;
  int *ie_error_x = nullptr;
  int *ie_error_y = nullptr;
  int *ie_error_z = nullptr;
#endif


  //------------------------------
  // Prepare color field parameters
  //------------------------------

  // ncolor: number of color fields

  int    ncolor  = field.groups()->size("color");

  // colorpt: the color 'array' (contains all color fields)
  enzo_float * colorpt = (enzo_float *) field.permanent();

  // coloff: offsets into the color array (for each color field)
  int * coloff   = (ncolor > 0) ? (new int [ncolor]) : NULL;
  int index_color = 0;
  for (int index_field = 0;
       index_field < field.field_count();
       index_field++) {
    std::string name = field.field_name(index_field);
    if (field.groups()->is_in(name,"color")) {
      coloff[index_color++]
	= (enzo_float *)(field.values(index_field)) - colorpt;
    }

  }

  // No subgrids, so colindex is NULL
  int *colindex         = NULL;

  /* Compute size (in enzo_floats) of the current grid. */

  int rank = cello::rank();

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

  int NumberOfSubgrids = 1;

  /* fix grid quantities so they are defined to at least 3 dims */

  for (int i = rank; i < 3; i++) {
    GridDimension[i]   = 1;
    GridStartIndex[i]  = 0;
    GridEndIndex[i]    = 0;
    //      GridGlobalStart[i] = 0;
  }

  /* allocate temporary space for solver (enough to fit 31 of the largest
     possible 2d slices plus 4*ncolor). */

  int tempsize = MAX(MAX(GridDimension[0]*GridDimension[1],
			 GridDimension[1]*GridDimension[2]),
		     GridDimension[2]*GridDimension[0]);

  enzo_float *temp = new enzo_float[tempsize*(32+ncolor*4)];

  /* create and fill in arrays which are easier for the solver to
     understand. */

  size = NumberOfSubgrids*3*(18+2*ncolor) + 1;
  int * array = new int[size];
  for (int i=0; i<size; i++) array[i] = 0;

  int * p = array;
  int *leftface  = p; p+=3;
  int *rightface = p; p+=3;
  int *istart    = p; p+=3;
  int *jstart    = p; p+=3;
  int *iend      = p; p+=3;
  int *jend      = p; p+=3;
  int *dindex    = p; p+=3*2; // 3 axes 2 faces
  int *Eindex    = p; p+=3*2;
  int *uindex    = p; p+=3*2;
  int *vindex    = p; p+=3*2;
  int *windex    = p; p+=3*2;
  int *geindex   = p; p+=3*2;

  // Offsets computed from the "standard" pointer to the start of each
  // flux data

  FluxData * flux_data = data()->flux_data();

  const int nx = mx - 2*gx;
  const int ny = my - 2*gy;
  const int nz = mz - 2*gz;

  enzo_float * standard = temp;

  // int l3[3] = {gx,gy,gz};
  // int u3[3] = {mx-gx,my-gy,mz-gz};
  int l3[3] = {gx,gy,gz};
  int u3[3] = {mx-gx-1,my-gy-1,mz-gz-1};
  const int nf = flux_data->num_fields();
  for (int i_f=0; i_f <nf; i_f++) {
    int * flux_index = 0;
    const int index_field = flux_data->index_field(i_f);
    const std::string field_name = field.field_name(index_field);

    if (field_name == "density")         flux_index = dindex;
    if (field_name == "velocity_x")      flux_index = uindex;
    if (field_name == "velocity_y")      flux_index = vindex;
    if (field_name == "velocity_z")      flux_index = windex;
    if (field_name == "total_energy")    flux_index = Eindex;
    if (field_name == "internal_energy") flux_index = geindex;

    for (int axis=0; axis<rank; axis++) {
      leftface[axis] = l3[axis];
      rightface[axis] = u3[axis];

      const int axis_i = (axis == 0) ? 1 : 0;
      const int axis_j = (axis == 2) ? 1 : 2;
      istart[axis] = l3[axis_i];
      jstart[axis] = l3[axis_j];
      iend[axis] = u3[axis_i];
      jend[axis] = u3[axis_j];
      for (int face=0; face<2; face++) {
        FaceFluxes * ff_b = flux_data->block_fluxes(axis,face,i_f);
        int dx,dy,dz;
        flux_index[axis*2+face] =
          ((enzo_float *)ff_b->flux_array(&dx,&dy,&dz).data()) - standard;
      }
    }
  }

  //==================================================
  //    enzo_float *standard = SubgridFluxes[0]->LeftFluxes[0][0];

  /* If using comoving coordinates, multiply dx by a(n+1/2).
     In one fell swoop, this recasts the equations solved by solver
     in comoving form (except for the expansion terms which are taken
     care of elsewhere). */

  enzo_float cosmo_a = 1.0, cosmo_dadt = 0.0;

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  if (comoving_coordinates) {
    cosmology->compute_expansion_factor(&cosmo_a, &cosmo_dadt, time+0.5*dt);
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
	CellWidthTemp[dim][i] = (cosmo_a*CellWidth[dim]);
    } else {
      for (int i=0; i<GridDimension[dim]; i++)
	CellWidthTemp[dim][i] = 1.0;
    }
  }

  /* call a Fortran routine to actually compute the hydro equations
     on this grid.
     Notice that it is hard-wired for three dimensions, but it does
     the right thing for < 3 dimensions. */
  /* note: Start/EndIndex are zero based */

  int gravity_on = (acceleration_x != NULL) ? 1 : 0;

  int iconsrec = 0;
  int iposrec = 0;

  int error = 0;

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
     &ncolor, colorpt, coloff, colindex,
     &error, ie_error_x,ie_error_y,ie_error_z,&num_ie_error
     );

  if (error != 0) {
    char buffer[256];
    snprintf (buffer,255,"Error %d in call to ppm_de block %s",error,name().c_str());
  }

#ifdef IE_ERROR_FIELD
  if (num_ie_error > 0) {
    CkPrintf ("DEBUG_IE_ERROR num_ie_error = %d\n",num_ie_error);

    enzo_float * error_ie  = (enzo_float*)field.values("internal_energy_error");
    for (int i=0; i<mx*my*mz; i++) {
      error_ie[i] = 0.0;
    }
    for (int k=0; k<num_ie_error; k++) {
      int i = ie_error_x[k] + mx*(ie_error_y[k] + my*ie_error_z[k]);
      CkPrintf ("DEBUG_IE_ERROR setting error_ie[%d (%d %d %d)]\n",
                i,ie_error_x[k],ie_error_y[k],ie_error_z[k]);
      error_ie[i] = 1.0;
    }
  }
#endif

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    delete [] CellWidthTemp[dim];
  }

  /* deallocate temporary space for solver */

  delete [] temp;
  if (rank < 2) delete [] velocity_y;
  if (rank < 3) delete [] velocity_z;

  delete [] array;

  delete [] coloff;
#ifdef IE_ERROR_FIELD
  delete [] ie_error_x;
  delete [] ie_error_y;
  delete [] ie_error_z;
#endif

  return ENZO_SUCCESS;

}
