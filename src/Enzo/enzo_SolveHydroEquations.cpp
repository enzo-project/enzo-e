// See LICENSE_ENZO file for license and copyright information

/// @file      enzo_SolveHydroEquations.cpp
/// @author    Greg Bryan
/// @date      November, 1994
/// @brief     Solve the hydro equations, saving subgrid fluxes

#include "cello.hpp"
#include "enzo.hpp"
#include <stdio.h>

#define NEW_FLUX
// #define DEBUG_NEW_FLUX
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

  //------------------------------
  // Prepare colour field parameters
  //------------------------------

  // ncolour: number of colour fields

  int    ncolour  = field.groups()->size("colour");

  // colourpt: the color 'array' (contains all color fields)
  enzo_float * colourpt = (enzo_float *) field.permanent();

  // coloff: offsets into the color array (for each color field)
  int * coloff   = (ncolour > 0) ? (new int [ncolour]) : NULL;
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

#ifdef NEW_FLUX
  //==================================================
  
  FluxData * flux_data = data()->flux_data();

  const int nx = mx - 2*gx;
  const int ny = my - 2*gy;
  const int nz = mz - 2*gz;
  auto field_names = field.groups()->group_list("conserved");
  const int nf = field_names.size();
  std::vector<int> field_list;
  field_list.resize(nf);
  for (int i=0; i<nf; i++) {
    field_list[i] = field.field_id(field_names[i]);
  }
  flux_data->allocate (nx,ny,nz,this->level(),dt,field_list);
  
  std::vector<double> flux_array[3][2][nf];

  for (int i_f=0; i_f <nf; i_f++) {
    const int index_field = field_list[i_f];
    for (int axis=0; axis<rank; axis++) {
      for (int face=0; face<2; face++) {
        int mx,my,mz;
        int dx,dy,dz;
        FaceFluxes * ff = flux_data->block_fluxes(axis,face,index_field);
        ff->get_dimensions(&mx,&my,&mz);
        flux_array[axis][face][i_f] = ff->flux_array(&dx,&dy,&dz);
      }
    }
  }

  //==================================================

  enzo_float * standard = temp;

  // int l3[3] = {gx,gy,gz};
  // int u3[3] = {mx-gx,my-gy,mz-gz};
  int l3[3] = {gx,gy,gz};
  int u3[3] = {mx-gx-1,my-gy-1,mz-gz-1};
  for (int i_f=0; i_f<nf; i_f++) {

    int * flux_index = 0;
    std::string field_name = field_names[i_f];
    
    if (field_name == "density") flux_index = dindex;
    if (field_name == "velocity_x") flux_index = uindex;
    if (field_name == "velocity_y") flux_index = vindex;
    if (field_name == "velocity_z") flux_index = windex;
    if (field_name == "total_energy") flux_index = Eindex;
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

      for (int face = 0; face < 2; face++) {
        flux_index[axis*2+face] = &flux_array[axis][face][i_f][0] - standard;
      }
    }
  }

  int n3[3] = {nx,ny,nz};
  for (int i_f=0; i_f<nf; i_f++) {
    for (int axis=0; axis<rank; axis++) {
      const int n1 = n3[(axis+1)%3];
      const int n2 = n3[(axis+2)%3];
      for (int face=0; face<2; face++) {
        auto & fa = flux_array[axis][face][i_f];
        double sumabs=0.0;
        for (int i1=0; i1<n1; i1++) {
          for (int i2=0; i2<n2; i2++) {
            const int i = i1 + n1*i2;
            sumabs += std::abs(fa[i]);
          }
        }
      }
    }
  }
  
#endif    

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

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    delete [] CellWidthTemp[dim];
  }

#ifdef NEW_FLUX  
#ifdef DEBUG_NEW_FLUX
  for (int i_f=0; i_f<nf; i_f++) {
    for (int axis=0; axis<rank; axis++) {
      const int n1 = n3[(axis+1)%3];
      const int n2 = n3[(axis+2)%3];
      for (int face=0; face<2; face++) {
        auto & fa = flux_array[axis][face][i_f];
        double sumabs=0.0;
        for (int i2=0; i2<n2; i2++) {
          for (int i1=0; i1<n1; i1++) {
            const int i = i1 + n1*i2;
            CkPrintf ("DEBUG_FLUX %s  %d %d %d  (%d %d) %g\n",
                      name().c_str(),axis,face,i_f,i1,i2,fa[i]);
            sumabs += std::abs(fa[i]);
          }
        }
      }
    }
  }
#endif
#endif
  
  /* deallocate temporary space for solver */

  delete [] temp;
  if (rank < 2) delete [] velocity_y;
  if (rank < 3) delete [] velocity_z;

  delete [] array;

  if (ncolour > 0) delete [] coloff;

  return ENZO_SUCCESS;

}
