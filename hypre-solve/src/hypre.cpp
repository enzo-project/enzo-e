//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Interface routines to HYPRE

/**
 * 
 * 
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 *
 */

#include <mpi.h>
#include <assert.h>
#include <string>
#include <vector>
#include <map>

#include "HYPRE_sstruct_ls.h"

#include "hdf5.h"

#include "hypre-solve.hpp"

//----------------------------------------------------------------------

const int debug        = 0;
const int trace        = 1;
const int trace_hypre  = 1;

//----------------------------------------------------------------------

#include "mpi.hpp"
#include "scalar.hpp"
#include "error.hpp"
#include "constants.hpp"
#include "point.hpp"
#include "faces.hpp"
#include "domain.hpp"
#include "grid.hpp"
#include "level.hpp"
#include "hierarchy.hpp"
#include "sphere.hpp"
#include "parameters.hpp"
#include "problem.hpp"
#include "hypre.hpp"
#include "error.hpp"

const Scalar matrix_scale = 1.0;  // 1.0:  1 1 1 -6 1 1 1

FILE *mpi_fp;
char mpi_file[20];

//======================================================================

void debug_print(Grid * grid)
{
  const int index = 2;
  if (debug) {
    int n0=grid->n(0);
    int n1=grid->n(1);
    int n2=grid->n(2);
    int i,i0,i1,i2;
    i0=0; i1=index; i2=index; i=Grid::index(i0,i1,i2,n0,n1,n2);
    printf ("X(0,%d,%d) = %g\n",index,index,(grid->values())[i]);
    i0=index; i1=0; i2=index; i=Grid::index(i0,i1,i2,n0,n1,n2);
    printf ("X(%d,0,%d) = %g\n",index,index,(grid->values())[i]);
    i0=index; i1=index; i2=0; i=Grid::index(i0,i1,i2,n0,n1,n2);
    printf ("X(%d,%d,0) = %g\n",index,index,(grid->values())[i]);
  } // debug
  return;
}
						     
//======================================================================
// PUBLIC MEMBER FUNCTIONS
//======================================================================

/// Hypre constructor

Hypre::Hypre (Parameters & parameters)
  : grid_(0),
    stencil_(0),
    graph_(0),
    A_(0),
    B_(0),
    X_(0),
    solver_(0),
    parameters_(parameters),
    resid_(-1.0),
    iter_(-1)
{
  if (trace_hypre) {
    sprintf (mpi_file,"hypre-solve.out.%d",pmpi->ip());
    mpi_fp = fopen (mpi_file,"w");
  } // trace_hypre
} // Hypre::Hypre

//----------------------------------------------------------------------

/// Initialize the Grid Hierarchy

/** Creates a hypre grid, with one part per level and one box per Grid
    patch object, for an AMR problem.  Sets grid box extents, grid
    part variables, and periodicity of the root-level grid part. */

void Hypre::init_hierarchy (Parameters & parameters,
			    Hierarchy  & hierarchy, 
			    Mpi        & mpi)
{

  int dim       = hierarchy.dimension();
  int num_parts = hierarchy.num_levels();

  // Create the hypre grid
  
  
  HYPRE_SStructGridCreate (MPI_COMM_WORLD, dim, num_parts, &grid_);
  if (trace_hypre) {
    fprintf (mpi_fp,"%s:%d %d HYPRE_SStructGridCreate (MPI_COMM_WORLD, %d, %d, %p)\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    dim, num_parts, &grid_);
    fflush(mpi_fp);
  } // trace_hypre

  ItHierarchyLevels itl (hierarchy);

  int part = 0;

  while (Level * level = itl++) {

    ItLevelGridsLocal itgl (*level);

    while (Grid * grid = itgl++) {

      int lower[3] = {grid->i_lower(0),grid->i_lower(1),grid->i_lower(2)};
      int upper[3] = {grid->i_upper(0),grid->i_upper(1),grid->i_upper(2)};

      // Set extents for boxes that comprise the hypre grid

      HYPRE_SStructGridSetExtents(grid_, part, lower, upper);
      if (trace_hypre) {fprintf (mpi_fp, "%s:%d %d HYPRE_SStructGridSetExtents(%p,%d, %d %d %d, %d %d %d);\n",
				__FILE__,__LINE__,pmpi->ip(),
				&grid_, part, 
				lower[0], lower[1], lower[2], 
				upper[0], upper[1], upper[2]
				);
	fflush(mpi_fp);
      } // trace_hypre
      
    } // while grid = itgl++

    // Create a single cell-centered variable for each grid part (level)

    HYPRE_SStructVariable variable_types[] = { HYPRE_SSTRUCT_VARIABLE_CELL };
    const int numvars = 1;

    HYPRE_SStructGridSetVariables(grid_, part, numvars, variable_types);
    if (trace_hypre) {
      fprintf (mpi_fp, "%s:%d %d HYPRE_SStructGridSetVariables(%p,%d,%d,%d);\n",
	      __FILE__,__LINE__,pmpi->ip(),
	      &grid_, part, numvars, variable_types[0]
	      );
      fflush(mpi_fp);
    } // trace_hypre

    // Set grid part to be periodic, with periodicity determined by the root
    // level size, current level, and refinement factor (ASSUMED TO BE 2)

    const int r = 2;
    int periodicity[3];

    // Determine periodicity of Level

    for (int i=0; i<3; i++) {
      periodicity[i] = hierarchy.level(0).zones(i);
      for (int k=0; k < part; k++) periodicity[i] *= r;
    }

    if (parameters.value("boundary") == "dirichlet") {
      periodicity[0] = 0;
      periodicity[1] = 0;
      periodicity[2] = 0;
    } 

    if (debug) printf ("%s:%d  Level = %d Periodicity = (%d,%d,%d)\n",
		       __FILE__,__LINE__, part,
		       periodicity[0],periodicity[1],periodicity[2]);

    HYPRE_SStructGridSetPeriodic (grid_, part, periodicity);
    if (trace_hypre) {
      fprintf (mpi_fp, "%s:%d %d HYPRE_SStructGridSetPeriodic (%p,%d,%d %d %d);\n",
	      __FILE__,__LINE__,pmpi->ip(),
	      &grid_, part, periodicity[0], periodicity[1], periodicity[2]
	      );
      fflush(mpi_fp);
    } // trace_hypre

    ++ part;
  } // while level = itl++

  // When finished, assemble the hypre grid

  HYPRE_SStructGridAssemble (grid_);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructGridAssemble (%p);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &grid_
	    );
    fflush(mpi_fp);
  } // trace_hypre
  
} // Hypre::init_hierarchy()

//----------------------------------------------------------------------

/// Initialize the discretization stencils.  

/** Creates and initializes a stencil object.  Supports 1, 2, or 3
    dimensional stencils. */

void Hypre::init_stencil (Hierarchy & hierarchy)

{

  int dim = hierarchy.dimension();

  HYPRE_SStructStencilCreate (dim,dim*2+1,&stencil_);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructStencilCreate (%d,%d,%p);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    dim,dim*2+1,&stencil_
	    );
    fflush(mpi_fp);
  } // trace_hypre

  int entries[][3] = { {  0, 0, 0 },     // center
		       {  1, 0, 0 },     // X+
		       { -1, 0, 0 },     // X-
		       {  0, 1, 0 },     // Y+
		       {  0,-1, 0 },     // Y-
		       {  0, 0, 1 },     // Z+
		       {  0, 0,-1 } };   // Z-

  if (dim >= 1) HYPRE_SStructStencilSetEntry (stencil_, 0, entries[0], 0);
  if (dim >= 1) HYPRE_SStructStencilSetEntry (stencil_, 1, entries[1], 0);
  if (dim >= 1) HYPRE_SStructStencilSetEntry (stencil_, 2, entries[2], 0);
  if (dim >= 2) HYPRE_SStructStencilSetEntry (stencil_, 3, entries[3], 0);
  if (dim >= 2) HYPRE_SStructStencilSetEntry (stencil_, 4, entries[4], 0);
  if (dim >= 3) HYPRE_SStructStencilSetEntry (stencil_, 5, entries[5], 0);
  if (dim >= 3) HYPRE_SStructStencilSetEntry (stencil_, 6, entries[6], 0);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructStencilSetEntry (%p,%d,%d %d %d,0)\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &stencil_, 0, entries[0][0], entries[0][1], entries[0][2]
	    );
    fflush(mpi_fp);
  } // trace_hypre
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructStencilSetEntry (%p,%d,%d %d %d,0)\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &stencil_, 1, entries[1][0], entries[1][1], entries[1][2]
	    );
    fflush(mpi_fp);
  }  // trace_hypre
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructStencilSetEntry (%p,%d,%d %d %d,0)\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &stencil_, 2, entries[2][0],entries[2][1],entries[2][2]
	    );
    fflush(mpi_fp);
  }  // trace_hypre
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructStencilSetEntry (%p,%d,%d %d %d,0)\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &stencil_, 3, entries[3][0],entries[3][1],entries[3][2]
	    );
    fflush(mpi_fp);
  }  // trace_hypre
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructStencilSetEntry (%p,%d,%d %d %d,0)\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &stencil_, 4, entries[4][0],entries[4][1],entries[4][2]
	    );
    fflush(mpi_fp);
  }  // trace_hypre
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructStencilSetEntry (%p,%d,%d %d %d,0)\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &stencil_, 5, entries[5][0],entries[5][1],entries[5][2]
	    );
    fflush(mpi_fp);
  }  // trace_hypre
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructStencilSetEntry (%p,%d,%d %d %d,0)\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &stencil_, 6, entries[6][0],entries[6][1],entries[6][2]
	    );
    fflush(mpi_fp);
  }  // trace_hypre

} // Hypre::init_stencil()

//----------------------------------------------------------------------

/// Initialize the graph.

/** Creates a graph containing the matrix non-zero structure.  Graph
    edges include both those for non-zeros from the stencil within
    each part (level), and non-zeros for graph entries connecting
    linked parts.

    Setting up the matrix nonzero structure is performed using the
    following steps:

    1. Define stencil connections within each level

    2. Define connections for unknowns adjacent to coarse unknowns

    3. Define connections for unknowns adjacent to fine unknowns 

    The matrix nonzero structure is generally nonsymmetric.  Only step
    1 is required for unigrid problems.

*/

void Hypre::init_graph (Hierarchy & hierarchy)

{
  // Create the hypre graph object

  HYPRE_SStructGraphCreate (MPI_COMM_WORLD, grid_, &graph_);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructGraphCreate (MPI_COMM_WORLD, %p,%p);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &grid_, &graph_
	    );
    fflush(mpi_fp);
  } // trace_hypre

  //  HYPRE_SStructGraphSetObjectType (graph_, HYPRE_SSTRUCT);
  //  if (trace_hypre) {
  //    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructGraphSetObjectType (%p,%d);\n",
  //	    __FILE__,__LINE__,pmpi->ip(),
  //	    &graph_, HYPRE_SSTRUCT
  //	    );
  //    fflush(mpi_fp);
  //  }  // trace_hypre

  ItHierarchyLevels itl (hierarchy);

  int part = 0;

  while (Level * level = itl++) {

    // 1. Define stencil connections within each level

    HYPRE_SStructGraphSetStencil (graph_, part, 0, stencil_);
    if (trace_hypre) {
      fprintf (mpi_fp, "%s:%d %d HYPRE_SStructGraphSetStencil (%p,%d,0,%p);\n",
	      __FILE__,__LINE__,pmpi->ip(),
	      &graph_, part, &stencil_
	      );
      fflush(mpi_fp);
    }  // trace_hypre

    ++ part;

    // 2. Define matrix nonzero structure connecting grids in level
    // with next-coarser.

    if (level->index() > 0) {
      ItLevelGridsAll itag (*level);

      while (Grid * grid = itag++) {

	// Define nonstencil entries for the grid

	init_graph_nonstencil_(*grid);

      } // while grid = itag++
    } // if level > 0

    ItLevelGridsAll itag (*level);

    while (Grid * grid = itag++) {

      // Clear the nonstencil entry counter for subsequent matrix
      // nonstencil entries

      int dim = hierarchy.dimension();
      grid->init_counter(dim*2+1);
    } // while grid = itag++
  } // while level = itl++

  // Assemble the graph

  HYPRE_SStructGraphAssemble (graph_);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructGraphAssemble (%p);\n",
	    __FILE__,__LINE__,pmpi->ip(),&graph_
	    );
    fflush(mpi_fp);
  } // trace_hypre

} // Hypre::init_graph()


//----------------------------------------------------------------------

/// Initialize the matrix A and right-hand-side vector b

/** Creates a matrix with a given non-zero structure, and sets nonzero
    values.

    Setting up the matrix elements is done with the following
    steps:

    1. Set stencil values

    2. Clean up stencil connections between parts and under overlapped grids

    3. Set values for unknowns between parent and children

    The matrix is generally nonsymmetric.  Only step 1 is required for
    unigrid problems.

*/

void Hypre::init_linear (Parameters          & parameters,
			 Hierarchy           & hierarchy,
			 std::vector<Point *>  points,
			 std::vector<Sphere *> spheres)

{
  // Create the hypre matrix A_, solution X_, and right-hand side B_ objects

  HYPRE_SStructMatrixCreate (MPI_COMM_WORLD, graph_, &A_);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructMatrixCreate (%p,%p);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &graph_, &A_
	    );
    fflush(mpi_fp);
  } // trace_hypre
  HYPRE_SStructVectorCreate (MPI_COMM_WORLD, grid_,  &X_);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructVectorCreate (%p,%p);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &grid_,  &X_
	    );
    fflush(mpi_fp);
  } // trace_hypre
  HYPRE_SStructVectorCreate (MPI_COMM_WORLD, grid_,  &B_);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructVectorCreate (%p,%p);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &grid_,  &B_
	    );
    fflush(mpi_fp);
  } // trace_hypre

  // Set the object types

  HYPRE_SStructMatrixSetObjectType (A_,HYPRE_SSTRUCT);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructMatrixSetObjectType (%p,%d);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &A_,HYPRE_SSTRUCT
	    );
    fflush(mpi_fp);
  } // trace_hypre
  HYPRE_SStructVectorSetObjectType (X_,HYPRE_SSTRUCT);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructVectorSetObjectType (%p,%d);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &X_,HYPRE_SSTRUCT
	    );
    fflush(mpi_fp);
  } // trace_hypre
  HYPRE_SStructVectorSetObjectType (B_,HYPRE_SSTRUCT);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructVectorSetObjectType (%p,%d);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &B_,HYPRE_SSTRUCT
	    );
    fflush(mpi_fp);
  } // trace_hypre

  // Initialize the hypre matrix and vector objects

  HYPRE_SStructMatrixInitialize (A_);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructMatrixInitialize (%p);\n",
	    __FILE__,__LINE__,pmpi->ip(),&A_
	    );
    fflush(mpi_fp);
  } // trace_hypre
  HYPRE_SStructVectorInitialize (X_);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructVectorInitialize (%p);\n",
	    __FILE__,__LINE__,pmpi->ip(),&X_
	    );
    fflush(mpi_fp);
  } // trace_hypre
  HYPRE_SStructVectorInitialize (B_);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructVectorInitialize (%p);\n",
	    __FILE__,__LINE__,pmpi->ip(),&B_
	    );
    fflush(mpi_fp);
  } // trace_hypre

  //--------------------------------------------------
  // Initialize the matrix A_
  //--------------------------------------------------

  ItHierarchyLevels itl (hierarchy);

  while (Level * level = itl++) {

    ItLevelGridsLocal itlg (*level);
    
    // 1. Set stencil values within level

    while (Grid * grid = itlg++) {

      init_matrix_stencil_(*grid);
    } // while grid = itlg++

    if (level->index() > 0) {

      //      // 2. Clean up stencil connections between levels

      //      init_matrix_clear_(level->index());

      // 3. Set matrix values between levels

      // WARNING: POSSIBLE SCALING ISSUE.  Below we loop over all grids;
      // however, we only need too loop over parent-child pairs such
      // that either child or parent is local to this MPI process.
 
      ItLevelGridsAll itag (*level);

      while (Grid * grid = itag++) {

 	init_matrix_nonstencil_(*grid);

      } // while grid = itag++
    } // while level > 0
  } // while level = itl++
  
  for (int part = 1; part < hierarchy.num_levels(); part++) {

    // 2. Clean up stencil connections between levels

    init_matrix_clear_(part);

  } // for part

  //--------------------------------------------------
  // Initialize B_ according to density
  //--------------------------------------------------

  Scalar local_shift_b_sum = 0.0;

  local_shift_b_sum += init_vector_points_  (hierarchy,points);
  local_shift_b_sum += init_vector_spheres_ (hierarchy,spheres);

  Scalar shift_b_sum = 0.0;

  MPI_Allreduce (&local_shift_b_sum, &shift_b_sum, 1, 
		 MPI_SCALAR, MPI_SUM, MPI_COMM_WORLD);

  // Shift B to zero out the null space if problem is periodic

  if ( parameters.value("boundary") == "periodic" ) {

    // Compute the shift

    int part = 0;
    long long shift_b_count = 0;
    ItHierarchyLevels itl (hierarchy);
    while (Level * level = itl++) {
      ItLevelGridsAll itg (*level);
      while (Grid * grid = itg++) {
	shift_b_count += grid->num_unknowns();
      } // while grid = itg++
    } // while level = itl++

    Scalar shift_b_amount = - shift_b_sum / shift_b_count;
    if (debug) printf ("Periodic shift = %g\n",shift_b_amount);

    // Perform the shift
    
    part = 0;
    while (Level * level = itl++) {
      ItLevelGridsLocal itg (*level);
      while (Grid * grid = itg++) {
	int lower[3] = { grid->i_lower(0), grid->i_lower(1), grid->i_lower(2) };
	int upper[3] = { grid->i_upper(0), grid->i_upper(1), grid->i_upper(2) };
	Scalar * values = new Scalar[grid->num_unknowns()];
	for (int i=0; i<grid->num_unknowns(); i++) values[i] = shift_b_amount;
	HYPRE_SStructVectorAddToBoxValues (B_,part,lower,upper,0,values);
	if (trace_hypre) {
	  fprintf (mpi_fp, "%s:%d %d HYPRE_SStructVectorAddToBoxValues (%p,%d, %d %d %d, %d %d %d, 0, %g...);\n",
		  __FILE__,__LINE__,pmpi->ip(),
		  &B_,part,
		  lower[0],lower[1],lower[2],
		  upper[0],upper[1],upper[2],0,values[0]
		  );
	  fflush(mpi_fp);
	} // trace_hypre

	delete [] values;
      } // while grid = itg++


      ++part;
    } // while level = itl++

    if (debug) printf ("%s:%d shift (count,sum,amount) = (%lld,%g,%g)\n",
		       __FILE__,__LINE__,shift_b_count,shift_b_sum,
		       shift_b_amount);
  } // if periodic

  // Assemble the matrix and vectors

  HYPRE_SStructMatrixAssemble (A_);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructMatrixAssemble (%p);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &A_
	    );
    fflush(mpi_fp);
  } // trace_hypre
  HYPRE_SStructVectorAssemble (B_);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructVectorAssemble (%p);\n",
	    __FILE__,__LINE__,pmpi->ip(),&B_
	    );
    fflush(mpi_fp);
  } // trace_hypre
  HYPRE_SStructVectorAssemble (X_);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructVectorAssemble (%p);\n",
	    __FILE__,__LINE__,pmpi->ip(),&X_
	    );
    fflush(mpi_fp);
  } // trace_hypre

  // Write the vector to a file for debugging

  if (parameters.value("dump_a") == "true") HYPRE_SStructMatrixPrint ("A",A_,0);
  if (parameters.value("dump_b") == "true") HYPRE_SStructVectorPrint ("B",B_,0);
  if (parameters.value("dump_x") == "true") HYPRE_SStructVectorPrint ("X0",X_,0);

} // Hypre::init_linear()

//----------------------------------------------------------------------

/// Initialize and solve the linear solver

void Hypre::solve (Parameters & parameters,
		   Hierarchy & hierarchy)

{
  std::string solver = parameters.value("solver");
  int         levels = hierarchy.num_levels();

  int    itmax  = 0;
  double restol = 0.0;

  // Check solver parameters
  std::string sitmax  = parameters.value("solver_itmax");
  std::string srestol = parameters.value("solver_restol");

  // If not defined, then define them
  if (sitmax == "")  parameters.add_parameter ("solver_itmax","200");
  if (srestol == "") parameters.add_parameter ("solver_restol","1e-6");

  // Set local variables
  itmax  = atoi(sitmax.c_str());
  restol = atof(srestol.c_str());

  if        (solver == "fac"  && levels > 1) {

    solve_fac_(hierarchy,itmax,restol);

  } else if (solver == "bicgstab") {

    solve_bicgstab_(hierarchy,itmax,restol);

  } else if (solver == "pfmg" && levels == 1) {

    solve_pfmg_(hierarchy,itmax,restol);

  } else {
    char error_message[100];
    sprintf (error_message, "Hypre::solve called with illegal combination of "
	     "solver %s on %d levels", solver.c_str(),levels);
    ERROR(error_message);
  }
  
  if (parameters.value("dump_x") == "true") HYPRE_SStructVectorPrint ("X",X_,1);

} // Hypre::solve()

//----------------------------------------------------------------------

/// Evaluate the success of the solve

void Hypre::evaluate (Hierarchy & hierarchy)

{

  char filename[80];

  ItHierarchyGridsLocal itg(hierarchy);
  while (Grid * grid = itg++) {
    int level = grid->level();
    int igg3[3][2];
    grid->indices(igg3);
    int lower[3],upper[3];
    lower[0] = igg3[0][0];
    lower[1] = igg3[1][0];
    lower[2] = igg3[2][0];
    upper[0] = igg3[0][1] - 1;
    upper[1] = igg3[1][1] - 1;
    upper[2] = igg3[2][1] - 1;
    int n3[3];
    n3[0] = grid->n(0);
    n3[1] = grid->n(1);
    n3[2] = grid->n(2);

    grid->allocate();

    HYPRE_SStructVectorGetBoxValues (X_,level,lower,upper,0,grid->values());  
    if (trace_hypre) {
      fprintf (mpi_fp, "%s:%d %d HYPRE_SStructVectorGetBoxValues (%p,%d, [%d,%d,%d], [%d,%d,%d],0,%g...%g);\n",
	       __FILE__,__LINE__,pmpi->ip(),
	       X_,level,
	       lower[0],lower[1],lower[2],
	       upper[0],upper[1],upper[2],
	       (grid->values())[0],(grid->values())[grid->n()-1]
	       );
      fflush(mpi_fp);
    } // trace_hypre

    debug_print(grid);
    if (grid->is_local()) sprintf (filename,"X.%d",grid->id());
    grid->write(filename);
    
    HYPRE_SStructVectorGetBoxValues (B_,level,lower,upper,0,grid->values());  
    if (trace_hypre) {
      fprintf (mpi_fp, "%s:%d %d HYPRE_SStructVectorGetBoxValues (%p,%d, [%d,%d,%d], [%d,%d,%d],0,%g...%g);\n",
	       __FILE__,__LINE__,pmpi->ip(),
	       B_,level,
	       lower[0],lower[1],lower[2],
	       upper[0],upper[1],upper[2],
	       (grid->values())[0],(grid->values())[grid->n()-1]
	       );
      fflush(mpi_fp);
    } // trace_hypre

    sprintf (filename,"B.%d",grid->id());
    grid->write(filename);

  } // grid = itg++

} // Hypre::evaluate()

//======================================================================
// PRIVATE MEMBER FUNCTIONS
//======================================================================

/// init_nonstencil_() is called twice: first by
/// init_graph_nonstencil_(), with phase == "graph", to set nonstencil
/// graph entries, and again by init_matrix_nonstencil_(), with phase
/// == "matrix", to set nonstencil matrix entries.

void Hypre::init_nonstencil_ (Grid & grid, std::string phase)
{
  _TRACE_;

  int id = grid.id();
  char filename[10];
  sprintf (filename,"grid.%03d",id);
  grid.write(filename);

  if (phase != "graph" && phase != "matrix") {
    char error_message[80];
    strcpy (error_message,"init_matrix_nonstencil_ called with phase = ");
    strcat (error_message,phase.c_str());
    ERROR(error_message);
  } // if phase unexpected

  const int r = 2; // WARNING: hard-coded refinement factor r = 2

  int face,axis,ig1,ig2;

  // global grid index limits

  int ig3[3][2];
  grid.indices(ig3);
  
  // Loop over each face zone in the grid, adding non-stencil graph
  // entries wherever a zone is adjacent to a coarse zone.  Both
  // fine-to-coarse and coarse-to-fine entries are added.

  int level_fine   = grid.level();
  int level_coarse = grid.level() - 1;

  bool discret_const = parameters_.value("discret") == "constant";

  assert (level_coarse >= 0);

  for (axis=0; axis<3; axis++) {

    // j0:    axis normal to face
    // j1,j2: axes within face

    int j0 = axis;
    int j1 = (axis+1)%3;
    int j2 = (axis+2)%3;

    // n0     grid size normal to face
    // n1,n2: grid size within face

    int n0 = ig3[j0][1] - ig3[j0][0];
    int n1 = ig3[j1][1] - ig3[j1][0];
    int n2 = ig3[j2][1] - ig3[j2][0];

    double h0 = grid.h(j0);
    double h1 = grid.h(j1);
    double h2 = grid.h(j2);

    double H0 = r*grid.h(j0);
    double H1 = r*grid.h(j1);
    double H2 = r*grid.h(j2);

    // ig3[][] should be divisible by r**level.  Just test r here.

    bool l0 = (ig3[j1][0]/r)*r == ig3[j1][0];
    bool l1 = (ig3[j1][1]/r)*r == ig3[j1][1];

    if (!l0) printf ("ig3[%d][0] = %d\n",j1,ig3[j1][0]);
    assert (l0);
    if (!l1) printf ("ig3[%d][1] = %d\n",j1,ig3[j1][1]);
    assert (l1);

    for (face=0; face<2; face++) {

      // Loop over face zones that are aligned with coarse zones (hence "+= r")

      for (ig1=0; ig1<n1; ig1 += r) {
	for (ig2=0; ig2<n2; ig2 += r) {

	  Grid * adjacent   = grid.faces().adjacent(axis,face,ig1,ig2);

	  Faces::Label & fz = grid.faces().label(axis,face,ig1,ig2);

	  // Add graph entries iff grid or adjacent grid is local, and if
	  // adjacent grid (if it exists) is in the next-coarser level

	  bool is_local = 
	    (adjacent != NULL) &&  (adjacent->is_local() || grid.is_local());

	  if (is_local && fz == Faces::_coarse_) {

	    // (fine) grid global indices

	    int igg3[3]; 

	    igg3[j0] = ig3[j0][0] + face*(n0 - r);
	    igg3[j1] = ig3[j1][0] + ig1;
	    igg3[j2] = ig3[j2][0] + ig2;

	    // (coarse) adjacent global indices

	    int ign3[3]; 

	    ign3[j0] = (igg3[j0]) / r  + (face*r-1);
	    ign3[j1] = (igg3[j1]) / r;
	    ign3[j2] = (igg3[j2]) / r;

	    //--------------------------------------------------
	    // GRAPH ENTRY: FINE-TO-COARSE 
	    //--------------------------------------------------

	    if (discret_const) {

	      //--------------------------------------------------
	      // (*) CONSTANT
	      //     Scale        = 2/3
	      //     Coefficients = 1
	      //--------------------------------------------------

	      int diggs[][3] = {{face*(r-1),0,0},
				{0,1,0},
				{0,0,1},
				{0,-1,0},
				{-face*(r-1),0,-1}};

	      if (grid.is_local()) {

		if (debug) {
		  printf ("ip=%d %s:%d fine-coarse %d - %d\n",
			  pmpi->ip(),__FILE__,__LINE__,grid.id(),adjacent->id());
		} // debug
		if (phase == "graph") {

		  int k=0;

		  igg3[j0]+=diggs[k][0];
		  igg3[j1]+=diggs[k][1];
		  igg3[j2]+=diggs[k][2];
		  for (int k=1; k<5; k++) {
		    HYPRE_SStructGraphAddEntries 
		      (graph_, level_fine, igg3, 0, level_coarse, ign3, 0);
		    if (trace_hypre) {
		      fprintf (mpi_fp, "%s:%d %d HYPRE_SStructGraphAddEntries (%p,%d, [%d,%d,%d] 0, %d, [%d,%d,%d],0);\n",
			      __FILE__,__LINE__,pmpi->ip(),
			      &graph_, 
			       level_fine, 
			       igg3[0], igg3[1], igg3[2],
			       level_coarse,
			       ign3[0], ign3[1], ign3[2]
			      );
		      fflush(mpi_fp);
		    } // trace_hypre
		    igg3[j0]+=diggs[k][0];
		    igg3[j1]+=diggs[k][1];
		    igg3[j2]+=diggs[k][2];
		  } // for k = 1:4

		} // phase == "graph"
		else if (phase == "matrix") {

		  // fine->coarse off-diagonal

		  double val_h = h1*h2/h0;
		  double val_s = 2. / 3.;

		  int entry;
		  double val_a;
		  double val;

		  int k=0;

		  igg3[j0]+=diggs[k][0];
		  igg3[j1]+=diggs[k][1];
		  igg3[j2]+=diggs[k][2];

		  for (int k=1; k<5; k++) {

		    val_a = 1.0; // DIFFUSION COEFFICIENT GOES HERE

		    // Update off-diagonal

		    entry = grid.counter(igg3)++;
		    val   = matrix_scale * val_h * val_s * val_a;
		    HYPRE_SStructMatrixAddToValues 
		      (A_, level_fine, igg3, 0, 1, &entry, &val);
		    if (trace_hypre) {
		      fprintf (mpi_fp, "%s:%d %d HYPRE_SStructMatrixAddToValues (%p, %d, 0, 1, %d, %g);\n",
			      __FILE__,__LINE__,pmpi->ip(),
			      &A_, 
			      level_fine,
			       //			      igg3[0], igg3[1], igg3[2],
			      entry, val
			      );
		      fflush(mpi_fp);
		    } // trace_hypre

		    // Update diagonal

		    entry = 0;
		    val = -val;
		    HYPRE_SStructMatrixAddToValues 
		      (A_, level_fine, igg3, 0, 1, &entry, &val);
		    if (trace_hypre) {
		      fprintf (mpi_fp, "%s:%d %d HYPRE_SStructMatrixAddToValues (%p, %d, 0, 1, %d, %g\n",
			      __FILE__,__LINE__,pmpi->ip(),
			      &A_, 
			      level_fine, 
			       //			      igg3[0], igg3[1], igg3[2],
			      entry, val
			      );
		      fflush(mpi_fp);
		    } // trace_hypre

		    igg3[j0]+=diggs[k][0];
		    igg3[j1]+=diggs[k][1];
		    igg3[j2]+=diggs[k][2];
		  } // for k = 1,4
		} // if phase == "matrix"
	      } // if grid.is_local()
	    } // if discret_const
	    else {
	      char error_message[80];
	      strcpy (error_message,"Unknown parameter discret = ");
	      strcat (error_message,parameters_.value("discret").c_str());
	      ERROR(error_message);
	    } // if discret unexpected

	    //--------------------------------------------------
	    // GRAPH ENTRY: COARSE-TO-FINE
	    //--------------------------------------------------

	    _TRACE_;
	    if (adjacent->is_local()) {
	      _TRACE_;
	      if (phase == "graph") {

		int diggs[][3] = {{1,0,0},
				  {0,1,0},
				  {-1,0,0},
				  {0,0,1},
				  {1,0,0},
				  {0,-1,0},
				  {-1,0,0},
				  {0,0,-1}};

		for (int k=0; k<8; k++) {
		  // DIES HERE IN r266 FOR e.g. N16.P211.L2.O1.S0.cg
		  // grid=3  adjacent=0
		  //  k=0 j=(0 1 2)
		  //  hypre.cpp:1031 0 HYPRE_SStructGraphAddEntries (0x7fff7354db68,0, [8,0,0] 0, 1, [14,0,0] 0);
		  // Corresponds to coarse zone on processor 0 connecting
		  // to fine zone on processor 1 where fine grid parent
		  // is adjacent to coarse grid
		  //
		  //   zones        grids      procs
		  // +-------+   +-------+   +-------+   
		  // |   |   |   | 0 | 1 |   | 0 | 1 |
		  // +---+   |   +---+   |   +---+   |
		  // | |X|X  |   |2|3|   |   |2|3|   |
		  // +---+---+   +---+---+   +---+---+
		  // 
		  // Problem is in hierarchy.cpp.  If either grid 3 or grid 1
		  // were on processor 0, then it would work.  Neither are,
		  // so grid 3 sets face cells that are neither flagged as
		  // neighbors or coarse neighbors as parent.  Since parent 0
		  // is on processor 0, adjacent->is_local() is true,
		  // but grid 3 and "adjacent" (mislabled as parent 0)
		  // are not really neighbors, and adding graph entries
		  // goes out-of-bounds for parent grid.

		  if (trace_hypre) {
		    fprintf (mpi_fp, "grid=%d  adjacent=%d\n",grid.id(),adjacent->id());
		    fprintf (mpi_fp, "k=%d j=(%d %d %d)\n",k,j0,j1,j2);
		    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructGraphAddEntries (%p,%d, [%d,%d,%d] 0, %d, [%d,%d,%d] 0);\n",
			     __FILE__,__LINE__,pmpi->ip(),
			    &graph_,
			    level_coarse, 
			     ign3[0], ign3[1], ign3[2],
			     level_fine,
			       igg3[0], igg3[1], igg3[2]
			    );
		    fprintf (mpi_fp, "axis=%d face=%d\n",axis,face);
		    fflush(mpi_fp);
		  } // trace_hypre
		  HYPRE_SStructGraphAddEntries 
		    (graph_, level_coarse, ign3, 0, level_fine, igg3, 0);
		  if (trace_hypre) {
		    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructGraphAddEntries (%p,%d, [%d,%d,%d] 0, %d, [%d,%d,%d] 0);\n",
			    __FILE__,__LINE__,pmpi->ip(),
			    &graph_,
			    level_coarse, 
			     ign3[0], ign3[1], ign3[2],
			     level_fine,
			       igg3[0], igg3[1], igg3[2]
			    );
		    fflush(mpi_fp);
		  } // trace_hypre
		  igg3[0] += diggs[k][0];
		  igg3[1] += diggs[k][1];
		  igg3[2] += diggs[k][2];
		} // for k=0,7

	      } // if phase == "graph"
	      else if (phase == "matrix") {

		double val_h = H1*H2/H0;
		double val_s = 1.;
		double val_a = 1.0; // DIFFUSION COEFFICIENT GOES HERE
		double val   = matrix_scale * val_h * val_s * val_a;
		int    entry;
		double value;

		// coarse->fine off-diagonal

		for (int i=0; i<8; i++) {
		  entry = adjacent->counter(ign3)++;
		  value = (1./8.) * val;
		  HYPRE_SStructMatrixAddToValues 
		    (A_, level_coarse,ign3, 0, 1, &entry, &value);
		  if (trace_hypre) {
		    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructMatrixAddToValues (%p,%d, [%d %d %d], 0, 1, %d %g);\n",
			    __FILE__,__LINE__,pmpi->ip(),
			    &A_, level_coarse,
			     ign3[0],ign3[1],ign3[2], 
			     entry, value
			    );
		    fflush(mpi_fp);
		  } // trace_hypre
		} // for i=0,7

		// coarse->coarse diagonal

		double val_diag = -val;
		entry = 0;

		//	      _TEMPORARY_;
		//	      val_diag*=0.00;
		// HYPRE_SStructMatrixAddToValues 
		// (A_, level_coarse, ign3, 0, 1, &entry, &val_diag);

	      } // if phase == "matrix"
	    } // if adjacent.is_local()
	  } // if is_local && fz == Faces::_coarse_
	} // for ig2
      } // for ig1
    } // for face
  } // for axis
} // Hypre::init_nonstencil_()

//------------------------------------------------------------------------

/// Set matrix stencil values for the grid interior

void Hypre::init_matrix_stencil_ (Grid & grid)

{
  int low[3]     = { grid.i_lower(0), grid.i_lower(1), grid.i_lower(2) };
  int up[3]      = { grid.i_upper(0), grid.i_upper(1), grid.i_upper(2) };
  int n          = grid.num_unknowns();
  int entries[7] = { 0,1,2,3,4,5,6 };
  double h3[3]   = {grid.h(0),grid.h(1),grid.h(2)};
  int    n3[3]   = {grid.n(0),grid.n(1),grid.n(2)};

  double h120 = h3[1]*h3[2] / h3[0];
  double h201 = h3[2]*h3[0] / h3[1];
  double h012 = h3[0]*h3[1] / h3[2];

  double hhh = h3[0]*h3[1]*h3[2];

  double * v0;         // Diagonal elements
  double * v1[3][2];   // Off-diagonal elements

  // Allocate storage

  v0 = new double [n];
  for (int axis=0; axis<3; axis++) {
    for (int face=0; face<2; face++) {
      v1[axis][face] = new double [n];
    } // for face
  } // for axis

  //-----------------------------------------------------------
  // Set stencil for all unknowns, ignoring boundary conditions
  //-----------------------------------------------------------

  int i0,i1,i2,i;
  for (i2 = 0; i2 < n3[2]; i2++) {

    // DIFFUSION COEFFICIENTS HERE

    double azp = 1.0;
    double azm = 1.0;

    for (i1 = 0; i1 < n3[1]; i1++) {

      // DIFFUSION COEFFICIENTS HERE

      double ayp = 1.0;
      double aym = 1.0;

      for (i0 = 0; i0 < n3[0]; i0++) {

	// DIFFUSION COEFFICIENTS HERE

	double axp = 1.0;
	double axm = 1.0;

	i = Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);

	v1[0][0][i] = matrix_scale * h120 * axm;
	v1[0][1][i] = matrix_scale * h120 * axp;
	v1[1][0][i] = matrix_scale * h201 * aym;
	v1[1][1][i] = matrix_scale * h201 * ayp;
	v1[2][0][i] = matrix_scale * h012 * azm;
	v1[2][1][i] = matrix_scale * h012 * azp;

	v0[i] = -( v1[0][0][i] + v1[0][1][i] +
		   v1[1][0][i] + v1[1][1][i] + 
		   v1[2][0][i] + v1[2][1][i] );

      } // for i0
    } // for i1
  } // for i2

  if (debug) {
    Scalar sum3[3];
    i0=0;
    sum3[0]=0.0;
    for (i1=0; i1<n3[1]; i1++) {
      for (i2=0; i2<n3[2]; i2++) {
	i=Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	sum3[0] += v1[0][0][i];
      } // for i2
    } // for i1
    i1=0;
    sum3[1]=0.0;
    for (i0=0; i0<n3[0]; i0++) {
      for (i2=0; i2<n3[2]; i2++) {
	i=Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	sum3[1] += v1[1][0][i];
      }// for i2
    } // for i0
    i2=0;
    sum3[2]=0.0;
    for (i0=0; i0<n3[0]; i0++) {
      for (i1=0; i1<n3[1]; i1++) {
	i=Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	sum3[2] += v1[2][0][i];
      } // for i1
    } // for i0
  } // if debug
  //-----------------------------------------------------------
  // Adjust stencil at grid boundaries
  //-----------------------------------------------------------

  Faces & faces = grid.faces();

  int axis,face;

  int level = grid.level();
  int debug3[3];
  debug3[0]=0;
  debug3[1]=0;
  debug3[2]=0;
  // X faces
  for (face=0; face<2; face++) {
    axis = 0;
    i0 = (face==0) ? 0 : n3[0]-1;
    for (i1=0; i1<n3[1]; i1++) {
      for (i2=0; i2<n3[2]; i2++) {
	int f = faces.label (axis,face,i1,i2);
	if (f != Faces::_boundary_ && f != Faces::_neighbor_) {
	  // Clear off-diagonal element
	  // DIFFUSION COEFFICIENTS HERE
	  double axp = 1.0;
	  double axm = 1.0;
	  double ax = (face==0) ? axm : axp;
	  i = Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  Scalar a = matrix_scale * h120 * ax;
	  v1[axis][face][i] -= a;
	  v0[i]             += a;
	} // if f neither boundary nor neighbor
      } // for i2
    } // for i1
    axis = 1;
    i1 = (face==0) ? 0 : n3[1]-1;
    for (i2=0; i2<n3[2]; i2++) {
      for (i0=0; i0<n3[0]; i0++) {
	int f = faces.label (axis,face,i2,i0);
	// Clear off-diagonal element
	if (f != Faces::_boundary_ && f != Faces::_neighbor_) {
	  i  = Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  // DIFFUSION COEFFICIENTS HERE
	  double ayp = 1.0;
	  double aym = 1.0;
	  double ay = (face==0) ? aym : ayp;
	  Scalar a = matrix_scale * h201 * ay;
	  v1[axis][face][i] -= a;
	  v0[i]             += a;
	} // if f neither boundary nor neighbor
      } // for i0
    } // for i2
    // Z faces
    axis = 2;
    i2 = (face==0) ? 0 : n3[2]-1;
    for (i0=0; i0<n3[0]; i0++) {
      for (i1=0; i1<n3[1]; i1++) {
	int f = faces.label (axis,face,i0,i1);
	if (f != Faces::_boundary_ && f != Faces::_neighbor_) {
	  // Clear off-diagonal element
	  // DIFFUSION COEFFICIENTS HERE
	  double azp = 1.0;
	  double azm = 1.0;
	  double az = (face==0) ? azm : azp;
	  Scalar a = matrix_scale * h012 * az;
	  i = Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);
	  v1[axis][face][i] -= a;
	  v0[i]             += a;
	} // if f neither boundary nor neighbor
      } // for i1
    } // for i0
  } // for face

  if (debug && level==1) {
    i = Grid::index(0,0,0,n3[0],n3[1],n3[2]);
    printf ("DEBUG %s:%d 000 %g %g %g  %g  %g %g %g\n",
			__FILE__,__LINE__,
			v0[i],
			v1[0][1][i],v1[0][0][i],
			v1[1][1][i],v1[1][0][i],
			v1[2][1][i],v1[2][0][i]);

    i = Grid::index(1,1,1,n3[0],n3[1],n3[2]);
    printf ("DEBUG %s:%d 111 %g %g %g  %g  %g %g %g\n",
			__FILE__,__LINE__,
			v0[i],
			v1[0][1][i],v1[0][0][i],
			v1[1][1][i],v1[1][0][i],
			v1[2][1][i],v1[2][0][i]);

  } // if debug && level == 1
  HYPRE_SStructMatrixSetBoxValues (A_,level,low,up,0,1,&entries[0],v0);
  HYPRE_SStructMatrixSetBoxValues (A_,level,low,up,0,1,&entries[1],v1[0][1]);
  HYPRE_SStructMatrixSetBoxValues (A_,level,low,up,0,1,&entries[2],v1[0][0]);
  HYPRE_SStructMatrixSetBoxValues (A_,level,low,up,0,1,&entries[3],v1[1][1]);
  HYPRE_SStructMatrixSetBoxValues (A_,level,low,up,0,1,&entries[4],v1[1][0]);
  HYPRE_SStructMatrixSetBoxValues (A_,level,low,up,0,1,&entries[5],v1[2][1]);
  HYPRE_SStructMatrixSetBoxValues (A_,level,low,up,0,1,&entries[6],v1[2][0]);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructMatrixSetBoxValues (%p,%d,%d %d %d, %d %d %d, 0,1,%d,%g);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &A_,level,
	    low[0], low[1], low[2],
	    up[0],up[1],up[2],entries[0],v0);
    fflush(mpi_fp);
  } // trace_hypre
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructMatrixSetBoxValues (%p,%d,%d %d %d, %d %d %d, 0,1,%d,%g);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &A_,level,
	    low[0], low[1], low[2],
	    up[0],up[1],up[2],entries[1],v1[0][1]);
    fflush(mpi_fp);
  } // trace_hypre
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructMatrixSetBoxValues (%p,%d,%d %d %d, %d %d %d, 0,1,%d,%g);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &A_,level,
	    low[0], low[1], low[2],
	    up[0],up[1],up[2],entries[2],v1[0][0]);
    fflush(mpi_fp);
  } // trace_hypre
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructMatrixSetBoxValues (%p,%d,%d %d %d, %d %d %d, 0,1,%d,%g);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &A_,level,
	    low[0], low[1], low[2],
	    up[0],up[1],up[2],entries[3],v1[1][1]);
    fflush(mpi_fp);
  } // trace_hypre
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructMatrixSetBoxValues (%p,%d,%d %d %d, %d %d %d, 0,1,%d,%g);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &A_,level,
	    low[0], low[1], low[2],
	    up[0],up[1],up[2],entries[4],v1[1][0]);
    fflush(mpi_fp);
  } // trace_hypre
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructMatrixSetBoxValues (%p,%d,%d %d %d, %d %d %d, 0,1,%d,%g);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &A_,level,
	    low[0], low[1], low[2],
	    up[0],up[1],up[2],entries[5],v1[2][1]);
    fflush(mpi_fp);
  } // trace_hypre
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructMatrixSetBoxValues (%p,%d,%d %d %d, %d %d %d, 0,1,%d,%g);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &A_,level,
	    low[0], low[1], low[2],
	    up[0],up[1],up[2],entries[6],v1[2][0]);
    fflush(mpi_fp);
  } // trace_hypre
  delete [] v0;
  for (axis=0; axis<3; axis++) {
    for (face=0; face<2; face++) {
      delete [] v1[axis][face];
    } // for face
  } // for axis

} // Hypre::init_matrix_stencil_()

//------------------------------------------------------------------------

/// Clean up stencil connections between parts for FAC solver

void Hypre::init_matrix_clear_ (int part)
{
  // WARNING: hard-coding refinement factor of 2
  int r_factors[3] = {2,2,2}; 
  //  if (part > 0) {

    // Clear stencil values from coarse to fine part

  
  printf ("HYPRE_SStructFACZeroCFSten\n");
  HYPRE_SStructFACZeroCFSten (A_,grid_, part, r_factors);
    if (trace_hypre) {
      fprintf (mpi_fp, "%s:%d %d HYPRE_SStructFACZeroCFSten (%p,%p,%d,%d %d %d);\n",
	      __FILE__,__LINE__,pmpi->ip(),
	      &A_,&grid_, part, r_factors[0],r_factors[1],r_factors[2]
	      );
      fflush(mpi_fp);
    } // trace_hypre

    // Clear stencil values from fine to coarse part

    printf ("HYPRE_SStructFACZeroFCSten\n");
    HYPRE_SStructFACZeroFCSten (A_,grid_, part);
    if (trace_hypre) {
      fprintf (mpi_fp, "%s:%d %d HYPRE_SStructFACZeroFCSten (%p,%p,%d);\n",
	      __FILE__,__LINE__,pmpi->ip(),
	      &A_,&grid_, part
	      );
      fflush(mpi_fp);
    } // trace_hypre

    // Set overlapped areas of part with identity

    if (part > 0) {

      printf ("HYPRE_SStructFACZeroAMRMatrixData\n");
      HYPRE_SStructFACZeroAMRMatrixData (A_, part-1, r_factors);
    if (trace_hypre) {
      fprintf (mpi_fp, "%s:%d %d HYPRE_SStructFACZeroAMRMatrixData (%p,%d,%d %d %d);\n",
	      __FILE__,__LINE__,pmpi->ip(),
	      &A_, part-1, r_factors[0], r_factors[1], r_factors[2]
	      );
      fflush(mpi_fp);
    } // trace_hypre
    } // part > 0

    // Need to clear under rhs also
    //   HYPRE_SStructFACZeroAMRVectorData(B_, plevels, prefinements);

} // Hypre::init_matrix_clear_()

//------------------------------------------------------------------------

/// Add contributions from point sources to right-hand side B

Scalar Hypre::init_vector_points_ (Hierarchy            & hierarchy,
				   std::vector<Point *> & points)

{

  const Scalar scaling0 = -4.0*Constants::G()*Constants::pi();

  Scalar shift_b_sum = 0.0;

  int i;
  for (i=0; i<int(points.size()); i++) {
    Point & point      = *points[i];
    Grid & grid        = hierarchy.grid(point.igrid());
    if (grid.is_local()) {

      Scalar cell_volume = grid.h(0) * grid.h(1) * grid.h(2);
      Scalar density     = point.mass() / cell_volume;
      Scalar value       = scaling0 * density;

      // Add contribution of the point to the right-hand side vector

      int index[3];
      Scalar lower[3],upper[3];
      grid.x_lower(lower[0],lower[1],lower[2]);
      grid.x_upper(upper[0],upper[1],upper[2]);
      for (int k=0; k<3; k++) {
	Scalar ap = point.x(k)      - lower[k];
	Scalar ag = upper[k] - lower[k];
	int    ig = grid.num_unknowns(k);
	int    i0 = grid.i_lower(k);
	index[k] = int (ap/ag*ig) + i0;
      } // for k=0,2
      if (index[0] < grid.i_lower(0) || index[0] > grid.i_upper(0) ||
	  index[1] < grid.i_lower(1) || index[1] > grid.i_upper(1) ||
	  index[2] < grid.i_lower(2) || index[2] > grid.i_upper(2)) {
	printf ("WARNING: Point apparently not in grid: \n");
	printf ("WARNING:    Point: (%g,%g,%g)\n",point.x(0),point.x(1),point.x(2));
	printf ("WARNING:    Grid:  (%g,%g,%g) - (%g,%g,%g)\n",
		lower[0],lower[1],lower[2],upper[0],upper[1],upper[2]);
      } // if index not in grid
      if (debug) {
	point.print();
	grid.print();
	printf ("Point index  = %d %d %d)\n",index[0],index[1],index[2]);
	printf ("Cell size    = %g %g %g\n",grid.h(0),grid.h(1),grid.h(2));
	printf ("Cell volume  = %g\n",cell_volume);
	printf ("Cell density = %g\n",density);
	printf ("RHS contribution = %g\n",value);
      } // if debug
    
      shift_b_sum += value;

      HYPRE_SStructVectorAddToValues (B_, grid.level(), index, 0, &value);
      if (trace_hypre) {
	fprintf (mpi_fp, "%s:%d %d HYPRE_SStructVectorAddToValues (%p,%d, [%d %d %d], 0, %g);\n",
		__FILE__,__LINE__,pmpi->ip(),
		&B_, grid.level(), index[0],index[1],index[2],value
		);
	fflush(mpi_fp);
      } // if trace_hypre
    } // if grid.is_local()
  } // for i=0 to # points
  return shift_b_sum;
} // Hypre::init_vector_points_()

//------------------------------------------------------------------------

/// Add contributions from sphere sources to right-hand side B

Scalar Hypre::init_vector_spheres_ (Hierarchy             & hierarchy,
				    std::vector<Sphere *> & spheres)

{

  if (spheres.size() > 0) {  
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    NOT_IMPLEMENTED("Contribution of sphere mass to right-hand side");
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  }

  return 0.0;
} // Hypre::init_vector_spheres_()

//------------------------------------------------------------------------

/// Initialize the PFMG hypre solver

void Hypre::solve_pfmg_ (Hierarchy & hierarchy, int itmax, double restol)

{

  // Create and initialize the solver

  HYPRE_SStructSysPFMGCreate    (MPI_COMM_WORLD, &solver_);

  // stopping criteria

  if (itmax != 0 )   HYPRE_SStructSysPFMGSetMaxIter(solver_,itmax);
  if (restol != 0.0) HYPRE_SStructSysPFMGSetTol(solver_,    restol);

  HYPRE_SStructSysPFMGSetLogging(solver_, 1);
  HYPRE_SStructSysPFMGSetup     (solver_,A_,B_,X_);

  // Solve the linear system

  HYPRE_SStructSysPFMGSolve     (solver_,A_,B_,X_);

  // Write out some diagnostic info about the solve

  HYPRE_SStructSysPFMGGetNumIterations (solver_,&iter_);
  HYPRE_SStructSysPFMGGetFinalRelativeResidualNorm (solver_,&resid_);

  printf ("HYPRE_SStructSysPFMGSolve num iterations: %d\n",iter_);
  printf ("HYPRE_SStructSysPFMGSolve final relative residual norm: %g\n",resid_);

  // Delete the solver

  HYPRE_SStructSysPFMGDestroy (solver_);
  solver_ = 0;

} // Hypre::solve_pfmg_()

//------------------------------------------------------------------------

/// Initialize the FAC hypre solver

void Hypre::solve_fac_ (Hierarchy & hierarchy, int itmax, double restol)

{
  int i;

  const int r = 2; // WARNING: hard-coded refinement factor r = 2

  // Create the solver

  HYPRE_SStructFACCreate(MPI_COMM_WORLD, &solver_);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructFACCreate(MPI_COMM_WORLD, %p);\n",
	    __FILE__,__LINE__,pmpi->ip(),&solver_
	    );
    fflush(mpi_fp);
  } // trace_hypre

  // Initialize parts

  int num_parts = hierarchy.num_levels();
  HYPRE_SStructFACSetMaxLevels(solver_,  num_parts);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructFACSetMaxLevels(%p,%d);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &solver_,  num_parts
	    );
    fflush(mpi_fp);
  } // trace_hypre
  int *parts  = new int [num_parts];
  for (i=0; i<num_parts; i++) parts[i] = i;
  HYPRE_SStructFACSetPLevels(solver_, num_parts, parts);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructFACSetPLevels(%p,%d, %d %d);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &solver_, num_parts, parts[0],parts[1]
	    );
    fflush(mpi_fp);
  } // trace_hypre

  // Initialize refinement factors

  typedef int int3[3];
  int3 *refinements = new int3 [num_parts];
  
  for (i=0; i<num_parts; i++) {
    refinements[i][0] = r;
    refinements[i][1] = r;
    refinements[i][2] = r;
  } // for i=0 to num_parts-1

  HYPRE_SStructFACSetPRefinements(solver_, num_parts, refinements);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructFACSetPRefinements(%p,%d,[%d %d %d],[%d %d %d]);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &solver_, num_parts, 
	     refinements[0][0], 
	     refinements[0][1], 
	     refinements[0][2],
	     refinements[1][0], 
	     refinements[1][1], 
	     refinements[1][2]
	    );
    fflush(mpi_fp);
  } // trace_hypre


  // solver parameters

  int npre   = 2;
  int npost  = 2;
  int csolve = 1;
  int relax  = 2;

  HYPRE_SStructFACSetNumPreRelax(solver_,      npre);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructFACSetNumPreRelax(%p,%d);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &solver_,npre
	    );
    fflush(mpi_fp);
  } // trace_hypre
  HYPRE_SStructFACSetNumPostRelax(solver_,     npost);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructFACSetNumPostRelax(%p,%d);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &solver_,     npost
	    );
    fflush(mpi_fp);
  } // trace_hypre
  HYPRE_SStructFACSetCoarseSolverType(solver_, csolve);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructFACSetCoarseSolverType(%p,%d);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &solver_, csolve
	    );
    fflush(mpi_fp);
  } // trace_hypre
  HYPRE_SStructFACSetRelaxType(solver_,        relax);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructFACSetRelaxType(%p,%d);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &solver_,        relax
	    );
    fflush(mpi_fp);
  } // trace_hypre

  // stopping criteria

  if (itmax != 0 ) {
    HYPRE_SStructFACSetMaxIter(solver_,itmax);
    if (trace_hypre) {
      fprintf (mpi_fp, "%s:%d %d HYPRE_SStructFACSetMaxIter(%p,%d);\n",
	      __FILE__,__LINE__,pmpi->ip(),
	      &solver_,itmax
	      );
      fflush(mpi_fp);
    } // trace_hypre
  } // if itmax != 0
  if (restol != 0.0) {
    HYPRE_SStructFACSetTol(solver_,    restol);
    if (trace_hypre) {
      fprintf (mpi_fp, "%s:%d %d HYPRE_SStructFACSetTol(%p,%g);\n",
	      __FILE__,__LINE__,pmpi->ip(),
	      &solver_,    restol
	      );
      fflush(mpi_fp);
    } // trace_hypre
  } // if restol != 0.0

  // output amount

  HYPRE_SStructFACSetLogging(solver_, 1);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructFACSetLogging(%p,1);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &solver_
	    );
    fflush(mpi_fp);
  } // trace_hypre

  // prepare for solve

  // DIES HERE 2008-03-07
  HYPRE_SStructFACSetup2(solver_, A_, B_, X_);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructFACSetup2(%p,%p,%p,%p);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &solver_, &A_, &B_, &X_
	    );
    fflush(mpi_fp);
  } // trace_hypre

  // Solve the linear system

  HYPRE_SStructFACSolve3(solver_, A_, B_, X_);
  if (trace_hypre) {
    fprintf (mpi_fp, "%s:%d %d HYPRE_SStructFACSolve3(%p,%p,%p,%p);\n",
	    __FILE__,__LINE__,pmpi->ip(),
	    &solver_, &A_, &B_, &X_
	    );
    fflush(mpi_fp);
  } // trace_hypre

  // Write out some diagnostic info about the solve

  HYPRE_SStructFACGetNumIterations(solver_, &iter_);
  HYPRE_SStructFACGetFinalRelativeResidualNorm(solver_, &resid_);

  printf ("HYPRE_SStructFACSolve3 num iterations: %d\n",iter_);
  printf ("HYPRE_SStructFACSolve3 final relative residual norm: %g\n",resid_);

  // Delete the solver

  HYPRE_SStructFACDestroy2(solver_);
  solver_ = 0;

  // Delete local dynamic storage
  delete [] parts;
  delete [] refinements;
} // Hypre::solve_fac_()

//------------------------------------------------------------------------

/// Initialize the BICGSTAB hypre solver

void Hypre::solve_bicgstab_ (Hierarchy & hierarchy, int itmax, double restol)

{
  _TRACE_;

  // Create the solver

  HYPRE_SStructBiCGSTABCreate(MPI_COMM_WORLD, &solver_);

  _TRACE_;

  // stopping criteria

  if (itmax != 0 )   HYPRE_SStructBiCGSTABSetMaxIter(solver_,itmax);
  if (restol != 0.0) HYPRE_SStructBiCGSTABSetTol(solver_,    restol);

  // output amount

  HYPRE_SStructBiCGSTABSetLogging(solver_, 1);

  // Initialize the solver

  HYPRE_SStructBiCGSTABSetup(solver_, A_, B_, X_);

  // Solve the linear system

  HYPRE_SStructBiCGSTABSolve(solver_, A_, B_, X_);

  // Write out some diagnostic info about the solve

  HYPRE_SStructBiCGSTABGetNumIterations(solver_, &iter_);
  HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm(solver_, &resid_);

  if (debug) printf ("HYPRE_SStructBiCGSTABSolve3 num iterations: %d\n",iter_);
  if (debug) printf ("HYPRE_SStructBiCGSTABSolve3 final relative residual norm: %g\n",resid_);


  // Delete the solver

  HYPRE_SStructBiCGSTABDestroy(solver_);
  solver_ = 0;

} // Hypre::solve_bicgstab_()

