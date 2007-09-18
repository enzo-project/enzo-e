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

#include "hypre-solve.hpp"

#include "mpi.hpp"
#include "scalar.hpp"
#include "constants.hpp"
#include "point.hpp"
#include "faces.hpp"
#include "grid.hpp"
#include "level.hpp"
#include "domain.hpp"
#include "hierarchy.hpp"
#include "sphere.hpp"
#include "parameters.hpp"
#include "problem.hpp"
#include "hypre.hpp"

const int debug  = 0;
const int trace  = 0;

//======================================================================
// PUBLIC MEMBER FUNCTIONS
//======================================================================

/// Hypre constructor

Hypre::Hypre ()
  : grid_(0),
    graph_(0),
    stencil_(0),
    A_(0),
    B_(0),
    X_(0),
    solver_(0)
{
  
}

//----------------------------------------------------------------------

/// Initialize the Grid Hierarchy

/** Creates a hypre grid, with one part per level and one box per Grid
    patch object, for an AMR problem.  Sets grid box extents, grid
    part variables, and periodicity of the root-level grid part. */

void Hypre::init_hierarchy (Parameters & parameters,
			    Hierarchy  & hierarchy, 
			    Mpi        & mpi)

{

  Grid::set_mpi (mpi);
  
  int dim       = hierarchy.dimension();
  int num_parts = hierarchy.num_levels();

  // Create the hypre grid
  
  // *******************************************************************
  HYPRE_SStructGridCreate (MPI_COMM_WORLD, dim, num_parts, &grid_);
  // *******************************************************************

  ItHierarchyLevels itl (hierarchy);

  int part = 0;

  while (Level * level = itl++) {

    ItLevelGridsLocal itg (*level);

    while (Grid * grid = itg++) {

      int lower[3] = {grid->i_lower(0),grid->i_lower(1),grid->i_lower(2)};
      int upper[3] = {grid->i_upper(0),grid->i_upper(1),grid->i_upper(2)};

      // Set extents for boxes that comprise the hypre grid

      // *******************************************************************
      HYPRE_SStructGridSetExtents(grid_, part, lower, upper);
      // *******************************************************************
      
    }

    // Create a single cell-centered variable for each grid part (level)

    HYPRE_SStructVariable variable_types[] = { HYPRE_SSTRUCT_VARIABLE_CELL };
    const int numvars = 1;

    // *******************************************************************
    HYPRE_SStructGridSetVariables(grid_, part, numvars, variable_types);
    // *******************************************************************

    // Set grid part to be periodic, with periodicity determined by the root
    // level size, current level, and refinement factor (assumed to be 2)

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

    // *******************************************************************
    HYPRE_SStructGridSetPeriodic (grid_, part, periodicity);
    // *******************************************************************

    ++ part;
  }

  // When finished, assemble the hypre grid

  // *******************************************************************
  HYPRE_SStructGridAssemble (grid_);
  // *******************************************************************
  
}

//----------------------------------------------------------------------

/// Initialize the discretization stencils.  

/** Creates and initializes a stencil object.  Supports 1, 2, or 3
    dimensional stencils. */

void Hypre::init_stencil (Hierarchy & hierarchy)

{

  int dim = hierarchy.dimension();

  // *******************************************************************
  HYPRE_SStructStencilCreate (dim,dim*2+1,&stencil_);
  // *******************************************************************

  int entries[][3] = { {  0, 0, 0 },
		       {  1, 0, 0 },
		       { -1, 0, 0 },
		       {  0, 1, 0 },
		       {  0,-1, 0 },
		       {  0, 0, 1 },
		       {  0, 0,-1 } };

  // *******************************************************************
  HYPRE_SStructStencilSetEntry (stencil_, 0, entries[0], 0);
  // *******************************************************************

  if (dim >= 1) {
    // *******************************************************************
    HYPRE_SStructStencilSetEntry (stencil_, 1, entries[1], 0);
    HYPRE_SStructStencilSetEntry (stencil_, 2, entries[2], 0);
    // *******************************************************************
  }
  if (dim >= 2) {
    // *******************************************************************
    HYPRE_SStructStencilSetEntry (stencil_, 3, entries[3], 0);
    HYPRE_SStructStencilSetEntry (stencil_, 4, entries[4], 0);
    // *******************************************************************
  }
  if (dim >= 3) {
    // *******************************************************************
    HYPRE_SStructStencilSetEntry (stencil_, 5, entries[5], 0);
    HYPRE_SStructStencilSetEntry (stencil_, 6, entries[6], 0);
    // *******************************************************************
  }

}

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

  // *******************************************************************
  HYPRE_SStructGraphCreate (MPI_COMM_WORLD, grid_, &graph_);
  // *******************************************************************

  ItHierarchyLevels itl (hierarchy);

  int part = 0;

  while (Level * level = itl++) {

    // 1. Define stencil connections within each level

    // *******************************************************************
    HYPRE_SStructGraphSetStencil (graph_, part, 0, stencil_);
    // *******************************************************************

    // WARNING: POSSIBLE SCALING ISSUE.
    //
    // Below we loop over all grids; however, we only need
    // too loop over parent-child pairs such either child or parent is
    // local to this MPI process.  Could be done with two loops,
    // looping over local parents, then local children, but need to be
    // careful not to add the same entries twice.  Could process local
    // parents first since they're unique, and only process local
    // children if their parent is remote.

    ItLevelGridsAll itag (*level);

    while (Grid * grid = itag++) {

      // 2. Define connections between grid patches in adjacent parts

      init_graph_nonstencil_(*grid);

    }
    ++ part;
  }

  // Assemble the graph

  // *******************************************************************
  HYPRE_SStructGraphAssemble (graph_);
  // *******************************************************************

}


//----------------------------------------------------------------------

/// Initialize the right-hand-side vector b

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

  // *******************************************************************
  HYPRE_SStructMatrixCreate (MPI_COMM_WORLD, graph_, &A_);
  HYPRE_SStructVectorCreate (MPI_COMM_WORLD, grid_, &X_);
  HYPRE_SStructVectorCreate (MPI_COMM_WORLD, grid_, &B_);
  // *******************************************************************

  // Set the matrix type

  // *******************************************************************
  HYPRE_SStructMatrixSetObjectType (A_,HYPRE_SSTRUCT);
  // *******************************************************************

  // Initialize the hypre matrix and vector objects

  // *******************************************************************
  HYPRE_SStructMatrixInitialize (A_);
  HYPRE_SStructVectorInitialize (X_);
  HYPRE_SStructVectorInitialize (B_);
  // *******************************************************************
 
  ItHierarchyLevels itl (hierarchy);

  int num_parts = hierarchy.num_levels();

  int part = 0;
  while (Level * level = itl++) {

    ItLevelGridsLocal itlg (*level);

    while (Grid * grid = itlg++) {

      // 1. Set stencil values

      init_matrix_stencil_(*grid);
    }

    // 2. Clean up stencil connections between parts

    init_matrix_clear_(*level);

    // WARNING: POSSIBLE SCALING ISSUE.
    //
    // Below we loop over all grids; however, we only need
    // too loop over parent-child pairs such either child or parent is
    // local to this MPI process.  Could be done with two loops,
    // looping over local parents, then local children, but need to be
    // careful not to add the same entries twice.  Could process local
    // parents first since they're unique, and only process local
    // children if their parent is remote.
 
    ItLevelGridsAll itag (*level);

    while (Grid * grid = itag++) {

      // 3. Set values for unknowns between parent and children

      init_matrix_nonstencil_(*grid);

    }
    ++ part;
  }

  // Initialize B_ according to density

  //    scaling0 = scaling factor for root-level, assuming matrix
  //    coefficients [-1,-1,-1,6,-1,-1,-1]


  Scalar local_shift_b_sum = 0.0;

  local_shift_b_sum += init_vector_points_  (hierarchy,points);
  local_shift_b_sum += init_vector_spheres_ (hierarchy,spheres);

  Scalar shift_b_sum = 0.0;

  MPI_Allreduce (&local_shift_b_sum, &shift_b_sum, 1, MPI_SCALAR, MPI_SUM, MPI_COMM_WORLD);

  if (debug) printf ("b_sum (local,global) = (%g,%g)\n",local_shift_b_sum,shift_b_sum);


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
      }
    }

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
	delete [] values;
      }
      ++part;
    }
    if (debug) printf ("%s:%d shift (count,sum,amount) = (%lld,%g,%g)\n",
		       __FILE__,__LINE__,shift_b_count,shift_b_sum,shift_b_amount);
  
  }

  // Assemble the matrix vectors

  // *******************************************************************
  HYPRE_SStructMatrixAssemble (A_);
  HYPRE_SStructVectorAssemble (B_);
  HYPRE_SStructVectorAssemble (X_);
  // *******************************************************************

  // Write the vector to a file for debugging

  // *******************************************************************
  if (parameters.value("dump_a") == "true") HYPRE_SStructMatrixPrint ("A",A_,1);
  if (parameters.value("dump_b") == "true") HYPRE_SStructVectorPrint ("B",B_,1);  
  // *******************************************************************

}

//----------------------------------------------------------------------

/// Initialize and solve the linear solver

void Hypre::solve (Parameters & parameters,
		   Hierarchy & hierarchy)

{
  if (hierarchy.num_levels() > 1) {
    solve_fac_(hierarchy);
  } else {
    solve_pfmg_(hierarchy);
  }
  
  // *******************************************************************
  if (parameters.value("dump_x") == "true") HYPRE_SStructVectorPrint ("X",X_,1);
  // *******************************************************************

}

//----------------------------------------------------------------------

/// Evaluate the success of the solve

void Hypre::evaluate (Hierarchy & hierarchy)

{
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  NOT_IMPLEMENTED("Hypre::evaluate()");
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
}

//======================================================================
// PRIVATE MEMBER FUNCTIONS
//======================================================================

/// Add nonzeros connecting grid with all children grid patches

void Hypre::init_graph_nonstencil_ (Grid & grid)

{

  NOT_IMPLEMENTED("Hypre::init_graph_coarse_()");

  ItGridChildren itchild (grid);

  int ilower_grid[3];

  grid.i_lower(ilower_grid);
  
  printf ("DEBUG %s:%d  grid i_lower (%d %d %d)\n",
	  __FILE__,__LINE__,
	  ilower_grid[0],ilower_grid[1],ilower_grid[2]);

  while (Grid * child = itchild++) {

    _TRACE_;

    if (grid.is_local() || child->is_local()) {
      int ilower_child[3];

      child->i_lower(ilower_child);

      printf ("DEBUG %s:%d  child i_lower (%d %d %d)\n",
	      __FILE__,__LINE__,
	      ilower_child[0],ilower_child[1],ilower_child[2]);

      // Add the matrix entries between fine and coarse grids iff
      // one of the grids is MPI-local

      // *******************************************************************
      //      HYPRE_SStructGraphAddEntries (graph_,
      //				    grid->level(),
      //				    grid->INDEX
      //				    neighbor.level(),
      //				    neighbor.INDEX,
      //				    0);
      // *******************************************************************
    }
  }
}

//------------------------------------------------------------------------

/// Set matrix stencil values for the grid

void Hypre::init_matrix_stencil_ (Grid & grid)

{
  _TRACE_;

  int part         = grid.level();
  int lower[3]     = { grid.i_lower(0), grid.i_lower(1), grid.i_lower(2) };
  int upper[3]     = { grid.i_upper(0), grid.i_upper(1), grid.i_upper(2) };
  int var          = 0;
  int count        = grid.num_unknowns();
  int entries[7]   = { 0,1,2,3,4,5,6 };

  double * values0  = new double [count];
  double * values1  = new double [count];

  for (int i=0; i<count; i++) {
    values0[i] = 6;
    values1[i] = -1;
  }

  // *******************************************************************
  HYPRE_SStructMatrixSetBoxValues (A_,part,lower,upper,
				   var,1,&entries[0],values0);
  // *******************************************************************
  for (int stencil = 1; stencil < 7; stencil ++) {
    // *******************************************************************
    HYPRE_SStructMatrixSetBoxValues (A_,part,lower,upper,
				     var, 1, &entries[stencil], values1);
    // *******************************************************************
  }

 delete [] values1;
 delete [] values0;

}

//------------------------------------------------------------------------

/// Clean up stencil connections between parts

void Hypre::init_matrix_clear_ (Level & level)
{
  // WARNING: hard-coding refinement factor of 2
  int r_factors[3] = {2,2,2}; 
  int part = level.index();
  if (part > 0) {

    // Clear stencil values from coarse to fine part

    HYPRE_SStructFACZeroCFSten (A_,grid_, part, r_factors);

    // Clear stencil values from fine to coarse part

    HYPRE_SStructFACZeroFCSten (A_,grid_, part);

    // Set overlapped areas of part with identity

    HYPRE_SStructFACZeroAMRMatrixData (A_, part-1, r_factors);
    // Need to clear under rhs also
    //   HYPRE_SStructFACZeroAMRVectorData(b, plevels, prefinements);
  }
}

//------------------------------------------------------------------------

/// Define matrix elements connecting parent and children grid patches

void Hypre::init_matrix_nonstencil_ (Grid & grid)

{

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  NOT_IMPLEMENTED("Hypre::init_matrix_fine_()");
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ItGridChildren itchild (grid);

  int ilower_grid[3];

  grid.i_lower(ilower_grid);
  
  printf ("DEBUG %s:%d  grid i_lower (%d %d %d)\n",
	  __FILE__,__LINE__,
	  ilower_grid[0],ilower_grid[1],ilower_grid[2]);

  while (Grid * child = itchild++) {

    _TRACE_;

    if (grid.is_local() || child->is_local()) {

      int ilower_child[3];

      child->i_lower(ilower_child);

      printf ("DEBUG %s:%d  child i_lower (%d %d %d)\n",
	      __FILE__,__LINE__,
	      ilower_child[0],ilower_child[1],ilower_child[2]);
      // *******************************************************************
      //	HYPRE_SStructMatrixAddEntries (grid->hypre_matrix(),
      //				      grid->level(),
      //				      grid->INDEX
      //				      neighbor.level(),
      //				      neighbor.INDEX,
      //				      variable);
      // *******************************************************************
    }
  }
}

//------------------------------------------------------------------------

/// Add contributions from point sources to right-hand side B

Scalar Hypre::init_vector_points_ (Hierarchy            & hierarchy,
				   std::vector<Point *> & points)

{

  Scalar scaling0 = -4.0*Constants::G()*Constants::pi();

  Scalar shift_b_sum = 0.0;

  int i;
  for (i=0; i<points.size(); i++) {
    Point & point      = *points[i];
    Grid & grid        = hierarchy.grid(point.igrid());
    if (grid.is_local()) {

      Scalar cell_volume = grid.h(0) * grid.h(1) * grid.h(2);
      Scalar density     = point.mass() / cell_volume;
      Scalar value       = scaling0 * density;

      // Add contribution of the point to the right-hand side vector

      int index[3];
      for (int k=0; k<3; k++) {
	Scalar ap = point.x(k)      - grid.x_lower(k);
	Scalar ag = grid.x_upper(k) - grid.x_lower(k);
	int    ig = grid.num_unknowns(k);
	int    i0 = grid.i_lower(k);
	index[k] = int (ap/ag*ig) + i0;
     
      }
      if (debug) {
	point.print();
	grid.print();
	printf ("Point index  = %d %d %d)\n",index[0],index[1],index[2]);
	printf ("Cell size    = %g %g %g\n",grid.h(0),grid.h(1),grid.h(2));
	printf ("Cell volume  = %g\n",cell_volume);
	printf ("Cell density = %g\n",density);
	printf ("RHS contribution = %g\n",value);
      }

      int part = grid.level();
      int var = 0;
    
      shift_b_sum += value;

      // *******************************************************************
      HYPRE_SStructVectorAddToValues (B_, part, index, var, &value);
      // *******************************************************************
    }
  }
  return shift_b_sum;
}

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
}

//------------------------------------------------------------------------

/// Initialize the PFMG hypre solver

void Hypre::solve_pfmg_ (Hierarchy & hierarchy)

{

  // Create and initialize the solver

  HYPRE_SStructSysPFMGCreate    (MPI_COMM_WORLD, &solver_);
  HYPRE_SStructSysPFMGSetLogging(solver_, 1);
  HYPRE_SStructSysPFMGSetup     (solver_,A_,B_,X_);

  // Solve the linear system

  HYPRE_SStructSysPFMGSolve     (solver_,A_,B_,X_);

  // Write out some diagnostic info about the solve

  int num_iterations;
  HYPRE_SStructSysPFMGGetNumIterations (solver_,&num_iterations);
  if (debug) printf ("HYPRE_SStructSysPFMGSolve num iterations: %d\n",num_iterations);

  double residual;
  HYPRE_SStructSysPFMGGetFinalRelativeResidualNorm (solver_,&residual);
  if (debug) printf ("HYPRE_SStructSysPFMGSolve final relative residual norm: %g\n",residual);

  // Delete the solver

  HYPRE_SStructSysPFMGDestroy (solver_);
  solver_ = 0;

  
}

//------------------------------------------------------------------------

/// Initialize the FAC hypre solver

void Hypre::solve_fac_ (Hierarchy & hierarchy)

{
  int i;

  // Create the solver

  HYPRE_SStructFACCreate(MPI_COMM_WORLD, &solver_);

  // Initialize parts

  int num_parts = hierarchy.num_levels();
  int *parts  = new int [num_parts];
  for (i=0; i<num_parts; i++) parts[i] = i;

  HYPRE_SStructFACSetMaxLevels(solver_,  num_parts);
  HYPRE_SStructFACSetPLevels(solver_, num_parts, parts);

  // Initialize refinement factors

  typedef int int3[3];
  int3 *refinements = new int3 [num_parts];
  
  for (i=0; i<num_parts; i++) {
    refinements[i][0] = 2;
    refinements[i][1] = 2;
    refinements[i][2] = 2;
  }

  HYPRE_SStructFACSetPRefinements(solver_, num_parts, refinements);

  // solver parameters

  int npre   = 2;
  int npost  = 2;
  int csolve = 2;
  int relax  = 2;

  HYPRE_SStructFACSetNumPreRelax(solver_,      npre);
  HYPRE_SStructFACSetNumPostRelax(solver_,     npost);
  HYPRE_SStructFACSetCoarseSolverType(solver_, csolve);
  HYPRE_SStructFACSetRelaxType(solver_,        relax);

  // stopping criteria

  int itmax   = 20;
  double rtol = 1e-6;

  HYPRE_SStructFACSetRelChange(solver_, 0);
  HYPRE_SStructFACSetMaxIter(solver_,    itmax);
  HYPRE_SStructFACSetTol(solver_,        rtol);

  // output amount

  HYPRE_SStructFACSetLogging(solver_, 1);

  // prepare for solve

  HYPRE_SStructFACSetup2(solver_, A_, B_, X_);

  // Solve the linear system

  HYPRE_SStructFACSolve3(solver_, A_, B_, X_);

  // Write out some diagnostic info about the solve

  int num_iterations;
  HYPRE_SStructFACGetNumIterations(solver_, &num_iterations);
  double residual;
  HYPRE_SStructFACGetFinalRelativeResidualNorm(solver_, &residual);

  // Delete the solver

  HYPRE_SStructFACDestroy2(solver_);
  solver_ = 0;

  // Delete local dynamic storage
  delete [] parts;
  delete [] refinements;
}

