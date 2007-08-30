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
#include "hierarchy.hpp"
#include "sphere.hpp"
#include "domain.hpp"
#include "problem.hpp"
#include "hypre.hpp"

const int debug         = 0;

const int dump_matrix   = 0;
const int dump_vector   = 0;
const int dump_solution = 0;

const int trace         = 0;

const int DIRICHLET     = 0;

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

void Hypre::init_hierarchy (Hierarchy & hierarchy, 
			    Mpi       & mpi)

{

  Grid::set_mpi (mpi);
  
  int dim    = hierarchy.dimension();
  int levels = hierarchy.num_levels();

  // Create the hypre grid
  
  // *******************************************************************
  HYPRE_SStructGridCreate (MPI_COMM_WORLD, dim, levels, &grid_);
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

    // Determine periodicity of root Level

    for (int i=0; i<3; i++) {
      // periodicity of root level
      periodicity[i] = hierarchy.level(0).zones(i);
      // adjust for periodicity of given level
      for (int k=0; k < part; k++) {
	periodicity[i] *= r;
      }
    }

    if (DIRICHLET) {
      periodicity[0] = 0;
      periodicity[1] = 0;
      periodicity[2] = 0;
    }

    if (debug) printf ("%s:%d Periodicity = (%d,%d,%d)\n",__FILE__,__LINE__,
		       periodicity[0],periodicity[1],periodicity[2]);

    // *******************************************************************
    HYPRE_SStructGridSetPeriodic (grid_, 0, periodicity);
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

    Setting up the matrix nonzero structure is done with the following
    steps:

     - Define the stencil for all interior grid elements

     - Handle coarse unknowns adjacent to fine unknowns in child
     - Handle fine unknowns adjacent to coarse unknowns in parent

     - Handle coarse unknowns adjacent to fine unknowns in neighbor's child
     - Handle fine unknowns adjacent to coarse unknowns in parent's neighbor
     
    The matrix nonzero structure is generally nonsymmetric.

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

    // Define stencil connections within each part

    // *******************************************************************
    HYPRE_SStructGraphSetStencil (graph_, part, 0, stencil_);
    // *******************************************************************

    ItLevelGridsLocal itg (*level);

    while (Grid * grid = itg++) {

      // Define connections for unknowns adjacent to fine unknowns in children

      init_graph_children_(*grid);

      // Define connections for unknowns adjacent to coarse unknowns in parent

      init_graph_parent_(hierarchy,*grid);

      // Define connections for unknowns adjacent to fine unknowns in neighbor's child

      init_graph_neighbors_children_(*grid);

      // Define connections for unknowns adjacent to coarse unknowns
      // in parent's neighbor

      init_graph_parents_neighbor_(*grid);

    }
    ++ part;
  }

  // Assemble the graph

  // *******************************************************************
  HYPRE_SStructGraphAssemble (graph_);
  // *******************************************************************

}

//----------------------------------------------------------------------

/// Initialize the matrix A

/** Creates a matrix with a given non-zero structure, and sets nonzero
    values.

    Setting up the matrix elements is done with the following
    steps:

     - Set stencil values
     - Clear stencil values in overlapped grids

     - Set coarse unknowns adjacent to fine unknowns in child
     - Set fine unknowns adjacent to coarse unknowns in parent

     - Set coarse unknowns adjacent to fine unknowns in neighbor's child
     - Set fine unknowns adjacent to coarse unknowns in parent's neighbor
     
    The matrix is generally nonsymmetric.

*/

void Hypre::init_matrix (Hierarchy & hierarchy)

{
  // Create the hypre matrix object

  // *******************************************************************
  HYPRE_SStructMatrixCreate (MPI_COMM_WORLD, graph_, &A_);
  // *******************************************************************

  // Set the matrix type

  // *******************************************************************
  HYPRE_SStructMatrixSetObjectType (A_,HYPRE_SSTRUCT);
  // *******************************************************************

  // Initialize the hypre matrix object

  // *******************************************************************
  HYPRE_SStructMatrixInitialize (A_);
  // *******************************************************************

  ItHierarchyLevels itl (hierarchy);

  int part = 0;

  while (Level * level = itl++) {

    ItLevelGridsLocal itg (*level);

    while (Grid * grid = itg++) {

      // Set matrix stencil values for each grid

      init_matrix_stencil_(*grid);

      // Clear matrix stencil values under child grids

      init_matrix_clear_(*grid);

      // Set matrix values for unknowns adjacent to fine unknowns in children

      init_matrix_children_(*grid);

      // Set matrix values for unknowns adjacent to coarse unknowns in parent

      init_matrix_parent_(hierarchy,*grid);

      // Set matrix values for unknowns adjacent to fine unknowns in neighbor's child

      init_matrix_neighbors_children_(*grid);

      // Set matrix values for unknowns adjacent to coarse unknowns in
      // parent's neighbor

      init_matrix_parents_neighbor_(*grid);

    }
    ++ part;
  }

  // Assemble the matrix

  // *******************************************************************
  HYPRE_SStructMatrixAssemble (A_);
  // *******************************************************************

  // Write the matrix to a file for debugging

  if (dump_matrix) {

    // *******************************************************************
    HYPRE_SStructMatrixPrint ("A",A_,1);
    // *******************************************************************

  }

}

//----------------------------------------------------------------------

/// Initialize the right-hand-side vector b

void Hypre::init_vectors (Hierarchy & hierarchy,
			  std::vector<Point *> points,
			  std::vector<Sphere *> spheres)

{
  // Create the hypre solution x and right-hand side b vector object

  // *******************************************************************
  HYPRE_SStructVectorCreate (MPI_COMM_WORLD, grid_, &B_);
  // *******************************************************************

  // *******************************************************************
  HYPRE_SStructVectorCreate (MPI_COMM_WORLD, grid_, &X_);
  // *******************************************************************

  // Initialize the hypre vector objects

  // *******************************************************************
  HYPRE_SStructVectorInitialize (B_);
  // *******************************************************************

  // *******************************************************************
  HYPRE_SStructVectorInitialize (X_);
  // *******************************************************************


  // Clear vectors

  if (debug) printf ("%s:%d Clear hypre vectors here: assumed done at initialization\n",
		     __FILE__,__LINE__);

  // Initialize B_ according to density

  //    scaling0 = scaling factor for root-level, assuming matrix
  //    coefficients [-1,-1,-1,6,-1,-1,-1]


  Scalar G  = Constants::G();
  Scalar pi = Constants::pi();
  Scalar scaling0 = -4.0*G*pi;

  Scalar shift_b_sum = 0.0;

  int i;
  for (i=0; i<points.size(); i++) {
    Point & point      = *points[i];
    Grid & grid        = hierarchy.grid(point.igrid());
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
      printf ("Index of point in grid (%d,%d,%d)\n",index[0],index[1],index[2]);
    }

    int part = grid.level();
    int var = 0;
    
    shift_b_sum += value;

    // *******************************************************************
    HYPRE_SStructVectorAddToValues (B_, part, index, var, &value);
    // *******************************************************************

   
  }

  if (spheres.size() > 0) {  
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    NOT_IMPLEMENTED("Contribution of sphere mass to right-hand side");
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  }

  // Shift B to zero out the null space if problem is periodic

  if (! DIRICHLET) {

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

    // Perform the shift
    
    part = 0;
    while (Level * level = itl++) {
      ItLevelGridsLocal itg (*level);
      while (Grid * grid = itg++) {
	int lower[3] = { grid->i_lower(0), grid->i_lower(1), grid->i_lower(2) };
	int upper[3] = { grid->i_upper(0), grid->i_upper(1), grid->i_upper(2) };
	Scalar * values = new Scalar[grid->num_unknowns()];
	for (i=0; i<grid->num_unknowns(); i++) values[i] = shift_b_amount;
	HYPRE_SStructVectorAddToBoxValues (B_,part,lower,upper,0,values);
	delete [] values;
      }
      ++part;
    }
    if (debug) printf ("%s:%d shift (count,sum,amount) = (%lld,%g,%g)\n",
		       __FILE__,__LINE__,shift_b_count,shift_b_sum,shift_b_amount);
  
  }
  // Assemble the vectors

  // *******************************************************************
  HYPRE_SStructVectorAssemble (B_);
  HYPRE_SStructVectorAssemble (X_);
  // *******************************************************************

  // Write the vector to a file for debugging

  if (dump_vector) {

    // *******************************************************************
    HYPRE_SStructVectorPrint ("B",B_,1);
    // *******************************************************************

  }
}

//----------------------------------------------------------------------

/// Initialize and solve the linear solver

void Hypre::solve (Hierarchy & hierarchy)

{
  if (hierarchy.num_levels() > 1) {
    solve_fac_(hierarchy);
  } else {
    solve_pfmg_(hierarchy);
  }
  if (dump_vector || dump_solution) {
    // *******************************************************************
    HYPRE_SStructVectorPrint ("X",X_,1);
    // *******************************************************************
  }

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

/// Handle coarse unknowns adjacent to fine unknowns in child

void Hypre::init_graph_children_ (Grid & grid)

{

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  NOT_IMPLEMENTED("init_graph_children_()");
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ItGridChildren itchild (grid);
  while (Grid * child = itchild++) {
    _TRACE_;
    // *******************************************************************
    //	HYPRE_SStructGraphAddEntries (grid->hypre_graph(),
    //				      grid->level(),
    //				      grid->INDEX
    //				      neighbor.level(),
    //				      neighbor.INDEX,
    //				      variable);
    // *******************************************************************
  }
}

//------------------------------------------------------------------------

/// Handle fine unknowns adjacent to coarse unknowns in parent

void Hypre::init_graph_parent_ (Hierarchy & hierarchy, Grid & grid)

{

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  NOT_IMPLEMENTED("Hypre::init_graph_parent_()");
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Grid * parent = hierarchy.parent(grid);
  _TRACE_;
  // *******************************************************************
  //	HYPRE_SStructGraphAddEntries (grid->hypre_graph(),
  //				      grid->level(),
  //				      grid->INDEX
  //				      neighbor.level(),
  //				      neighbor.INDEX,
  //				      variable);
  // *******************************************************************
}

//------------------------------------------------------------------------

/// Handle coarse unknowns adjacent to fine unknowns in neighbor's child

void Hypre::init_graph_neighbors_children_ (Grid & grid)

{
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  NOT_IMPLEMENTED("Hypre::init_graph_neighbors_children_()");
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  // *******************************************************************
  //	HYPRE_SStructGraphAddEntries (grid->hypre_graph(),
  //				      grid->level(),
  //				      grid->INDEX
  //				      neighbor.level(),
  //				      neighbor.INDEX,
  //				      variable);
  // *******************************************************************
}

//------------------------------------------------------------------------

/// Handle fine unknowns adjacent to coarse unknowns in parent's neighbor

void Hypre::init_graph_parents_neighbor_ (Grid & grid)

{
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  NOT_IMPLEMENTED("Hypre::init_graph_parents_neighbor_()");
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  // *******************************************************************
  //	HYPRE_SStructGraphAddEntries (grid->hypre_graph(),
  //				      grid->level(),
  //				      grid->INDEX
  //				      neighbor.level(),
  //				      neighbor.INDEX,
  //				      variable);
  // *******************************************************************
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

/// Clear matrix stencil values under child grids

void Hypre::init_matrix_clear_ (Grid & grid)

{
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  NOT_IMPLEMENTED("Hypre::init_matrix_clear_()");
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ItGridChildren itchildren (grid);
  
  while (Grid * child = itchildren++) {
    _TRACE_;
  }

  // *******************************************************************
  //	HYPRE_SStructMatrixAddEntries (grid->hypre_matrix(),
  //				      grid->level(),
  //				      grid->INDEX
  //				      neighbor.level(),
  //				      neighbor.INDEX,
  //				      variable);
  // *******************************************************************

}

//------------------------------------------------------------------------

/// Define connections for unknowns adjacent to fine unknowns in children

void Hypre::init_matrix_children_ (Grid & grid)

{

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  NOT_IMPLEMENTED("Hypre::init_matrix_children_()");
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ItGridChildren itchild (grid);
  while (Grid * child = itchild++) {

    _TRACE_;

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

//------------------------------------------------------------------------

/// Set matrix values for unknowns adjacent to coarse unknowns in parent

void Hypre::init_matrix_parent_ (Hierarchy & hierarchy, Grid & grid)

{

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  NOT_IMPLEMENTED("Hypre::init_matrix_parent_()");
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Grid * parent = hierarchy.parent(grid);
  _TRACE_;
  // *******************************************************************
  //	HYPRE_SStructMatrixAddEntries (grid->hypre_matrix(),
  //				      grid->level(),
  //				      grid->INDEX
  //				      neighbor.level(),
  //				      neighbor.INDEX,
  //				      variable);
  // *******************************************************************
}

//------------------------------------------------------------------------

/// Set matrix values for unknowns adjacent to fine unknowns in neighbor's child

void Hypre::init_matrix_neighbors_children_ (Grid & grid)

{

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  NOT_IMPLEMENTED("Hypre::init_matrix_neighbors_children_()");
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  // *******************************************************************
  //	HYPRE_SStructMatrixAddEntries (grid->hypre_matrix(),
  //				      grid->level(),
  //				      grid->INDEX
  //				      neighbor.level(),
  //				      neighbor.INDEX,
  //				      variable);
  // *******************************************************************
}

//------------------------------------------------------------------------

/// Set matrix values for unknowns adjacent to coarse unknowns in
/// parent's neighbor

void Hypre::init_matrix_parents_neighbor_ (Grid & grid)

{

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  NOT_IMPLEMENTED("Hypre::init_matrix_parents_neighbor_()");
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  // *******************************************************************
  //	HYPRE_SStructMatrixAddEntries (grid->hypre_matrix(),
  //				      grid->level(),
  //				      grid->INDEX
  //				      neighbor.level(),
  //				      neighbor.INDEX,
  //				      variable);
  // *******************************************************************
}

//------------------------------------------------------------------------

/// Initialize the PFMG hypre solver

void Hypre::solve_pfmg_ (Hierarchy & hierarchy)

{
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  NOT_IMPLEMENTED("Hypre::solve_pfmg_()");
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  HYPRE_SStructSysPFMGCreate    (MPI_COMM_WORLD, &solver_);
  HYPRE_SStructSysPFMGSetLogging(solver_, 1);
  HYPRE_SStructSysPFMGSetup     (solver_,A_,B_,X_);
  HYPRE_SStructSysPFMGSolve     (solver_,A_,B_,X_);

  int num_iterations;
  HYPRE_SStructSysPFMGGetNumIterations (solver_,&num_iterations);
  if (debug) printf ("HYPRE_SStructSysPFMGSolve num iterations: %d\n",num_iterations);

  double residual;
  HYPRE_SStructSysPFMGGetFinalRelativeResidualNorm (solver_,&residual);
  if (debug) printf ("HYPRE_SStructSysPFMGSolve final relative residual norm: %g\n",residual);

  HYPRE_SStructSysPFMGDestroy (solver_);
}

//------------------------------------------------------------------------

/// Initialize the FAC hypre solver

void Hypre::solve_fac_ (Hierarchy & hierarchy)

{
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  NOT_IMPLEMENTED("Hypre::init_solver_fac_()");
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
}

