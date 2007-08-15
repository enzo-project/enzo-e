
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

#include "mpi.hpp"
#include "scalar.hpp"
#include "point.hpp"
#include "discret.hpp"
#include "grid.hpp"
#include "level.hpp"
#include "hierarchy.hpp"
#include "sphere.hpp"
#include "domain.hpp"
#include "problem.hpp"
#include "hypre.hpp"

const int debug_hypre           = 0;

const int debug_hypre_hierarchy = 0;
const int debug_hypre_graph     = 0;

const int debug_discret = 0;

//======================================================================

/// Hypre constructor

Hypre::Hypre ()

{
  
}

//======================================================================

/// Initialize the Grid Hierarchy

/** Creates a hierarchy of hypre grids for an AMR problem.  */

void Hypre::init_hierarchy (Hierarchy & H, Mpi &mpi)

{

  Grid::set_mpi (mpi);

  int dim    = H.dimension();
  int levels = H.num_levels();

  for (int i=0; i<H.num_grids(); i++) {
    if (H.grid(i).is_local()) {

      // Create hypre grids

      init_hierarchy_create_grid_(H.grid(i),dim,levels);

      // Set hypre grid extents

      init_hierarchy_set_grid_extents_(H.grid(i));

      // Set the hypre grid variables

      init_hierarchy_set_grid_variables_(H.grid(i));

    }
  }

  // All grids must exist before assembling

  mpi.barrier(); // MAY BE UNNECESSARY?

  for (int i=0; i<H.num_grids(); i++) {

    if (H.grid(i).is_local()) {

      // Assemble grids

      init_hierarchy_assemble_grids_(H.grid(i));
    }
  }
  
}

//------------------------------------------------------------------------

/// Create hypre Grids.

void Hypre::init_hierarchy_create_grid_(
					Grid & grid,
					int dim,
					int levels
					)

{

  // Create the HYPRE grid structure

  if (debug_hypre_hierarchy) {
    printf ("DEBUG_HYPRE_HIERARCHY HYPRE_SStructGridCreate (\n");
    printf ("DEBUG_HYPRE_HIERARCHY      dim    = %d\n", dim);
    printf ("DEBUG_HYPRE_HIERARCHY      levels = %d\n", levels);
  }

  HYPRE_SStructGridCreate (MPI_COMM_WORLD, 
			   dim, 
			   levels, 
			   &grid.hypre_grid());

}

//------------------------------------------------------------------------

/// Set hypre Grid extents

void Hypre::init_hierarchy_set_grid_extents_(Grid & grid)

{

  // Set HYPRE grid extents

  int lower_extents[3] = { 0, 0, 0};
  int upper_extents[3] = {grid.num_unknowns(0)-1,
			  grid.num_unknowns(1)-1,
			  grid.num_unknowns(2)-1};
  int part = grid.level();

  if (debug_hypre_hierarchy) {
    printf ("DEBUG_HYPRE_HIERARCHY HYPRE_SStructGridSetExtents\n");
    printf ("DEBUG_HYPRE_HIERARCHY          part = %d\n",part);
    printf ("DEBUG_HYPRE_HIERARCHY lower_extents = %d %d %d\n",
	    lower_extents[0],lower_extents[1],lower_extents[2]);
    printf ("DEBUG_HYPRE_HIERARCHY upper_extents = %d %d %d\n",
	    upper_extents[0],upper_extents[1],upper_extents[2]);
  }

  HYPRE_SStructGridSetExtents(grid.hypre_grid(),
			      part,
			      lower_extents,
			      upper_extents);
}

//------------------------------------------------------------------------

/// Set hypre Grid variable type

void Hypre::init_hierarchy_set_grid_variables_(Grid & grid)

{

  HYPRE_SStructVariable variable_types[] = { HYPRE_SSTRUCT_VARIABLE_CELL };

  if (debug_hypre_hierarchy) {
    printf ("DEBUG_HYPRE_HIERARCHY HYPRE_SStructGridSetVariables(\n");
    printf ("DEBUG_HYPRE_HIERARCHY         grid.id() = %d\n",grid.id());
    printf ("DEBUG_HYPRE_HIERARCHY                     1 \n");
    printf ("DEBUG_HYPRE_HIERARCHY    variable_types = %d\n",variable_types[0]);
  }

  HYPRE_SStructGridSetVariables(grid.hypre_grid(),
				grid.level(),
				1,
				variable_types);
}

//------------------------------------------------------------------------

/// Assemble the hypre Grid hierarchy.

void Hypre::init_hierarchy_assemble_grids_(Grid & grid)
{

  if (debug_hypre_hierarchy) {
    printf ("DEBUG_HYPRE_HIERARCHY HYPRE_SStructGridAssemble()\n",
	    __FILE__,__LINE__);
  }

  HYPRE_SStructGridAssemble (grid.hypre_grid());
}


//======================================================================

/// Initialize the discretization stencils.  

/** Creates a stencil object and initializes entries.  Supports 1, 2, 
    or 3 dimensions. */

void Hypre::init_stencil (Hierarchy & hierarchy)

{

  int dim = hierarchy.dimension();

  HYPRE_SStructStencilCreate (dim,dim*2+1,&stencil_);

  int entries[][3] = { {  0, 0, 0 },
		       {  1, 0, 0 },
		       { -1, 0, 0 },
		       {  0, 1, 0 },
		       {  0,-1, 0 },
		       {  0, 0, 1 },
		       {  0, 0,-1 } };

  HYPRE_SStructStencilSetEntry (stencil_, 0, entries[0], 0);

  if (dim >= 1) HYPRE_SStructStencilSetEntry (stencil_, 1, entries[1], 0);
  if (dim >= 1) HYPRE_SStructStencilSetEntry (stencil_, 2, entries[2], 0);
  if (dim >= 2) HYPRE_SStructStencilSetEntry (stencil_, 3, entries[3], 0);
  if (dim >= 2) HYPRE_SStructStencilSetEntry (stencil_, 4, entries[4], 0);
  if (dim >= 3) HYPRE_SStructStencilSetEntry (stencil_, 5, entries[5], 0);
  if (dim >= 3) HYPRE_SStructStencilSetEntry (stencil_, 6, entries[6], 0);

}

//======================================================================

/// Initialize the graph.

/** Creates a graph containing the matrix non-zero structure.  Graph
    edges include both those for non-zeros from the stencil within
    each part (AMR hierarchy level), and for graph entries connecting
    linked parts.  

    Setting up the matrix elements is done using the following
    steps:

    1. Define the stencil for all interior grid elements

    2. Zero-out matrix elements covered by a refined grid

    3. Handle matrix elements connecting neighboring grids

     - Handle matrix elements defining the boundary conditions

     - Handle coarse unknowns adjacent to fine unknowns in child
     - Handle fine unknowns adjacent to coarse unknowns in parent

     - Handle coarse unknowns adjacent to fine unknowns in neighbor's child
     - Handle fine unknowns adjacent to coarse unknowns in parent's neighbor
     
     - Handle any remaining connections

    Note that the matrix generated is not symmetric.

*/

void Hypre::init_graph (Hierarchy & hierarchy)

{
  int i;

  const int variable = 0;

  for (i=0; i<hierarchy.num_grids(); i++) {

    Grid & grid = hierarchy.grid(i);

    if (grid.is_local()) {

      // Create the hypre graph for the grid

      if (debug_hypre_graph) printf ("DEBUG_HYPRE_GRAPH HYPRE_SStructGraphCreate()\n",
			       __FILE__,__LINE__);

      HYPRE_SStructGraphCreate (MPI_COMM_WORLD, 
				grid.hypre_grid(), 
				&grid.hypre_graph());

      // Set the stencil for the grid

      if (debug_hypre_graph) printf ("DEBUG_HYPRE_GRAPH HYPRE_SStructGraphSetStencil()\n",
			       __FILE__,__LINE__);

      HYPRE_SStructGraphSetStencil (grid.hypre_graph(),
				    grid.id(),
				    0,
				    stencil_);

      // Set any additional entries between grid and parent

      if (grid.level() > 0) {

	//  Loop over boundary zones
	//    Based on boundary zones types
	//	HYPRE_SStructGraphAddEntries (grid.hypre_graph(),
	//				      grid.level(),
	//				      grid.INDEX
	//				      neighbor.level(),
	//				      neighbor.INDEX,
	//				      variable);

      }

      // Assemble the graph

      if (debug_hypre_graph) printf ("DEBUG_HYPRE_GRAPH HYPRE_SStructGraphAssemble()\n",
			       __FILE__,__LINE__);
      HYPRE_SStructGraphAssemble (grid.hypre_graph());
    }
  }
}

//======================================================================

/// Initialize the matrix A

void Hypre::init_matrix (Hierarchy & hierarchy)

{
  printf ("Hypre::init_matrix() is not implemented yet\n"); 
}

//======================================================================

/// Initialize the right-hand-side vector b

void Hypre::init_rhs (Hierarchy & hierarchy)

{
  printf ("Hypre::init_rhs() is not implemented yet\n"); 
}

//======================================================================

/// Initialize the linear solver

void Hypre::init_solver (Hierarchy & hierarchy)

{
  printf ("Hypre::init_solver() is not implemented yet\n"); 
}

//======================================================================

/// Solve the linear system Ax = b

void Hypre::solve (Hierarchy & hierarchy)

{
  printf ("Hypre::solve() is not implemented yet\n"); 
}

//======================================================================

/// Evaluate the success of the solve

void Hypre::evaluate (Hierarchy & hierarchy)

{
  printf ("Hypre::evaluate() is not implemented yet\n"); 
}

