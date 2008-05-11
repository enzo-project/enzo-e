//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Minimalistic driver for hypre AMR problem 1: two nested grids on two procs

/**
 * 
 * 
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 *
 */

#include <stdio.h>
#include <mpi.h>
#include "HYPRE_sstruct_ls.h"

const bool debug = true;

main(int argc, char * argv[])
{

  //------------------------------------------------------------
  // Initialize MPI
  //------------------------------------------------------------

  int ip,np;

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &np);
  MPI_Comm_rank (MPI_COMM_WORLD, &ip);


  //------------------------------------------------------------
  // Parse command-line arguments
  //------------------------------------------------------------

  int N;
  if (argc == 2) {
    N = atoi(argv[1]);
  } else {
    if (ip==0) fprintf (stderr, "Usage: %s <size>\n",argv[0]);
    MPI_Finalize();
    exit(1);
  }

  //------------------------------------------------------------
  // Initialize hypre grid
  //------------------------------------------------------------

  HYPRE_SStructGrid grid;


  // Create the grid

  HYPRE_SStructGridCreate (MPI_COMM_WORLD, 3, 2, &grid);

  // Initialize grid part extents

  int fine_lower[3]   = {0,0,0};
  int fine_upper[3]   = {N-1,N-1,N-1};
  int coarse_lower[3] = {N/4,N/4,N/4};
  int coarse_upper[3] = {3*N/4-1,3*N/4-1,3*N/4-1};

  HYPRE_SStructGridSetExtents(grid, 0, fine_lower,   fine_upper);

  HYPRE_SStructGridSetExtents(grid, 1, coarse_lower, coarse_upper);

  if (debug) {
    printf ("Fine grid extents:   (%d %d %d) to (%d %d %d)\n",
	    fine_lower[0],fine_lower[1],fine_lower[2],
	    fine_upper[0],fine_upper[1],fine_upper[2]);

    printf ("Coarse grid extents: (%d %d %d) to (%d %d %d)\n",
	    coarse_lower[0],coarse_lower[1],coarse_lower[2],
	    coarse_upper[0],coarse_upper[1],coarse_upper[2]);
  }

  // Set grid variables

  HYPRE_SStructVariable variable_type[] = { HYPRE_SSTRUCT_VARIABLE_CELL };
  HYPRE_SStructGridSetVariables(grid, 0, 1, variable_type);
  HYPRE_SStructGridSetVariables(grid, 1, 1, variable_type);

  // Assemble the grid
  
  HYPRE_SStructGridAssemble (grid);


  //------------------------------------------------------------
  // Create hypre stencil
  //------------------------------------------------------------

  HYPRE_SStructStencil stencil;
  
  int entries[][3] = { {  0, 0, 0 },
		       {  1, 0, 0 },
		       { -1, 0, 0 },
		       {  0, 1, 0 },
		       {  0,-1, 0 },
		       {  0, 0, 1 },
		       {  0, 0,-1 } };

  HYPRE_SStructStencilCreate (3,7,&stencil);

  HYPRE_SStructStencilSetEntry (stencil, 0, entries[0], 0);
  HYPRE_SStructStencilSetEntry (stencil, 1, entries[1], 0);
  HYPRE_SStructStencilSetEntry (stencil, 2, entries[2], 0);
  HYPRE_SStructStencilSetEntry (stencil, 3, entries[3], 0);
  HYPRE_SStructStencilSetEntry (stencil, 4, entries[4], 0);
  HYPRE_SStructStencilSetEntry (stencil, 5, entries[5], 0);
  HYPRE_SStructStencilSetEntry (stencil, 6, entries[6], 0);

  //------------------------------------------------------------
  // Create hypre graph
  //------------------------------------------------------------

  HYPRE_SStructGraph   graph;

  // Create the graph

  HYPRE_SStructGraphCreate (MPI_COMM_WORLD, grid, &graph);

  // Initialize stencil part
  
  HYPRE_SStructGraphSetStencil (graph, 0, 0, stencil);
  HYPRE_SStructGraphSetStencil (graph, 1, 0, stencil);

  // Initialize nonstencil part

  
  //  HYPRE_SStructGraphAddEntries 
  //    (graph, level_fine, igg3, 0, level_coarse, ign3, 0);
  //  HYPRE_SStructGraphAssemble (graph_);

  //------------------------------------------------------------
  // Create hypre matrix and vectors
  //------------------------------------------------------------

  //  HYPRE_SStructMatrix  A_;       // hypre matrix
  //  HYPRE_SStructVector  B_;       // hypre vector right-hand side
  //  HYPRE_SStructVector  X_;       // hypre vector solution

  //  HYPRE_SStructMatrixCreate (MPI_COMM_WORLD, graph_, &A_);
  //  HYPRE_SStructVectorCreate (MPI_COMM_WORLD, grid_,  &X_);
  //  HYPRE_SStructMatrixSetObjectType (A_,HYPRE_SSTRUCT);
  //  HYPRE_SStructVectorSetObjectType (X_,HYPRE_SSTRUCT);
  //  HYPRE_SStructMatrixInitialize (A_);
  //  HYPRE_SStructVectorInitialize (X_);
  //	HYPRE_SStructVectorAddToBoxValues (B_,part,lower,upper,0,values);
  //  HYPRE_SStructMatrixAssemble (A_);
  //  HYPRE_SStructVectorAssemble (B_);
  //		    HYPRE_SStructMatrixAddToValues 
  //		      (A_, level_fine, igg3, 0, 1, &entry, &val);
  //  HYPRE_SStructMatrixSetBoxValues (A_,level,low,up,0,1,&entries[0],v0);

  //------------------------------------------------------------
  // Solve the linear system
  //------------------------------------------------------------

  //  HYPRE_SStructSolver  solver_;  // hypre solver

  //  HYPRE_SStructSysPFMGCreate    (MPI_COMM_WORLD, &solver_);
  //  HYPRE_SStructSysPFMGSetLogging(solver_, 1);
  //  HYPRE_SStructSysPFMGSetup     (solver_,A_,B_,X_);
  //  HYPRE_SStructSysPFMGSolve     (solver_,A_,B_,X_);
  //  HYPRE_SStructSysPFMGGetNumIterations (solver_,&iter_);
  //  HYPRE_SStructSysPFMGGetFinalRelativeResidualNorm (solver_,&resid_);

  //  printf ("HYPRE_SStructSysPFMGSolve num iterations: %d\n",iter_);
  //  printf ("HYPRE_SStructSysPFMGSolve final relative residual norm: %g\n",resid_);
  //  HYPRE_SStructSysPFMGDestroy (solver_);

  //------------------------------------------------------------
  // Exit
  //------------------------------------------------------------

  printf ("Finished!\n");
  MPI_Finalize();
}
