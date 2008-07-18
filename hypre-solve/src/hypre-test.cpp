//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Minimalistic driver for hypre AMR problem 1: two nested grids on two procs

/**
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 *
 * $Id$
 *
 */

#include <stdio.h>
#include <mpi.h>
#include <assert.h>

#include "HYPRE_sstruct_ls.h"

const double PI = 3.14159265358979323;
const double G  = 6.67428e-8;
const double MASS = 1e43;
const double BOX = 8e9;

const int num_parts   = 2; // Number of hypre "parts"
const int part_coarse = 0; // Hypre "part" id's
const int part_fine   = 1;

const int r = 2;           // refinement factor TESTED FOR 2 ONLY

int min (int i, int j) { return i<j ? i : j;};
int index(int i0, int i1, int i2, int N) {
  return i0 + N*(i1 + N*i2);
}

const bool debug = true;

#define TRACE { printf ("%s:%d TRACE\n",__FILE__,__LINE__); fflush(stdout); }
#define PTRACE { printf ("%d %s:%d TRACE\n",ip_this,__FILE__,__LINE__); fflush(stdout); }

main(int argc, char * argv[])
{

  //------------------------------------------------------------
  // Action: Initialize MPI
  //    Out: np (processor count)
  //    Out: ip_this (processor rank)
  //------------------------------------------------------------

  TRACE;

  int ip_this,np;

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &np);
  MPI_Comm_rank (MPI_COMM_WORLD, &ip_this);

  int ip_coarse   = 0;            // MPI rank for grids
  int ip_fine     = min(1,np-1); // either (0,0) (serial) or (0,1) (parallel)

  assert (ip_fine == 0 || ip_fine == 1);

  //------------------------------------------------------------
  // Action: Parse arguments and exit with message on error
  //    Out: N  (problem size)
  //------------------------------------------------------------

  PTRACE;

  int N;
  if (argc == 2) {
    N = atoi(argv[1]);
  } else {
    if (ip_this==0) fprintf (stderr, "Usage: %s <size>\n",argv[0]);
    MPI_Finalize();
    exit(1);
  }

  //------------------------------------------------------------
  // Action: create and initialize hypre 'grid'
  //    Out: grid  (a two-level hypre grid)
  //------------------------------------------------------------

  PTRACE;

  HYPRE_SStructGrid grid;


  // Create the grid

  HYPRE_SStructGridCreate (MPI_COMM_WORLD, 3, num_parts, &grid);

  // Initialize grid part extents

  int lower_coarse[3] = {    0,      0,      0};
  int upper_coarse[3] = {  N-1,    N-1,    N-1};
  int lower_fine[3]   = {  N/r,    N/r,    N/r};
  int upper_fine[3]   = {3*N/r-1,3*N/r-1,3*N/r-1};

  HYPRE_SStructGridSetExtents(grid, part_coarse, lower_coarse, upper_coarse);
  HYPRE_SStructGridSetExtents(grid, part_fine,   lower_fine,   upper_fine);

  if (debug) {
    printf ("Fine grid extents:   (%d %d %d) to (%d %d %d)\n",
	    lower_fine[0],lower_fine[1],lower_fine[2],
	    upper_fine[0],upper_fine[1],upper_fine[2]);

    printf ("Coarse grid extents: (%d %d %d) to (%d %d %d)\n",
	    lower_coarse[0],lower_coarse[1],lower_coarse[2],
	    upper_coarse[0],upper_coarse[1],upper_coarse[2]);
  }

  // Set grid variables

  HYPRE_SStructVariable variable_type[] = { HYPRE_SSTRUCT_VARIABLE_CELL };
  HYPRE_SStructGridSetVariables(grid, part_coarse, 1, variable_type);
  HYPRE_SStructGridSetVariables(grid, part_fine,   1, variable_type);

  // Assemble the grid
  
  HYPRE_SStructGridAssemble (grid);

  //------------------------------------------------------------
  // Action: create hypre stencil
  //    Out: stencil  (structured connections)
  //------------------------------------------------------------

  PTRACE;

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
  // Action: create hypre graph
  //    Out: graph  (unstructured connections)
  //------------------------------------------------------------

  PTRACE;

  HYPRE_SStructGraph   graph;

  // Create the graph

  HYPRE_SStructGraphCreate (MPI_COMM_WORLD, grid, &graph);

  // Initialize stencil part
  
  HYPRE_SStructGraphSetStencil (graph, part_coarse, 0, stencil);
  HYPRE_SStructGraphSetStencil (graph, part_fine,   0, stencil);

  // Initialize nonstencil part

  int axis,face;     // 0 <= axis <= 2;  0 <= face <= 1
  int ic0,ic1;       // coarse face indices
  int ind_fine[3];   // fine grid indices
  int ind_coarse[3]; // coarse grid indices

  for (axis=0; axis<3; axis++) {

    int j0 = axis;
    int j1 = (axis+1)%3;
    int j2 = (axis+2)%3;

    for (face=0; face<2; face++) {

      // loop over coarse face zones

      for (ic0=0; ic0<N/r; ic0++) {
	for (ic1=0; ic1<N/r; ic1++) {

	  // coarse zone index

	  ind_coarse[j0] =    (1-face) * (N/4-1) + (face) * (3*N/4);
	  ind_coarse[j1] = N/4 + ic0;
	  ind_coarse[j2] = N/4 + ic1;

	  // fine zone index 000

	  ind_fine[j0]   = r*((1-face) * (N/4)   + (face) * (3*N/4-1));
	  ind_fine[j1]   = r*ind_coarse[j1];
	  ind_fine[j2]   = r*ind_coarse[j2];

	  // FINE-TO-COARSE

	  // fine zone index 000
	  HYPRE_SStructGraphAddEntries (graph, 
					part_fine,   ind_fine, 0, 
					part_coarse, ind_coarse, 0);
	  ++ ind_fine[j1]; 	  

	  // fine zone index 010
	  HYPRE_SStructGraphAddEntries (graph, 
					part_fine,   ind_fine, 0, 
					part_coarse, ind_coarse, 0);
	  ++ ind_fine[j2];	  

	  // fine zone index 011
	  HYPRE_SStructGraphAddEntries (graph, 
					part_fine,   ind_fine, 0, 
					part_coarse, ind_coarse, 0);
	  -- ind_fine[j1];	  

	  // fine zone index 001
	  HYPRE_SStructGraphAddEntries (graph, 
					part_fine,   ind_fine, 0, 
					part_coarse, ind_coarse, 0);
	  -- ind_fine[j2];	  

	  // fine zone index 000

	  // COARSE-TO-FINE

	  HYPRE_SStructGraphAddEntries (graph, 
					part_coarse, ind_coarse, 0, 
					part_fine, ind_fine, 0);
	  ++ ind_fine[j0]; 	  

	  // fine zone index 100
	  HYPRE_SStructGraphAddEntries (graph, 
					part_coarse, ind_coarse, 0, 
					part_fine, ind_fine, 0);
	  ++ ind_fine[j1];	  

	  // fine zone index 110
	  HYPRE_SStructGraphAddEntries (graph, 
					part_coarse, ind_coarse, 0, 
					part_fine, ind_fine, 0);
	  -- ind_fine[j0];	  

	  // fine zone index 010
	  HYPRE_SStructGraphAddEntries (graph, 
					part_coarse, ind_coarse, 0, 
					part_fine, ind_fine, 0);
	  ++ ind_fine[j2];	  

	  // fine zone index 011
	  HYPRE_SStructGraphAddEntries (graph, 
					part_coarse, ind_coarse, 0, 
					part_fine, ind_fine, 0);
	  ++ ind_fine[j0]; 	  

	  // fine zone index 111
	  HYPRE_SStructGraphAddEntries (graph, 
					part_coarse, ind_coarse, 0, 
					part_fine, ind_fine, 0);
	  -- ind_fine[j1];	  

	  // fine zone index 101
	  HYPRE_SStructGraphAddEntries (graph, 
					part_coarse, ind_coarse, 0, 
					part_fine, ind_fine, 0);
	  -- ind_fine[j0];	  

	  // fine zone index 001
	  HYPRE_SStructGraphAddEntries (graph, 
					part_coarse, ind_coarse, 0, 
					part_fine, ind_fine, 0);
	  -- ind_fine[j2];	  
	  // [fine zone index 000]

	}
      }
    }
  }

  HYPRE_SStructGraphAssemble (graph);

  //------------------------------------------------------------
  // Create hypre matrix and vectors
  //------------------------------------------------------------

  PTRACE;

  HYPRE_SStructMatrix  A;       // hypre matrix
  HYPRE_SStructVector  X;       // hypre vector solution
  HYPRE_SStructVector  B;       // hypre vector right-hand side

  HYPRE_SStructMatrixCreate (MPI_COMM_WORLD, graph, &A);
  HYPRE_SStructVectorCreate (MPI_COMM_WORLD, grid,  &B);
  HYPRE_SStructVectorCreate (MPI_COMM_WORLD, grid,  &X);

  HYPRE_SStructMatrixSetObjectType (A,HYPRE_SSTRUCT);
  HYPRE_SStructVectorSetObjectType (X,HYPRE_SSTRUCT);
  HYPRE_SStructVectorSetObjectType (B,HYPRE_SSTRUCT);

  HYPRE_SStructMatrixInitialize (A);
  HYPRE_SStructVectorInitialize (X);
  HYPRE_SStructVectorInitialize (B);

  // Initialize B

  double coarse_h =       BOX / N;  // Coarse mesh width
  double fine_h   = 0.5 * BOX / N;  // Fine mesh width

  int ind3[3] = {N-1, N-1, N-1};
  double bval = -4.0 * G * PI * MASS / (fine_h*fine_h*fine_h);

  HYPRE_SStructVectorAddToValues (B, 1, ind3, 0, &bval);
  ++ ind3[0];	  // fine zone index 100
  HYPRE_SStructVectorAddToValues (B, 1, ind3, 0, &bval);
  ++ ind3[1];	  // fine zone index 110
  HYPRE_SStructVectorAddToValues (B, 1, ind3, 0, &bval);
  -- ind3[0];	  // fine zone index 010
  HYPRE_SStructVectorAddToValues (B, 1, ind3, 0, &bval);
  ++ ind3[2];	  // fine zone index 011
  HYPRE_SStructVectorAddToValues (B, 1, ind3, 0, &bval);
  ++ ind3[0]; 	  // fine zone index 111
  HYPRE_SStructVectorAddToValues (B, 1, ind3, 0, &bval);
  -- ind3[1];	  // fine zone index 101
  HYPRE_SStructVectorAddToValues (B, 1, ind3, 0, &bval);
  -- ind3[0];	  // fine zone index 001
  HYPRE_SStructVectorAddToValues (B, 1, ind3, 0, &bval);
  -- ind3[2];	  // fine zone index 000

  // Initialize A

  // Stencil entries

  int nums[7] = {0,1,2,3,4,5,6};

  double * v0 = new double [N*N*N];
  double * v1 = new double [N*N*N];

  int i;
  for (i=0; i<N*N*N; i++) {
    v0[i] = -6*coarse_h;
    v1[i] =  1*coarse_h;
  }

  HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				   0,1,&nums[0],v0);
  HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				   0,1,&nums[1],v1);
  HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				   0,1,&nums[2],v1);
  HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				   0,1,&nums[3],v1);
  HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				   0,1,&nums[4],v1);
  HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				   0,1,&nums[5],v1);
  HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				   0,1,&nums[6],v1);

  for (i=0; i<N*N*N; i++) {
    v0[i] = -6*fine_h;
    v1[i] =  1*fine_h;
  }

  HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				   0,1,&nums[0],v0);
  HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				   0,1,&nums[1],v1);
  HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				   0,1,&nums[2],v1);
  HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				   0,1,&nums[3],v1);
  HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				   0,1,&nums[4],v1);
  HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				   0,1,&nums[5],v1);
  HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				   0,1,&nums[6],v1);

  delete [] v0;
  delete [] v1;

  PTRACE;

  // Declare counts for zones
  int * count_coarse = new int [N*N*N];
  int * count_fine   = new int [N*N*N];

  PTRACE;
  printf ("DEBUG count_coarse = %p\n",count_coarse);
  printf ("DEBUG count_fine   = %p\n",count_fine);
  for (i=0; i<N*N*N; i++) {
    count_coarse[i] = 0;
    count_fine[i] = 7;
  }

  // Matrix graph entries

  for (axis=0; axis<3; axis++) {

    int j0 = axis;
    int j1 = (axis+1)%3;
    int j2 = (axis+2)%3;

    for (face=0; face<2; face++) {

      // loop over coarse face zones

      for (ic0=0; ic0<N/r; ic0++) {
	for (ic1=0; ic1<N/r; ic1++) {

	  // coarse zone index

	  ind_coarse[j0] =    (1-face) * (N/4-1) + (face) * (3*N/4);
	  ind_coarse[j1] = N/4 + ic0;
	  ind_coarse[j2] = N/4 + ic1;

	  // fine zone index 000

	  ind_fine[j0]   = r*((1-face) * (N/4)   + (face) * (3*N/4-1));
	  ind_fine[j1]   = r*ind_coarse[j1];
	  ind_fine[j2]   = r*ind_coarse[j2];

	  // FINE-TO-COARSE
	  
	  int    entry;
	  double val;
	  
	  val = fine_h;
	  int icount;
	  
	  icount = index(ind_fine[0],ind_fine[1],ind_fine[2],N);
	  icount -= index(N/r,N/r,N/r,N);
	  if (icount < 0 || icount > N*N*N) {
	    // Display error if icount out of range
	    PTRACE;
	    printf ("%d %d %d  %d\n",
		    ind_fine[0],ind_fine[1],ind_fine[2],
		    icount);
	  }
	  assert (icount >= 0);
	  assert (icount < N*N*N);
	  entry = count_fine[icount]++;
	  HYPRE_SStructMatrixAddToValues 
	    (A, part_fine, ind_fine, 0, part_coarse, &entry, &val);
	  ++ ind_fine[j1]; 	  

	  // fine zone index 010
	  icount = index(ind_fine[0],ind_fine[1],ind_fine[2],N);
	  icount -= index(N/r,N/r,N/r,N);
	  if (icount < 0 || icount > N*N*N) {
	    // Display error if icount out of range
	    PTRACE;
	    printf ("%d %d %d  %d\n",
		    ind_fine[0],ind_fine[1],ind_fine[2],
		    icount);
	  }
	  assert (icount >= 0);
	  assert (icount < N*N*N);
	  entry = count_fine[icount]++;
	  HYPRE_SStructMatrixAddToValues 
	    (A, part_fine, ind_fine, 0, part_coarse, &entry, &val);
	  ++ ind_fine[j2];	  

	  // fine zone index 011
	  icount = index(ind_fine[0],ind_fine[1],ind_fine[2],N);
	  icount -= index(N/r,N/r,N/r,N);
	  if (icount < 0 || icount > N*N*N) {
	    // Display error if icount out of range
	    PTRACE;
	    printf ("%d %d %d  %d\n",
		    ind_fine[0],ind_fine[1],ind_fine[2],
		    icount);
	  }
	  assert (icount >= 0);
	  assert (icount < N*N*N);
	  entry = count_fine[icount]++;
	  HYPRE_SStructMatrixAddToValues 
	    (A, part_fine, ind_fine, 0, part_coarse, &entry, &val);
	  -- ind_fine[j1];	  

	  // fine zone index 001
	  icount = index(ind_fine[0],ind_fine[1],ind_fine[2],N);
	  icount -= index(N/r,N/r,N/r,N);
	  if (icount < 0 || icount > N*N*N) {
	    // Display error if icount out of range
	    PTRACE;
	    printf ("%d %d %d  %d\n",
		    ind_fine[0],ind_fine[1],ind_fine[2],
		    icount);
	  }
	  assert (icount >= 0);
	  assert (icount < N*N*N);
	  entry = count_fine[icount]++;
	  HYPRE_SStructMatrixAddToValues 
	    (A, part_fine, ind_fine, 0, part_coarse, &entry, &val);
	  -- ind_fine[j2];	  

	  // fine zone index [000]

	  // COARSE-TO-FINE

	  HYPRE_SStructGraphAddEntries (graph, 
					part_coarse, ind_coarse, 0, 
					part_fine, ind_fine, 0);
	  ++ ind_fine[j0]; 	  

	  // fine zone index 100
	  HYPRE_SStructGraphAddEntries (graph, 
					part_coarse, ind_coarse, 0, 
					part_fine, ind_fine, 0);
	  ++ ind_fine[j1];	  

	  // fine zone index 110
	  HYPRE_SStructGraphAddEntries (graph, 
					part_coarse, ind_coarse, 0, 
					part_fine, ind_fine, 0);
	  -- ind_fine[j0];	  

	  // fine zone index 010
	  HYPRE_SStructGraphAddEntries (graph, 
					part_coarse, ind_coarse, 0, 
					part_fine, ind_fine, 0);
	  ++ ind_fine[j2];	  

	  // fine zone index 011
	  HYPRE_SStructGraphAddEntries (graph, 
					part_coarse, ind_coarse, 0, 
					part_fine, ind_fine, 0);
	  ++ ind_fine[j0]; 	  

	  // fine zone index 111
	  HYPRE_SStructGraphAddEntries (graph, 
					part_coarse, ind_coarse, 0, 
					part_fine, ind_fine, 0);
	  -- ind_fine[j1];	  

	  // fine zone index 101
	  HYPRE_SStructGraphAddEntries (graph, 
					part_coarse, ind_coarse, 0, 
					part_fine, ind_fine, 0);
	  -- ind_fine[j0];	  

	  // fine zone index 001
	  HYPRE_SStructGraphAddEntries (graph, 
					part_coarse, ind_coarse, 0, 
					part_fine, ind_fine, 0);
	  -- ind_fine[j2];	  
	  // fine zone index 000

	}
      }
      
    }
  }

  PTRACE;

  printf ("DEBUG count_coarse = %p\n",count_coarse);
  delete [] count_coarse;
  delete [] count_fine;

  PTRACE;

  // Assemble the matrix and vectors

  HYPRE_SStructMatrixAssemble (A);
  HYPRE_SStructVectorAssemble (X);
  HYPRE_SStructVectorAssemble (B);

  //  HYPRE_SStructMatrixAddToValues 
  //		      (A_, level_fine, igg3, 0, 1, &entry, &val);


  PTRACE;

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
