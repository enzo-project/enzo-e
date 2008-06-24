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
#include "HYPRE_sstruct_ls.h"

const double PI = 3.14159265358979323;
const double G  = 6.67428e-8;
const double MASS = 1e43;
const double BOX = 8e9;

const int COARSE = 0;
const int FINE   = 1;

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

  int coarse_lower[3]   = {0,0,0};
  int coarse_upper[3]   = {N-1,N-1,N-1};
  int fine_lower[3] = {N/2,N/2,N/2};
  int fine_upper[3] = {3*N/2-1,3*N/2-1,3*N/2-1};

  HYPRE_SStructGridSetExtents(grid, 1, fine_lower,   fine_upper);
  HYPRE_SStructGridSetExtents(grid, 0, coarse_lower, coarse_upper);

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

  int axis,face; // axis = 0:2;  face = 0:1
  int ic0,ic1;     // 
  int ifg3[3]; // fine grid indices
  int icg3[3]; // coarse grid indices

  for (axis=0; axis<3; axis++) {

    int j0 = axis;
    int j1 = (axis+1)%3;
    int j2 = (axis+2)%3;

    for (face=0; face<2; face++) {

      // loop over coarse face zones

      for (ic0=0; ic0<N/2; ic0++) {
	for (ic1=0; ic1<N/2; ic1++) {

	  // coarse zone index

	  icg3[j0] = (1-face) * (N/4-1) + (face) * (3*N/4);
	  icg3[j1] = N/4 + ic0;
	  icg3[j2] = N/4 + ic1;

	  // fine zone index 00

	  ifg3[j0] = 2 * (icg3[j0] + (1-face)*1 + face * (-1));
	  ifg3[j1] = 2 * icg3[j1];
	  ifg3[j2] = 2 * icg3[j2];

	  // fine-to-coarse
	  
	  HYPRE_SStructGraphAddEntries (graph, FINE, ifg3, 0, COARSE, icg3, 0);
	  ++ ifg3[j1]; 	  // fine zone index 010
	  HYPRE_SStructGraphAddEntries (graph, FINE, ifg3, 0, COARSE, icg3, 0);
	  ++ ifg3[j2];	  // fine zone index 011
	  HYPRE_SStructGraphAddEntries (graph, FINE, ifg3, 0, COARSE, icg3, 0);
	  -- ifg3[j1];	  // fine zone index 001
	  HYPRE_SStructGraphAddEntries (graph, FINE, ifg3, 0, COARSE, icg3, 0);
	  -- ifg3[j2];	  // fine zone index 000

	  // coarse-to-fine

	  HYPRE_SStructGraphAddEntries (graph, COARSE, icg3, 0, FINE, ifg3, 0);
	  ++ ifg3[j0]; 	  // fine zone index 100
	  HYPRE_SStructGraphAddEntries (graph, COARSE, icg3, 0, FINE, ifg3, 0);
	  ++ ifg3[j1];	  // fine zone index 110
	  HYPRE_SStructGraphAddEntries (graph, COARSE, icg3, 0, FINE, ifg3, 0);
	  -- ifg3[j0];	  // fine zone index 010
	  HYPRE_SStructGraphAddEntries (graph, COARSE, icg3, 0, FINE, ifg3, 0);
	  ++ ifg3[j2];	  // fine zone index 011
	  HYPRE_SStructGraphAddEntries (graph, COARSE, icg3, 0, FINE, ifg3, 0);
	  ++ ifg3[j0]; 	  // fine zone index 111
	  HYPRE_SStructGraphAddEntries (graph, COARSE, icg3, 0, FINE, ifg3, 0);
	  -- ifg3[j1];	  // fine zone index 101
	  HYPRE_SStructGraphAddEntries (graph, COARSE, icg3, 0, FINE, ifg3, 0);
	  -- ifg3[j0];	  // fine zone index 001
	  HYPRE_SStructGraphAddEntries (graph, COARSE, icg3, 0, FINE, ifg3, 0);
	  -- ifg3[j2];	  // fine zone index 000

	}
      }
      
    }
  }

  HYPRE_SStructGraphAssemble (graph);

  //------------------------------------------------------------
  // Create hypre matrix and vectors
  //------------------------------------------------------------

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

  double coarse_h =       BOX / N;
  double fine_h   = 0.5 * BOX / N;

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
    v1[i] = coarse_h;
  }

  HYPRE_SStructMatrixSetBoxValues (A,COARSE,coarse_lower,coarse_upper,
				   0,1,&nums[0],v0);
  HYPRE_SStructMatrixSetBoxValues (A,COARSE,coarse_lower,coarse_upper,
				   0,1,&nums[1],v1);
  HYPRE_SStructMatrixSetBoxValues (A,COARSE,coarse_lower,coarse_upper,
				   0,1,&nums[2],v1);
  HYPRE_SStructMatrixSetBoxValues (A,COARSE,coarse_lower,coarse_upper,
				   0,1,&nums[3],v1);
  HYPRE_SStructMatrixSetBoxValues (A,COARSE,coarse_lower,coarse_upper,
				   0,1,&nums[4],v1);
  HYPRE_SStructMatrixSetBoxValues (A,COARSE,coarse_lower,coarse_upper,
				   0,1,&nums[5],v1);
  HYPRE_SStructMatrixSetBoxValues (A,COARSE,coarse_lower,coarse_upper,
				   0,1,&nums[6],v1);

  for (i=0; i<N*N*N; i++) {
    v0[i] = -6*fine_h;
    v1[i] = fine_h;
  }

  HYPRE_SStructMatrixSetBoxValues (A,FINE,fine_lower,fine_upper,
				   0,1,&nums[0],v0);
  HYPRE_SStructMatrixSetBoxValues (A,FINE,fine_lower,fine_upper,
				   0,1,&nums[1],v1);
  HYPRE_SStructMatrixSetBoxValues (A,FINE,fine_lower,fine_upper,
				   0,1,&nums[2],v1);
  HYPRE_SStructMatrixSetBoxValues (A,FINE,fine_lower,fine_upper,
				   0,1,&nums[3],v1);
  HYPRE_SStructMatrixSetBoxValues (A,FINE,fine_lower,fine_upper,
				   0,1,&nums[4],v1);
  HYPRE_SStructMatrixSetBoxValues (A,FINE,fine_lower,fine_upper,
				   0,1,&nums[5],v1);
  HYPRE_SStructMatrixSetBoxValues (A,FINE,fine_lower,fine_upper,
				   0,1,&nums[6],v1);


  // Matrix graph entries

  for (axis=0; axis<3; axis++) {

    int j0 = axis;
    int j1 = (axis+1)%3;
    int j2 = (axis+2)%3;

    //
    //
    //  +-+-+-+-+-+-+-+-+
    //  | | | | | | | | |
    //  +-+-+-+-+-+-+-+-+
    //  |
    //  |
    //  |
    //
    //
    //
    //
    //
    //
    //

    for (face=0; face<2; face++) {

      // loop over coarse face zones

      for (ic0=0; ic0<N/2; ic0++) {
	for (ic1=0; ic1<N/2; ic1++) {

	  // coarse zone index

	  icg3[j0] = (1-face) * (N/4-1) + (face) * (3*N/4);
	  icg3[j1] = N/4 + ic0;
	  icg3[j2] = N/4 + ic1;

	  // fine zone index 00

	  ifg3[j0] = 2 * (icg3[j0] + (1-face)*1 + face * (-1));
	  ifg3[j1] = 2 * icg3[j1];
	  ifg3[j2] = 2 * icg3[j2];

	  printf ("coarse: %d %d %d   fine:%d %d %d\n",
		  icg3[0],icg3[1],icg3[2],
		  ifg3[0],ifg3[1],ifg3[2]);
	  // fine-to-coarse
	  
	  int * coarse_count = new int [N*N*N];
	  int * fine_count = new int [N*N*N];

	  for (i=0; i<N*N*N; i++) {
	    coarse_count[i] = 0;
	    fine_count[i] = 7;
	  }

	  int    entry;
	  double val;

	  
	  val = fine_h;
	  i = ifg3[0] + N*(ifg3[1] + N*(ifg3[2]));
	  entry = fine_count[i]++;
	  HYPRE_SStructMatrixAddToValues 
	    (A, FINE, ifg3, 0, COARSE, &entry, &val);
	  ++ ifg3[j1]; 	  // fine zone index 010
	  i = ifg3[0] + N*(ifg3[1] + N*(ifg3[2]));
	  entry = fine_count[i]++;
	  HYPRE_SStructMatrixAddToValues 
	    (A, FINE, ifg3, 0, COARSE, &entry, &val);
	  ++ ifg3[j2];	  // fine zone index 011
	  i = ifg3[0] + N*(ifg3[1] + N*(ifg3[2]));
	  entry = fine_count[i]++;
	  HYPRE_SStructMatrixAddToValues 
	    (A, FINE, ifg3, 0, COARSE, &entry, &val);
	  -- ifg3[j1];	  // fine zone index 001
	  i = ifg3[0] + N*(ifg3[1] + N*(ifg3[2]));
	  entry = fine_count[i]++;
	  HYPRE_SStructMatrixAddToValues 
	    (A, FINE, ifg3, 0, COARSE, &entry, &val);
	  -- ifg3[j2];	  // fine zone index 000

	  // coarse-to-fine

	  HYPRE_SStructGraphAddEntries (graph, COARSE, icg3, 0, FINE, ifg3, 0);
	  ++ ifg3[j0]; 	  // fine zone index 100
	  HYPRE_SStructGraphAddEntries (graph, COARSE, icg3, 0, FINE, ifg3, 0);
	  ++ ifg3[j1];	  // fine zone index 110
	  HYPRE_SStructGraphAddEntries (graph, COARSE, icg3, 0, FINE, ifg3, 0);
	  -- ifg3[j0];	  // fine zone index 010
	  HYPRE_SStructGraphAddEntries (graph, COARSE, icg3, 0, FINE, ifg3, 0);
	  ++ ifg3[j2];	  // fine zone index 011
	  HYPRE_SStructGraphAddEntries (graph, COARSE, icg3, 0, FINE, ifg3, 0);
	  ++ ifg3[j0]; 	  // fine zone index 111
	  HYPRE_SStructGraphAddEntries (graph, COARSE, icg3, 0, FINE, ifg3, 0);
	  -- ifg3[j1];	  // fine zone index 101
	  HYPRE_SStructGraphAddEntries (graph, COARSE, icg3, 0, FINE, ifg3, 0);
	  -- ifg3[j0];	  // fine zone index 001
	  HYPRE_SStructGraphAddEntries (graph, COARSE, icg3, 0, FINE, ifg3, 0);
	  -- ifg3[j2];	  // fine zone index 000

	}
      }
      
    }
  }

  // Assemble the matrix and vectors

  HYPRE_SStructMatrixAssemble (A);
  HYPRE_SStructVectorAssemble (X);
  HYPRE_SStructVectorAssemble (B);

  //  HYPRE_SStructMatrixAddToValues 
  //		      (A_, level_fine, igg3, 0, 1, &entry, &val);

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
