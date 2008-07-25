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
#include <math.h>
#include <assert.h>

#include "HYPRE_sstruct_ls.h"

//----------------------------------------------------------------------
// PARAMETERS
//----------------------------------------------------------------------

const double MASS           = 1e43;
const double BOX            = 8e9;
const int    itmax          = 50;
const double restol         = 1e-6;
const bool   use_fac_solver = false;
const bool   debug          = true;

//----------------------------------------------------------------------
// CONSTANTS
//----------------------------------------------------------------------

const double PI             = M_PI;
const double G              = 6.67428e-8;
const int    num_parts      = 2; // Number of hypre "parts"
const int    part_coarse    = 0; // Coarse part id
const int    part_fine      = 1; // Fine part id
const int    r              = 2; // refinement factor TESTED FOR 2 ONLY

//----------------------------------------------------------------------
// FUNCTIONS
//----------------------------------------------------------------------

int min (int i, int j) { return i<j ? i : j;};
int index(int i0, int i1, int i2, int N) {
  return i0 + N*(i1 + N*i2);
}

//----------------------------------------------------------------------
// MACROS
//----------------------------------------------------------------------

#define TRACE \
  { \
    printf ("%s:%d TRACE\n",__FILE__,__LINE__); fflush(stdout); \
    fflush(stdout); \
  }

#define PTRACE \
  { \
    printf ("%d %s:%d TRACE\n",mpi_rank,__FILE__,__LINE__); \
    fflush(stdout); \
  }

//======================================================================
// MAIN
//======================================================================

int main(int argc, char * argv[])
{

  //------------------------------------------------------------
  // INITIALIZE MPI
  //    Out: mpi_size      (number of MPI processors)
  //    Out: mpi_rank      (rank of this MPI process)
  //    Out: is_mpi_coarse (whether MPI process owns coarse grid)
  //    Out: is_mpi_fine   (whether MPI process owns fine grid)
  //------------------------------------------------------------

  TRACE;

  int mpi_rank,mpi_size;

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);

  // MPI rank for grids

  const int mpi_rank_coarse   = 0;            
  const int mpi_rank_fine     = min(1,mpi_size-1);
  assert (mpi_rank_fine == 0 || mpi_rank_fine == 1);

  // Boolean variables for whether MPI process owns coarse or fine grids

  const bool is_mpi_coarse = (mpi_rank == mpi_rank_coarse);
  const bool is_mpi_fine   = (mpi_rank == mpi_rank_fine);

  printf ("%d %d %d\n",mpi_rank, is_mpi_coarse, is_mpi_fine);

  //------------------------------------------------------------
  // PARSE AND CHECK ARGUMENTS
  //    Out: N  (problem size)
  //------------------------------------------------------------

  PTRACE;

  int N;
  if (argc == 2) {
    N = atoi(argv[1]);
  } else {
    if (mpi_rank==0) fprintf (stderr, "Usage: %s <size>\n",argv[0]);
    MPI_Finalize();
    exit(1);
  }

  //------------------------------------------------------------
  // CREATE AND INITIALIZE HYPRE 'GRID'
  //    Out: grid  (a two-level hypre grid)
  //------------------------------------------------------------

  PTRACE;

  HYPRE_SStructGrid grid;

  // Create the grid

  HYPRE_SStructGridCreate (MPI_COMM_WORLD, 3, num_parts, &grid);

  // Initialize coarse grid extents on coarse grid MPI process

  HYPRE_SStructVariable variable_type[] = { HYPRE_SSTRUCT_VARIABLE_CELL };

  int lower_coarse[3] = {    0,      0,      0};
  int upper_coarse[3] = {  N-1,    N-1,    N-1};

  if (is_mpi_coarse) {

    HYPRE_SStructGridSetExtents(grid, part_coarse, lower_coarse, upper_coarse);

    if (debug) {
      printf ("Coarse grid extents: (%d %d %d) to (%d %d %d)\n",
	      lower_coarse[0],lower_coarse[1],lower_coarse[2],
	      upper_coarse[0],upper_coarse[1],upper_coarse[2]);
    }

    HYPRE_SStructGridSetVariables(grid, part_fine,   1, variable_type);

  }

  // Initialize fine grid extents on fine grid MPI process

  int lower_fine[3]   = {  N/r,    N/r,    N/r};
  int upper_fine[3]   = {3*N/r-1,3*N/r-1,3*N/r-1};

  if (is_mpi_fine) {

    HYPRE_SStructGridSetExtents(grid, part_fine,   lower_fine,   upper_fine);

    if (debug) {
      printf ("Fine grid extents:   (%d %d %d) to (%d %d %d)\n",
	      lower_fine[0],lower_fine[1],lower_fine[2],
	      upper_fine[0],upper_fine[1],upper_fine[2]);
    }

    HYPRE_SStructGridSetVariables(grid, part_coarse, 1, variable_type);
  }

  // Assemble the grid
  
  HYPRE_SStructGridAssemble (grid);

  //------------------------------------------------------------
  // CREATE HYPRE STENCIL
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
  // CREATE HYPRE GRAPH
  //    Out: graph  (unstructured connections)
  //------------------------------------------------------------

  PTRACE;

  HYPRE_SStructGraph   graph;

  // Create the graph

  HYPRE_SStructGraphCreate (MPI_COMM_WORLD, grid, &graph);

  // Initialize stencil part
  
  if (is_mpi_coarse) {
    HYPRE_SStructGraphSetStencil (graph, part_coarse, 0, stencil);
  }
  if (is_mpi_fine) {
    HYPRE_SStructGraphSetStencil (graph, part_fine,   0, stencil);
  }

  // Initialize nonstencil part

  int axis,face;     // 0 <= axis <= 2;  0 <= face <= 1
  int ic0,ic1;       // coarse face indices
  int ind_fine[3];   // fine grid indices
  int ind_coarse[3]; // coarse grid indices

  //----------------------------------------
  // GRAPH ENTRIES: FINE-TO-COARSE
  //----------------------------------------

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
	}
      }
    }
  }

  //----------------------------------------
  // GRAPH ENTRIES: COARSE-TO-FINE
  //----------------------------------------

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

  // Create the hypre vector right-hand side B

  HYPRE_SStructVector  B;       
  HYPRE_SStructVectorCreate (MPI_COMM_WORLD, grid,  &B);
  HYPRE_SStructVectorSetObjectType (B,HYPRE_SSTRUCT);
  HYPRE_SStructVectorInitialize (B);

  double coarse_h =       BOX / N;  // Coarse mesh width
  double fine_h   = 0.5 * BOX / N;  // Fine mesh width

  int ind3[8][3] = {
    {N-1, N-1, N-1},
    {N-1, N-1, N},
    {N-1, N,   N-1},
    {N-1, N,   N},
    {N,   N-1, N-1},
    {N,   N-1, N},
    {N,   N,   N-1},
    {N,   N,   N}};

  double bval = -4.0 * G * PI * MASS / (fine_h*fine_h*fine_h);

  HYPRE_SStructVectorAddToValues (B, part_fine, ind3[0], 0, &bval);
  HYPRE_SStructVectorAddToValues (B, part_fine, ind3[1], 0, &bval);
  HYPRE_SStructVectorAddToValues (B, part_fine, ind3[2], 0, &bval);
  HYPRE_SStructVectorAddToValues (B, part_fine, ind3[3], 0, &bval);
  HYPRE_SStructVectorAddToValues (B, part_fine, ind3[4], 0, &bval);
  HYPRE_SStructVectorAddToValues (B, part_fine, ind3[5], 0, &bval);
  HYPRE_SStructVectorAddToValues (B, part_fine, ind3[6], 0, &bval);
  HYPRE_SStructVectorAddToValues (B, part_fine, ind3[7], 0, &bval);

  HYPRE_SStructVectorAssemble (B);

  // Create the hypre matrix A

  HYPRE_SStructMatrix  A;
  HYPRE_SStructMatrixCreate (MPI_COMM_WORLD, graph, &A);
  HYPRE_SStructMatrixSetObjectType (A,HYPRE_SSTRUCT);
  HYPRE_SStructMatrixInitialize (A);

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

  //--------------------------------------------------
  // MATRIX A FINE-TO-CORSE ENTRIES
  //--------------------------------------------------

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

	}
      }
    }
  }

  //--------------------------------------------------
  // MATRIX A COARSE-TO-FINE ENTRIES
  //--------------------------------------------------

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

  HYPRE_SStructMatrixAssemble (A);

  // Create the hypre vector solution X

  HYPRE_SStructVector  X;
  HYPRE_SStructVectorCreate (MPI_COMM_WORLD, grid,  &X);
  HYPRE_SStructVectorSetObjectType (X,HYPRE_SSTRUCT);
  HYPRE_SStructVectorInitialize (X);
  HYPRE_SStructVectorAssemble (X);

  // Assemble the matrix and vectors


  //  HYPRE_SStructMatrixAddToValues 
  //		      (A_, level_fine, igg3, 0, 1, &entry, &val);


  PTRACE;

  //------------------------------------------------------------
  // Solve the linear system
  //------------------------------------------------------------

  HYPRE_SStructSolver  solver;  // hypre solver
  int iter;
  double resid;

  if (use_fac_solver) {
    
    HYPRE_SStructFACCreate(MPI_COMM_WORLD, &solver);

    int num_parts = 2;
    HYPRE_SStructFACSetMaxLevels(solver,  num_parts);

    int parts[2] = {0,1};
    HYPRE_SStructFACSetPLevels(solver, num_parts, parts);

    int refinements[2][3] = {{2,2,2},{2,2,2}};
    HYPRE_SStructFACSetPRefinements(solver, num_parts, refinements);

    // solver parameters

    HYPRE_SStructFACSetNumPreRelax(solver,      2);
    HYPRE_SStructFACSetNumPostRelax(solver,     2);
    HYPRE_SStructFACSetCoarseSolverType(solver, 1);
    HYPRE_SStructFACSetRelaxType(solver,        2);

    HYPRE_SStructFACSetMaxIter(solver,itmax);
    HYPRE_SStructFACSetTol(solver,    restol);

    HYPRE_SStructFACSetLogging(solver, 1);

    HYPRE_SStructFACSetup2(solver, A, B, X);

    // SOLVE
    HYPRE_SStructFACSolve3(solver, A, B, X);

    HYPRE_SStructFACGetNumIterations(solver, &iter);
    HYPRE_SStructFACGetFinalRelativeResidualNorm(solver, &resid);

  } else {
    HYPRE_SStructBiCGSTABCreate(MPI_COMM_WORLD, &solver);
    HYPRE_SStructBiCGSTABSetLogging(solver, 1);
    HYPRE_SStructBiCGSTABSetup(solver, A, B, X);
    // SOLVE
    HYPRE_SStructBiCGSTABSolve(solver, A, B, X);
    HYPRE_SStructBiCGSTABGetNumIterations(solver, &iter);
    HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm(solver, &resid);
  }

  printf ("HYPRE_SStructSysPFMGSolve num iterations:               %d\n",iter);
  printf ("HYPRE_SStructSysPFMGSolve final relative residual norm: %g\n",resid);

  if (use_fac_solver) {
    HYPRE_SStructFACDestroy2(solver);
  } else {
    HYPRE_SStructBiCGSTABDestroy(solver);
  }

  //------------------------------------------------------------
  // Exit
  //------------------------------------------------------------

  printf ("Finished!\n");
  MPI_Finalize();
}

// Write Grid
//  if (fp == 0) fp = stdout;
//  fprintf (fp,"Grid\n"
//	  "   id             %d\n"
//	  "   parent id      %d\n"
//	  "   processor      %d\n"
//	  "   lower position "SCALAR_PRINTF SCALAR_PRINTF SCALAR_PRINTF"\n"
//	  "   upper position "SCALAR_PRINTF SCALAR_PRINTF SCALAR_PRINTF"\n"
//	  "   lower index    %d %d %d\n"
//	  "   zones          %d %d %d\n"
//	  "   level          %d\n",
//	  id_,id_parent_,ip_,
//	  xl_[0],xl_[1],xl_[2],
//	  xu_[0],xu_[1],xu_[2],
//	  il_[0],il_[1],il_[2],
//	   n_ [0],n_ [1],n_ [2],
//	   level_);
//  if (u_ && ! brief) {
//    for (int i0=0; i0<n_[0]; i0++) {
//      for (int i1=0; i1<n_[1]; i1++) {
//	for (int i2=0; i2<n_[2]; i2++) {
//	  int i = index(i0,i1,i2,n_[0],n_[1],n_[2]);
//	  fprintf (fp,"%d %d %d %g\n",i0,i1,i2,u_[i]);
//	}
//      }
//    }
//  }
