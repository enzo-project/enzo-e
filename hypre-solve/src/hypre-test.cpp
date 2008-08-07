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

const double MASS    = 1e43;
const double BOX     = 8e9;
const int    itmax   = 50;
const double restol  = 1e-6;
const bool   use_fac = false;
const bool   debug   = false;
const bool   trace   = true;

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
  if (trace) {							\
    printf ("%s:%d TRACE\n",__FILE__,__LINE__); fflush(stdout); \
    fflush(stdout); \
  }

#define PTRACE \
  if (trace) {						    \
    printf ("%d %s:%d TRACE\n",mpi_rank,__FILE__,__LINE__); \
    fflush(stdout); \
  }

#define ASSERT_BOUND(A,B,C) \
  { \
    if ((B) < (A) || (C) < (B) ) { \
      printf ("%s:%d  ASSERT_BOUND: %d %d %d\n", \
	      __FILE__,__LINE__,A,B,C); \
      assert (0); \
    } \
  }

//======================================================================
// MAIN
//======================================================================

int main(int argc, char * argv[])
{

  int i,i0,i1,i2;

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

  if (debug) printf ("mpi_rank=%d  is_mpi_coarse=%d is_mpi_fine=%d\n",
		     mpi_rank, is_mpi_coarse, is_mpi_fine);

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

    HYPRE_SStructGridSetVariables(grid, part_coarse, 1, variable_type);

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

    HYPRE_SStructGridSetVariables(grid, part_fine,   1, variable_type);
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

  PTRACE;
  if (is_mpi_coarse) {
    HYPRE_SStructGraphSetStencil (graph, part_coarse, 0, stencil);
  }
  if (is_mpi_fine) {
    HYPRE_SStructGraphSetStencil (graph, part_fine,   0, stencil);
  }
  PTRACE;

  // Initialize nonstencil part

  int axis,face;     // 0 <= axis <= 2;  0 <= face <= 1
  int ic0,ic1;       // coarse face indices
  int ind_fine[3];   // fine grid indices
  int ind_coarse[3]; // coarse grid indices

  //----------------------------------------
  // GRAPH ENTRIES: FINE-TO-COARSE
  //----------------------------------------

  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // @@@@@ THIS LOOP CRASHES IN PARALLEL Rev. 374 @@@@@
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

	   // ind_fine computed above corresponds to coarse zone
	  if (face == 1)  ++ ind_fine[j0];

	  // 000 ---------------------------------------------
	  HYPRE_SStructGraphAddEntries 
	    (graph, part_fine, ind_fine, 0, part_coarse, ind_coarse, 0);
	  ++ ind_fine[j1]; 	  
	  // 010 ---------------------------------------------
	  HYPRE_SStructGraphAddEntries
	    (graph, part_fine, ind_fine, 0, part_coarse, ind_coarse, 0);
	  ++ ind_fine[j2];	  
	  // 011 ---------------------------------------------
	  HYPRE_SStructGraphAddEntries
	    (graph, part_fine, ind_fine, 0, part_coarse, ind_coarse, 0);
	  -- ind_fine[j1];	  
	  // 001 ---------------------------------------------
	  HYPRE_SStructGraphAddEntries
	    (graph, part_fine, ind_fine, 0, part_coarse, ind_coarse, 0);
	  -- ind_fine[j2];	  
	  // 000 ---------------------------------------------
	}
      }
    }
  }

  PTRACE;

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

	  // 100 ---------------------------------------------
	  HYPRE_SStructGraphAddEntries 
	    (graph, part_coarse, ind_coarse, 0, part_fine, ind_fine, 0);
	  ++ ind_fine[j1];	  
	  // 110 ---------------------------------------------
	  HYPRE_SStructGraphAddEntries 
	    (graph, part_coarse, ind_coarse, 0, part_fine, ind_fine, 0);
	  -- ind_fine[j0];	  
	  // 010 ---------------------------------------------
	  HYPRE_SStructGraphAddEntries 
	    (graph, part_coarse, ind_coarse, 0, part_fine, ind_fine, 0);
	  ++ ind_fine[j2];	  
	  // 011 ---------------------------------------------
	  HYPRE_SStructGraphAddEntries 
	    (graph, part_coarse, ind_coarse, 0, part_fine, ind_fine, 0);
	  ++ ind_fine[j0]; 	  
	  // 111 ---------------------------------------------
	  HYPRE_SStructGraphAddEntries 
	    (graph, part_coarse, ind_coarse, 0, part_fine, ind_fine, 0);
	  -- ind_fine[j1];	  
	  // 101 ---------------------------------------------
	  HYPRE_SStructGraphAddEntries 
	    (graph, part_coarse, ind_coarse, 0, part_fine, ind_fine, 0);
	  -- ind_fine[j0];	  
	  // 001 ---------------------------------------------
	  HYPRE_SStructGraphAddEntries 
	    (graph, part_coarse, ind_coarse, 0, part_fine, ind_fine, 0);
	  -- ind_fine[j2];	  
	  // 000 ---------------------------------------------
	}
      }
    }
  }

  PTRACE;

  HYPRE_SStructGraphAssemble (graph);

  //============================================================
  // CREATE HYPRE MATRIX A AND VECTORS X,B
  //============================================================

  PTRACE;

  //--------------------------------------------------
  // Create the hypre vector right-hand side B
  //--------------------------------------------------

  HYPRE_SStructVector  B;       
  HYPRE_SStructVectorCreate (MPI_COMM_WORLD, grid,  &B);
  HYPRE_SStructVectorSetObjectType (B,HYPRE_SSTRUCT);
  HYPRE_SStructVectorInitialize (B);

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

  if (is_mpi_fine) {

    double bval = -4.0 * G * PI * MASS / (fine_h*fine_h*fine_h);

    HYPRE_SStructVectorAddToValues (B, part_fine, ind3[0], 0, &bval);
    HYPRE_SStructVectorAddToValues (B, part_fine, ind3[1], 0, &bval);
    HYPRE_SStructVectorAddToValues (B, part_fine, ind3[2], 0, &bval);
    HYPRE_SStructVectorAddToValues (B, part_fine, ind3[3], 0, &bval);
    HYPRE_SStructVectorAddToValues (B, part_fine, ind3[4], 0, &bval);
    HYPRE_SStructVectorAddToValues (B, part_fine, ind3[5], 0, &bval);
    HYPRE_SStructVectorAddToValues (B, part_fine, ind3[6], 0, &bval);
    HYPRE_SStructVectorAddToValues (B, part_fine, ind3[7], 0, &bval);
  }

  HYPRE_SStructVectorAssemble (B);

  // Create the hypre matrix A

  HYPRE_SStructMatrix  A;
  HYPRE_SStructMatrixCreate (MPI_COMM_WORLD, graph, &A);
  HYPRE_SStructMatrixSetObjectType (A,HYPRE_SSTRUCT);
  HYPRE_SStructMatrixInitialize (A);

  // Initialize A

  // Stencil entries

  int nums[7] = {0,1,2,3,4,5,6};

  double * v0  = new double [N*N*N];
  double * vxp = new double [N*N*N];
  double * vxm = new double [N*N*N];
  double * vyp = new double [N*N*N];
  double * vym = new double [N*N*N];
  double * vzp = new double [N*N*N];
  double * vzm = new double [N*N*N];

  double coarse_h =       BOX / N;  // Coarse mesh width

  if (is_mpi_coarse) {

    //--------------------------------------------------
    // MATRIX STENCIL VALUES: GRID INTERIOR
    //--------------------------------------------------

    for (i=0; i<N*N*N; i++) {
      v0[i]  = -6*coarse_h;
      vxp[i] =  1*coarse_h;
      vxm[i] =  1*coarse_h;
      vyp[i] =  1*coarse_h;
      vym[i] =  1*coarse_h;
      vzp[i] =  1*coarse_h;
      vzm[i] =  1*coarse_h;
    }

    //--------------------------------------------------
    // MATRIX STENCIL VALUES: GRID BOUNDARY
    //--------------------------------------------------

//     for (i0=0; i0<N; i0++) {
//       for (i1=0; i1<N; i1++) {
//  	vxp[index(N-1,i0,i1,N)] = 0.0;
//  	vxm[index(  0,i0,i1,N)] = 0.0;
//  	vyp[index(i0,N-1,i1,N)] = 0.0;
//  	vym[index(i0,  0,i1,N)] = 0.0;
//  	vzp[index(i0,i1,N-1,N)] = 0.0;
//  	vzm[index(i0,i1,  0,N)] = 0.0;
//       }
//     }

    //--------------------------------------------------
    // SET HYPRE MATRIX STENCIL VALUES
    //--------------------------------------------------

    HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				     0,1,&nums[0],v0);
    HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				     0,1,&nums[1],vxp);
    HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				     0,1,&nums[2],vxm);
    HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				     0,1,&nums[3],vyp);
    HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				     0,1,&nums[4],vym);
    HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				     0,1,&nums[5],vzp);
    HYPRE_SStructMatrixSetBoxValues (A,part_coarse,lower_coarse,upper_coarse,
				     0,1,&nums[6],vzm);
  }

  if (is_mpi_fine) {

    //--------------------------------------------------
    // MATRIX STENCIL VALUES: GRID INTERIOR
    //--------------------------------------------------

    for (i=0; i<N*N*N; i++) {
      v0[i]  = -6*fine_h;
      vxp[i] =  1*fine_h;
      vxm[i] =  1*fine_h;
      vyp[i] =  1*fine_h;
      vym[i] =  1*fine_h;
      vzp[i] =  1*fine_h;
      vzm[i] =  1*fine_h;
    }

    //--------------------------------------------------
    // MATRIX STENCIL VALUES: GRID BOUNDARY
    //--------------------------------------------------

    for (i0=0; i0<N; i0++) {
      for (i1=0; i1<N; i1++) {
	vxp[index(N-1,i0,i1,N)] -= fine_h;
	v0 [index(N-1,i0,i1,N)] += fine_h;
	vxm[index(  0,i0,i1,N)] -= fine_h;
	v0 [index(  0,i0,i1,N)] += fine_h;
	vyp[index(i0,N-1,i1,N)] -= fine_h;
	v0 [index(i0,N-1,i1,N)] += fine_h;
	vym[index(i0,  0,i1,N)] -= fine_h;
	v0 [index(i0,  0,i1,N)] += fine_h;
	vzp[index(i0,i1,N-1,N)] -= fine_h;
	v0 [index(i0,i1,N-1,N)] += fine_h;
	vzm[index(i0,i1,  0,N)] -= fine_h;
	v0 [index(i0,i1,  0,N)] += fine_h;
      }
    }

    //--------------------------------------------------
    // SET HYPRE MATRIX STENCIL VALUES
    //--------------------------------------------------

    HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				     0,1,&nums[0],v0);
    HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				     0,1,&nums[1],vxp);
    HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				     0,1,&nums[2],vxm);
    HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				     0,1,&nums[3],vyp);
    HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				     0,1,&nums[4],vym);
    HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				     0,1,&nums[5],vzp);
    HYPRE_SStructMatrixSetBoxValues (A,part_fine,lower_fine,upper_fine,
				     0,1,&nums[6],vzm);

    //--------------------------------------------------
    // MATRIX CLEAR STENCIL OVERLAP
    //--------------------------------------------------

    if (is_mpi_fine) {
      int refinements[3] = {2,2,2};

      HYPRE_SStructFACZeroCFSten (A,grid, part_fine, refinements);

      // This seems to do nothing:

      HYPRE_SStructFACZeroFCSten (A,grid, part_fine);

      // @@@ THIS SEEMS TO BREAK, BUT SEEMS TO BE IN hypre-solve @@@

      HYPRE_SStructFACZeroAMRMatrixData (A, part_coarse, refinements);
      
    }

  }

  delete [] v0;
  delete [] vxp;
  delete [] vxm;
  delete [] vyp;
  delete [] vym;
  delete [] vzp;
  delete [] vzm;

  PTRACE;

  //--------------------------------------------------
  // MATRIX NONSTENCIL ENTRIES
  //--------------------------------------------------

  // Declare counts for zones
  int * count_coarse = new int [N*N*N];
  int * count_fine   = new int [N*N*N];

  PTRACE;
  if (debug) {
    printf ("DEBUG count_coarse = %p\n",count_coarse);
    printf ("DEBUG count_fine   = %p\n",count_fine);
  }

  for (i=0; i<N*N*N; i++) {
    count_coarse[i] = 7;
    count_fine[i]   = 7; // stencils
  }

  //--------------------------------------------------
  // MATRIX ENTRIES: FINE-TO-CORSE 
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
  
	   // ind_fine computed above corresponds to coarse zone
	  if (face == 1)  ++ ind_fine[j0];

	  int    entry;
	  double value;

	  int icount;
	  int ishift = -index(N/r,N/r,N/r,N);
	  int num_entries = 1;
	  
	  double a = 2.0/3.0;

	  // 000 ---------------------------------------------
	  icount = ishift + index(ind_fine[0],ind_fine[1],ind_fine[2],N);
	  ASSERT_BOUND(0,icount,N*N*N);
	  // diagonal
	  entry = count_fine[icount]++;
	  value = a * fine_h;
	  HYPRE_SStructMatrixAddToValues 
	    (A, part_fine, ind_fine, 0, num_entries, &entry, &value);
	  // off-diagonal
	  entry = 0;
	  value = -value;
	  HYPRE_SStructMatrixAddToValues 
	    (A, part_fine, ind_fine, 0, num_entries, &entry, &value);
	  ++ ind_fine[j1]; 	  
	  // 010 ---------------------------------------------
	  icount = ishift + index(ind_fine[0],ind_fine[1],ind_fine[2],N);
	  ASSERT_BOUND(0,icount,N*N*N);
	  // diagonal
	  entry = count_fine[icount]++;
	  value = a * fine_h;
	  HYPRE_SStructMatrixAddToValues 
	    (A, part_fine, ind_fine, 0, num_entries, &entry, &value);
	  // off-diagonal
	  entry = 0;
	  value = -value;
	  HYPRE_SStructMatrixAddToValues 
	    (A, part_fine, ind_fine, 0, num_entries, &entry, &value);
	  ++ ind_fine[j2];	  
	  // 011 ---------------------------------------------
	  icount = ishift + index(ind_fine[0],ind_fine[1],ind_fine[2],N);
	  ASSERT_BOUND(0,icount,N*N*N);
	  // diagonal
	  entry = count_fine[icount]++;
	  value = a * fine_h;
	  HYPRE_SStructMatrixAddToValues 
	    (A, part_fine, ind_fine, 0, num_entries, &entry, &value);
	  // off-diagonal
	  entry = 0;
	  value = -value;
	  HYPRE_SStructMatrixAddToValues 
	    (A, part_fine, ind_fine, 0, num_entries, &entry, &value);
	  -- ind_fine[j1];	  
	  // 001 ---------------------------------------------
	  icount = ishift + index(ind_fine[0],ind_fine[1],ind_fine[2],N);
	  ASSERT_BOUND(0,icount,N*N*N);
	  // diagonal
	  entry = count_fine[icount]++;
	  value = a * fine_h;
	  HYPRE_SStructMatrixAddToValues 
	    (A, part_fine, ind_fine, 0, num_entries, &entry, &value);
	  // off-diagonal
	  entry = 0;
	  value = -value;
	  HYPRE_SStructMatrixAddToValues 
	    (A, part_fine, ind_fine, 0, num_entries, &entry, &value);
	  -- ind_fine[j2];	  
	  // 000 ---------------------------------------------

	}
      }
    }
  }

  //--------------------------------------------------
  // MATRIX ENTRIES: COARSE-TO-FINE
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
	  double value;

	  int icount = index(ind_coarse[0],ind_coarse[1],ind_coarse[2],N);
	  int num_entries = 1;

 	  for (int k=0; k<8; k++) {
 	    ASSERT_BOUND(0,icount,N*N*N);
 	    entry = count_coarse[icount]++;
 	    value = (1./8.) * coarse_h;
 	    HYPRE_SStructMatrixAddToValues 
 	      (A, part_coarse, ind_coarse, 0, num_entries, &entry, &value);
 	  }
	}
      }
    }
  }

  PTRACE;

  if (debug) printf ("DEBUG count_coarse = %p\n",count_coarse);
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
  //		      (A_, level_fine, igg3, 0, 1, &entry, &value);


  PTRACE;

  //------------------------------------------------------------
  // Solve the linear system
  //------------------------------------------------------------

  HYPRE_SStructSolver  solver;  // hypre solver
  int iter;
  double resid;

  if (use_fac) {
    
    //--------------------------------------------------
    // FAC SOLVER
    //--------------------------------------------------

    HYPRE_SStructFACCreate(MPI_COMM_WORLD, &solver);

    int num_parts = 2;
    HYPRE_SStructFACSetMaxLevels(solver,  num_parts);

    int parts[2] = {0,1};
    HYPRE_SStructFACSetPLevels(solver, num_parts, parts);

    int refinements[2][3] = {{2,2,2},{2,2,2}};
    HYPRE_SStructFACSetPRefinements(solver, num_parts, refinements);

    // solver parameters

    HYPRE_SStructFACSetNumPreRelax(solver,      4);
    HYPRE_SStructFACSetNumPostRelax(solver,     4);
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

    //--------------------------------------------------
    // BiCGSTAB SOLVER
    //--------------------------------------------------

    HYPRE_SStructBiCGSTABCreate(MPI_COMM_WORLD, &solver);
    HYPRE_SStructBiCGSTABSetLogging(solver, 1);
    HYPRE_SStructBiCGSTABSetup(solver, A, B, X);
    // SOLVE
    HYPRE_SStructBiCGSTABSolve(solver, A, B, X);
    HYPRE_SStructBiCGSTABGetNumIterations(solver, &iter);
    HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm(solver, &resid);
  }

  printf ("solver:     %s\n", use_fac ? "FAC" : "BiCGSTAB");
  printf ("iterations: %d\n",iter);
  printf ("residual:   %g\n",resid);

  if (use_fac) {
    HYPRE_SStructFACDestroy2(solver);
  } else {
    HYPRE_SStructBiCGSTABDestroy(solver);
  }

  //------------------------------------------------------------
  // WRITE MATRIX
  //------------------------------------------------------------

  HYPRE_SStructMatrixPrint ("A",A,1);

  //------------------------------------------------------------
  // WRITE SOLUTION
  //------------------------------------------------------------

  FILE *fp;

  //-------------------------------------------------------------
  // WRITE COARSE GRID
  //-------------------------------------------------------------

  fp = fopen("X.0","w");
  fprintf (fp,"Grid\n");
  fprintf (fp,"   id             0\n");
  fprintf (fp,"   parent id      -1\n");
  fprintf (fp,"   processor      %d\n",mpi_rank);
  fprintf (fp,"   lower position %d %d %d\n",
	   lower_coarse[0],lower_coarse[1],lower_coarse[2]);
  fprintf (fp,"   upper position %d %d %d\n",
	   upper_coarse[0],upper_coarse[1],upper_coarse[2]);
  fprintf (fp,"   lower index    0 0 0\n");
  fprintf (fp,"   zones          %d %d %d\n",N,N,N);
  fprintf (fp,"   level          0\n");

  double * x_coarse = new double [N*N*N];
  HYPRE_SStructVectorGetBoxValues (X,part_coarse,lower_coarse,upper_coarse,0,x_coarse);  
  for (i0=0; i0<N; i0++) {
    for (i1=0; i1<N; i1++) {
      for (i2=0; i2<N; i2++) {
	i = i0 + N*(i1 + N*i2);
	fprintf (fp,"%d %d %d %g\n",i0,i1,i2,x_coarse[i]);
      }
    }
  }
  delete [] x_coarse;

  fclose(fp);

  //-------------------------------------------------------------
  // WRITE FINE GRID
  //-------------------------------------------------------------

  fp = fopen("X.1","w");
  fprintf (fp,"Grid\n");
  fprintf (fp,"   id             1\n");
  fprintf (fp,"   parent id      0\n");
  fprintf (fp,"   processor      %d\n",mpi_rank);
  fprintf (fp,"   lower position %d %d %d\n",
	   lower_fine[0],lower_fine[1],lower_fine[2]);
  fprintf (fp,"   upper position %d %d %d\n",
	   upper_fine[0],upper_fine[1],upper_fine[2]);
  fprintf (fp,"   lower index    %d %d %d\n",
	   lower_fine[0],lower_fine[1],lower_fine[2]);
  fprintf (fp,"   zones          %d %d %d\n",N,N,N);
  fprintf (fp,"   level          1\n");

  double * x_fine = new double [N*N*N];
  HYPRE_SStructVectorGetBoxValues (X,part_fine,lower_fine,upper_fine,0,x_fine);  
  for (int i0=0; i0<N; i0++) {
    for (int i1=0; i1<N; i1++) {
      for (int i2=0; i2<N; i2++) {
	i = i0 + N*(i1 + N*i2);
	fprintf (fp,"%d %d %d %g\n",i0,i1,i2,x_fine[i]);
      }
    }
  }
  delete [] x_fine;

  fclose(fp);

  //------------------------------------------------------------
  // EXIT
  //------------------------------------------------------------

  printf ("Finished!\n");
  MPI_Finalize();
}

