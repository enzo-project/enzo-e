#include <mpi.h>
#include <stdio.h>
#include "HYPRE_sstruct_ls.h"

#define TRACE printf ("%d %s:%d\n",ip,__FILE__,__LINE__); fflush(stdout);

main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);

  int ip;
  MPI_Comm_rank (MPI_COMM_WORLD, &ip);

  //-----------------------------------------------
  // Create HYPRE grids

  HYPRE_SStructGrid g;

  TRACE;
  HYPRE_SStructGridCreate (MPI_COMM_WORLD, 3,2,&g);

  //-----------------------------------------------
  // Set HYPRE grid extents

  int g_lower[] = { 0, 0, 0};
  int g_upper[] = {15,15,15};

  TRACE;
  HYPRE_SStructGridSetExtents(g,ip,g_lower,g_upper);

  //-----------------------------------------------
  // Define HYPRE grid variables

  HYPRE_SStructVariable variable_types[] = { HYPRE_SSTRUCT_VARIABLE_CELL };

  TRACE;
  HYPRE_SStructGridSetVariables(g,ip,1,variable_types);

  //-----------------------------------------------
  // Set HYPRE grid neighbors

  int index_map[] = {0,1,2};
  if (ip==0) {
    int g1_exts[][3] = {{16,0,0},{16,15,15}};
    int n1_exts[][3] = {{ 0,0,0},{ 0,15,15}};
    TRACE;
    HYPRE_SStructGridSetNeighborBox(g,
				    0,g1_exts[0],g1_exts[1],
				    1,n1_exts[0],n1_exts[1],
				    index_map);
  } else {
    int g2_exts[][3] = {{-1,0,0},{-1,15,15}};
    int n2_exts[][3] = {{15,0,0},{15,15,15}};
    TRACE;
    HYPRE_SStructGridSetNeighborBox(g,
				    1,g2_exts[0],g2_exts[1],
				    0,n2_exts[0],n2_exts[1],
				    index_map);
  }

  //-----------------------------------------------
  // Assemble HYPRE grids

  TRACE;
  HYPRE_SStructGridAssemble (g);  // SIGSEGV

  TRACE;
  MPI_Barrier(MPI_COMM_WORLD);
  TRACE;
}
