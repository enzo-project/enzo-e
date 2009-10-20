#include <stdio.h>
#include <mpi.h>

#include "error.hpp"
#include "parameters.hpp"
#include "performance.hpp"

//----------------------------------------------------------------------
//
// PURPOSE
// 
// This is a parallel Jacobi-like method test problem for evaluating
// different MPI communication strategies on regular grids
//
//----------------------------------------------------------------------
//
// PARAMETERS
//
// mpi_type
//    
//     B    Blocking send/receive
//     I    Non-blocking send/receive
//     BB   Buffered blocking send/receive
//     BI   Buffered non-blocking send/receive
//     GF   One-sided get with fence synchronization
//     G4   One-sided get with start / complete / post / wait synchronization
//     GL   One-sided get with lock synchronization
//     PF   One-sided get with fence synchronization
//     P4   One-sided get with start / complete / post / wait synchronization
//     PL   One-sided get with lock synchronization
//
// data_type          int
//
//      1   In-place per-block
//      2   In-place multi-block
//      3   Full block copy per-block
//      4   Ghost-only copy per-block
//      5   Ghost-only copy multi-block
//
// grid_size    [nx, ny, nz]
//
// block_size   [bx, by, bz]
//
// proc_size    [px, py, pz]
//
// ghost_depth  [gx, gy, gz]
//
// task_order
//
//      1   natural
//      2   hilbert
//
// levels           [cores, cpus, nodes] (optional: for task ordering)
//
//----------------------------------------------------------------------


int main (int argc, char ** argv)
{
  //--------------------------------------------------
  // Initialize MPI
  //--------------------------------------------------

  int np,ip;
  MPI_Init(&argc,&argv);
  MPI_Comm_size (MPI_COMM_WORLD, &np);
  MPI_Comm_rank (MPI_COMM_WORLD, &ip);

  //--------------------------------------------------
  // Read parameter file
  //--------------------------------------------------

  if (argc != 2) {
    if (ip==0) {
      fprintf (stderr,"Usage: %s <filename>\n",argv[0]);
    }
    MPI_Abort(MPI_COMM_WORLD,1);
    exit(1);
  }

  Parameters parameters;

  //--------------------------------------------------
  // Read in input file
  //--------------------------------------------------

  FILE * fp = fopen (argv[1],"r");

  if (fp == NULL) {
    if (ip==0) {
      fprintf (stderr,
	       "%s:%d Error: filename %s not found\n",
	       __FILE__,__LINE__,argv[1]);
    }
    MPI_Abort(MPI_COMM_WORLD,1);
    exit(1);
  }

  parameters.read(fp);

  //--------------------------------------------------
  // Declare parameters
  //--------------------------------------------------

  std::string mpi_type;  // B I BB BI GF G4 GL PF P4 PL

  int         data_type; // 1-5
  int         task_order; // 1: natural 2: hilbert

  int         nx,ny,nz;  // grid size
  int         bx,by,bz;  // block size
  int         px,py,pz;  // processor grid size
  int         tx,ty,tz;  // task grid size

  int         gx,gy,gz;  // ghost zone depth

  int         levels[3] = {0};

  //--------------------------------------------------
  // Read parameters
  //--------------------------------------------------

  parameters.set_group ("Mpi_array");
  mpi_type   = parameters.value_string  ("mpi_type","B");
  data_type  = parameters.value_integer ("data_type",1);
  task_order = parameters.value_integer ("task_order",1);
  nx         = parameters.list_value_integer (0,"grid_size",16);
  ny         = parameters.list_value_integer (1,"grid_size",16);
  nz         = parameters.list_value_integer (2,"grid_size",16);
  bx         = parameters.list_value_integer (0,"block_size",8);
  by         = parameters.list_value_integer (1,"block_size",8);
  bz         = parameters.list_value_integer (2,"block_size",8);
  px         = parameters.list_value_integer (0,"proc_size",2);
  py         = parameters.list_value_integer (1,"proc_size",2);
  pz         = parameters.list_value_integer (2,"proc_size",1);
  gx         = parameters.list_value_integer (0,"ghost_depth",1);
  gy         = parameters.list_value_integer (1,"ghost_depth",1);
  gz         = parameters.list_value_integer (2,"ghost_depth",1);
  
  tx = nx / bx;
  ty = ny / by;
  tz = nz / bz;

  //--------------------------------------------------
  // Write parameters
  //--------------------------------------------------
  
    printf ("mpi_type   = %s\n",mpi_type.c_str());
    printf ("data_type  = %d\n",data_type);
    printf ("task_order = %d\n",task_order);
    printf ("nx,ny,nz   = [%d %d %d]\n",nx,ny,nz);
    printf ("bx,by,bz   = [%d %d %d]\n",bx,by,bz);
    printf ("px,py,pz   = [%d %d %d]\n",px,py,pz);
    printf ("tx,ty,tz   = [%d %d %d]\n",tx,ty,tz);
    printf ("gx,gy,gz   = [%d %d %d]\n",gx,gy,gz);
    printf ("levels     = [%d %d %d]\n",levels[0],levels[1],levels[2]);

    MPI_Barrier(MPI_COMM_WORLD);
  fflush(stdout);
  MPI_Finalize();
}
