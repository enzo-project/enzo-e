#include <stdio.h>
#include <mpi.h>

#include "cello.h"

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
//      1   B    Blocking send/receive
//      2   I    Non-blocking send/receive
//      3   BB   Buffered blocking send/receive
//      4   BI   Buffered non-blocking send/receive
//      5   GF   One-sided get with fence synchronization
//      6   G4   One-sided get with start / complete / post / wait
//      7   GL   One-sided get with lock synchronization
//      8   PF   One-sided get with fence synchronization
//      9   P4   One-sided get with start / complete / post / wait
//     10   PL   One-sided get with lock synchronization
//
// data_type          int
//
//      1   In-place per-block
//      3   Copy per-block
//      2   In-place multi-block
//      4   Ghost-copy per-block
//      5   Ghost-copy multi-block
//
// grid_size    [nx, ny, nz]
//
// block_size   [bx, by, bz]
//
// proc_size    [px, py, pz]
//
// ghost_depth  [gx, gy, gz]
//
// levels           [cores, cpus, nodes] (optional: for task ordering)
//
//----------------------------------------------------------------------

void check_range (int ip,
		  std::string name,
		  int var,
		  int min,
		  int max)
{	  
  if ( ! (min <= var && var <= max)) {
    if (ip==0) {
      fprintf (stderr,"%s = %d is out of range [%d,%d]\n",
	       name.c_str(),var,min,max);
    }
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}


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

  enum enum_mpi_type {
    mpi_type_unknown,
    mpi_type_B,
    mpi_type_I,
    mpi_type_BB,
    mpi_type_BI,
    mpi_type_GF,
    mpi_type_G4,
    mpi_type_GL,
    mpi_type_PF,
    mpi_type_P4,
    mpi_type_PL,
    mpi_type_maximum = mpi_type_PL
  };
  int mpi_type;

  enum enum_data_type {
    data_type_unknown,
    data_type_alias,
    data_type_copy,
    data_type_alias_packed,
    data_type_copy_ghost,
    data_type_copy_ghost_packed,
    data_type_maximum = data_type_copy_ghost_packed
  };
  int  data_type;

  int         nx,ny,nz;  // grid size
  int         bx,by,bz;  // block size
  int         px,py,pz;  // processor grid size

  int         gx,gy,gz;  // ghost zone depth

  int         levels[3] = {0};

  //--------------------------------------------------
  // Read parameters
  //--------------------------------------------------

  parameters.set_group ("Mpi_array");
  mpi_type   = parameters.value_integer ("mpi_type",0);
  data_type  = parameters.value_integer ("data_type",0);
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
  
  //--------------------------------------------------
  // Check parameters
  //--------------------------------------------------

  check_range (ip,"mpi_type",mpi_type,1,mpi_type_maximum);
  check_range (ip,"data_type",data_type,1,data_type_maximum);
  check_range (ip,"nx",nx,4,256);
  check_range (ip,"ny",ny,4,256);
  check_range (ip,"nz",nz,4,256);
  check_range (ip,"bx",bx,4,128);
  check_range (ip,"by",by,4,128);
  check_range (ip,"bz",bz,4,128);
  check_range (ip,"px",px,1,128);
  check_range (ip,"py",py,1,128);
  check_range (ip,"pz",pz,1,128);
  check_range (ip,"gx",gx,1,8);
  check_range (ip,"gy",gy,1,8);
  check_range (ip,"gz",gz,1,8);

  // Check processor count

  if ((px * py * pz) != np) {
    if (ip==0) {
      fprintf (stderr,
	       "%s:%d Error: px*py*pz [%d %d %d] != np [%d]\n",
	       __FILE__,__LINE__,px,py,pz,np);
    }
    MPI_Abort(MPI_COMM_WORLD,1);
    exit(1);
  }

  // Check number of processors * processor block size = problem size

  int qx = nx / px; // array size per processor 
  int qy = ny / py;
  int qz = nz / pz;
  
  if (px*qx != nx || py*qy != ny || pz*qz != nz) {
    if (ip==0) {
      fprintf (stderr,
	       "%s:%d Error: [nx ny nz] = [%d %d %d] not divisible by [px py pz] = [%d %d %d]\n",
	       __FILE__,__LINE__,nx,ny,nz,px,py,pz);
    }
    MPI_Abort(MPI_COMM_WORLD,1);
    exit(1);
  }

  // Check tasks per processor * task size = number of processors

  int tx = nx / bx; // number of tasks in domain
  int ty = ny / by;
  int tz = nz / bz;

  int tpx = tx / px; // Number of tasks per processor
  int tpy = ty / py;
  int tpz = tz / pz;

  if (tx*tpx != px || ty*tpy != py || tz*tpz != pz) {
    if (ip==0) {
      fprintf (stderr,
	       "%s:%d Error: [px py pz] = [%d %d %d] not divisible by [tx ty tz] = [%d %d %d]\n",
	       __FILE__,__LINE__,px,py,pz,tx,ty,tz);
    }
    MPI_Abort(MPI_COMM_WORLD,1);
    exit(1);
  }

  //--------------------------------------------------
  // Write parameters
  //--------------------------------------------------
  
  if (ip==0) {
    printf ("mpi_type   = %d\n",mpi_type);
    printf ("data_type  = %d\n",data_type);
    printf ("nx,ny,nz   = [%d %d %d]\n",nx,ny,nz);
    printf ("bx,by,bz   = [%d %d %d]\n",bx,by,bz);
    printf ("px,py,pz   = [%d %d %d]\n",px,py,pz);
    printf ("tx,ty,tz   = [%d %d %d]\n",tx,ty,tz);
    printf ("gx,gy,gz   = [%d %d %d]\n",gx,gy,gz);
    printf ("levels     = [%d %d %d]\n",levels[0],levels[1],levels[2]);
  }


  //--------------------------------------------------
  // Initialize task arrays
  //--------------------------------------------------
  
  // Allocate processor-local array

  int ndx = qx + 2*gx;
  int ndy = qy + 2*gy;
  int ndz = qz + 2*gz;
  Scalar * array = new Scalar [ ndx * ndy * ndz ];

  // Allocate task array pointers

  Scalar ** task_blocks = new Scalar * [ tpx * tpy * tpz ];

  int go = gx + ndx * (gy + ndy * gz);

  for (int itpz = 0; itpz < tpz; itpz ++) {
    for (int itpy = 0; itpy < tpy; itpy ++) {
      for (int itpx = 0; itpx < tpx; itpx ++) {
	int itp = itpx + tpx * (itpy + tpy * itpz);
	int ia = (itpx*bx) + ndx * ( (itpy*by) + ndy * (itpz*bz) );
	task_blocks [ itp ] = & array [ ia ];
      }
    }
  }

  //--------------------------------------------------
  // Run test
  //--------------------------------------------------


  //--------------------------------------------------
  // Exit
  //--------------------------------------------------

  MPI_Barrier(MPI_COMM_WORLD);
  fflush(stdout);
  MPI_Finalize();
}
