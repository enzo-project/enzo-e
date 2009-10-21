const bool trace = false;

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

//----------------------------------------------------------------------
// INCLUDES
//----------------------------------------------------------------------

#include <stdio.h>
#include <mpi.h>

#include "cello.h"

#include "error.hpp"
#include "parameters.hpp"
#include "performance.hpp"

//----------------------------------------------------------------------
// ENUMERATIONS
//----------------------------------------------------------------------

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

  enum enum_data_type {
    data_type_unknown,
    data_type_alias,
    data_type_copy,
    data_type_alias_packed,
    data_type_copy_ghost,
    data_type_copy_ghost_packed,
    data_type_maximum = data_type_copy_ghost_packed
  };

//----------------------------------------------------------------------
// FUNCTION PROTOTYPES
//----------------------------------------------------------------------

void  init_block 
(
 Scalar * task_block,
 int bx, int by, int bz,
 int gx, int gy, int gz
 );

void compute_B 
(
 enum_data_type data_type,
 Scalar ** task_blocks,
 int sx, int sy, int sz,
 int bx, int by, int bz,
 int gx, int gy, int gz
 );

void update_block 
(
 Scalar * task_block,
 int bx, int by, int bz,
 int gx, int gy, int gz
 );


void check_range 
(
 int ip, 
 std::string name, 
 int var, 
 int min, 
 int max
 );

void ghost_exchange 
(
 Scalar * task_block,
 enum_data_type data_type,
 int bx, int by, int bz,
 int gx, int gy, int gz
 );

void ghost_exchange_packed 
(
 Scalar ** task_blocks,
 enum_data_type data_type,
 int sx, int sy, int sz,
 int bx, int by, int bz,
 int gx, int gy, int gz
 );

//======================================================================
// MAIN
//======================================================================

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

  TRACE;

  Parameters parameters;

  //--------------------------------------------------
  // Read in input file
  //--------------------------------------------------

  TRACE;
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

  TRACE;
  parameters.read(fp);

  //--------------------------------------------------
  // Declare parameters
  //--------------------------------------------------

  int mpi_type;

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
  
  TRACE;
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

  TRACE;
  if ((px * py * pz) != np) {
    if (ip==0) {
      fprintf (stderr,
	       "%s:%d Error: px*py*pz [%d %d %d] != np [%d]\n",
	       __FILE__,__LINE__,px,py,pz,np);
    }
    MPI_Abort(MPI_COMM_WORLD,1);
    exit(1);
  }
  TRACE;

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
  TRACE;

  int tx = nx / bx; // number of tasks in domain
  int ty = ny / by;
  int tz = nz / bz;


  TRACE;
  int sx = tx / px; // Number of tasks per processor
  int sy = ty / py;
  int sz = tz / pz;



  TRACE;
  if (px*sx != tx || py*sy != ty || pz*sz != tz) {
    if (ip==0) {
      fprintf (stderr,
	       "%s:%d Error: [tx ty tz] = [%d %d %d] not divisible by [px py pz] = [%d %d %d]\n",
	       __FILE__,__LINE__,tx,ty,tz,px,py,pz);
      fflush(stderr);
    }
    MPI_Abort(MPI_COMM_WORLD,1);
    exit(1);
  }

  TRACE;
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

  TRACE;

  //--------------------------------------------------
  // Initialize task arrays
  //--------------------------------------------------
  
  // Allocate processor-local array

  int ndx = qx + 2*gx;
  int ndy = qy + 2*gy;
  int ndz = qz + 2*gz;
  Scalar * array = new Scalar [ ndx * ndy * ndz ];

  // Allocate task array pointers

  Scalar ** task_blocks = new Scalar * [ sx * sy * sz ];

  int ga = gx + ndx * (gy + ndy * gz);

  for (int isz = 0; isz < sz; isz ++) {
    for (int isy = 0; isy < sy; isy ++) {
      for (int isx = 0; isx < sx; isx ++) {
	int is = isx + sx * (isy + sy * isz);
	int ia = (isx*bx) + ndx * ( (isy*by) + ndy * (isz*bz) );
	task_blocks [ is ] = & array [ ia + ga ];
      }
    }
  }

  TRACE;
  //--------------------------------------------------
  // Run test
  //--------------------------------------------------

  switch (mpi_type) {

  case mpi_type_B:

    for (int is=0; is<sx*sy*sz; is++) {
      init_block (task_blocks[is],bx,by,bz,gx,gy,gz);
    }

    compute_B ((enum_data_type)data_type,
	       task_blocks,
	       sx,sy,sz,
	       bx,by,bz,
	       gx,gy,gz);
    
    break;

  default:

    if (ip==0) {
      fprintf (stderr,
	       "%s:%d Error: unknown mpi_type %d\n\n",
	       __FILE__,__LINE__,mpi_type);
    }
    MPI_Abort (MPI_COMM_WORLD,1);

    break;

  }

  //--------------------------------------------------
  // Exit
  //--------------------------------------------------

  MPI_Barrier(MPI_COMM_WORLD);
  fflush(stdout);
  TRACE;
  MPI_Finalize();
}

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

//----------------------------------------------------------------------

void  init_block (Scalar * task_block,
		  int bx, int by, int bz,
		  int gx, int gy, int gz)
{
  // Compute block size with ghosts

  int nx = bx + 2*gx;
  int ny = by + 2*gy;
  int nz = bz + 2*gz;

  // Initialize blocks including ghosts

  for (int isz = -gz; isz < bz + gz; isz++) {
    for (int isy = -gy; isy < by + gy; isy++) {
      for (int isx = -gx; isx < bx + gx; isx++) {
	int is = isx + nx * (isy + ny* isz);
	task_block [is] = isx + isy + isz;
      }
    }
  }
}

//----------------------------------------------------------------------

void compute_B (enum_data_type data_type,
		Scalar ** task_blocks,
		int sx, int sy, int sz,
		int bx, int by, int bz,
		int gx, int gy, int gz)
{

  // packed ghost zone exchange

  ghost_exchange_packed (task_blocks,
			 data_type,
			 sx,sy,sz,
			 bx,by,bz,
			 gx,gy,gz);

  // loop over blocks (isx,isy,isz)

  for (int isz = 0; isz < sz; isz++) {
    for (int isy = 0; isy < sy; isy++) {
      for (int isx = 0; isx < sx; isx++) {

	int is = isx + sx * (isy + isy * isz);

	// non-packed ghost zone exchange

	ghost_exchange (task_blocks[is],data_type,bx,by,bz,gx,gy,gz);


	// update the block

	update_block (task_blocks[is],bx,by,bz,gx,gy,gz);

      }
    }
  }
}

//----------------------------------------------------------------------

void update_block (Scalar * task_block,
		   int bx, int by, int bz,
		   int gx, int gy, int gz)
{
  // Compute block size with ghosts

  int nx = bx + 2*gx;
  int ny = by + 2*gy;
  int nz = bz + 2*gz;

  // Update blocks including ghosts

  int dx = 1;
  int dy = bx;
  int dz = bx*by;
  Scalar d = 1.0 / 6.0;

  for (int isz = 0; isz < bz; isz++) {
    for (int isy = 0; isy < by; isy++) {
      for (int isx = 0; isx < bx; isx++) {
	int is = isx + nx * (isy + ny* isz);
	task_block [is] = 
	  d * (task_block [is + dx] + task_block [is - dx] +
	       task_block [is + dy] + task_block [is - dy] +
	       task_block [is + dz] + task_block [is - dz]);
      }
    }
  }

}

void ghost_exchange 
(
 Scalar * task_block,
 enum_data_type data_type,
 int bx, int by, int bz,
 int gx, int gy, int gz
 )
{

  switch (data_type) {
  case data_type_alias:
    break;
  case data_type_copy:
    break;
  case data_type_copy_ghost:
    break;
  case data_type_alias_packed:
  case data_type_copy_ghost_packed:
    // Skip
    break;
  case data_type_unknown:
  default:
    // Error
    break;
  }
}

void ghost_exchange_packed 
(
 Scalar ** task_blocks,
 enum_data_type data_type,
 int sx, int sy, int sz,
 int bx, int by, int bz,
 int gx, int gy, int gz)
{
  switch (data_type) {
  case data_type_alias_packed:
    break;
  case data_type_copy_ghost_packed:
    break;
  case data_type_alias:
  case data_type_copy:
  case data_type_copy_ghost:
    // Skip
    break;
  case data_type_unknown:
  default:
    // Error
    break;
  }
}

