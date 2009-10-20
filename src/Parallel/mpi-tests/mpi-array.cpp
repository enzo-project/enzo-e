#include <stdio.h>
#include <mpi.h>

#include "error.hpp"
#include "performance.hpp"

//----------------------------------------------------------------------
//
// PURPOSE
// 
// This is a Jacobi method test problem for evaluating different MPI
// communication strategies on regular arrays
//
//----------------------------------------------------------------------
//
// COMMUNICATION METHOD VARIATIONS
//
// B    Blocking send/receive
// I    Non-blocking send/receive
// BB   Buffered blocking send/receive
// BI   Buffered non-blocking send/receive
// GF   One-sided get with fence synchronization
// G4   One-sided get with start / complete / post / wait synchronization
// GL   One-sided get with lock synchronization
// PF   One-sided get with fence synchronization
// P4   One-sided get with start / complete / post / wait synchronization
// PL   One-sided get with lock synchronization
// 
//----------------------------------------------------------------------
//
// DATA LAYOUT VARIATIONS
//
// 1 in-place per-block
// 2 in-place multi-block
// 3 full block copy per-block
// 4 ghost-only copy per-block
// 5 ghost-only copy multi-block
//
//----------------------------------------------------------------------
//
// PARAMETERS
//
//  parameter    example
// ---------    ---------
// domain size    nx,ny,nz   [1024,1024,1024]
// block size     bx,by,bz   [   8,   8,   8]
// task grid      tx,ty,tz   [ 128, 128, 128]
// processor grid px,py,pz   [   8,   8,   8]
//
// ghost depth    1,2,3
// task ordering  natural         ipx + npx*(ipy + npy*ipz)
//                hilbert
//                hierarchical-natural
//                hierarchical-hilbert
//
//----------------------------------------------------------------------
//
// ISSUES
// 
//   * performance
//   * resiliency
//
//----------------------------------------------------------------------



int main (int argc, char ** argv)
{
  // Initialize MPI

  int np,ip;
  MPI_Init(&argc,&argv);
  MPI_Comm_size (MPI_COMM_WORLD, &np);
  MPI_Comm_rank (MPI_COMM_WORLD, &ip);

  MPI_Status status;

  if (ip==0) printf ("num_sends = %d\n",num_sends);
  MPI_Barrier(MPI_COMM_WORLD);
  if (ip == 0) {
    Timer timer;
    for (int count = 0; count < 3; count ++) {
    for (int ip2=1; ip2 < np; ip2++) {
      timer.start();
      for (int i=0; i<num_sends; i++) {
	MPI_Send(&value, 1, MPI_FLOAT, ip2, ip2*num_sends+i, MPI_COMM_WORLD);
	MPI_Recv(&value, 1, MPI_FLOAT, ip2, ip2*num_sends+i, MPI_COMM_WORLD, &status);
      }
      timer.stop();
      printf ("%d %f\n",ip2,timer.value());
      timer.clear();
    }
    }

  } else {
    for (int count = 0; count < 3; count ++) {
    for (int i=0; i<num_sends; i++) {
      MPI_Recv(&value, 1, MPI_FLOAT, 0, ip*num_sends+i, MPI_COMM_WORLD, &status);
      MPI_Send(&value, 1, MPI_FLOAT, 0, ip*num_sends+i, MPI_COMM_WORLD);
    }
    }
  }
  fflush(stdout);
  MPI_Finalize();
}
