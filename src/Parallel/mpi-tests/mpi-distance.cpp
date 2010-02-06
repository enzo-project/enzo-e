/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @bug       
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */

#include <stdio.h>
#include <mpi.h>

#include "error.hpp"
#include "performance.hpp"

// Estimate the latency from the root process to all other processes

int main (int argc, char ** argv)
{
  MPI_Status status;
  float value = 2;
  int np,ip;
  int num_sends = 1;
  if (argc > 1) {
    num_sends = atoi(argv[1]);
  }
    
  MPI_Init(&argc,&argv);
  MPI_Comm_size (MPI_COMM_WORLD, &np);
  MPI_Comm_rank (MPI_COMM_WORLD, &ip);
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
