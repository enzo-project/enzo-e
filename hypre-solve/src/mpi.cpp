
/// Mpi class for basic communication routines and data

/**
 * 
 * 
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 *
 */

#include <stdio.h>
#include <mpi.h>
#include "mpi.hpp"

const int debug = 0;

//----------------------------------------------------------------------

Mpi::Mpi ()
  : comm_(0),
    np_(1),
    ip_(0)
{
}

//----------------------------------------------------------------------

Mpi::~Mpi ()
{
  // MPI_Finalize ();
}

//----------------------------------------------------------------------

Mpi::Mpi (int *argc, char ***argv)
{
  MPI_Init (argc,argv);

  comm_ = MPI_COMM_WORLD;

  MPI_Comm_size (comm_, &np_);
  MPI_Comm_rank (comm_, &ip_);

  if (debug) printf ("%s:%d  np=%d  ip=%d\n",__FILE__,__LINE__,np_,ip_);

}

//----------------------------------------------------------------------
