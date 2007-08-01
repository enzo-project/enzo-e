//======================================================================
//
//        File: mpi.cpp
//
//     Summary: Communication routines
//
// Description: Basic communication routines and data
//
//      Author: James Bordner <jobordner@ucsd.edu>
//
//        Date: 2007-04-10
//
//======================================================================

#include <stdio.h>
#include <mpi.h>
#include "mpi.hpp"

const int debug = 0;

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
