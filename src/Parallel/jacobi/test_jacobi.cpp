// $Id: test_Parallel.cpp 1694 2010-08-04 05:51:33Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_combine.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Aug  9 10:14:56 PDT 2010
/// @brief    Prototype test program for combining charm++, MPI, and OMP


// #include <stdio.h>
// #include <stdlib.h>

#include "jacobi.hpp"

// ------------------------------ INCLUDES ------------------------------
#ifdef CONFIG_USE_MPI
#   include <mpi.h>
#endif

#ifdef CONFIG_USE_CHARM
#   include "jacobi.decl.h"
#endif

#include "jacobi_Block.hpp"

// ------------------------------ MAIN ------------------------------

#ifdef CONFIG_USE_CHARM
CProxy_Main mainProxy;

class Main : public CBase_Main
{
public:
  Main(CkArgMsg* main) 
#else
    int main(int argc, char ** argv)
#endif

  {

  // Initialize default comm rank and size

  int ip = 0;
  int np = 1;

  // Determine actual comm rank and size

#ifdef CONFIG_USE_MPI
  MPI_Init (&ARGC,&ARGV);
  MPI_Comm_rank (MPI_COMM_WORLD,&ip);
  MPI_Comm_size (MPI_COMM_WORLD,&np);
#endif

#ifdef CONFIG_USE_CHARM
#endif

  // Initialize problem parameters

  if (ARGC != 3) {
    if (ip==0) {
      PRINTF ("\nUsage: %s N M\n\n",ARGV[0]);
      PRINTF ("   N: problem size = N*N*N\n");
      PRINTF ("   M: block size   = M*M*M\n\n");
    }
#ifdef CONFIG_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
#ifdef CONFIG_USE_CHARM
    CkExit();
#endif
  }    

  PRINTF ("(ip,np) = (%d,%d)\n",ip,np);

  int n = atoi(ARGV[1]);  // Problem size = n*n*n
  int m = atoi(ARGV[2]);  // Block size = m*m*m
  int nb = n/m;           // Number of blocks
  int n3 = n*n*n;
  int m3 = m*m*m;
  int nb3 = nb*nb*nb;

  PRINTF ("n=%d\n",n);

  // Create data blocks

  Block * block[nb];
  for (int i=0; i<nb; i++) {
    double xm = 0.0;
    double xp = 1.0;
    double ym = 0.0;
    double yp = 1.0;
    double zm = 0.0;
    double zp = 1.0;
    block[i] = new Block(m, xm,xp,ym,yp,zm,zp);
  }

  for (int i=0; i<nb; i++) {
    delete block[i];
  }

#ifdef CONFIG_USE_MPI
  MPI_Finalize();
#endif
#ifdef CONFIG_USE_CHARM
    CkExit();
#endif

  };

#ifdef CONFIG_USE_CHARM
};
#   include "jacobi.def.h"
#endif
