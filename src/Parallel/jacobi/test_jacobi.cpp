// $Id: test_Parallel.cpp 1694 2010-08-04 05:51:33Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_combine.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Aug  9 10:14:56 PDT 2010
/// @brief    Prototype test program for combining charm++, MPI, and OMP


#include <stdio.h>
#include <stdlib.h>

#include "parallel.def"

#include "jacobi.hpp"

// ------------------------------ MAIN ------------------------------

CProxy_Main mainProxy;

class Main : public CBase_Main {

public:
	   
  Main(CkArgMsg* main)


  {

    if (PARALLEL_ARGC != 3) {

      PARALLEL_PRINTF ("\nUsage: %s N M\n\n",PARALLEL_ARGV[0]);
      PARALLEL_PRINTF ("   N: problem size = N*N*N\n");
      PARALLEL_PRINTF ("   M: block size   = M*M*M\n\n");
    
      PARALLEL_EXIT;
    }    
 
    int n = atoi(PARALLEL_ARGV[1]);  // Problem size = n*n*n
    int m = atoi(PARALLEL_ARGV[2]);  // Block size = m*m*m
    int nb = n/m;           // Number of blocks per dimension

    PARALLEL_INIT;

    CProxy_Patch patch = CProxy_Patch::ckNew(nb,nb,nb);

    patch.advance(m);


    //    patch.print();


    //    PARALLEL_EXIT;

  };
};

#include "test_jacobi.def.h"
