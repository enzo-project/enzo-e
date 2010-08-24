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
 
    // Problem size = n*n*n
    int problem_size = atoi(PARALLEL_ARGV[1]);  

    // Block size = m*m*m
    int block_size   = atoi(PARALLEL_ARGV[2]);  

    // Number of blocks per dimension
    int block_count  = problem_size/block_size;

    PARALLEL_INIT;

    CProxy_Patch patch = 
      CProxy_Patch::ckNew(block_count,block_count,block_count);

    patch.p_evolve(block_count,block_size);


    //    patch.print();


    //    PARALLEL_EXIT;
    CkPrintf ("MAIN EXIT\n");

  };
};

#include "test_jacobi.def.h"
