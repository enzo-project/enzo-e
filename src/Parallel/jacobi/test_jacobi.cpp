// $Id: test_Parallel.cpp 1694 2010-08-04 05:51:33Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_combine.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Aug  9 10:14:56 PDT 2010
/// @brief    Prototype test program for combining charm++, MPI, and OMP


#include <stdio.h>
#include <stdlib.h>

#include "counter.hpp"
#include "jacobi.hpp"

// ------------------------------ MAIN ------------------------------

CProxy_Main mainProxy;

class Main : public CBase_Main {

private:

  // command line parameters
  int problem_size_;
  int block_size_;
  int iteration_max_;

  int iteration_;
  Counter * blocks_;
  int block_count_;
  CProxy_Patch patches_;

public:
	   
  Main(CkArgMsg* main)
    
  {
    if (main->argc != 1 + 3) {

      CkPrintf ("\nUsage: %s problem_size  block_size  iterations\n\n",PARALLEL_ARGV[0]);
    
      CkExit();
    }    
 
    problem_size_  = atoi(main->argv[1]);
    block_size_    = atoi(main->argv[2]);
    iteration_max_ = atoi(main->argv[3]);

    iteration_     = 0;
    block_count_   = problem_size_/block_size_;
    blocks_        = new Counter (block_count_*block_count_*block_count_);

    patches_       = 
      CProxy_Patch::ckNew(block_count_,block_size_,thisProxy,
			  block_count_,block_count_,block_count_);

    CkPrintf ("Evolve(%d)\n",iteration_);
    patches_.p_evolve();


    CkPrintf ("MAIN EXIT\n");

  };

  void p_next()
  {
    if (blocks_->wait()) {
      if (++iteration_ >= iteration_max_) {
	CkPrintf ("End computation\n");
	CkExitAfterQuiescence();
      } else {
	CkPrintf ("Evolve(%d)\n",iteration_);
	patches_.p_evolve();
      }
    }
    
  }
};

#include "test_jacobi.def.h"
