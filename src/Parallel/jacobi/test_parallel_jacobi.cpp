// See LICENSE_CELLO file for license and copyright information

/// @file     test_parallel_jacobi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Aug  9 10:14:56 PDT 2010
/// @brief    Prototype test program for CHARM++


#include <stdio.h>
#include <stdlib.h>

#include "parallel_jacobi.hpp"

// ------------------------------ MAIN ------------------------------

CProxy_Main proxy_main;

class Main : public CBase_Main {

private:

  // command line options
  double norm_;

  int problem_size_;
  int block_size_;
  int iteration_max_;

  int iteration_;
  int num_blocks_;
  jacobi::Counter * blocks_;
  CProxy_CharmPatch patches_;

public:
	   
  Main(CkArgMsg* main)
    
  {
    if (main->argc != 1 + 3) {

      CkPrintf ("\nUsage: %s problem_size  block_size  iterations\n\n",PARALLEL_ARGV[0]);
    
      CkExit();
    }    
 
    norm_ = 0.0;
    problem_size_  = atoi(main->argv[1]);
    block_size_    = atoi(main->argv[2]);
    iteration_max_ = atoi(main->argv[3]);

    iteration_     = 0;
    num_blocks_   = problem_size_/block_size_;
    blocks_        = new jacobi::Counter (num_blocks_*num_blocks_*num_blocks_);

    patches_       = 
      CProxy_CharmPatch::ckNew
      (num_blocks_,block_size_,iteration_max_, thisProxy, 
       num_blocks_,num_blocks_,num_blocks_);

    CkPrintf ("Evolve(%d)\n",iteration_);
    patches_.p_evolve();


    CkPrintf ("MAIN EXIT\n");

  };

  void p_next(int ix, int iy, int iz, double s2)
  {
    norm_ += s2;
    if (blocks_->remaining() == 0) {
      if (++iteration_ >= iteration_max_) {
	CkPrintf ("End computation %g time = %g\n",norm_,CkWallTimer());
	CkExitAfterQuiescence();
      } else {
	CkPrintf ("Evolve(%d) %g\n",iteration_,norm_);
	norm_ = 0.0;
	patches_.p_evolve();
      }
    }
  }
    
};

using namespace jacobi;
#include "test_parallel_jacobi.def.h"
