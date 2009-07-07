//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Main driver to test for memory leaks in hypre-solve

/**
 * @file      test-mem.cpp
 * @brief     Main driver to test for memory leaks in hypre-solve
 * @author    James Bordner
 * @bug       none
 *
 * $Id: hypre-solve.cpp 686 2009-06-23 20:27:24Z bordner $
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

#include <map>
#include <string>
#include <vector>

#include "HYPRE_sstruct_ls.h"

#include "newgrav-hypre-solve.h"
#include "jbMem.h"

//----------------------------------------------------------------------

const int debug    = 0;
const int trace    = 0;

//----------------------------------------------------------------------

#include "newgrav-scalar.h"
#include "newgrav-constants.h"
#include "newgrav-error.h"
#include "newgrav-performance.h"
#include "newgrav-point.h"
#include "newgrav-faces.h"
#include "newgrav-mpi.h"
#include "newgrav-domain.h"
#include "newgrav-grid.h"
#include "newgrav-level.h"
#include "newgrav-hierarchy.h"
#include "newgrav-parameters.h"
#include "newgrav-problem.h"
#include "newgrav-hypre.h"

//======================================================================

#define PRINT_MEM \
  printf ("%s:%d jbMem bytes = %lld\n", \
  __FILE__,__LINE__,jbMem::bytes_);

//======================================================================

const char * pass_string = " [ Pass ]";
const char * fail_string = " [ FAIL ]";

//======================================================================
// BEGIN MAIN
//======================================================================

int main(int argc, char **argv)
{

  // --------------------------------------------------
  // MPI initialization
  // --------------------------------------------------

  pmpi = new Mpi (&argc,&argv);

  const int num_grids = 10;
  int     id = 0;
  int     id_parent = -1;
  int     ip = 0; 
  Scalar xl[3] = {0.0, 0.0, 0.0};
  Scalar xu[3] = {1.0, 1.0, 1.0};
  int    il[3] = {0,0,0};
  int    n[3] = {32,32,32};
  int bytes_begin = jbMem::bytes_;
  for (int i=0; i<num_grids; i++) {

    Grid * grid = new Grid (id, id_parent, ip, xl, xu, il, n);

    delete grid;

  }
  int bytes_end = jbMem::bytes_;

  // grid memory test passes if bytes allocated are same before and after
  bool grid_pass = bytes_begin == bytes_end;

  // Set grid result string to be pass or fail
  std::string grid_result = grid_pass ? pass_string : fail_string;
  
  // Print grid resultA
  printf ("Grid memory: %s\n",grid_result.c_str());
    

}


