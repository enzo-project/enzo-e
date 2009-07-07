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

void test_grid();
void test_level();
void test_problem();
void test_hierarchy();
void print_result (std::string class_name, int bytes_begin, int bytes_end);

//======================================================================
// BEGIN MAIN
//======================================================================

int main(int argc, char **argv)
{

  // --------------------------------------------------
  // MPI initialization
  // --------------------------------------------------

  pmpi = new Mpi (&argc,&argv);

  int bytes_begin;
  int bytes_end;

  // Test Grid class

  MPI_Barrier(MPI_COMM_WORLD);
  bytes_begin = jbMem::bytes_;
  test_grid();
  bytes_end = jbMem::bytes_;
  print_result ("Grid",bytes_begin,bytes_end);
  MPI_Barrier(MPI_COMM_WORLD);

  // Test Level class

  MPI_Barrier(MPI_COMM_WORLD);
  bytes_begin = jbMem::bytes_;
  test_level();
  bytes_end = jbMem::bytes_;
  print_result ("Level",bytes_begin,bytes_end);
  MPI_Barrier(MPI_COMM_WORLD);

  // Test Problem class

  MPI_Barrier(MPI_COMM_WORLD);
  bytes_begin = jbMem::bytes_;
  test_problem();
  bytes_end = jbMem::bytes_;
  print_result ("Problem",bytes_begin,bytes_end);
  MPI_Barrier(MPI_COMM_WORLD);

  // Test Hierarchy class

  MPI_Barrier(MPI_COMM_WORLD);
  bytes_begin = jbMem::bytes_;
  test_hierarchy();
  bytes_end = jbMem::bytes_;
  print_result ("Hierarchy",bytes_begin,bytes_end);
  MPI_Barrier(MPI_COMM_WORLD);

  delete pmpi;
  MPI_Finalize();

}

void test_grid()
{
  const int num_grids = 10;
  int     id = 0;
  int     id_parent = -1;
  int     ip = 0; 
  Scalar xl[3] = {0.0, 0.0, 0.0};
  Scalar xu[3] = {1.0, 1.0, 1.0};
  int    il[3] = {0,0,0};
  int    n[3] = {32,32,32};

  for (int i=0; i<num_grids; i++) {

    Grid * grid = new Grid (id, id_parent, ip, xl, xu, il, n);

    delete grid;

  }

}

void test_level()
{
  int     id = 0;
  int     id_parent = -1;
  int     ip = 0; 
  Scalar xl[3] = {0.0, 0.0, 0.0};
  Scalar xu[3] = {1.0, 1.0, 1.0};
  int    il[3] = {0,0,0};
  int    n[3] = {32,32,32};

  // Create grid
  Grid * grid = new Grid (id, id_parent, ip, xl, xu, il, n);

  // Create level
  Level * level = new Level(0);

  // Insert grid in level
  level->insert_grid (*grid);

  // Delete grids in level
  ItLevelGridsAll itg (*level);
  while (Grid * grid = itg++) {
    delete grid;
  }

  // Delete level
  delete level;

}

void test_problem()
{
  Problem * problem = new Problem;
  problem->read ("in.test-mem");
  delete problem;
}

void test_hierarchy()
{
  int     id = 0;
  int     id_parent = -1;
  int     ip = 0; 
  Scalar xl[3] = {0.0, 0.0, 0.0};
  Scalar xu[3] = {1.0, 1.0, 1.0};
  int    il[3] = {0,0,0};
  int    n[3] = {32,32,32};

  Hierarchy * hierarchy = new Hierarchy();

  Scalar dl[3] = {-1.0,-1.0,-1.0};
  Scalar du[3] = {1.0,1.0,1.0};
  Domain * domain = new Domain (3,dl,du);
  Grid * grid = new Grid (id, id_parent, ip, xl, xu, il, n);

  hierarchy->insert_grid(grid);

  hierarchy->initialize(*domain,*pmpi, true);

  delete hierarchy;
  delete domain;

}

void print_result (std::string class_name, int bytes_begin, int bytes_end)
{

  // memory test passes if bytes allocated are same before and after
  bool pass = bytes_begin == bytes_end;

  // Set grid result string to be pass or fail
  std::string result = pass ? pass_string : fail_string;
  
  // Print result
  printf ("%s memory: %s %d %d\n",class_name.c_str(),result.c_str(),
	  bytes_begin,bytes_end);
}
