
/// Main driver for hypre AMR test solver

/**
 * 
 * 
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 *
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

#include "scalar.hpp"
#include "point.hpp"
#include "sphere.hpp"
#include "discret.hpp"
#include "mpi.hpp"
#include "grid.hpp"
#include "level.hpp"
#include "hierarchy.hpp"
#include "problem.hpp"
#include "hypre.hpp"

const int debug = 0;

//======================================================================
// BEGIN MAIN
//======================================================================

int main(int argc, char **argv)
{

  // Determine executable name

  std::string exec_name (argv[0]);
  int size = exec_name.rfind("/");
  exec_name.replace(0,size+1,"");

  int np,ip;

  // --------------------------------------------------
  // MPI initialization
  // --------------------------------------------------

  if (debug) printf ("DEBUG %s:%d\n",__FILE__,__LINE__);

  Mpi mpi (&argc,&argv);

  if (argc==2) {

    std::string filename (argv[1]);

    // --------------------------------------------------
    // Problem initialization
    // --------------------------------------------------

    // create a new problem and read it in

    if (debug) printf ("DEBUG %s:%d\n",__FILE__,__LINE__);

    Problem problem;

    problem.read  (filename);

    // --------------------------------------------------
    // Initialize the grid hierarchy
    // --------------------------------------------------

    if (debug) printf ("DEBUG %s:%d\n",__FILE__,__LINE__);

    Hypre hypre;

    hypre.init_hierarchy (problem.hierarchy(),mpi);

    if (debug) printf ("DEBUG %s:%d\n",__FILE__,__LINE__);

    // --------------------------------------------------
    // Initialize the stencils
    // --------------------------------------------------

    //    hypre.init_stencil (problem.hierarchy());

    // --------------------------------------------------
    // Initialize the graph
    // --------------------------------------------------

    //    hypre.init_graph (problem.hierarchy());

    // --------------------------------------------------
    // Initialize the matrix A
    // --------------------------------------------------

    //    hypre.init_matrix (problem.hierarchy());

    // --------------------------------------------------
    // Initialize the right-hand-side vector b
    // --------------------------------------------------

    //    hypre.init_rhs (problem.hierarchy());

    // --------------------------------------------------
    // Solve the linear system Ax = b
    // --------------------------------------------------

    // hypre.solve (problem.hierarchy());

    // --------------------------------------------------
    // Evaluate the solution
    // --------------------------------------------------

    //    hypre.evaluate (problem.hierarchy());

    // --------------------------------------------------
    // MPI Finalize
    // --------------------------------------------------

    mpi.barrier();

    if (mpi.ip() == 0) { printf ("End %s\n",exec_name.c_str()); fflush(stdout); }


    MPI_Finalize ();

  } else {

    printf ("Usage: %s <input-file>\n",argv[0]);

  }


}

