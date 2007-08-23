
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

#include "hypre-solve.hpp"

#include "scalar.hpp"
#include "point.hpp"
#include "sphere.hpp"
#include "faces.hpp"
#include "mpi.hpp"
#include "grid.hpp"
#include "level.hpp"
#include "hierarchy.hpp"
#include "domain.hpp"
#include "problem.hpp"
#include "hypre.hpp"

const int debug = 1;
const int trace = 0;

//======================================================================
// BEGIN MAIN
//======================================================================

int main(int argc, char **argv)
{

  // Determine executable name

  std::string exec_name (argv[0]);
  int size = exec_name.rfind("/");
  exec_name.replace(0,size+1,"");

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

    _TRACE_;

    Problem problem;
    _TRACE_;

    problem.read  (filename);

    _TRACE_;
    
    _TRACE_;
    Hierarchy & hierarchy = problem.hierarchy();

    // determine interconnections between grids

    hierarchy.init_grids();

    // categorize grid boundary zones according to levels of adjacent
    // or containing grids

    hierarchy.init_faces();
  
    if (debug) problem.print ();

    // --------------------------------------------------
    // Initialize the grid hierarchy
    // --------------------------------------------------

    _TRACE_;

    Hypre hypre;

    hypre.init_hierarchy (hierarchy,mpi);

    // --------------------------------------------------
    // Initialize the stencils
    // --------------------------------------------------
    
    _TRACE_;

    hypre.init_stencil (hierarchy);

    // --------------------------------------------------
    // Initialize the graph
    // --------------------------------------------------

    _TRACE_;
    hypre.init_graph (hierarchy);

    // --------------------------------------------------
    // Initialize the matrix A
    // --------------------------------------------------

    //    _TRACE_;
    //    hypre.init_matrix (hierarchy);

    // --------------------------------------------------
    // Initialize the right-hand-side vector b
    // --------------------------------------------------

    //    _TRACE_;
    //    hypre.init_rhs (hierarchy);

    // --------------------------------------------------
    // Solve the linear system Ax = b
    // --------------------------------------------------

    //    _TRACE_;
    //    hypre.solve (hierarchy);

    // --------------------------------------------------
    // Evaluate the solution
    // --------------------------------------------------

    _TRACE_;
    hypre.evaluate (hierarchy);

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

