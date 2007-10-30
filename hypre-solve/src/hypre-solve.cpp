//345678901234567890123456789012345678901234567890123456789012345678901234567890

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
#include "performance.hpp"
#include "point.hpp"
#include "sphere.hpp"
#include "faces.hpp"
#include "mpi.hpp"
#include "domain.hpp"
#include "grid.hpp"
#include "level.hpp"
#include "hierarchy.hpp"
#include "parameters.hpp"
#include "problem.hpp"
#include "hypre.hpp"

const int debug    = 0;
const int trace    = 0;

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

  if (debug) printf ("DEBUG %s:%d mpi (ip,np) = (%d,%d)\n",
		     __FILE__,__LINE__,mpi.ip(),mpi.np());

  Grid::set_mpi (mpi);

  // --------------------------------------------------
  // JBPERF initialization
  // --------------------------------------------------


  // ***************
  JBPERF_BEGIN("HS");
  JBPERF_START("0-total");
  // ***************

  if (argc==2) {

    std::string filename (argv[1]);

    // --------------------------------------------------
    // Problem initialization
    // --------------------------------------------------

    // create a new problem and read it in

    Problem problem;

    // ***************
    JBPERF_START("1-problem");
    // ***************

    problem.read  (filename);

    // ***************
    JBPERF_STOP("1-problem");
    // ***************

    // --------------------------------------------------
    // Initialize the hierarchy
    // --------------------------------------------------

    Hierarchy & hierarchy = problem.hierarchy();
    Grid::set_domain (problem.domain());
    Level::set_domain (problem.domain());

    // determine interconnections between grids

    // ***************
    JBPERF_START("2-grids");
    // ***************

    hierarchy.initialize(problem.domain(), mpi);

    // ***************
    JBPERF_STOP("2-grids");
    // ***************

    if (debug) problem.print ();

    // WARNING: problem-size returns the size of the root-grid, not the entire hierarchy

    JBPERF_GLOBAL("problem-size",problem.hierarchy().num_unknowns0());

    // --------------------------------------------------
    // Initialize hypre
    // --------------------------------------------------

    Hypre hypre;

    // ***************
    JBPERF_START("4-hypre-init");
    // ***************

    hypre.init_hierarchy (problem.parameters(),hierarchy,mpi);

    // Initialize the stencils
    
    hypre.init_stencil (hierarchy);

    // Initialize the graph

    hypre.init_graph (hierarchy);

    // Initialize the linear system

    hypre.init_linear (problem.parameters(),
		       hierarchy,
		       problem.points(),
		       problem.spheres());

    // ***************
    JBPERF_STOP("4-hypre-init");
    // ***************

    // --------------------------------------------------
    // Solve the linear system Ax = b
    // --------------------------------------------------

    // ***************
    JBPERF_START("5-hypre-solve");
    // ***************

    hypre.solve (problem.parameters(),hierarchy);

    // ***************
    JBPERF_STOP("5-hypre-solve");
    // ***************

    // --------------------------------------------------
    // Evaluate the solution
    // --------------------------------------------------

    hypre.evaluate (hierarchy);

    // --------------------------------------------------
    // jbPerf Finalize
    // --------------------------------------------------

    // ****************
    JBPERF_STOP("0-total");
    JBPERF_END("HS");
    // ***************

    // --------------------------------------------------
    // MPI Finalize
    // --------------------------------------------------

    mpi.barrier();


    if (mpi.ip() == 0) { printf ("End %s\n",exec_name.c_str()); fflush(stdout); }

    MPI_Finalize ();

    // MPI_Finalize() called by Mpi destructor

  } else {

    printf ("Usage: %s <input-file>\n",argv[0]);

  }

}


