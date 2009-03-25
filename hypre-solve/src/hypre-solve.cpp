//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Main driver for hypre AMR test solver

/**
 * @file      hypre-solve.cpp
 * @brief     General test problem driver
 * @author    James Bordner
 * @bug       none
 *
 * $Id$
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

//----------------------------------------------------------------------

const int debug    = 0;
const int trace    = 0;
const int barrier  = 0;

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

  pmpi = new Mpi (&argc,&argv);

  if (debug) printf ("DEBUG %s:%d mpi (ip,np) = (%d,%d)\n",
		     __FILE__,__LINE__,pmpi->ip(),pmpi->np());

  Grid::set_mpi (*pmpi);

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
    _TRACE_;

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
    _TRACE_;

    Hierarchy & hierarchy = problem.hierarchy();
    Grid::set_domain (problem.domain());
    Level::set_domain (problem.domain());

    // determine interconnections between grids

    // ***************
    JBPERF_START("2-grids");
    // ***************

    bool is_periodic = problem.parameters().value("boundary") == "periodic";
    hierarchy.initialize(problem.domain(), *pmpi,is_periodic);

    // ***************
    JBPERF_STOP("2-grids");
    // ***************

    if (debug) problem.print ();

    // WARNING: problem-size returns the size of the root-grid, not the entire hierarchy

    JBPERF_GLOBAL("problem-size",problem.hierarchy().num_unknowns0());
    JBPERF_GLOBAL("problem-levels",problem.hierarchy().num_levels());

    // --------------------------------------------------
    // Initialize hypre
    // --------------------------------------------------
    _TRACE_;

    Hypre hypre (hierarchy,problem.parameters());

    // ***************
    JBPERF_START("4-hypre-init");
    // ***************

    if (barrier) pmpi->barrier();
    hypre.init_hierarchy (*pmpi);

    // Initialize the stencils
    
    if (barrier) pmpi->barrier();
    hypre.init_stencil ();

    // Initialize the graph

    if (barrier) pmpi->barrier();
    hypre.init_graph ();

    // Initialize the elements of matrix A and vector B

    if (barrier) pmpi->barrier();
    hypre.init_elements (problem.points());

    // ***************
    JBPERF_STOP("4-hypre-init");
    // ***************

    // --------------------------------------------------
    // Solve the linear system A X = B
    // --------------------------------------------------
    _TRACE_;

    // ***************
    JBPERF_START("5-hypre-solve");
    // ***************

    if (barrier) pmpi->barrier();
    hypre.solve ();

    // ***************
    JBPERF_STOP("5-hypre-solve");
    // ***************

    // --------------------------------------------------
    // Evaluate the solution
    // --------------------------------------------------
    _TRACE_;

    if (barrier) pmpi->barrier();
    hypre.evaluate ();

    // --------------------------------------------------
    // jbPerf Finalize
    // --------------------------------------------------
    _TRACE_;

    // ****************
    JBPERF_STOP("0-total");
    JBPERF_END("HS");
    // ***************

    // --------------------------------------------------
    // MPI Finalize
    // --------------------------------------------------
    _TRACE_;

    //    pmpi->barrier();

    if (barrier) pmpi->barrier();

    if (pmpi->is_root()) { 
      bool success = true;
      // Residual too high
      int itmax     = atoi(problem.parameters().value("solver_itmax").c_str());
      double restol = atof(problem.parameters().value("solver_restol").c_str());
      if (hypre.residual() > restol) {
	printf ("Diverged %s: %g > %g\n",exec_name.c_str(),
		hypre.residual(),restol); fflush(stdout); 
	success = false;
      }
      // Iterations reached limit
      if (hypre.iterations() >= itmax) {
	printf ("Stalled %s: %d >= %d\n",exec_name.c_str(),
		hypre.iterations(),itmax); fflush(stdout); 
	success = false;
      }
      // Appears to have completed successfully
      if (success) {
	printf ("Success %s\n",exec_name.c_str()); fflush(stdout); 
      }
    }

    MPI_Finalize();
    delete pmpi;

  } else {

    printf ("Usage: %s <input-file>\n",argv[0]);

  }

}


