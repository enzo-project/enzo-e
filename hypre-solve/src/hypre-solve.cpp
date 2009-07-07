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


  if (argc==2) {

    // --------------------------------------------------
    // JBPERF initialization
    // --------------------------------------------------

    jbPerf.begin("EL");

    std::string filename (argv[1]);

    // --------------------------------------------------
    // Problem initialization
    // --------------------------------------------------

    // create a new problem and read it in

    jbPerf.start("problem");
    Problem problem;
    jbPerf.stop("problem");

    jbPerf.start("problem-read");
    problem.read  (filename);
    jbPerf.stop("problem-read");

    // --------------------------------------------------
    // Initialize the hierarchy
    // --------------------------------------------------

    Hierarchy & hierarchy = problem.hierarchy();
    Grid::set_domain (problem.domain());  // only needed by geomview viz
    Level::set_domain (problem.domain()); // only needed by geomview viz

    // determine interconnections between grids

    bool is_periodic = problem.parameters().value("boundary") == "periodic";
    jbPerf.start("hierarchy-initialize");
    hierarchy.initialize(problem.domain(), *pmpi,is_periodic);
    jbPerf.stop("hierarchy-initialize");

    if (debug) problem.print ();

    // WARNING: problem-size returns the size of the root-grid, not the entire hierarchy

    // --------------------------------------------------
    // Initialize hypre
    // --------------------------------------------------

    Hypre hypre (hierarchy,problem.parameters());


    jbPerf.start("hypre-init-hierarchy");
    hypre.init_hierarchy (*pmpi);
    jbPerf.stop("hypre-init-hierarchy");

    // Initialize the stencils
    
    jbPerf.start("hypre-init-stencil");
    hypre.init_stencil ();
    jbPerf.stop("hypre-init-stencil");

    // Initialize the graph

    jbPerf.start("hypre-init-graph");
    hypre.init_graph ();
    jbPerf.stop("hypre-init-graph");

    // Initialize the elements of matrix A and vector B

    jbPerf.start("hypre-init-elements");
    hypre.init_elements (problem.points());
    jbPerf.stop("hypre-init-elements");

    // --------------------------------------------------
    // Solve the linear system A X = B
    // --------------------------------------------------

    jbPerf.start("hypre-solve");
    hypre.solve ();
    jbPerf.stop("hypre-solve");

    // --------------------------------------------------
    // Evaluate the solution
    // --------------------------------------------------

    jbPerf.start("hypre-evaluate");
    hypre.evaluate ();
    jbPerf.stop("hypre-evaluate");

    jbPerf.start("hypre-delete");
    // --------------------------------------------------
    // jbPerf Finalize
    // --------------------------------------------------

    // --------------------------------------------------
    // MPI Finalize
    // --------------------------------------------------


    MPI_Finalize();
    delete pmpi;

    jbPerf.end("EL");

  } else {

    printf ("Usage: %s <input-file>\n",argv[0]);

  }

}


