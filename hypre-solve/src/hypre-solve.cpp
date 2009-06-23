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


  if (argc==2) {

    std::string filename (argv[1]);

    // --------------------------------------------------
    // Problem initialization
    // --------------------------------------------------

    // create a new problem and read it in

    Problem problem;

    problem.read  (filename);

    // --------------------------------------------------
    // Initialize the hierarchy
    // --------------------------------------------------

    Hierarchy & hierarchy = problem.hierarchy();
    Grid::set_domain (problem.domain());  // only needed by geomview viz
    Level::set_domain (problem.domain()); // only needed by geomview viz

    // determine interconnections between grids

    bool is_periodic = problem.parameters().value("boundary") == "periodic";
    hierarchy.initialize(problem.domain(), *pmpi,is_periodic);

    if (debug) problem.print ();

    // WARNING: problem-size returns the size of the root-grid, not the entire hierarchy

    // --------------------------------------------------
    // Initialize hypre
    // --------------------------------------------------

    Hypre hypre (hierarchy,problem.parameters());


    hypre.init_hierarchy (*pmpi);

    // Initialize the stencils
    
    hypre.init_stencil ();

    // Initialize the graph

    hypre.init_graph ();

    // Initialize the elements of matrix A and vector B

    hypre.init_elements (problem.points());

    // --------------------------------------------------
    // Solve the linear system A X = B
    // --------------------------------------------------

    hypre.solve ();

    // --------------------------------------------------
    // Evaluate the solution
    // --------------------------------------------------

    hypre.evaluate ();

    // --------------------------------------------------
    // jbPerf Finalize
    // --------------------------------------------------

    // --------------------------------------------------
    // MPI Finalize
    // --------------------------------------------------


    MPI_Finalize();
    delete pmpi;

  } else {

    printf ("Usage: %s <input-file>\n",argv[0]);

  }

}


