//======================================================================
//
//        File: hypre-grav.C
//
//     Summary: Main driver for hypre AMR gravity solver
//
// Description:
//
//      Author: James Bordner <jobordner@ucsd.edu>
//
//        Date: 2007-03-26
//----------------------------------------------------------------------
//1345678901234567890123456789012345678901234567890123456789012345678901234567890
//======================================================================

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string>
#include <vector>
#include <mpi.h>

#include "HYPRE_sstruct_ls.h"

#include "scalar.hpp"
#include "point.hpp"
#include "sphere.hpp"
#include "discret.hpp"
#include "grid.hpp"
#include "level.hpp"
#include "hierarchy.hpp"
#include "mpi.hpp"
#include "problem.hpp"
#include "hypre.hpp"

const int debug = 1;

//======================================================================
// BEGIN MAIN
//======================================================================

int main(int argc, char **argv)
{

  int np,ip;

  // --------------------------------------------------
  // MPI initialization
  // --------------------------------------------------
  
  Mpi mpi (&argc,&argv);

  if (argc==2) {

    std::string filename (argv[1]);

    // --------------------------------------------------
    // Problem initialization
    // --------------------------------------------------

    // create a new problem and read it in

    Problem problem;

    problem.read  (filename);

    // Test printing, writing, and dumping problem to a file

    problem.print ();
    problem.write ();

    FILE *fp = fopen ("problem.out","w");
    problem.write (fp);
    fclose(fp);

    // --------------------------------------------------
    // Initialize the grid hierarchy
    // --------------------------------------------------

    Hypre hypre (mpi);

    if (debug) printf ("DEBUG ================================================\n");
    if (debug) printf ("DEBUG %s:%d\n",__FILE__,__LINE__);
    if (debug) printf ("DEBUG ================================================\n");

    hypre.init_hierarchy (problem.hierarchy());

    //    int i;
    //    for (i=0; i<Hierarchy.num_grids(); i++) {
    //      printf ("Initializing hypre grid %p\n",&Grid::grid(i));
    //      hypre.init_grid (Grid::grid(i));
    //    }

    // --------------------------------------------------
    // Initialize the stencils
    // --------------------------------------------------

    if (debug) printf ("DEBUG ================================================\n");
    if (debug) printf ("DEBUG %s:%d\n",__FILE__,__LINE__);
    if (debug) printf ("DEBUG ================================================\n");

    hypre.init_stencil (problem.hierarchy());

    // --------------------------------------------------
    // Initialize the graph
    // --------------------------------------------------

    if (debug) printf ("DEBUG ================================================\n");
    if (debug) printf ("DEBUG %s:%d\n",__FILE__,__LINE__);
    if (debug) printf ("DEBUG ================================================\n");

    hypre.init_graph (problem.hierarchy());

    // --------------------------------------------------
    // Initialize the matrix A
    // --------------------------------------------------

    // --------------------------------------------------
    // Initialize the right-hand-side vector b
    // --------------------------------------------------

    // --------------------------------------------------
    // Solve the linear system Ax = b
    // --------------------------------------------------

    // --------------------------------------------------
    // Evaluate the solution
    // --------------------------------------------------

    // --------------------------------------------------
    // MPI Finalize
    // --------------------------------------------------

    if (debug) printf ("DEBUG ================================================\n");
    if (debug) printf ("DEBUG %s:%d\n",__FILE__,__LINE__);
    if (debug) printf ("DEBUG ================================================\n");

    MPI_Finalize ();

  } else {

    printf ("Usage: %s <input-file>\n",argv[0]);

  }

  if (debug) printf ("DEBUG ================================================\n");
  if (debug) printf ("DEBUG %s:%d\n",__FILE__,__LINE__);
  if (debug) printf ("DEBUG ================================================\n");

}

