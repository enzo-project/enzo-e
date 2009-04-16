//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Problem class source file

/**
 * 
 * @file      problem.cpp
 * @brief     Implementation of the Problem class
 * @author    James Bordner
 * @bug       none
 *
 * $Id$
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <map>
#include <string>
#include <vector>

#include "HYPRE_sstruct_ls.h"

#include "newgrav-hypre-solve.h"


//----------------------------------------------------------------------

const int debug = 0;

//----------------------------------------------------------------------

#include "newgrav-scalar.h"
#include "newgrav-constants.h"
#include "newgrav-error.h"
#include "newgrav-parameters.h"
#include "newgrav-point.h"
#include "newgrav-faces.h"
#include "newgrav-mpi.h"
#include "newgrav-domain.h"
#include "newgrav-grid.h"
#include "newgrav-level.h"
#include "newgrav-hierarchy.h"
#include "newgrav-problem.h"

const int trace = 0;

//======================================================================

Problem::Problem () throw ()
{
  //
}

//----------------------------------------------------------------------

Problem::~Problem () throw ()
{
  //
}

//----------------------------------------------------------------------

Problem::Problem (const Problem & p) throw ()
{
  domain_    = p.domain_;
  points_    = p.points_;
  hierarchy_ = p.hierarchy_;
}

//----------------------------------------------------------------------

Problem & Problem::operator = (const Problem & p) throw ()
{
  domain_    = p.domain_;
  points_    = p.points_;
  hierarchy_ = p.hierarchy_;
  return *this;
}

//----------------------------------------------------------------------

void Problem::read (std::string filename) throw ()
{

  bool boundary_defined = false;

  parameters_.read(filename);
  ItParameters itp (parameters_);

  while (itp++) {

    std::string key   = itp.key();
    std::string value = itp.value();


    if (key == "dimension") {

      int d = atoi(value.c_str());
      hierarchy_.set_dim(d);
      Point::set_dim(d);

      // Domain

    } else if (key == "domain") {

      domain_.input (value);

      // Grid ...

    } else if (key == "grid") {

      hierarchy_.insert_grid(new Grid(value));

    } else if (key == "point") {

      points_.push_back(new Point(value));      

    } else if (key == "boundary") {
      boundary_defined = true;
    }
  }

  if (! boundary_defined ) {
    fprintf (stderr,"Input parameter 'boundary' must be defined\n");
    MPI_Finalize();
    exit(1);
  }
}

//----------------------------------------------------------------------

void Problem::print () throw ()
{
  int i;
  domain_.print();
  hierarchy_.print();
  for (i=0; i<num_points(); i++)  point(i).print();
}

//----------------------------------------------------------------------

void Problem::write (FILE *fp) throw ()
{
  if (fp == 0) fp = stdout;
  int i;
  domain_.write(fp);
  hierarchy_.write(fp);
  for (i=0; i<num_points(); i++)  point(i).write(fp);
}
