
/// Problem class source file

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

#include <map>
#include <string>
#include <vector>

#include "HYPRE_sstruct_ls.h"

#include "hypre-solve.hpp"

#include "scalar.hpp"
#include "parameters.hpp"
#include "point.hpp"
#include "faces.hpp"
#include "mpi.hpp"
#include "grid.hpp"
#include "level.hpp"
#include "domain.hpp"
#include "hierarchy.hpp"
#include "sphere.hpp"
#include "problem.hpp"

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
  spheres_   = p.spheres_;
  points_    = p.points_;
  hierarchy_ = p.hierarchy_;
}

//----------------------------------------------------------------------

Problem & Problem::operator = (const Problem & p) throw ()
{
  domain_    = p.domain_;
  spheres_   = p.spheres_;
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
      Sphere::set_dim(d);

      // Domain

    } else if (key == "domain") {

      domain_.read (value);

      // Grid ...

    } else if (key == "grid") {

      hierarchy_.insert_grid(new Grid(value));

      // Sphere ...

    } else if (key == "sphere") {

      spheres_.push_back(new Sphere(value));      

      // Point ...

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
  for (i=0; i<num_spheres(); i++) sphere(i).print();
  for (i=0; i<num_points(); i++)  point(i).print();
}

//----------------------------------------------------------------------

void Problem::write (FILE *fp) throw ()
{
  if (fp == 0) fp = stdout;
  int i;
  domain_.write(fp);
  hierarchy_.write(fp);
  for (i=0; i<num_spheres(); i++) sphere(i).write(fp);
  for (i=0; i<num_points(); i++)  point(i).write(fp);
}
