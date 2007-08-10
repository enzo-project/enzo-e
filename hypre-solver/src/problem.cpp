
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

#include "scalar.hpp"
#include "point.hpp"
#include "discret.hpp"
#include "mpi.hpp"
#include "grid.hpp"
#include "level.hpp"
#include "hierarchy.hpp"
#include "sphere.hpp"
#include "problem.hpp"

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
  spheres_   = p.spheres_;
  points_    = p.points_;
  hierarchy_ = p.hierarchy_;
}

//----------------------------------------------------------------------

Problem & Problem::operator = (const Problem & p) throw ()
{
  spheres_   = p.spheres_;
  points_    = p.points_;
  hierarchy_ = p.hierarchy_;
  return *this;
}

//----------------------------------------------------------------------

void Problem::read (std::string filename) throw ()
{
  FILE *fp = fopen(filename.c_str(),"r");
  char buffer[BUFFER_LENGTH];
  int i;

  // Clear the buffer
  for (i=0; i<BUFFER_LENGTH; i++) buffer[i]=0;

  while (readline_ (fp,buffer,BUFFER_LENGTH)) {

    char obj[BUFFER_LENGTH];
    for (i=0; i<BUFFER_LENGTH; i++) obj[i]=0;
    sscanf(buffer,"%s",obj);
    const char * args = buffer + strlen(obj) + 1;

    // dimension <dim>

    if (strcmp(obj,"dimension")==0) {
      int d = atoi(args);
      hierarchy_.set_dim(d);
      Point::set_dim(d);
      Sphere::set_dim(d);

    // Grid ...

    } else if (strcmp(obj,"grid")==0) {

      hierarchy_.insert_grid(new Grid(args));

    // Sphere ...

    } else if (strcmp(obj,"sphere")==0) {

      spheres_.push_back(new Sphere(args));      

    // Point ...

    } else if (strcmp(obj,"point")==0) {

      points_.push_back(new Point(args));      

    }

    // Clear the buffer
    for (i=0; i<BUFFER_LENGTH; i++) buffer[i]=0;
  }

  // All grids are inserted into the hierarchy: call hierarchy.init_levels()

  hierarchy_.init_levels();

}

//----------------------------------------------------------------------

void Problem::print () throw ()
{
  int i;
  hierarchy().print();
  for (i=0; i<num_spheres(); i++) sphere(i).print();
  for (i=0; i<num_points(); i++)  point(i).print();
}

//----------------------------------------------------------------------

void Problem::write (FILE *fp) throw ()
{
  if (fp == 0) fp = stdout;
  int i;
  hierarchy().write(fp);
  for (i=0; i<num_spheres(); i++) sphere(i).write(fp);
  for (i=0; i<num_points(); i++)  point(i).write(fp);
}

//======================================================================

int Problem::readline_ (FILE* fp, char * buffer, int n) throw()
{
  int i=0;
  int c;
  buffer[i] = c = fgetc(fp);
  while (c != EOF && c != '\n' && i < n-1) {
    ++i;
    buffer[i] = c = fgetc(fp);
  }
  if (buffer[i] == '\n') buffer[i] = '\0';
  if (i == n-1) {
    fprintf (stderr,"Line too long: %s\n",buffer);
    exit(1);
  }
  return (c != EOF);
}

