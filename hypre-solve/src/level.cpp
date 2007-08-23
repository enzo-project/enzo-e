
/// Level class source file

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
#include <string>
#include <vector>

#include "HYPRE_sstruct_ls.h"

#include "hypre-solve.hpp"

#include "scalar.hpp"
#include "faces.hpp"
#include "mpi.hpp"
#include "grid.hpp"
#include "level.hpp"

//----------------------------------------------------------------------

const int debug = 1;

//----------------------------------------------------------------------

Level::Level (int n) throw ()
  : n_ (n)
{
  // Initialize extents
  for (int i=0; i<3; i++) {
    il_[i] = INT_MAX;
    iu_[i] = INT_MIN;
  }
  grids0_.push_back (0);
}
	  
//----------------------------------------------------------------------

Level::~Level () throw ()
{
}

//======================================================================

void Level::insert_grid (Grid & grid) throw ()
{
  // Append Grid to list of all Grids in the Level

  grids0_[grids0_.size()-1] = & grid;
  grids0_.push_back (0);

  // Update Level extents given Grid indices

  for (int i=0; i<3; i++) {
    il_[i] = MIN(il_[i],grid.i_lower(i));
    iu_[i] = MAX(iu_[i],grid.i_upper(i));
  }
}

//======================================================================

void Level::print () throw ()
{
  printf ("Level %d\n",n_);
  printf ("   i_lower = (%d,%d,%d)\n",il_[0],il_[1],il_[2]);
  printf ("   i_upper = (%d,%d,%d)\n",iu_[0],iu_[1],iu_[2]);
  for (int i=0; i<num_grids(); i++) {
    grid(i).print();
  }
}

//----------------------------------------------------------------------

void Level::write (FILE *fp) throw ()
{
  if (fp == 0) fp = stdout;

  for (int i=0; i<num_grids(); i++) {
    fprintf (fp,"Level %d\n",n_);
    grid(i).write(fp);
  }
}

