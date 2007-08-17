
/// Faces class source file

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
#include <math.h>
#include <string>
#include <vector>

#include "HYPRE_sstruct_ls.h"

#include "scalar.hpp"
#include "faces.hpp"
#include "mpi.hpp"
#include "grid.hpp"

//----------------------------------------------------------------------

const int debug = 0;

const int Faces::_unknown_  = -1000;
const int Faces::_boundary_ = -2000;

//----------------------------------------------------------------------

Faces::Faces (int n[3]) throw ()
{
  alloc_(n);   // NOTE: inefficient for dimension < 3
}
	  
//----------------------------------------------------------------------

Faces::~Faces () throw ()
{
  dealloc_();
}

//--------------------------------------------------------------------
// PUBLIC MEMBER FUNCTIONS
//--------------------------------------------------------------------

void Faces::print() throw()
{
   printf ("Faces::debug()\n");
   printf ("   n1_ = (%d,%d,%d)\n",n1_[0],n1_[1],n1_[2]);
   printf ("   n2_ = (%d,%d,%d)\n",n2_[0],n2_[1],n2_[2]);
   printf ("    n_  = (%d,%d,%d)\n",n_[0],n_[1],n_[2]);
   int axis,face,i,j;
   for (axis=0; axis<3; axis++) {
     for (face=0; face<2; face++) {
       printf ("   Axis %d   Face %d\n",axis,face);
       for (i=0; i<n1_[axis]; i++) {
	 for (j=0; j<n2_[axis]; j++) {
	   int index = i+n1_[axis]*j;
	   int cell = neighbor_cell_[axis][face][index];
	   if (cell == _unknown_)       printf ("??");
	   else if (cell == _boundary_) printf ("BB");
	   else printf ("%2d",cell);
	 }
	 printf ("\n");
       }
     }
   }


}

//======================================================================
// PRIVATE MEMBER FUNCTIONS
//======================================================================

void Faces::alloc_ (int n[3]) throw ()
//
// Allocate and initialize storage for neighbor_cell_[][]
//
{
  int i;

  int N = n[0]*n[1]*n[2];
  
  for (int axis=0; axis<3; axis++) {
    // Determine boundary face sizes
    n1_[axis] = n[(axis+1)%3];
    n2_[axis] = n[(axis+2)%3];
    n_ [axis] = N / n[axis];
    if (debug) printf ("DEBUG %s:%d axis %d n1=%d n2=%d n=%d\n",
		       __FILE__,__LINE__,axis,n1_[axis],n2_[axis],n_[axis]);

    // Allocate and initialize boundary faces

    for (int face=0; face<2; face++) {
      neighbor_cell_[axis][face] = new int [n_[axis]];
      for (i=0; i<n_[axis]; i++) neighbor_cell_[axis][face][i] = _unknown_;
    }
  }
}

//----------------------------------------------------------------------

void Faces::dealloc_ () throw ()
{
  for (int axis=0; axis<3; axis++) {
    for (int face=0; face<2; face++) {
      delete [] neighbor_cell_[axis][face];
      neighbor_cell_[axis][face] = 0;
      
    }
  }
}
