
/// Discret class source file

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
#include "discret.hpp"
#include "mpi.hpp"
#include "grid.hpp"

//----------------------------------------------------------------------

const int debug = 1;

//----------------------------------------------------------------------

Discret::Discret (int n[3]) throw ()
{
  alloc_(n);   // NOTE: inefficient for dimension < 3
}
	  
//----------------------------------------------------------------------

Discret::~Discret () throw ()
{
  dealloc_();
}

//--------------------------------------------------------------------
// PUBLIC MEMBER FUNCTIONS
//--------------------------------------------------------------------

void Discret::print() throw()
{
  return;
   printf ("Discret::debug()\n");
   printf ("   n1_ = (%d,%d,%d)\n",n1_[0],n1_[1],n1_[2]);
   printf ("   n2_ = (%d,%d,%d)\n",n2_[0],n2_[1],n2_[2]);
   printf ("    n_  = (%d,%d,%d)\n",n_[0],n_[1],n_[2]);
   int axis,face,i,j;
   for (axis=0; axis<3; axis++) {
     for (face=0; face<2; face++) {
       printf ("   Axis %d   Face %d\n",axis,face);
       for (i=0; i<n1_[axis]; i++) {
	 for (j=0; j<n2_[axis]; j++) {
	   char c;
	   int index = i+n1_[axis]*j;
	   if (neighbor_cell_[axis][face][index]==_same_)    c='=';
	   if (neighbor_cell_[axis][face][index]==_coarse_)  c='>';
	   if (neighbor_cell_[axis][face][index]==_fine_)    c='<';
	   if (neighbor_cell_[axis][face][index]==_unknown_) c='?';
	   printf ("%c ",c);
	 }
	 printf ("\n");
       }
     }
   }


}

//======================================================================
// PRIVATE MEMBER FUNCTIONS
//======================================================================

void Discret::alloc_ (int n[3]) throw ()
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
    for (int face=0; face<2; face++) {
      // Allocate boundary face
      neighbor_cell_[axis][face] = new enum_neighbor_cell [n_[axis]];
      // Clear boundary face
      for (i=0; i<n_[axis]; i++) neighbor_cell_[axis][face][i] = _unknown_;
    }
  }
}

//----------------------------------------------------------------------

void Discret::dealloc_ () throw ()
{
  for (int axis=0; axis<3; axis++) {
    for (int face=0; face<2; face++) {
      delete [] neighbor_cell_[axis][face];
      neighbor_cell_[axis][face] = 0;
      
    }
  }
}
