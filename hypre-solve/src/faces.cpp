//345678901234567890123456789012345678901234567890123456789012345678901234567890

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
#include <mpi.h>

// #include "HYPRE_sstruct_ls.h"

const int debug = 0;

#include "hypre-solve.hpp"

#include "scalar.hpp"
#include "faces.hpp"
#include "mpi.hpp"
#include "domain.hpp"
#include "grid.hpp"

//----------------------------------------------------------------------

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// LabelName[] VALUES SHOULD MATCH Label enum ENTRIES IN faces.hpp
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

const char * Faces::LabelName[] = {
    "unknown",
    "boundary",
    "coarse",
    "fine",
    "neighbor",
    "covered",
    "adjacent_covered",
    "error"
  };

const int _num_types_(Faces::_last_ - Faces::_first_ + 1);

//----------------------------------------------------------------------

Faces::Faces (int *n) throw ()
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

int & Faces::entry_ (int axis, 
		     int face, 
		     int i, int j, 
		     int * entry[3][2]) throw ()
{
  // Remap axis to minimal one if face zone lies on edge or corner
  if (axis==1) {
    if (j==0 || j==n_[0]-1) {
      // Axis 1 to axis 0
      int face0=face;
      int i0=i;
      axis = 0;
      face = (j==0) ? 0 : 1;
      i = (face0==0) ? 0 : (n_[1]-1);
      j = i0;
    }
  } else if (axis==2) {
    if (i==0 || i==n_[0]-1) {
      // Axis 2 to axis 0
      int face0=face;
      axis = 0;
      face = (i==0) ? 0 : 1;
      i = j;
      j = (face0==0) ? 0 : (n_[2]-1);
    } else if (j==0 || j==n_[1]-1) {
      // Axis 2 to axis 1
      int face0=face;
      int i0=i;
      axis = 1;
      face = (j==0) ? 0 : 1;
      i = (face0==0) ? 0 : (n_[2]-1);
      j = i0;
    }
  }
  return entry[axis][face][i+n1_[axis]*j];
}

//----------------------------------------------------------------------

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
	   Label cell = label_[axis][face][index];
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

/// Allocate and initialize storage for label_[][] and adjacent_[][]
void Faces::alloc_ (int *n) throw ()
{
  int N = n[0]*n[1]*n[2];
  
  for (int axis=0; axis<3; axis++) {

    int i;

    // Determine face zone sizes

    n1_[axis] = n[(axis+1)%3];
    n2_[axis] = n[(axis+2)%3];
    n_ [axis] = N / n[axis];

    // Allocate and clear face zone labels and grid neighbors

    for (int face=0; face<2; face++) {

      // allocate and clear face zone labels

      label_[axis][face] = new Label [n_[axis]];
      for (i=0; i<n_[axis]; i++) label_[axis][face][i] = _unknown_;

      // allocate and clear face zone nonzero entries

      adjacent_[axis][face] = new pGrid [n_[axis]];
      for (i=0; i<n_[axis]; i++) adjacent_[axis][face][i] = NULL;
    }
  }
}

//----------------------------------------------------------------------

void Faces::dealloc_ () throw ()
{
  for (int axis=0; axis<3; axis++) {
    for (int face=0; face<2; face++) {

      // deallocate and clear face zone labels

      delete [] label_[axis][face];
      label_[axis][face] = NULL;

      // deallocate and clear grid neighbors

      delete [] adjacent_[axis][face];
      adjacent_[axis][face] = NULL;
      
    }
  }
}
