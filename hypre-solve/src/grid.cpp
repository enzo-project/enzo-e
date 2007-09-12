
/// Grid class source file

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

#include "hypre-solve.hpp"

#include "scalar.hpp"
#include "faces.hpp"
#include "mpi.hpp"
#include "grid.hpp"

//======================================================================

const int debug       = 0;
const int debug_input = 0;

Mpi Grid::mpi_;

//======================================================================

Grid::Grid (std::string parms) throw ()
  : level_ (-1)

{
  // Initialize 0-sentinels in arrays

  neighbors0_.push_back (0);
  children0_.push_back (0);

  // Define a grid given text parameters, typically from a file

  read (parms);

  faces_ = new Faces(n_);
}
	  
//======================================================================

Grid::~Grid () throw ()
{
}

//======================================================================

void Grid::print () throw ()
{
  printf ("Grid\n"
	  "   id             %d\n"
	  "   parent id      %d\n"
	  "   processor      %d\n"
	  "   lower position "SCALAR_PRINTF SCALAR_PRINTF SCALAR_PRINTF"\n"
	  "   upper position "SCALAR_PRINTF SCALAR_PRINTF SCALAR_PRINTF"\n"
	  "   lower index    %d %d %d\n"
	  "   zones          %d %d %d\n"
	  "   level          %d\n",
	  id_,id_parent_,ip_,
	  xl_[0],xl_[1],xl_[2],
	  xu_[0],xu_[1],xu_[2],
	  il_[0],il_[1],il_[2],
	  n_ [0],n_ [1],n_ [2],
	  level_);
}

//======================================================================

void Grid::write (FILE *fp) throw ()
{
  if (fp == 0) fp = stdout;

  fprintf (fp,"grid "
	   "%d %d %d "
	   SCALAR_PRINTF SCALAR_PRINTF SCALAR_PRINTF
	   SCALAR_PRINTF SCALAR_PRINTF SCALAR_PRINTF
	   "%d %d %d\n"
	   "%d %d %d\n",
	   id_,id_parent_,ip_,
	   xl_[0],xl_[1],xl_[2],
	   xu_[0],xu_[1],xu_[2],
	   il_[0],il_[1],il_[2],
	   n_ [0],n_ [1],n_ [2]);
}

//======================================================================

void Grid::read (std::string parms) throw ()
{

  sscanf (parms.c_str(),
	  "%d%d%d" 
	  SCALAR_SCANF SCALAR_SCANF SCALAR_SCANF
	  SCALAR_SCANF SCALAR_SCANF SCALAR_SCANF
	  "%d%d%d"
	  "%d%d%d",
	  &id_, &id_parent_, &ip_,
	  &xl_[0],&xl_[1],&xl_[2],
	  &xu_[0],&xu_[1],&xu_[2],
	  &il_[0],&il_[1],&il_[2],
	  &n_[0],&n_[1],&n_[2]);
  if (debug_input) print();

}

//======================================================================

bool Grid::is_adjacent (Grid & g2) throw ()
{
  Grid & g1 = *this;

  // Must be in same level to be adjacent

  if (g1.level_ != g2.level_) return false;

  // hh is a tolerance to avoid problems with
  // comparing floating point numbers.  It is
  // taken to be 1/2 the mesh width.
  double hh; 
  // Assume they are adjacent
  bool far = false;  
  for (int i=0; i<3; i++) {
    hh = 0.5 * (g1.xu_[i]-g1.xl_[i])/g1.n_[i]; // hh should be same for g2
    far = far || (g1.xu_[i] < (g2.xl_[i] - hh));
    far = far || (g2.xu_[i] < (g1.xl_[i] - hh));
  }
  return ! far;
}

//======================================================================

/// Return indices of zones adjacent to neighboring grid.

/** Return lower and upper coordinates of cells adjacent to the input
    neighboring grid.  If the input grid is not a neighbor, return
    false. */

bool Grid::find_neighbor_indices (Grid & neighbor, 
				  int *gl, int *gu)
{
  //  +------+
  //  |......|
  //  |......+---+    Return coordinates of unknowns X in grid G's local
  //  |..N..X|...|    coordinates (gl[] and gu[]).  Return false if
  //  |.....X|.G.|    there are none (e.g. if grids are not neighbors, 
  //  +------+...|    or only adjacent at corners or edges not faces.
  //         |...| 
  //         +---+

  Grid & g1 = *this;
  Grid & g2 = neighbor;

  // Check grid spacing

  int i;
  Scalar h1[3],h2[3];
  for (i=0; i<3; i++) {
    h1[i] = ( g1.x_upper(i) - g1.x_lower(i) ) / g1.num_unknowns(i);
    h2[i] = ( g2.x_upper(i) - g2.x_lower(i) ) / g2.num_unknowns(i);
    // Check that mesh sizes are close.  Can fail on bad input
    if (debug) printf ("DEBUG %s:%d h1[%d] = %g  h2 = %g\n",
		       __FILE__,__LINE__,i,h1[i],h2[i]);
    assert (fabs(h1[i] - h2[i])/h1[i] < 0.01); 
  }

  // Find which grid planes are close

  int face[3];   // -1 g1.lower adjacent to g2.upper
                 //  0 not adjacent along face
                 // +1 g1.upper adjacent to g2.lower

  for (i=0; i<3; i++) {
    face[i] = 0;  // Assume not close
    if (fabs(g1.x_lower(i) - g2.x_upper(i)) < 0.5*h1[i]) face[i] = -1;
    if (fabs(g1.x_upper(i) - g2.x_lower(i)) < 0.5*h1[i]) face[i] =  1;
    if (debug) printf ("DEBUG %s:%d face[%d] = %d\n",
		       __FILE__,__LINE__,i,face[i]);
  }

  // Return if grids are not adjacent along a single face 
  // (then must be adjacent along an edge, a corner, or not neighbors)
  
  int face_count = abs(face[0]) + abs(face[1]) + abs(face[2]);
  if (face_count != 1 ) return false;

  // Find extent of 'intersecting grid'

  Scalar lower[3],upper[3];

  for (i=0; i<3; i++) {
    lower[i] = MAX(g1.x_lower(i),g2.x_lower(i));
    upper[i] = MIN(g1.x_upper(i),g2.x_upper(i));
    if (debug) printf ("DEBUG %s:%d lower[%d] = %g  upper = %g\n",
		       __FILE__,__LINE__,i,lower[i],upper[i]);
  }

  if (debug) printf ("g1.num_unknowns = (%d %d %d)\n",
		     g1.num_unknowns(0),
		     g1.num_unknowns(1),
		     g1.num_unknowns(2));
  if (debug) printf ("g2.num_unknowns = (%d %d %d)\n",
		     g2.num_unknowns(0),
		     g2.num_unknowns(1),
		     g2.num_unknowns(2));
  for (i=0; i<3; i++) {
    if (face[i] == -1) {
      gl[i] = gu[i] = -1;
    } else if (face[i] == 0) {
      gl[i] = int ( (lower[i] - g1.x_lower(i)) / h1[i] );
      gu[i] = int ( (upper[i] - g1.x_lower(i)) / h1[i] - 1);
    } else if (face[i] == 1) {
      gl[i] = gu[i] = g1.num_unknowns(i);
    }
  }
  if (debug)  printf ("DEBUG %s:%d gl[] = (%d %d %d)\n",
		      __FILE__,__LINE__,gl[0],gl[1],gl[2]);
  if (debug)  printf ("DEBUG %s:%d gu[] = (%d %d %d)\n",
		      __FILE__,__LINE__,gu[0],gu[1],gu[2]);

  return true;
}
//======================================================================
// PRIVATE MEMBER FUNCTIONS
//======================================================================

