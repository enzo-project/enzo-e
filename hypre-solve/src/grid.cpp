
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
#include "domain.hpp"
#include "grid.hpp"

//======================================================================

const int debug       = 1;
const int debug_input = 0;

const int trace = 1;

Mpi    Grid::mpi_;
Domain Grid::domain_;

//======================================================================

Grid::Grid (std::string parms) throw ()
  : level_ (-1)

{
  // Initialize 0-sentinels in arrays

  neighbors0_.push_back (0);
  children0_.push_back (0);

  // Define a grid given text parameters, typically from a file

  read (parms);

  if (is_local()) {
    faces_ = new Faces(n_);
  } else {
    // Note: faces_ may be needed for grids adjacent to local grids
    //       these will be allocated as needed
    faces_ = 0;
  }
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

void Grid::geomview_grid (FILE *fpr, bool full) throw ()
{

  if (full) {
    fprintf (fpr,"VECT\n");
    fprintf (fpr,"6 18 2\n");
    fprintf (fpr,"1 1 8 3 3 2\n");
    fprintf (fpr,"1 0 1 0 0 0\n");

    // Print points at domain boundaries to provide geomview with bounding box

    Scalar dl[3],du[3];
    Grid::domain_.lower(dl[0],dl[1],dl[2]);
    Grid::domain_.upper(du[0],du[1],du[2]);
    fprintf (fpr,"%g %g %g\n",dl[0],dl[1],dl[2]);
    fprintf (fpr,"%g %g %g\n",du[0],du[1],du[2]);

  }
  fprintf (fpr,"%g %g %g\n",xl_[0],xl_[1],xl_[2]);
  fprintf (fpr,"%g %g %g\n",xu_[0],xl_[1],xl_[2]);
  fprintf (fpr,"%g %g %g\n",xu_[0],xu_[1],xl_[2]);
  fprintf (fpr,"%g %g %g\n",xl_[0],xu_[1],xl_[2]);
  fprintf (fpr,"%g %g %g\n",xl_[0],xl_[1],xl_[2]);
  fprintf (fpr,"%g %g %g\n",xl_[0],xl_[1],xu_[2]);
  fprintf (fpr,"%g %g %g\n",xu_[0],xl_[1],xu_[2]);
  fprintf (fpr,"%g %g %g\n",xu_[0],xl_[1],xl_[2]);
  fprintf (fpr,"%g %g %g\n",xl_[0],xu_[1],xl_[2]);
  fprintf (fpr,"%g %g %g\n",xl_[0],xu_[1],xu_[2]);
  fprintf (fpr,"%g %g %g\n",xl_[0],xl_[1],xu_[2]);
  fprintf (fpr,"%g %g %g\n",xu_[0],xu_[1],xl_[2]);
  fprintf (fpr,"%g %g %g\n",xu_[0],xu_[1],xu_[2]);
  fprintf (fpr,"%g %g %g\n",xu_[0],xl_[1],xu_[2]);
  fprintf (fpr,"%g %g %g\n",xl_[0],xu_[1],xu_[2]);
  fprintf (fpr,"%g %g %g\n",xu_[0],xu_[1],xu_[2]);

  if (full) {
    fprintf (fpr,"1 1 1 1\n");
    fprintf (fpr,"1 1 1 0\n");
  }
}

//======================================================================

void Grid::geomview_face (FILE *fpr, bool full) throw ()
{

  if (debug) printf ("DEBUG %s:%d grid ip = %d mpi ip = %d grid id = %d\n",
		     __FILE__,__LINE__,ip(),mpi_.ip(),id());
  if (debug) {
    int ip;
    MPI_Comm_rank(MPI_COMM_WORLD,&ip);
    printf ("DEBUG %s:%d mpi rank = %d\n",__FILE__,__LINE__,ip);
  }
  if (! is_local()) return;

  float bcolor[] = {1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0};
  float rcolor[] = {1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0};
  float gcolor[] = {1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
  float acolor[] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};

  if (full) {
    fprintf (fpr,"CQUAD\n");
    // Print points at domain boundaries to provide geomview with bounding box

    Scalar dl[3],du[3];
    Grid::domain_.lower(dl[0],dl[1],dl[2]);
    Grid::domain_.upper(du[0],du[1],du[2]);
    fprintf (fpr,
	     "%g %g %g 0 0 0 1 "
	     "%g %g %g 0 0 0 1 "
	     "%g %g %g 0 0 0 1 "
	     "%g %g %g 0 0 0 1\n",
	     dl[0],dl[1],dl[2],
	     dl[0],dl[1],dl[2],
	     dl[0],dl[1],dl[2],
	     dl[0],dl[1],dl[2]);
    fprintf (fpr,
	     "%g %g %g 0 0 0 1 "
	     "%g %g %g 0 0 0 1 "
	     "%g %g %g 0 0 0 1 "
	     "%g %g %g 0 0 0 1\n",
	     du[0],du[1],du[2],
	     du[0],du[1],du[2],
	     du[0],du[1],du[2],
	     du[0],du[1],du[2]);
  }

  if (debug) printf ("DEBUG %s:%d grid id = %d\n",__FILE__,__LINE__,id());

  int axis,face,i1,i2;
  int j0,j1,j2;
  Scalar d0,d1,d2;
  int i;

  Scalar zc[3];   // zc = zone center
  Scalar fz[3];  // fz = offset of face center from zone center
  Scalar fc[3];  // fc = face centers
  Scalar vf[3];  // vf = offset of face vertices from face center
  Scalar vc[4][3]; // vf = offset of face vertices from zone center

  const Scalar svf = 0.4;  // svf = scaling of offset of vertices from face center
  //                            (0.5 = full side)
  Scalar sfz       = 0.0;  // sfz = scaling of offset of face center from zone center
  //                            (0.5 = actual face location; set dynamically belew)

  for (axis=0; axis<3; axis++) {
    for (face = 0; face<2; face++) {
      for (i1=0; i1<faces().n1(axis); i1++) {
	for (i2=0; i2<faces().n2(axis); i2++) {
	  // If face zone has the given label, print it
	  // determine index position (j0,j1,j2) of the face zone

	  if (face==0) sfz = -svf;
	  if (face==1) sfz = +svf;

	  if (axis==0) {
	    if (face==0) j0 = 0;
	    if (face==1) j0 = n_[0] - 1;
	    j1 = i1;
	    j2 = i2;
	    vc[0][0] = + sfz*h(0);  vc[0][1] = + svf*h(1);  vc[0][2] = + svf*h(2);
	    vc[1][0] = + sfz*h(0);  vc[1][1] = + svf*h(1);  vc[1][2] = - svf*h(2);
	    vc[2][0] = + sfz*h(0);  vc[2][1] = - svf*h(1);  vc[2][2] = - svf*h(2);
	    vc[3][0] = + sfz*h(0);  vc[3][1] = - svf*h(1);  vc[3][2] = + svf*h(2);
	  } else if (axis==1) {
	    if (face==0) j1 = 0;
	    if (face==1) j1 = n_[1] - 1;
	    j2 = i1;
	    j0 = i2;
	    vf[0] = svf*h(0);
	    vf[1] = 0;
	    vf[2] = svf*h(2);
	    fz[0] = 0;
	    fz[1] = h(1);
	    fz[2] = 0;
	    vc[0][0] = + svf*h(0);  vc[0][1] = + sfz*h(1);  vc[0][2] = + svf*h(2);
	    vc[1][0] = + svf*h(0);  vc[1][1] = + sfz*h(1);  vc[1][2] = - svf*h(2);
	    vc[2][0] = - svf*h(0);  vc[2][1] = + sfz*h(1);  vc[2][2] = - svf*h(2);
	    vc[3][0] = - svf*h(0);  vc[3][1] = + sfz*h(1);  vc[3][2] = + svf*h(2);
	  } else if (axis==2) {
	    if (face==0) j2 = 0;
	    if (face==1) j2 = n_[2] - 1;
	    j0 = i1;
	    j1 = i2;
	    vf[0] = svf*h(0);
	    vf[1] = svf*h(1);
	    vf[2] = 0;
	    fz[0] = 0;
	    fz[1] = 0;
	    fz[2] = h(2);
	    vc[0][0] = + svf*h(0);  vc[0][1] = + svf*h(1);  vc[0][2] = + sfz*h(2);
	    vc[1][0] = - svf*h(0);  vc[1][1] = + svf*h(1);  vc[1][2] = + sfz*h(2);
	    vc[2][0] = - svf*h(0);  vc[2][1] = - svf*h(1);  vc[2][2] = + sfz*h(2);
	    vc[3][0] = + svf*h(0);  vc[3][1] = - svf*h(1);  vc[3][2] = + sfz*h(2);
	  }

	  // Determine cell center
	  zone(j0,j1,j2,zc[0],zc[1],zc[2]);

	  int label = faces().label(axis,face,i1,i2);

	  for (i=0; i<4; i++) {
	    fprintf (fpr,"%g %g %g  %f %f %f %f  ",
		     zc[0]+vc[i][0], zc[1]+vc[i][1], zc[2]+vc[i][2],
		     bcolor[label], rcolor[label], gcolor[label], acolor[label]);
	  }
	  fprintf (fpr,"\n");
	}
      }
    }
  }
}

//--------------------------------------------------------------------

/// Write the Face data to the given open file in geomview format

// void Faces::geomview_face (char fileprefix[]) throw ()
// {

//   // For each Label type, write all face-zones with that label to a geomview vect file

//   for (Label label = _first_;
//        label <= _last_;
//        label = Label(label + 1))
//     {
//       // Open
//       std::string filename = (std::string)(fileprefix) + "-" + LabelName[label] + ".vect";
//       if (debug) printf ("DEBUG %s:%d %s\n",__FILE__,__LINE__,filename.c_str());


//     }
//   NOT_IMPLEMENTED("Faces::geomview()");
// }

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
    h1[i] = ( g1.xu_[i] - g1.xl_[i] ) / g1.num_unknowns(i);
    h2[i] = ( g2.xu_[i] - g2.xl_[i] ) / g2.num_unknowns(i);
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
    if (fabs(g1.xl_[i] - g2.xu_[i]) < 0.5*h1[i]) face[i] = -1;
    if (fabs(g1.xu_[i] - g2.xl_[i]) < 0.5*h1[i]) face[i] =  1;
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
    lower[i] = MAX(g1.xl_[i],g2.xl_[i]);
    upper[i] = MIN(g1.xu_[i],g2.xu_[i]);
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
      gl[i] = int ( (lower[i] - g1.xl_[i]) / h1[i] );
      gu[i] = int ( (upper[i] - g1.xl_[i]) / h1[i] - 1);
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

