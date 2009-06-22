//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Grid class source file

/**
 * @file      grid.cpp
 * @brief     Implemtation of the Grid class
 * @author    James Bordner
 * @bug       none
 *
 * $Id$
 */

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>
#include <mpi.h>

// #include "HYPRE_sstruct_ls.h"

#include "newgrav-hypre-solve.h"

//======================================================================

const int debug       = 0;
const int debug_input = 0;
const int trace       = 0;

//======================================================================

#include "newgrav-error.h"
#include "newgrav-scalar.h"
#include "newgrav-constants.h"
#include "newgrav-faces.h"
#include "newgrav-mpi.h"
#include "newgrav-domain.h"
#include "newgrav-grid.h"

//======================================================================

Mpi    Grid::mpi_;
Domain Grid::domain_;

//======================================================================

Grid::Grid (std::string parms) throw ()
  : level_ (-1),
    u_(0),
    counters_ (0)

{
  // Initialize 0-sentinels in arrays

  neighbors0_.push_back (0);
  children0_.push_back (0);

  // Define a grid given text parameters, typically from a file

  input (parms);

  // Allocate Faces was here.

  printf ("%s:%d n_ = %d %d %d\n",__FILE__,__LINE__,n_[0],n_[1],n_[2]);

  faces_ = new Faces(n_);

  // Allocate counters_ here.

  counters_ = new int [n_[0]*n_[1]*n_[2]];

}

//======================================================================

Grid::Grid (int     id, 
	    int     id_parent, 
	    int     ip, 
	    Scalar *xl,
	    Scalar *xu,
	    int    *il,
	    int    *n) throw ()
{
 // Initialize 0-sentinels in arrays

  neighbors0_.push_back (0);
  children0_.push_back (0);

  // Define a grid given text parameters, typically from a file

  input (id,id_parent,ip,xl,xu,il,n);

  // Allocate Faces was here.

  printf ("%s:%d n_ = %d %d %d\n",__FILE__,__LINE__,n_[0],n_[1],n_[2]);

  faces_ = new Faces(n_);

  // Allocate counters_ here.

  counters_ = new int [n_[0]*n_[1]*n_[2]];

}

//======================================================================

Grid::Grid (FILE *fp) throw ()
  : id_(-1),
    id_parent_(-1),
    ip_(-1),
    faces_(NULL),
    level_(-1),
    u_(NULL),
    counters_(NULL)
{
  this->read(fp);
}
//======================================================================

Grid::~Grid () throw ()
{
  delete [] u_;
  delete faces_;
  delete [] counters_;
}

//======================================================================

void Grid::print () throw ()
{
  this->write(stdout,true);
}

//======================================================================

void Grid::write (FILE *fp, bool brief) throw ()
{
  if (fp == 0) fp = stdout;
  fprintf (fp,"Grid\n"
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
  if (u_ && ! brief) {
    for (int i0=0; i0<n_[0]; i0++) {
      for (int i1=0; i1<n_[1]; i1++) {
	for (int i2=0; i2<n_[2]; i2++) {
	  int i = index(i0,i1,i2,n_[0],n_[1],n_[2]);
	  fprintf (fp,"%d %d %d %22.15e\n",i0,i1,i2,u_[i]);
	}
      }
    }
  }
}

//======================================================================

void Grid::read (FILE *fp, bool brief) throw ()
{
  if (fp == 0) fp = stdin;
  fscanf (fp,"Grid");
  fscanf (fp,"   id             %d",&id_);
  if (debug) {printf ("DEBUG %s:%d %d\n",__FILE__,__LINE__,id_); fflush(stdout); }
  fscanf (fp,"   parent id      %d",&id_parent_);
  if (debug) {printf ("DEBUG %s:%d %d\n",__FILE__,__LINE__,id_parent_); fflush(stdout); }
  fscanf (fp,"   processor      %d",&ip_);
  if (debug) {printf ("DEBUG %s:%d %d\n",__FILE__,__LINE__,ip_); fflush(stdout); }
  fscanf (fp,"   lower position "SCALAR_SCANF SCALAR_SCANF SCALAR_SCANF,
	  &xl_[0],&xl_[1],&xl_[2]);
  if (debug) {printf ("DEBUG %s:%d %g %g %g\n",__FILE__,__LINE__,xl_[0],xl_[1],xl_[2]); fflush(stdout); }
  fscanf (fp,"   upper position "SCALAR_SCANF SCALAR_SCANF SCALAR_SCANF,
	  &xu_[0],&xu_[1],&xu_[2]);
  fscanf (fp,"   lower index    %d %d %d",
	  &il_[0],&il_[1],&il_[2]);
  fscanf (fp,"   zones          %d %d %d",
	  &n_ [0],&n_ [1],&n_ [2]);
  fscanf (fp,"   level          %d",	 &level_);

  this->allocate();
  if (u_ && ! brief) {
    int i0,i1,i2;
    Scalar u;
    int status;
    while ((status = fscanf(fp,"%d%d%d"SCALAR_SCANF, &i0,&i1,&i2,&u)) != EOF) {
      int i = index(i0,i1,i2,n_[0],n_[1],n_[2]);
      u_[i] = u;
    }
  }
}

//======================================================================

Scalar * Grid::values () throw ()
{
  assert (u_);
  return u_;
}

//======================================================================

void Grid::allocate () throw ()
{
  if (u_) deallocate();
  u_ = new Scalar [n_[0]*n_[1]*n_[2]];
}

//======================================================================

void Grid::deallocate () throw ()
{
  delete [] u_;
  u_ = 0;
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

void Grid::geomview_face (FILE         *fpr, 
			  bool          full) throw ()
{
  int num_types = Faces::_last_ - Faces::_first_ + 1;
  Faces::Label * types = new Faces::Label [num_types];
  for (int i=0; i<num_types; i++) types[i] = Faces::Label(i);
  geomview_face_type (fpr,types,num_types,full);
  delete [] types;
}

//======================================================================

void Grid::geomview_face_type (FILE         *fpr, 
			       Faces::Label *types, 
			       int           num_types, 
			       bool          full) throw ()
{

  if (debug) printf ("DEBUG %s:%d grid ip = %d mpi ip = %d grid id = %d\n",
		     __FILE__,__LINE__,ip(),mpi_.ip(),id());
  if (debug) {
    int ip;
    MPI_Comm_rank(MPI_COMM_WORLD,&ip);
    printf ("DEBUG %s:%d mpi rank = %d\n",__FILE__,__LINE__,ip);
  }
  if (! is_local()) return;

  // Default acolor = 0.0 to make them invisible

  // unknown            000   black
  // boundary           100   red
  // coarse             010   green
  // fine               001   blue
  // neighbor           011   cyan
  // covered            101   magenta
  // adjacent_covered   110   yellow
  // error              111   white

  float bcolor[] = {0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
  float rcolor[] = {0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0};
  float gcolor[] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0};
  float acolor[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  // Make requested types visible
  int i;
  for (i=0; i<num_types; i++) acolor[types[i]] = 1.0;

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

  Scalar zc[3];   // zc = zone center
  Scalar fz[3];  // fz = offset of face center from zone center
  Scalar vf[3];  // vf = offset of face vertices from face center
  Scalar vc[4][3]; // vf = offset of face vertices from zone center

  const Scalar svf = 0.2;  // svf = scaling of offset of vertices from face center
  //                            (0.5 = full side)
  Scalar sfz;              // sfz = scaling of offset of face center from zone center
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
	    vc[0][0] = +sfz*h(0); vc[0][1] = +svf*h(1); vc[0][2] = +svf*h(2);
	    vc[1][0] = +sfz*h(0); vc[1][1] = +svf*h(1); vc[1][2] = -svf*h(2);
	    vc[2][0] = +sfz*h(0); vc[2][1] = -svf*h(1); vc[2][2] = -svf*h(2);
	    vc[3][0] = +sfz*h(0); vc[3][1] = -svf*h(1); vc[3][2] = +svf*h(2);
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
	    vc[0][0] = +svf*h(0); vc[0][1] = +sfz*h(1); vc[0][2] = +svf*h(2);
	    vc[1][0] = +svf*h(0); vc[1][1] = +sfz*h(1); vc[1][2] = -svf*h(2);
	    vc[2][0] = -svf*h(0); vc[2][1] = +sfz*h(1); vc[2][2] = -svf*h(2);
	    vc[3][0] = -svf*h(0); vc[3][1] = +sfz*h(1); vc[3][2] = +svf*h(2);
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
	    vc[0][0] = +svf*h(0); vc[0][1] = +svf*h(1); vc[0][2] = +sfz*h(2);
	    vc[1][0] = -svf*h(0); vc[1][1] = +svf*h(1); vc[1][2] = +sfz*h(2);
	    vc[2][0] = -svf*h(0); vc[2][1] = -svf*h(1); vc[2][2] = +sfz*h(2);
	    vc[3][0] = +svf*h(0); vc[3][1] = -svf*h(1); vc[3][2] = +sfz*h(2);
	  }

	  // Determine cell center
	  zone(j0,j1,j2,zc[0],zc[1],zc[2]);

	  int label = faces().label(axis,face,i1,i2);
	  if (label < Faces::_first_ || label > Faces::_last_) {
	    label = Faces::_error_;
	  }

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

void Grid::input (std::string parms) throw ()
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

void Grid::input (int     id, 
		  int     id_parent, 
		  int     ip, 
		  Scalar *xl,
		  Scalar *xu,
		  int    *il,
		  int    *n) throw ()
{

  id_        = id;
  id_parent_ = id_parent;
  ip_        = ip;
  xl_[0]     = xl[0];
  xl_[1]     = xl[1];
  xl_[2]     = xl[2];
  xu_[0]     = xu[0];
  xu_[1]     = xu[1];
  xu_[2]     = xu[2];
  il_[0]     = il[0];
  il_[1]     = il[1];
  il_[2]     = il[2];
  n_[0]      = n[0];
  n_[1]      = n[1];
  n_[2]      = n[2];

}

//----------------------------------------------------------------------

bool Grid::is_adjacent (Grid & g2, Scalar period[3]) throw ()
{
  Grid & g1 = *this;

  // Must be in same level to be adjacent

  if (g1.level_ != g2.level_) return false;

  // Precompute tolerance as 0.5 * cell size 
  
  double hh[3];
  for (int i=0; i<3; i++) {
    hh[i] = 0.5 * (g1.xu_[i]-g1.xl_[i])/g1.n_[i];
  }

  // Assume g1 and g2 are adjacent

  bool far = false;

  // Handle the case where we're testing the grid against itself
  // (can be adjacent if domain is periodic)

  // If grids g1 and g2 are the same grid...

  if (g1.id() == g2.id()) {
    
    for (int i=0; i<3; i++) {
      if (fabs(g1.xl_[i] + period[i] - g1.xu_[i]) < hh[i]) {
	return true;
      }
    }
    return false;
  } else {

    for (int i=0; i<3; i++) {

      // Define tolerance hh to be 0.5 * cell width

      if (period[i] == 0) {

	// If not periodic ...

	// If g1's upper edge is below g2's lower edge, then they are far

	far = far || (g1.xu_[i] < g2.xl_[i] - hh[i]);

	// If g2's upper edge is below g1's lower edge, then they are far

	far = far || (g2.xu_[i] < g1.xl_[i] - hh[i]);

      } else {

	// If periodic ...

	// If g1's upper edge is below g2's lower edge, then they are far

	far = far || ( (g1.xu_[i] < g2.xl_[i] - hh[i]) &&
		       (g2.xu_[i] < g1.xl_[i] - hh[i] + period[i]));

	// If g2's upper edge is below g1's lower edge, then they are far

	far = far || ( (g2.xu_[i] < g1.xl_[i] - hh[i]) &&
		       (g1.xu_[i] < g2.xl_[i] - hh[i] + period[i]));

      }
    }

    // If g1 and g2 are not far, then they are adjacent

    //  if (! far) printf ("%s:%d grids %d and %d are adjacent\n",
    //		     __FILE__,__LINE__,g1.id(),g2.id());

    return ! far;
  }

}

//----------------------------------------------------------------------

bool Grid::contains (Grid * grid) throw ()
{
  // Assume initially that the other grid is contained in this one
  bool is_inside = true;

  // The other grid is not contained in this grid if it is in a coarser level
  if (grid->level() < this->level()) is_inside = false;

  double hh[3];
  for (int i=0; i<3; i++) {
    hh[i] = 0.5 * (grid->xu_[i]-grid->xl_[i])/grid->n_[i];
  }

  // grid is contained in (*this) if its projection along each axis is contained in this's
  for (int i=0; i<3; i++) {
    // Lower point is greater
    is_inside = is_inside && ((this->xl_[i]) < (grid->xl_[i] + hh[i]));
    // Upper point is less
    is_inside = is_inside && ((grid->xu_[i]) < (this->xu_[i] + hh[i]));
  }
  return is_inside;
}

//----------------------------------------------------------------------

/// Determine the axis, face, and range of indices of zones adjacent 
/// to the neighboring grid.  Returns false if the neighbor is not
/// actually a neighbor.  Assumes grids are in the same level.

bool Grid::neighbor_shared_face (Grid & neighbor, 
				 int & axis, int & face, 
				 int & il0, int & il1, 
				 int & iu0, int & iu1,
				 int iperiod[3],
				 int & count) throw ()
{

  Grid & grid = *this;

  // Get grid index bounds

  int index_grid[3][2];
  grid.indices(index_grid);

  // Get neighbor index bounds

  int index_neighbor[3][2];
  neighbor.indices(index_neighbor);

  // Find matching face, if any

  bool found_face = false;
  int  iaxis      = -1;
  int  iface      = -1;
  int  num        = 0;

  for (axis=0; axis<3; axis++) {
    for (face=0; face<2; face++) {

      bool b1 = index_grid[axis][face] == index_neighbor[axis][1-face];
      bool b2 = 
      	(index_grid[axis][face]   + (1-face)*iperiod[axis]) == 
      	(index_neighbor[axis][1-face] +     face*iperiod[axis]);

      if ((b1 || b2) && ! found_face) {
	if (num == count) {
	  found_face = true;
	  iaxis=axis;
	  iface=face;
	  count++;
	} else {
	  num++;
	}
      }
    }
  }

  axis=iaxis;
  face=iface;

  // Exit if face isn't found

  if (! found_face) return false;

  // face axes

  int j0=(axis+1)%3;
  int j1=(axis+2)%3;

  // Compute local indices intersection from global indices of each grid

  il0 = MAX(index_grid[j0][0],index_neighbor[j0][0]) - index_grid[j0][0];
  iu0 = MIN(index_grid[j0][1],index_neighbor[j0][1]) - index_grid[j0][0] - 1;

  il1 = MAX(index_grid[j1][0],index_neighbor[j1][0]) - index_grid[j1][0];
  iu1 = MIN(index_grid[j1][1],index_neighbor[j1][1]) - index_grid[j1][0] - 1;

  return true;

}

//----------------------------------------------------------------------

/// Determine the axis, face, and range of indices of zones adjacent 
/// to the adjacent grid in the next-coarser level.  Returns false if
/// the grid is not actually adjacent.  Assumes the adjacent grid is
/// in the next-coarser level, and is not the parent grid.

bool Grid::coarse_shared_face (Grid & coarse, 
			       int & axis, int & face, 
			       int & il0, int & il1, 
			       int & iu0, int & iu1,
			       int iperiod[3],
			       int & count) throw ()
{

  Grid & grid = *this;

  // Get grid index bounds

  int index_grid[3][2];
  grid.indices(index_grid);

  // Get coarse index bounds

  int index_coarse[3][2];
  coarse.indices(index_coarse);

  // Multiply index_coarse by r_factor_ to match index_grid

  for (axis=0; axis<3; axis++) {
    for (face=0; face<2; face++) {
      index_coarse[axis][face] *= r_factor_;
    }
  }

  // Find count'th matching face, if there is one

  bool found_face = false;
  int  iaxis      = -1;
  int  iface      = -1;
  int  num        = 0;

  for (axis=0; axis<3; axis++) {
    for (face=0; face<2; face++) {

      bool b1 = index_grid[axis][face] == index_coarse[axis][1-face];
      bool b2 = 
	(  index_grid[axis][face]   + (1-face)*iperiod[axis]) == 
	(index_coarse[axis][1-face] +     face*iperiod[axis]);

      if ((b1 || b2) && ! found_face) {
	if (num == count) {
	  found_face = true;
	  iaxis=axis;
	  iface=face;
	  count++;
	} else {
	  num++;
	}
      }
    }
  }

  axis=iaxis;
  face=iface;

  // Exit if face isn't found

  if (! found_face) return false;

  // face axes

  int j0=(axis+1)%3;
  int j1=(axis+2)%3;

  // Compute local indices intersection from global indices of each grid

  il0 = MAX(index_grid[j0][0],index_coarse[j0][0]) - index_grid[j0][0];
  iu0 = MIN(index_grid[j0][1],index_coarse[j0][1]) - index_grid[j0][0] - 1;

  il1 = MAX(index_grid[j1][0],index_coarse[j1][0]) - index_grid[j1][0];
  iu1 = MIN(index_grid[j1][1],index_coarse[j1][1]) - index_grid[j1][0] - 1;

  return true;
}

//----------------------------------------------------------------------

/// Determine the "count"th axis (indexing from 0), face and
/// corresponding range of grid indices of zones adjacent to
/// the containing parent grid, and increment "count".  Returns true
/// if the returned values are valid, or false if there is no
/// "count"th face.

bool Grid::parent_shared_face (Grid & parent, 
			       int & axis, int & face, 
			       int & il0, int & il1, 
			       int & iu0, int & iu1,
			       int & count) throw ()
{

  Grid & grid = *this;

  // Get grid index bounds

  int index_grid[3][2];
  grid.indices(index_grid);

  // Get parent index bounds

  int index_parent[3][2];
  parent.indices(index_parent);

  // Divide index_grid by r_factor_ to match index_parent

  for (axis=0; axis<3; axis++) {
    for (face=0; face<2; face++) {
      index_grid[axis][face] /= r_factor_;
    }
  }

  // Find count'th matching face, if there is one

  bool found_face = false;
  int  iaxis      = -1;
  int  iface      = -1;
  int  num        = 0;

  for (axis = 0; axis < 3; axis++) {
    for (face = 0; face < 2; face++) {
      if (index_grid[axis][face] == index_parent[axis][face] && !found_face) {
	if (num == count) {
	  found_face = true;
	  iaxis = axis;
	  iface = face;
	  count++;
	} else {
	  num++;
	}
      }
    }
  }

  axis = iaxis;
  face = iface;

  // Exit if face isn't found

  if (! found_face) return false;

  // face axes

  int j0=(axis+1)%3;
  int j1=(axis+2)%3;

  // Compute local indices intersection from global indices of each grid
  // Divide by r so that indices correspond to coarse grid

  il0 = MAX(index_grid[j0][0],index_parent[j0][0]) - index_parent[j0][0];
  iu0 = MIN(index_grid[j0][1],index_parent[j0][1]) - index_parent[j0][0] - 1;

  il1 = MAX(index_grid[j1][0],index_parent[j1][0]) - index_parent[j1][0];
  iu1 = MIN(index_grid[j1][1],index_parent[j1][1]) - index_parent[j1][0] - 1;

  return true;

}

//----------------------------------------------------------------------

/// Determine the "count"th axis (indexing from 0), face and
/// corresponding range of coarse-grid indices of zones adjacent to
/// the interior of the parent grid, and increment "count".  Returns true
/// if the returned values are valid, or false if there is no
/// "count"th face.   Indices are relative to the grid.

bool Grid::parent_interior_face (Grid & parent, 
				 int & axis, int & face, 
				 int & il0, int & il1, 
				 int & iu0, int & iu1,
				 int & count) throw ()
{

  Grid & grid = *this;

  // Get grid index bounds

  int index_grid[3][2];
  grid.indices(index_grid);

  // Get parent index bounds

  int index_parent[3][2];
  parent.indices(index_parent);

  // Multiply index_parent by r_factor_ to match index_grid

  for (axis=0; axis<3; axis++) {
    for (face=0; face<2; face++) {
      index_parent[axis][face] *= r_factor_;
    }
  }

  // Find count'th matching face, if there is one

  bool found_face = false;
  int  iaxis      = -1;
  int  iface      = -1;
  int  num        = 0;

  for (axis = 0; axis < 3; axis++) {
    for (face = 0; face < 2; face++) {
      if (index_grid[axis][face] != index_parent[axis][face] && !found_face) {
	if (num == count) {
	  found_face = true;
	  iaxis = axis;
	  iface = face;
	  count++;
	} else {
	  num++;
	}
      }
    }
  }

  axis = iaxis;
  face = iface;

  // Exit if face isn't found

  if (! found_face) return false;

  // face axes

  int j0=(axis+1)%3;
  int j1=(axis+2)%3;

  // Compute local indices intersection from global indices of each grid

  il0 = MAX(index_grid[j0][0],index_parent[j0][0]) - index_grid[j0][0];
  iu0 = MIN(index_grid[j0][1],index_parent[j0][1]) - index_grid[j0][0] - 1;

  il1 = MAX(index_grid[j1][0],index_parent[j1][0]) - index_grid[j1][0];
  iu1 = MIN(index_grid[j1][1],index_parent[j1][1]) - index_grid[j1][0] - 1;

  return true;

}

//----------------------------------------------------------------------


