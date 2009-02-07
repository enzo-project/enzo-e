//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Generate a PNG projections of a set of grid files

/**
 * 
 * @file      hypre-png.cpp
 * @brief     Create png projections of a collection of hypre-solve grid files
 * @author    James Bordner
 * @bug       none
 *
 * $Id$
 *
 */

#include <stdio.h>
#include <mpi.h>
#include <assert.h>

#include <string>
#include <vector>

const int debug = 0;
const int trace = 0;

#include "scalar.hpp"
#include "mpi.hpp"
#include "faces.hpp"
#include "domain.hpp"
#include "grid.hpp"
#include "hypre-solve.hpp"
#include "pngwriter.h"
#include "colormap.h"

#define MAX_LEVELS 5
#define TRACE if (trace) printf ("%s:%d TRACE \n",__FILE__,__LINE__);

void images_allocate(double *images3[3],
		     int     n3[3],
		     int     ni3[3],
		     double *rgb33[3][3],
		     int     top_level,
		     int     scale);
void images_grid_add(Grid   * grid,
		     double * images3[3],
		     int      n3[3],
		     int      ni3[3],
		     int      nl3[3],
		     int      nu3[3],
		     int      top_level,
		     int      scale);
void images_colormap (double *images3[3],
		      int     ni3[3],
		      int     rgbmap[3][256],
		      double *rgb33[3][3]);
void images_generate (int     ni3[3],
		      double *rgb33[3][3]);


//======================================================================
// BEGIN MAIN
//======================================================================

int main(int argc, char **argv)
{
  if (argc <= 3) {
    printf ("Usage: %s <colormap-id> <scaling> <grid-file-1> [<grid-file-2> ...]\n",
	    argv[0]);
	  
    exit (1);
  }

  pmpi = new Mpi;

  MPI_Init(&argc,&argv);
  int rgbmap[3][256];

  // Argument 1. Colormap number

  int iarg = 1;
  int i_colormap = atoi(argv[iarg++]);

  for (int line=0; line<256; line++) {
    rgbmap[0][line] = colormap[i_colormap][line][0];
    rgbmap[1][line] = colormap[i_colormap][line][1];
    rgbmap[2][line] = colormap[i_colormap][line][2];
  }

  // Argument 2. scaling

  int scale = atoi(argv[iarg++]);

  // Read in grids, and store them in vectors, one vector per level

  TRACE;

  typedef std::vector<Grid *> gridlist;
  gridlist levels[MAX_LEVELS];
  for (; iarg<argc; iarg++) {
    printf ("Reading %s\n",argv[iarg]);
    FILE * fp = fopen(argv[iarg],"r");
  TRACE;
    Grid * grid = new Grid(fp);
  TRACE;
    int level = grid->level();
  TRACE;
    if (level < MAX_LEVELS) levels[level].push_back(grid);
  TRACE;
    fclose(fp);
  }

  int l,i;
  unsigned j;

  // Determine top_level

  TRACE;

  int top_level=0;
  for (l=MAX_LEVELS-1; l>=0; l--) {
    if (levels[l].size() > 0) {
      top_level = l;
      break;
    }
  }

  // Determine png array sizes

  TRACE;

  int nl3[3],nu3[3],n3[3];  // Index ranges and size
  bool is_first = true;
  for (l=0; l<MAX_LEVELS; l++) {
    for (j=0; j<levels[l].size(); j++) {
      int ind[3][2];
      levels[l][j]->indices(ind);
      int kg = int(pow(2,top_level - levels[l][j]->level())); 
      printf ("Grid %d: (%d %d %d) to (%d %d %d)\n",
	      levels[l][j]->id(),
	      kg*ind[0][0],kg*ind[1][0],kg*ind[2][0],
	      kg*ind[0][1],kg*ind[1][1],kg*ind[2][1]);
      if (is_first) {
	for(i=0; i<3; i++) {
	  nl3[i] = kg*ind[i][0];
	  nu3[i] = kg*ind[i][1];
	}
      } else {
	for(i=0; i<3; i++) {
	  nl3[i] = MIN(nl3[i],kg*ind[i][0]);
	  nu3[i] = MAX(nu3[i],kg*ind[i][1]);
	}
      }
      printf ("%d %d %d  %d %d %d\n",
	      nl3[0],nl3[1],nl3[2],
	      nu3[0],nu3[1],nu3[2]);
      is_first = false;
    }
  }
  // Determine n3[]

  TRACE;

  for (i=0; i<3; i++) {
    n3[i] = nu3[i] - nl3[i];
  }

  // Allocate images

  TRACE;

  double *images3[3];

  int ni3[3];
  double *rgb33[3][3];
  images_allocate(images3,n3,ni3,rgb33,top_level,scale);
  
  // Assemble images

  TRACE;

  // Loop over levels, coarse to fine

  for (l=0; l<=top_level; l++) {

    // Loop over grids in level

    for (j=0; j<levels[l].size(); j++) {

      Grid * grid = levels[l][j];

      // Overlay grid data in image

      images_grid_add(grid,images3,n3,ni3,nl3,nu3,top_level,scale);

    }
  }

  TRACE;

  images_colormap (images3,ni3,rgbmap,rgb33);

  TRACE;

  images_generate (ni3,rgb33);

  TRACE;
}

//----------------------------------------------------------------------

void images_allocate(double *images3[3],
		     int     n3[3],
		     int     ni3[3],
		     double *rgb33[3][3],
		     int     top_level,
		     int     scale)
{
  int k0 = scale;
  for (int axis=0; axis<3; axis++) {
    ni3[axis] = scale*n3[axis];
    int n1 = k0*n3[(axis+1)%3];
    int n2 = k0*n3[(axis+2)%3];
    printf ("n1 n2 = %d %d\n",n1,n2);
    images3[axis] = new double [n1*n2];
    rgb33[axis][0] = new double [n1*n2]; // red
    rgb33[axis][1] = new double [n1*n2]; // green
    rgb33[axis][2] = new double [n1*n2]; // blue
    for (int i=0; i<n1*n2; i++) {
      images3[axis][i] = 0;
      rgb33[axis][0][i] = 0;
      rgb33[axis][1][i] = 0;
      rgb33[axis][2][i] = 0;
    }
  }
}

//----------------------------------------------------------------------

void images_grid_add(Grid   * grid,
		     double * images3[3],
		     int      n3[3],
		     int      ni3[3],
		     int      nl3[3],
		     int      nu3[3],
		     int      top_level,
		     int      scale)
{

  // Add grid projection to images3

  // k*k*k is size of a grid cell in image3

  // scale: image zones / fine level zone
  // kg: fine level zones / this grid zones
  int kg = int(pow(2,top_level - grid->level())); 

  double *u = grid->values();
  int gl3[3];
  grid->i_lower(gl3[0],gl3[1],gl3[2]);

  printf ("grid lower (%d %d %d)\n",gl3[0],gl3[1],gl3[2]);
  printf ("image lower (%d %d %d)\n",nl3[0],nl3[1],nl3[2]);
  printf ("image upper (%d %d %d)\n",nu3[0],nu3[1],nu3[2]);
  printf ("image size (%d %d %d)\n",ni3[0],ni3[1],ni3[2]);

  int ng3[3];
  ng3[0] = grid->n(0);
  ng3[1] = grid->n(1);
  ng3[2] = grid->n(2);
  // i0,i1,i2     grid local index 0 to ng3[*]
  // ii0,ii1,ii2  grid global index nl3 to nu3
  for (int i0=0; i0<ng3[0]; i0++) {
    int ii0 = scale*(i0 + gl3[0])*kg;
    for (int i1=0; i1<ng3[1]; i1++) {
      int ii1 = scale*(i1 + gl3[1])*kg;
      for (int i2=0; i2<ng3[2]; i2++) {
	int ii2 = scale*(i2 + gl3[2])*kg;

	// 0 <= i* < n3[*]
	// nl3[*] <= ii* <= nu3[*]
	
	int ig = i0 + ng3[0]*(i1 + ng3[1]*i2);
	int ix = ii1 + ni3[1]*ii2;
	int iy = ii2 + ni3[2]*ii0;
	int iz = ii0 + ni3[0]*ii1;
	for (int ik0=0; ik0<scale*kg; ik0++) {
	  for (int ik1=0; ik1<scale*kg; ik1++) {
	    for (int ik2=0; ik2<scale*kg; ik2++) {
	      int ikx = ik1 + ni3[1]*ik2;
	      int iky = ik2 + ni3[2]*ik0;
	      int ikz = ik0 + ni3[0]*ik1;
	      assert (ix+ikx >= 0);
	      assert (iy+iky >= 0);
	      assert (iz+ikz >= 0);
	      assert (ix+ikx < ni3[1]*ni3[2]);
	      assert (iy+iky < ni3[2]*ni3[0]);
	      assert (iz+ikz < ni3[0]*ni3[1]);
	      images3[0][ix+ikx] += u[ig];
	      images3[1][iy+iky] += u[ig];
	      images3[2][iz+ikz] += u[ig];
	    }
	  }
	}
      }
    }
  }
}

//----------------------------------------------------------------------

void images_colormap (double *images3[3],
		      int     ni3[3],
		      int     rgbmap[3][256],
		      double *rgb33[3][3])
{
  int axis;
  double rmin3[3],rmax3[3];
  // Determine min and max for each image
  for (axis=0; axis<3; axis++) {
    rmin3[axis] = images3[axis][0];
    rmax3[axis] = images3[axis][0];
    int n1 = ni3[(axis+1)%3];
    int n2 = ni3[(axis+2)%3];
    for (int i=0; i<n1*n2; i++) {
      rmin3[axis] = MIN(rmin3[axis],images3[axis][i]);
      rmax3[axis] = MAX(rmax3[axis],images3[axis][i]);
    }
    int index;
    if (rmax3[axis] > rmin3[axis]) {
      for (int i=0; i<n1*n2; i++) {
	index = int(255*(images3[axis][i]-rmin3[axis])/(rmax3[axis]-rmin3[axis]));
	index = MAX(index,0);
	index = MIN(index,255);

	rgb33[axis][0][i] = rgbmap[0][index] / 255.0;
	rgb33[axis][1][i] = rgbmap[1][index] / 255.0;
	rgb33[axis][2][i] = rgbmap[2][index] / 255.0;
      }
    } else {
      index = 127;
      for (int i=0; i<n1*n2; i++) {
	rgb33[axis][0][i] = rgbmap[0][index] / 255.0;
	rgb33[axis][1][i] = rgbmap[1][index] / 255.0;
	rgb33[axis][2][i] = rgbmap[2][index] / 255.0;
      }
    }
  }
}

//----------------------------------------------------------------------

void images_generate (int     ni3[3],
		      double  *rgb33[3][3])
{

  char filename[40];

  for (int axis=0; axis<3; axis++) {
    sprintf (filename,"project-%d.png",axis);
    int n1 = ni3[(axis+1)%3];
    int n2 = ni3[(axis+2)%3];
    pngwriter png (n1,n2,0,filename);
    for (int i1=0; i1 < n1; i1++) {
      for (int i2=0; i2 < n2; i2++) {
	int i = i1 + n1*i2;
	png.plot(i1+1,i2+1,
		 rgb33[axis][0][i],
		 rgb33[axis][1][i],
		 rgb33[axis][2][i]);
      }
    }
    png.close();
  }

}

//----------------------------------------------------------------------

void images_generate_shift (int     ni3[3],
			    double  *rgb33[3][3])
{

  char filename[40];

  for (int axis=0; axis<3; axis++) {
    sprintf (filename,"project-shift-%d.png",axis);
    int n1 = ni3[(axis+1)%3];
    int n2 = ni3[(axis+2)%3];
    pngwriter png (n1,n2,0,filename);
    for (int i1=0; i1 < n1; i1++) {
      int k1 = (i1 + n1/2) % n1;
      for (int i2=0; i2 < n2; i2++) {
	int k2 = (i2 + n2/2) % n2;
	int i = k1 + n1*k2;
	png.plot(i1+1,i2+1,
		 rgb33[axis][0][i],
		 rgb33[axis][1][i],
		 rgb33[axis][2][i]);
      }
    }
    png.close();
  }

}

