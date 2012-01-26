// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_test.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-12
/// @brief    Functions for Mesh test programs

#include "mesh.hpp"

#include <assert.h>

//----------------------------------------------------------------------

Tree * test_tree_22()
{

  // +-------+---+-+-+  
  // |       |   +-+-+  
  // |       +-+-+-+-+  
  // |       +-+-+   |  
  // +-+-+---+-+-+---+  
  // +-+-+   |   |   |  
  // +-+-+-+-+---+---|  
  // |   +-+-+   |   |  
  // +---+-+-+-------+
  //
  // level 0:   level 1:     level 2:
  // +---+     +---+---+     +---+---+---+---+
  // | 0 |     | 3 | 4 |     |   |   | 23| 24|
  // +---+     +---+---+     +---+---+---+---+
  //           | 1 | 2 |     |   |   | 21| 22|
  //           +---+---+     +---+---+---+---+
  //                         | 7 | 8 | 19| 20|
  //                         +---+---+---+---+
  //                         | 5 | 6 | 17| 18|
  //                         +---+---+---+---+
  // level 3:
  // +---+---+---+---+---+---+---+---+
  // |   |   |   |   |   |   | 31| 32|
  // +---+---+---+---+---+---+---+---+
  // |   |   |   |   |   |   | 29| 30|
  // +---+---+---+---+---+---+---+---+
  // |   |   |   |   | 27| 28|   |   |
  // +---+---+---+---+---+---+---+---+
  // |   |   |   |   | 25| 26|   |   |
  // +---+---+---+---+---+---+---+---+
  // | 15| 16|   |   |   |   |   |   |
  // +---+---+---+---+---+---+---+---+
  // | 13| 14|   |   |   |   |   |   |
  // +---+---+---+---+---+---+---+---+
  // |   |   | 11| 12|   |   |   |   |
  // +---+---+---+---+---+---+---+---+
  // |   |   | 9 | 10|   |   |   |   |
  // +---+---+---+---+---+---+---+---+



  // Refined nodes: root, 0 1 3 01 02 30 33
  // (0 1 2 4 6 7 13 16)

  int d = 2;
  int r = 2;

  int * array = new int [33];
  for (int i=0; i<33; i++) array[i] = i;
  int * a = array;

  Tree * tree = new Tree (d,r);

  NodeTrace node_trace (tree->root_node());
  
  Node * node ;

  tree->root_node()->set_data(((void *) a++));

  // root
  tree->refine_node (node_trace);

  node = node_trace.node();
  node->child(0)->set_data(((void *) a++));
  node->child(1)->set_data(((void *) a++));
  node->child(2)->set_data(((void *) a++));
  node->child(3)->set_data(((void *) a++));

  // 0
  node_trace.push(0);
  tree->refine_node (node_trace);
  node = node_trace.node();
  
  node->child(0)->set_data(((void *) a++));
  node->child(1)->set_data(((void *) a++));
  node->child(2)->set_data(((void *) a++));
  node->child(3)->set_data(((void *) a++));

  // 01
  node_trace.push(1);
  tree->refine_node (node_trace);
  node = node_trace.node();
  
  node->child(0)->set_data(((void *) a++));
  node->child(1)->set_data(((void *) a++));
  node->child(2)->set_data(((void *) a++));
  node->child(3)->set_data(((void *) a++));

  // 02
  node_trace.pop();
  node_trace.push(2);
  tree->refine_node (node_trace);
  node = node_trace.node();
  
  node->child(0)->set_data(((void *) a++));
  node->child(1)->set_data(((void *) a++));
  node->child(2)->set_data(((void *) a++));
  node->child(3)->set_data(((void *) a++));

  // 1
  node_trace.pop();
  node_trace.pop();
  node_trace.push(1);
  tree->refine_node (node_trace);
  node = node_trace.node();
  
  node->child(0)->set_data(((void *) a++));
  node->child(1)->set_data(((void *) a++));
  node->child(2)->set_data(((void *) a++));
  node->child(3)->set_data(((void *) a++));

  // 3
  node_trace.pop();
  node_trace.push(3);
  tree->refine_node (node_trace);
  node = node_trace.node();
  
  node->child(0)->set_data(((void *) a++));
  node->child(1)->set_data(((void *) a++));
  node->child(2)->set_data(((void *) a++));
  node->child(3)->set_data(((void *) a++));

  // 30
  node_trace.push(0);
  tree->refine_node (node_trace);
  node = node_trace.node();

  node->child(0)->set_data(((void *) a++));
  node->child(1)->set_data(((void *) a++));
  node->child(2)->set_data(((void *) a++));
  node->child(3)->set_data(((void *) a++));

  // 33
  node_trace.pop();
  node_trace.push(3);
  tree->refine_node (node_trace);
  node = node_trace.node();
  
  node->child(0)->set_data(((void *) a++));
  node->child(1)->set_data(((void *) a++));
  node->child(2)->set_data(((void *) a++));
  node->child(3)->set_data(((void *) a++));


  return tree;
}

Tree * test_tree_32()
{

  // 3D version of the following;
  // +-------+---+-+-+
  // |   |   |   |   |
  // +---+-+-+-+-+---+
  // |   +-+-+-+-+   |
  // +---+-+-+-+-+---+
  // |   +-+-+-+-+   |
  // +---+-+-+-+-+---+
  // |   |   |   |   |
  // +---+-+-+-------+
  //
  // Refined nodes: root, 0 07 1 16 2 25 3 34 4 43 5 52 6 61 7 70

  int d = 3;
  int r = 2;

  Tree * tree = new Tree (d,r);

  NodeTrace node_trace  (tree->root_node());

  // root
  tree->refine_node (node_trace);

  // 0
  node_trace.push(0);
  tree->refine_node (node_trace);
  // 07
  node_trace.push(7);
  tree->refine_node (node_trace);
  // 1
  node_trace.pop();
  node_trace.pop();
  node_trace.push(1);
  tree->refine_node (node_trace);
  // 16
  node_trace.push(6);
  tree->refine_node (node_trace);
  // 2
  node_trace.pop();
  node_trace.pop();
  node_trace.push(2);
  tree->refine_node (node_trace);
  // 25
  node_trace.push(5);
  tree->refine_node (node_trace);
  // 3
  node_trace.pop();
  node_trace.pop();
  node_trace.push(3);
  tree->refine_node (node_trace);
  // 34
  node_trace.push(4);
  tree->refine_node (node_trace);
  // 4
  node_trace.pop();
  node_trace.pop();
  node_trace.push(4);
  tree->refine_node (node_trace);
  // 43
  node_trace.push(3);
  tree->refine_node (node_trace);
  // 5
  node_trace.pop();
  node_trace.pop();
  node_trace.push(5);
  tree->refine_node (node_trace);
  // 52
  node_trace.push(2);
  tree->refine_node (node_trace);
  // 6
  node_trace.pop();
  node_trace.pop();
  node_trace.push(6);
  tree->refine_node (node_trace);
  // 61
  node_trace.push(1);
  tree->refine_node (node_trace);
  // 7
  node_trace.pop();
  node_trace.pop();
  node_trace.push(7);
  tree->refine_node (node_trace);
  // 70
  node_trace.push(0);
  tree->refine_node (node_trace);

  return tree;
}

//----------------------------------------------------------------------

int * create_levels_from_image (const char * pngfile, 
			  int * nx, int * ny, int max_levels)
// return an array of integer values between [0 and max_levels)
// corresponding to the grayscale values of the input png file,
// indexed by ix + nx*iy
{

  pngwriter png;

  png.readfromfile(pngfile);

  (*nx) = png.getwidth();
  (*ny) = png.getheight();

  int size = (*nx)*(*ny);

  int * level_array = new int [size];

  for (int i=0; i<size; i++) level_array[i] = 0;


  for (int iy=0; iy<*ny; iy++) {
    for (int ix=0; ix<*nx; ix++) {
      int pixel = png.read(ix+1,iy+1);
      int i = ix + (*nx)*iy;
      level_array[i] = max_levels * 1.0*pixel / 256;
    }
  }

  png.close();

  return level_array;
}

//----------------------------------------------------------------------

void create_tree_from_levels 
(
 Tree * tree, 
 int * levels, 
 int nx, int ny, int nz=1
)
{

  Timer timer;
  timer.start();

  int r = tree->refinement();

  // Determine size of finest tree level
  //  int n=1;
  //  for (int i=0; i<max_level; i++) n*=r;

  //  n = nx;
  // Create tree
  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
	int i = ix + nx*(iy + ny*iz);
	int level = levels[i];

	double x = 1.0*ix / nx;
	double y = 1.0*iy / ny;
	double z = 1.0*iz / nz;
	NodeTrace node_trace (tree->root_node());
	while (--level > 0) {
	  if (node_trace.node()->is_leaf()) {
	    tree->refine_node (node_trace);
	  }
	  int irx = x * r;
	  int iry = y * r;
	  int irz = z * r;
	  int ir = irx + r*(iry + r*irz);

	  node_trace.push(ir);

	  x *= r; while (x >= 1.0) x -= 1.0;
	  y *= r; while (y >= 1.0) y -= 1.0;
	  z *= r; while (z >= 1.0) z -= 1.0;
	  
	}
      }
    }
  }
  printf ("create_tree_from_levels() time = %f s\n",timer.value());
}
//----------------------------------------------------------------------

inline void rotate
(double *xr, double *yr, double *zr,
 double x, double y, double z,
 double theta, double phi, double psi,
 bool ortho)
{
  double cph = cos(phi);
  double cps = cos(psi);
  double cth = cos(theta);
  double sph = sin(phi);
  double sps = sin(psi);
  double sth = sin(theta);

  x -= 0.5;
  y -= 0.5;
  z -= 0.5;

  double r11 = cps*cth*cph-sps*sph;
  double r12 = cps*cth*sph+sps*cph;
  double r13 = -cps*sth;

  double r21 = -sps*cth*cph-cps*sph;
  double r22 = -sps*cth*sph+cps*cph;
  double r23 = sps*sth;

  double r31 = sth*cph;
  double r32 = sth*sph;
  double r33 = cth;
    
  if (xr) (*xr) = r11*x + r12*y + r13*z;
  if (yr) (*yr) = r21*x + r22*y + r23*z;
  if (zr) (*zr) = r31*x + r32*y + r33*z;

  if (!ortho && zr) {
    const double distance     = 3.0;
    const double focal_length = 3.0;
    double scale = focal_length / (distance - (*zr));
    if (xr) (*xr) *= scale;
    if (yr) (*yr) *= scale;
    if (zr) (*zr) *= scale;
  }

  if (xr) (*xr) += 0.5;
  if (yr) (*yr) += 0.5;
  if (zr) (*zr) += 0.5;

}

void create_image_from_tree (Tree * tree, const char * filename, 
			     int nx, int ny,
			     int level_lower=0, int level_upper=1000,
			     double theta=0.0, double phi=0.0, double psi=0.0,
			     double scale=1.0, bool ortho=true,
			     int falloff=0)
/// @brief Generate a PNG image of a tree
///
/// @param tree is the 2D or 3D tree from which to generate the image
/// @param filename is the name of the image file, including extension
/// @param nx,ny are the size of the image
/// @param thx,thy,thz   are the rotation angles around x, y, and z axes
/// @param scale is the scaling factor of the image
/// @param ortho determines whether to generate an orthographic or
/// perspective projection
/// @param falloff refers to color blending with respect to levels: (level) ^ (-falloff)
{

  pngwriter png (nx+1,ny+1,0,filename);

  int i;

    // determine color

    //    magenta   blue   cyan   green   yellow   red

    //    r  1.0                            1.0    1.0 
    //    g                 0.5    1.0      1.0    0
    //    b  1.0      1     0.5                    0
    
  const int num_colors = 7;
  const double rc[] = {1, 0, 0, 0, 1, 1, 1};
  const double gc[] = {0, 0, 1, 1, 1, 0, 1};
  const double bc[] = {1, 1, 1, 0, 0, 0, 1};
  int num_levels = tree->max_level() + 1;
  double * ra = new double [num_levels];
  double * ga = new double [num_levels];
  double * ba = new double [num_levels];
  for (int i=0; i<num_levels; i++) {
    double ci = double(i)/num_levels;
    int ic = ci * (num_colors-1);
    double a = ci *(num_colors-1) - ic;
    ra[i] = (1-a)*rc[ic] + a*rc[ic+1];
    ga[i] = (1-a)*gc[ic] + a*gc[ic+1];
    ba[i] = (1-a)*bc[ic] + a*bc[ic+1];
  }

  ItNode it_node (tree,level_lower,level_upper);
  int r = tree->refinement();
  int d = tree->dimension();
  double rinv = 1.0/r;
  int count = 0;
  while ((it_node.next_leaf())) {
    count++;
    const NodeTrace * node_trace  = it_node.node_trace();
    double xmin = 0.0; double xmax = 1.0;
    double ymin = 0.0; double ymax = 1.0;
    double zmin = 0.0; double zmax = 1.0;
    double h = 1.0 / r;
    int level = node_trace->level();
    // determine node boundaries scaled by [0:1,0:1]
    for (i=1; i<=level; i++) {
      int index_curr = node_trace->index_level(i);
      int kx,ky,kz;
      tree->index(index_curr,&kx,&ky,&kz);
      xmin += h*kx;
      ymin += h*ky;
      zmin += h*kz;
      xmax = xmin+h;
      ymax = ymin+h;
      zmax = zmin+h;
      h*=rinv;
    }
    if (d==2) {zmin=0; zmax=0;}
    // Compute the rotated points
    double x000,y000,z000;
    double x001,y001,z001;
    double x010,y010,z010;
    double x011,y011,z011;
    double x100,y100,z100;
    double x101,y101,z101;
    double x110,y110,z110;
    double x111,y111,z111;

    rotate(&x000,&y000,&z000, xmin,ymin,zmin, theta,phi,psi,ortho);
    rotate(&x001,&y001,&z001, xmin,ymin,zmax, theta,phi,psi,ortho);
    rotate(&x010,&y010,&z010, xmin,ymax,zmin, theta,phi,psi,ortho);
    rotate(&x011,&y011,&z011, xmin,ymax,zmax, theta,phi,psi,ortho);
    rotate(&x100,&y100,&z100, xmax,ymin,zmin, theta,phi,psi,ortho);
    rotate(&x101,&y101,&z101, xmax,ymin,zmax, theta,phi,psi,ortho);
    rotate(&x110,&y110,&z110, xmax,ymax,zmin, theta,phi,psi,ortho);
    rotate(&x111,&y111,&z111, xmax,ymax,zmax, theta,phi,psi,ortho);

    // scale points
    
    double amin = MIN (nx,ny);

    x000 = amin * ((x000-0.5)*scale + 0.5) + 0.5*(nx - amin);
    x001 = amin * ((x001-0.5)*scale + 0.5) + 0.5*(nx - amin);
    x010 = amin * ((x010-0.5)*scale + 0.5) + 0.5*(nx - amin);
    x011 = amin * ((x011-0.5)*scale + 0.5) + 0.5*(nx - amin);
    x100 = amin * ((x100-0.5)*scale + 0.5) + 0.5*(nx - amin);
    x101 = amin * ((x101-0.5)*scale + 0.5) + 0.5*(nx - amin);
    x110 = amin * ((x110-0.5)*scale + 0.5) + 0.5*(nx - amin);
    x111 = amin * ((x111-0.5)*scale + 0.5) + 0.5*(nx - amin);

    y000 = amin * ((y000-0.5)*scale + 0.5) + 0.5*(ny - amin);
    y001 = amin * ((y001-0.5)*scale + 0.5) + 0.5*(ny - amin);
    y010 = amin * ((y010-0.5)*scale + 0.5) + 0.5*(ny - amin);
    y011 = amin * ((y011-0.5)*scale + 0.5) + 0.5*(ny - amin);
    y100 = amin * ((y100-0.5)*scale + 0.5) + 0.5*(ny - amin);
    y101 = amin * ((y101-0.5)*scale + 0.5) + 0.5*(ny - amin);
    y110 = amin * ((y110-0.5)*scale + 0.5) + 0.5*(ny - amin);
    y111 = amin * ((y111-0.5)*scale + 0.5) + 0.5*(ny - amin);

    // plot cube

    // compute "opacity"

    double o = pow(1.0*(node_trace->level()+1) / (num_levels), falloff);

    double o000 = o*0.5*(1+z000);
    double o001 = o*0.5*(1+z001);
    double o010 = o*0.5*(1+z010);
    double o011 = o*0.5*(1+z011);
    double o100 = o*0.5*(1+z100);
    double o101 = o*0.5*(1+z101);
    double o110 = o*0.5*(1+z110);
    double o111 = o*0.5*(1+z111);
    int k = node_trace->level();

    png.line_blend(int(x000),int(y000),int(x001),int(y001),o000,ra[k],ga[k],ba[k]);
    png.line_blend(int(x010),int(y010),int(x011),int(y011),o010,ra[k],ga[k],ba[k]);
    png.line_blend(int(x100),int(y100),int(x101),int(y101),o100,ra[k],ga[k],ba[k]);
    png.line_blend(int(x110),int(y110),int(x111),int(y111),o110,ra[k],ga[k],ba[k]);

    png.line_blend(int(x000),int(y000),int(x010),int(y010),o000,ra[k],ga[k],ba[k]);
    png.line_blend(int(x001),int(y001),int(x011),int(y011),o001,ra[k],ga[k],ba[k]);
    png.line_blend(int(x100),int(y100),int(x110),int(y110),o100,ra[k],ga[k],ba[k]);
    png.line_blend(int(x101),int(y101),int(x111),int(y111),o101,ra[k],ga[k],ba[k]);

    png.line_blend(int(x000),int(y000),int(x100),int(y100),o000,ra[k],ga[k],ba[k]);
    png.line_blend(int(x001),int(y001),int(x101),int(y101),o001,ra[k],ga[k],ba[k]);
    png.line_blend(int(x010),int(y010),int(x110),int(y110),o010,ra[k],ga[k],ba[k]);
    png.line_blend(int(x011),int(y011),int(x111),int(y111),o011,ra[k],ga[k],ba[k]);

  }
  delete [] ra;
  delete [] ga;
  delete [] ba;
  png.close();
}

