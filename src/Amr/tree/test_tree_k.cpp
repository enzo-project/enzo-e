/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @bug       
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */

#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>

#include <string>

#include "cello.h"

#include "error.hpp"
#include "memory.hpp"
#include "monitor.hpp"
#include "disk.hpp"
#include "amr_node2k.hpp"
#include "amr_node3k.hpp"
#include "amr_tree2k.hpp"
#include "amr_tree3k.hpp"

#define index(ix,iy,iz,n) ((ix) + (n)*((iy) + (n)*(iz)))

const bool debug    = false;
const bool geomview = false;

const int  cell_size = 1;
const int  line_width = 1;
const int  gray_threshold = 127;
const int  max_level = 6;
const int sphere_size = 128;

#include "image.h"

//----------------------------------------------------------------------

int * create_level_array (int * n0, int * ny, int max_levels);
int * create_level_array3 (int * n3, int max_levels);
int * create_sphere (int n3, int max_levels);
void write_image(std::string filename, float * image, int nx, int ny, int nz=0);

void create_tree ( int * level_array, int nx, int ny, int nz, int k,  int d, 
		   std::string name, bool full_nodes );

//----------------------------------------------------------------------

int main(int argc, char ** argv)
{
  // read in the gimp image into level

  int nx,ny,nz;
  int * level_array = create_level_array(&nx,&ny,max_level);
  nz = 1;

  //--------------------------------------------------
  // 2D tests
  //--------------------------------------------------

  int k,d;
  create_tree (level_array, nx, ny, nz, k=2, d=2, "tree2-2-f0", false);
  create_tree (level_array, nx, ny, nz, k=2, d=2, "tree2-2-f1", true);
  create_tree (level_array, nx, ny, nz, k=4, d=2, "tree2-4-f1", true);
  create_tree (level_array, nx, ny, nz, k=4, d=2, "tree2-4-f0", false);
  create_tree (level_array, nx, ny, nz, k=8, d=2, "tree2-8-f1", true);
  create_tree (level_array, nx, ny, nz, k=8, d=2, "tree2-8-f0", false);
//   create_tree (level_array, nx, ny, nz, k=16, d=2, "tree2-16-f1", true);
//   create_tree (level_array, nx, ny, nz, k=16, d=2, "tree2-16-f0", false);

  delete [] level_array;

  //--------------------------------------------------
  // 3D tests
  //--------------------------------------------------

  int n3;
  level_array = create_sphere(n3 = sphere_size,max_level);
  // level_array = create_level_array3(&n3,max_level);

  create_tree (level_array, n3, n3, n3, k=2, d=3, "tree3-2-f1", true);
  create_tree (level_array, n3, n3, n3, k=2, d=3, "tree3-2-f0", false);
  create_tree (level_array, n3, n3, n3, k=4, d=3, "tree3-4-f1", true);
  create_tree (level_array, n3, n3, n3, k=4, d=3, "tree3-4-f0", false);
  create_tree (level_array, n3, n3, n3, k=8, d=3, "tree3-8-f1", true);
  create_tree (level_array, n3, n3, n3, k=8, d=3, "tree3-8-f0", false);
//   create_tree (level_array, n3, n3, n3, k=16, d=3, "tree3-16-f1", true);
//   create_tree (level_array, n3, n3, n3, k=16, d=3, "tree3-16-f0", false);

  delete [] level_array;

  Memory::print();

}

//----------------------------------------------------------------------

// read in the gimp-generated image data into a level array
// values are set to [0:max_levels)

int * create_level_array (int * nx, int * ny, int max_levels)
{

  int size = (width > height) ? width: height;

  int * level_array = new int [size*size];

  for (int i=0; i<size*size; i++) level_array[i] = 0;

  *nx = size;
  *ny = size;

  int pixel[3];
  char * data = header_data;

  int ix0 = (size - width) / 2;
  int iy0 = (size - height) / 2;

  int max = 0;
  for (size_t iy=0; iy<height; iy++) {
    for (size_t ix=0; ix<width; ix++) {
      HEADER_PIXEL(data,pixel);
      int i = (iy+iy0) + size*(ix+ix0);
      float r = 1.0*pixel[0]/256;
      float g = 1.0*pixel[1]/256;
      float b = 1.0*pixel[2]/256;
      level_array[i] = max_levels * (r + g + b) / 3;
      if (level_array[i] > max) max = level_array[i];
    }
  }
  return level_array;
}

//----------------------------------------------------------------------
// read in the gimp-generated image data into a level array
// values are set to [0:max_levels)

int * create_level_array3 (int * n3, int max_levels)
{

  if (width != height) {
    printf ("%s:%d width = %d  height = %d\n",__FILE__,__LINE__,width, height);
    exit(1);
  }

  *n3 = width;
  int n = *n3;

  int * level_array = new int [n*n*n];

  for (int i=0; i<n*n*n; i++) level_array[i] = 0;

  int pixel[3];
  char * data = header_data;
 
  float r = 0.125;  // width of the 2D image in the 3D cube
  int nxm = n*(1.0-r)/2;
  int nxp = n*(1.0+r)/2;

  for (int iz=0; iz<n; iz++) {
    for (int iy=0; iy<n; iy++) {
      HEADER_PIXEL(data,pixel);
      for (int ix=0; ix<nxm; ix++) {
	level_array[index(iz,iy,ix,n)] = 0;
      }
      for (int ix=nxm; ix<nxp; ix++) {
	float r = 1.0*pixel[0]/256;
	float g = 1.0*pixel[1]/256;
	float b = 1.0*pixel[2]/256;
	int value = max_levels * (r + g + b) / 3;
	level_array[index(iz,iy,ix,n)] = value;
      }
      for (int ix=nxp; ix<n; ix++) {
	level_array[index(iz,iy,ix,n)] = 0;
      }
    }
  }
  return level_array;
}

//----------------------------------------------------------------------

int * create_sphere (int n3, int max_levels)
{

  int * level_array = new int [n3*n3*n3];

  for (int i=0; i<n3*n3*n3; i++) level_array[i] = 0;

  const double R = 0.3;  // radius
  double R2 = R*R;
  
  double x,y,z;

  for (int iz=0; iz<n3/2; iz++) {
    z = double(iz) / n3 - 0.5;
    for (int iy=0; iy<n3/2; iy++) {
      y = double(iy) / n3 - 0.5;
      for (int ix=0; ix<n3/2; ix++) {
	x = double(ix) / n3 - 0.5;
	double r2 = x*x + y*y + z*z;
	double v = r2 < R2 ? max_levels : 0;

	level_array[index(     ix,     iy,     iz,n3)] = v;
	level_array[index(n3-ix-1,     iy,     iz,n3)] = v;
	level_array[index(     ix,n3-iy-1,     iz,n3)] = v;
	level_array[index(n3-ix-1,n3-iy-1,     iz,n3)] = v;
	level_array[index(     ix,     iy,n3-iz-1,n3)] = v;
	level_array[index(n3-ix-1,     iy,n3-iz-1,n3)] = v;
	level_array[index(     ix,n3-iy-1,n3-iz-1,n3)] = v;
	level_array[index(n3-ix-1,n3-iy-1,n3-iz-1,n3)] = v;
      }
    }
  }
  return level_array;
}

//----------------------------------------------------------------------

void write_image(std::string filename, float * image, int nx, int ny, int nz)
{
  if (nx > 8192 || ny > 8192 || nz > 8192) {
    printf ("%s:%d (nx,ny,nz) = (%d,%d,%d)\n",__FILE__,__LINE__,nx,ny,nz);
    exit(1);
  }
  Hdf5 hdf5;
  hdf5.file_open((filename+".hdf5").c_str(),"w");
  printf ("write_image %d %d %d\n",nx,ny,nz);
  hdf5.dataset_open_write ("tree_image",nx,ny,nz);
  hdf5.write(image);
  hdf5.dataset_close ();
  hdf5.file_close();

  Monitor monitor;
  float min=image[0];
  float max=image[0];
  for (int i=0; i<nx*ny*nz; i++) {
    if (min > image[i]) min = image[i];
    if (max < image[i]) max = image[i];
  }
  //  int color_map[] = {0,0,0,1,1,1};
  //  monitor.plot_png ((filename+".png").c_str(),image,nx,ny,min,max,color_map, 2);

}

//----------------------------------------------------------------------

void create_tree 
(
 int * level_array, 
 int nx, int ny, int nz,
 int k,  int d, 
 std::string name, 
 bool full_nodes
 )
{

  float * image;

  TreeK * tree = 0;
  if (d==2) tree = new Tree2K(k);
  if (d==3) tree = new Tree3K(k);

  int num_nodes;
  float mem_per_node;

  Memory::reset();

  printf ("--------------------------------------------------\n");
  printf ("k=%d d=%d full=%d\n",k,d,full_nodes);
  printf ("--------------------------------------------------\n");

  //--------------------------------------------------
  // Refine the tree
  //--------------------------------------------------

  printf ("\nINITIAL TREE\n");

  Memory::set_active(true);
  tree->refine(level_array,nx,ny,nz,max_level,full_nodes);
  Memory::print();
  Memory::set_active(false);


  // Determine image size
  int image_size = cell_size + 2*line_width;
  for (int i=0; i<tree->levels(); i++) {
    image_size = 2*image_size - line_width;
  }

  if (d==2) {
    image = tree->create_image(image_size,line_width);
    write_image(name + "-0",image,image_size,image_size,1);
  } else {
    image = tree->create_image(image_size,line_width,0);
    write_image(name + "-0-x",image,image_size,image_size,1);
    image = tree->create_image(image_size,line_width,1);
    write_image(name + "-0-y",image,image_size,image_size,1);
    image = tree->create_image(image_size,line_width,2);
    write_image(name + "-0-z",image,image_size,image_size,1);
  }
  if (geomview) tree->geomview(name + "-0.gv");

  num_nodes = (d==2) ? Node2K::num_nodes() : Node3K::num_nodes();
  mem_per_node = (float) Memory::bytes(0) / num_nodes;
  printf ("nodes      = %d\n",num_nodes);
  printf ("levels     = %d\n",tree->levels());
  printf ("bytes/node = %g\n",mem_per_node);

  delete [] image;

  //--------------------------------------------------
  // Balance the tree
  //--------------------------------------------------

  printf ("\nBALANCED TREE\n");

  Memory::set_active(true);
  tree->balance(full_nodes);
  Memory::print();
  Memory::set_active(false);


  image = tree->create_image(image_size,line_width);
  write_image(name + "-1",image,image_size,image_size,1);
  if (geomview) tree->geomview(name + "-1.gv");

  num_nodes = (d==2) ? Node2K::num_nodes() : Node3K::num_nodes();
  mem_per_node = (float) Memory::bytes(0) / num_nodes;
  printf ("nodes      = %d\n",num_nodes);
  printf ("levels     = %d\n",tree->levels());
  printf ("bytes/node = %g\n",mem_per_node);

  delete image;


  //--------------------------------------------------
  // Optimize the tree
  //--------------------------------------------------

  printf ("\nOPTIMIZED TREE\n");

  Memory::set_active(true);
  tree->optimize();
  Memory::print();
  Memory::set_active(false);

  image = tree->create_image(image_size,line_width);
  write_image(name + "-2",image,image_size,image_size,1);
  if (geomview) tree->geomview(name + "-2.gv");

  num_nodes = (d==2) ? Node2K::num_nodes() : Node3K::num_nodes();
  mem_per_node = (float) Memory::bytes(0) / num_nodes;
  printf ("nodes      = %d\n",num_nodes);
  printf ("levels     = %d\n",tree->levels());
  printf ("bytes/node = %g\n",mem_per_node);

  delete image;

  delete tree;

}
