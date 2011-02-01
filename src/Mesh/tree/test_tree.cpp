// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_tree.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-28
/// @brief    Test program for Tree4 and Tree16 classes

#include "test.hpp"

#include "mesh.hpp"

#include "image.h"

/// @brief    Width in pixels of finest grid patches in image
const int  cell_size = 1;

/// @brief    Width of lines in image
const int  line_width = 1;

/// @brief    Bound on number of levels in tree, including root level
const int  max_level = 10;


//----------------------------------------------------------------------

int * create_level_array (int * n0, int * n1, int max_levels)
/// @brief    Read the gimp-generated image data in "image.h" into a level array with values 0 <= x < max_levels
/// @param    n0        Output width parameter (equals n1)
/// @param    n1        Output height parameter (equals n0)
/// @param    max_levels Bound on number of levels in tree, including root level
{
  // width and height are defined in image.h
  int size = (width > height) ? width: height;

  int * level_array = new int [size*size];

  for (int i=0; i<size*size; i++) level_array[i] = false;

  *n0 = size;
  *n1 = size;

  char * data = header_data;

  // If image is not square, find offset to center it in level_array
  int ix0 = (size - width) / 2;
  int iy0 = (size - height) / 2;

  for (size_t iy=0; iy<height; iy++) {
    for (size_t ix=0; ix<width; ix++) {
      int pixel[3];
      HEADER_PIXEL(data,pixel);
      int i = (iy+iy0) + size*(ix+ix0);
      float r = 1.0*pixel[0]/256;
      float g = 1.0*pixel[1]/256;
      float b = 1.0*pixel[2]/256;
      level_array[i] = max_levels * (r + g + b) / 3;
    }
  }
  return level_array;
}

//----------------------------------------------------------------------

void write_image(float * image, int nx, int ny, int k, int full, int step)
/// @brief    Write an image to an hdf5 file
/// @param    image     Image array to write to file
/// @param    nx        Image width
/// @param    ny        Image height
/// @param    k         Refinement factor
/// @param    full      Whether all tree nodes are fully refined
/// @param    step      Step in the refinement process

{
  FileHdf5 hdf5;
  char filename[40];
  sprintf (filename,"tree%d-f%d-%d-old.hdf5",k,full,step);
  hdf5.file_open(filename,"w");
  hdf5.dataset_open_write ("tree_image",nx,ny,1);
  hdf5.write(image);
  hdf5.dataset_close ();
  hdf5.file_close();
  delete [] image;
}

//----------------------------------------------------------------------

void write_stats (int num_nodes, int num_levels)
{
  Memory * memory = Memory::instance();
  printf ("nodes  = %d\n",num_nodes);
  printf ("levels = %d\n",num_levels);
  printf ("bytes/node = %g\n",
	  (float) memory->bytes(0) / num_nodes);
}

//----------------------------------------------------------------------

int image_size(int cell_size, int line_width, int num_levels)
{
  int n = cell_size + 2*line_width;
  for (int i=0; i<num_levels; i++) {
    n = 2*n - line_width;
  }
  return n;
}

void memory_start()
{
  //  Memory::set_active(true);
}

void memory_stop()
{
  //  Memory::print();
  //  Memory::set_active(false);
}

void memory_clear()
{
  // leads to negative numbers
  // Memory::reset();
}

void write_header (int k, int d, int full)
{
  printf ("--------------------------------------------------\n");
  printf ("Refinement factor = %d\n",k);
  printf ("Dimension         = %d\n",d);
  printf ("Full nodes        = %d\n",full);
  printf ("--------------------------------------------------\n");
}

//----------------------------------------------------------------------

int main(int argc, char ** argv)
{
  int n0,n1,n;

  // read in the gimp image into level_array

  int * level_array = create_level_array(&n0,&n1,max_level);

  //=================================================================

  Tree4 * tree4;
  
  for (int full_nodes = 0; full_nodes <= 1; full_nodes++) {

    tree4      = new Tree4;

    write_header (2,2,full_nodes);

    printf ("\nINITIAL\n");
    tree4->refine(level_array,n0,n1,max_level,full_nodes);
    write_stats (tree4->num_nodes(), tree4->levels());
    n = image_size(cell_size,line_width, tree4->levels());
    write_image(tree4->create_image(n,line_width),n,n,2,full_nodes,0);

    printf ("\nBALANCED\n");
    tree4->balance(full_nodes);
    write_stats (tree4->num_nodes(), tree4->levels());
    write_image(tree4->create_image(n,line_width),n,n,2,full_nodes,1);

    printf ("\nCOALESCED\n");
    tree4->optimize();
    write_stats (tree4->num_nodes(), tree4->levels());
    write_image(tree4->create_image(n,line_width),n,n,2,full_nodes,2);

    delete tree4;
  }

  //=================================================================

  Tree16 * tree16;

  for (int full_nodes = 0; full_nodes <= 1; full_nodes++) {

    tree16     = new Tree16;

    write_header (4,2,full_nodes);

    printf ("\nINITIAL\n");
    tree16->refine(level_array,n0,n1,max_level,full_nodes);
    write_stats (tree16->num_nodes(), tree16->levels());
    n = image_size(cell_size,line_width, tree16->levels());
    write_image(tree16->create_image(n,line_width),n,n,4,full_nodes,0);

    printf ("\nBALANCED\n");
    tree16->balance(full_nodes);
    write_stats (tree16->num_nodes(), tree16->levels());
    write_image(tree16->create_image(n,line_width),n,n,4,full_nodes,1);

    printf ("\nCOALESCED\n");
    tree16->optimize();
    write_stats (tree16->num_nodes(), tree16->levels());
    write_image(tree16->create_image(n,line_width),n,n,4,full_nodes,2);

    delete tree16;
  }

}
