#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>

#include <string>

#include "cello.h"

#include "array.hpp"
#include "disk.hpp"
#include "node4.h"
#include "tree4.h"

const bool debug = false;
const int  cell_size = 1;
const int  line_width = 1;
const int  gray_threshold = 127;
const int  max_level = 12;

#include "image.h"

// read in the gimp-generated image data into a level array
// values are set to [0:max_levels)

int * create_level_array (int * n0, int * n1, int max_levels)
{
  int size = (width > height) ? width: height;

  int * level_array = new int [size*size];
  for (int i=0; i<size*size; i++) level_array[i] = false;

  *n0 = size;
  *n1 = size;

  int pixel[3];
  char * data = header_data;

  int ix0 = (size - width) / 2;
  int iy0 = (size - height) / 2;

  int max = 0;
  for (int iy=0; iy<height; iy++) {
    for (int ix=0; ix<width; ix++) {
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

void write_image(float * image, int nx, int ny, const char * filename)
{
  Hdf5 hdf5;
  hdf5.file_open(filename,"w");
  ArraySerial tree_array (image,nx,ny,1);
  hdf5.dataset_open ("tree_image",tree_array);
  hdf5.write(tree_array);
  hdf5.dataset_close ();
  hdf5.file_close();
}

main()
{
  int n0,n1;

  // read in the gimp image into level

  int * level_array = create_level_array(&n0,&n1,max_level);

  // Create a tree with full nodes and refine to the image

  Tree4 * tree = new Tree4;
  bool full_nodes;
  tree->refine(level_array,n0,n1,max_level,full_nodes=true);
  printf ("Levels = %d\n",tree->levels());

  // Determine image size
  int n = cell_size + 2*line_width;
  for (int i=0; i<tree->levels(); i++) {
    n = 2*n - line_width;
  }


  float * image;

  printf ("\n");
  image = tree->create_image(n,line_width);
  printf ("nodes  = %d\n",Node4::num_nodes());
  write_image(image,n,n,"tree-0.hdf5");
  delete image;

  // Normalize the tree

  printf ("\n");
  tree->normalize(full_nodes);
  printf ("nodes  = %d\n",Node4::num_nodes());
  image = tree->create_image(n,line_width);
  write_image(image,n,n,"tree-1.hdf5");
  delete image;

  // Optimize the tree

  printf ("\n");
  tree->optimize(full_nodes);
  printf ("nodes  = %d\n",Node4::num_nodes());
  image = tree->create_image(n,line_width);
  write_image(image,n,n,"tree-2.hdf5");
  delete image;

  // Create a new tree and refine to the image with non-full nodes

  delete tree;

  tree = new Tree4;
  tree->refine(level_array,n0,n1,max_level,full_nodes=false);
  printf ("Levels = %d\n",tree->levels());

  // Determine image size
  n = cell_size + 2*line_width;
  for (int i=0; i<tree->levels(); i++) {
    n = 2*n - line_width;
  }

  printf ("\n");
  image = tree->create_image(n,line_width);
  printf ("nodes  = %d\n",Node4::num_nodes());
  write_image(image,n,n,"tree-3.hdf5");
  delete image;

  // Normalize the tree

  printf ("\n");
  tree->normalize(full_nodes);
  printf ("nodes  = %d\n",Node4::num_nodes());
  image = tree->create_image(n,line_width);
  write_image(image,n,n,"tree-4.hdf5");
  delete image;

  // Optimize the tree

  printf ("\n");
  tree->optimize(full_nodes);
  printf ("nodes  = %d\n",Node4::num_nodes());
  image = tree->create_image(n,line_width);
  write_image(image,n,n,"tree-5.hdf5");
  delete image;
  
  delete tree;



}
