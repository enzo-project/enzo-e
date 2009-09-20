#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>

#include <string>

#include "cello.h"

#include "array.hpp"
#include "disk.hpp"
#include "node4.h"
#include "tree4.h"
#include "node16.h"
#include "tree16.h"

const bool debug = false;
const int  cell_size = 1;
const int  line_width = 1;
const int  gray_threshold = 127;
const int  max_level = 8;

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

  // 1 

  Tree4 * tree4 = new Tree4;
  bool full_nodes;
  tree4->refine(level_array,n0,n1,max_level,full_nodes=true);

  // Determine image size
  int n = cell_size + 2*line_width;
  for (int i=0; i<tree4->levels(); i++) {
    n = 2*n - line_width;
  }


  float * image;

  printf ("\n");
  printf ("INITIAL FULL TREE 4\n");
  image = tree4->create_image(n,line_width);
  printf ("nodes  = %d\n",Node4::num_nodes());
  printf ("levels = %d\n",tree4->levels());
  write_image(image,n,n,"tree4-0.hdf5");
  delete image;

  // Normalize the tree

  printf ("\n");
  tree4->normalize(full_nodes);
  printf ("NORMALIZED FULL TREE 4\n");
  printf ("nodes  = %d\n",Node4::num_nodes());
  printf ("levels = %d\n",tree4->levels());
  image = tree4->create_image(n,line_width);
  write_image(image,n,n,"tree4-1.hdf5");
  delete image;

  // Optimize the tree

  printf ("\n");
  tree4->optimize(full_nodes);
  printf ("OPTIMIZED FULL TREE 4\n");
  printf ("nodes  = %d\n",Node4::num_nodes());
  printf ("levels = %d\n",tree4->levels());
  image = tree4->create_image(n,line_width);
  write_image(image,n,n,"tree4-2.hdf5");
  delete image;

  delete tree4;

  // Create a new tree and refine to the image with non-full nodes

  // 2

  tree4 = new Tree4;
  tree4->refine(level_array,n0,n1,max_level,full_nodes=false);

  // Determine image size
  n = cell_size + 2*line_width;
  for (int i=0; i<tree4->levels(); i++) {
    n = 2*n - line_width;
  }

  printf ("\n");
  printf ("INITIAL NON-FULL TREE 4\n");
  image = tree4->create_image(n,line_width);
  printf ("nodes  = %d\n",Node4::num_nodes());
  printf ("levels = %d\n",tree4->levels());
  write_image(image,n,n,"tree4-3.hdf5");
  delete image;

  // Normalize the tree

  printf ("\n");
  printf ("NORMALIZED NON-FULL TREE 4\n");
  tree4->normalize(full_nodes);
  printf ("nodes  = %d\n",Node4::num_nodes());
  printf ("levels = %d\n",tree4->levels());
  image = tree4->create_image(n,line_width);
  write_image(image,n,n,"tree4-4.hdf5");
  delete image;

  // Optimize the tree

  printf ("\n");
  printf ("OPTIMIZED NON-FULL TREE 4\n");
  tree4->optimize(full_nodes);
  printf ("nodes  = %d\n",Node4::num_nodes());
  printf ("levels = %d\n",tree4->levels());
  image = tree4->create_image(n,line_width);
  write_image(image,n,n,"tree4-5.hdf5");
  delete image;
  
  delete tree4;

  // Create a tree with full nodes and refine to the image

  // 3

  Tree16 * tree16 = new Tree16;
  full_nodes;
  tree16->refine(level_array,n0,n1,max_level,full_nodes=true);

  // Determine image size
  n = cell_size + 2*line_width;
  for (int i=0; i<tree16->levels(); i++) {
    n = 2*n - line_width;
  }


  printf ("\n");
  printf ("INITIAL FULL TREE 16\n");
  image = tree16->create_image(n,line_width);
  printf ("TREE 16\n");
  printf ("nodes  = %d\n",Node16::num_nodes());
  printf ("levels = %d\n",tree16->levels());
  write_image(image,n,n,"tree16-0.hdf5");
  delete image;

  // Normalize the tree

  printf ("\n");
  tree16->normalize(full_nodes);
  printf ("NORMALIZED FULL TREE 16\n");
  printf ("nodes  = %d\n",Node16::num_nodes());
  printf ("levels = %d\n",tree16->levels());
  image = tree16->create_image(n,line_width);
  write_image(image,n,n,"tree16-1.hdf5");
  delete image;

  // Optimize the tree

  printf ("\n");
  tree16->optimize(full_nodes);
  printf ("OPTIMIZED FULL TREE 16\n");
  printf ("nodes  = %d\n",Node16::num_nodes());
  printf ("levels = %d\n",tree16->levels());
  image = tree16->create_image(n,line_width);
  write_image(image,n,n,"tree16-2.hdf5");
  delete image;

  delete tree16;

  // Create a new tree and refine to the image with non-full nodes

  // 4

  tree16 = new Tree16;
  tree16->refine(level_array,n0,n1,max_level,full_nodes=false);

  // Determine image size
  n = cell_size + 2*line_width;
  for (int i=0; i<tree16->levels(); i++) {
    n = 2*n - line_width;
  }

  printf ("\n");
  image = tree16->create_image(n,line_width);
  printf ("INITIAL NON-FULL TREE 16\n");
  printf ("nodes  = %d\n",Node16::num_nodes());
  printf ("levels = %d\n",tree16->levels());
  write_image(image,n,n,"tree16-3.hdf5");
  delete image;

  // Normalize the tree

  printf ("\n");
  tree16->normalize(full_nodes);
  printf ("NORMALIZED NON-FULL TREE 16\n");
  printf ("nodes  = %d\n",Node16::num_nodes());
  printf ("levels = %d\n",tree16->levels());
  image = tree16->create_image(n,line_width);
  write_image(image,n,n,"tree16-4.hdf5");
  delete image;

  // Optimize the tree

  printf ("\n");
  tree16->optimize(full_nodes);
  printf ("OPTIMIZED NON-FULL TREE 16\n");
  printf ("nodes  = %d\n",Node16::num_nodes());
  printf ("levels = %d\n",tree16->levels());
  image = tree16->create_image(n,line_width);
  write_image(image,n,n,"tree16-5.hdf5");
  delete image;
  
  delete tree16;



}
