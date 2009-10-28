#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>

#include <string>

#include "cello.h"

#include "error.hpp"
#include "memory.hpp"
#include "disk.hpp"
#include "amr_node_k.hpp"
#include "amr_tree_k.hpp"

const bool debug = false;
const int  cell_size = 1;
const int  line_width = 1;
const int  gray_threshold = 127;
const int  max_level = 10;

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

void write_image(float * image, int nx, int ny, const char * filename)
{
  Hdf5 hdf5;
  hdf5.file_open(filename,"w");
  hdf5.dataset_open_write ("tree_image",nx,ny,1);
  hdf5.write(image);
  hdf5.dataset_close ();
  hdf5.file_close();
}

int main(int argc, char ** argv)
{
  int n0,n1;

  // read in the gimp image into level

  int * level_array = create_level_array(&n0,&n1,max_level);

  // Create a tree with full nodes and refine to the image

  // 1 

  Tree_k * tree4 = new Tree_k(2);

  bool full_nodes;

  Memory::reset();

  Memory::set_active(true);
  tree4->refine(level_array,n0,n1,max_level,full_nodes=true);
  Memory::print();
  Memory::set_active(false);
  printf ("Bytes / node = %g\n",(float)Memory::bytes(0) / Node_k::num_nodes());

  // Determine image size
  int n = cell_size + 2*line_width;
  for (int i=0; i<tree4->levels(); i++) {
    n = 2*n - line_width;
  }


  float * image;

  printf ("\n");
  printf ("INITIAL FULL TREE 4\n");
  image = tree4->create_image(n,line_width);
  printf ("nodes  = %d\n",Node_k::num_nodes());
  printf ("levels = %d\n",tree4->levels());
  write_image(image,n,n,"tree4-0.hdf5");
  delete image;

  // Balance the tree


  printf ("\n");

  Memory::set_active(true);
  tree4->balance(full_nodes);
  Memory::print();
  Memory::set_active(false);
  printf ("Bytes / node = %g\n",(float)Memory::bytes(0) / Node_k::num_nodes());

  printf ("BALANCED FULL TREE 4\n");
  printf ("nodes  = %d\n",Node_k::num_nodes());
  printf ("levels = %d\n",tree4->levels());
  image = tree4->create_image(n,line_width);
  write_image(image,n,n,"tree4-1.hdf5");
  delete image;

  // Optimize the tree

  printf ("\n");
  Memory::set_active(true);
  tree4->optimize();
  Memory::print();
  Memory::set_active(false);
  printf ("Bytes / node = %g\n",(float)Memory::bytes(0) / Node_k::num_nodes());

  printf ("OPTIMIZED FULL TREE 4\n");
  printf ("nodes  = %d\n",Node_k::num_nodes());
  printf ("levels = %d\n",tree4->levels());
  image = tree4->create_image(n,line_width);
  write_image(image,n,n,"tree4-2.hdf5");
  delete image;

  delete tree4;

  // Create a new tree and refine to the image with non-full nodes

  // 2

  tree4 = new Tree_k(2);

  Memory::set_active(true);
  tree4->refine(level_array,n0,n1,max_level,full_nodes=false);
  Memory::print();
  Memory::set_active(false);
  printf ("Bytes / node = %g\n",(float)Memory::bytes(0) / Node_k::num_nodes());

  Memory::print();

  // Determine image size
  n = cell_size + 2*line_width;
  for (int i=0; i<tree4->levels(); i++) {
    n = 2*n - line_width;
  }

  printf ("\n");
  printf ("INITIAL NON-FULL TREE 4\n");
  image = tree4->create_image(n,line_width);
  printf ("nodes  = %d\n",Node_k::num_nodes());
  printf ("levels = %d\n",tree4->levels());
  write_image(image,n,n,"tree4-3.hdf5");
  delete image;

  // Balance the tree

  printf ("\n");
  printf ("BALANCED NON-FULL TREE 4\n");
  Memory::set_active(true);
  tree4->balance(full_nodes);
  Memory::print();
  Memory::set_active(false);
  printf ("Bytes / node = %g\n",(float)Memory::bytes(0) / Node_k::num_nodes());

  Memory::print();

  printf ("nodes  = %d\n",Node_k::num_nodes());
  printf ("levels = %d\n",tree4->levels());
  image = tree4->create_image(n,line_width);
  write_image(image,n,n,"tree4-4.hdf5");
  delete image;

  // Optimize the tree

  printf ("\n");
  printf ("OPTIMIZED NON-FULL TREE 4\n");

  Memory::set_active(true);
  tree4->optimize();
  Memory::print();
  Memory::set_active(false);
  printf ("Bytes / node = %g\n",(float)Memory::bytes(0) / Node_k::num_nodes());

  printf ("nodes  = %d\n",Node_k::num_nodes());
  printf ("levels = %d\n",tree4->levels());
  image = tree4->create_image(n,line_width);
  write_image(image,n,n,"tree4-5.hdf5");
  delete image;
  
  delete tree4;

  // Create a tree with full nodes and refine to the image

  // 3

  Tree_k * tree16 = new Tree_k(4);
  Memory::set_active(true);
  tree16->refine(level_array,n0,n1,max_level,full_nodes=true);
  Memory::print();
  Memory::set_active(false);
  printf ("Bytes / node = %g\n",(float)Memory::bytes(0) / Node_k::num_nodes());

  // Determine image size
  n = cell_size + 2*line_width;
  for (int i=0; i<tree16->levels(); i++) {
    n = 2*n - line_width;
  }


  printf ("\n");
  printf ("INITIAL FULL TREE 16\n");
  image = tree16->create_image(n,line_width);
  printf ("TREE 16\n");
  printf ("nodes  = %d\n",Node_k::num_nodes());
  printf ("levels = %d\n",tree16->levels());
  write_image(image,n,n,"tree16-0.hdf5");
  delete image;

  // Balance the tree

  printf ("\n");

  Memory::set_active(true);
  tree16->balance(full_nodes);
  Memory::print();
  Memory::set_active(false);
  printf ("Bytes / node = %g\n",(float)Memory::bytes(0) / Node_k::num_nodes());

  printf ("BALANCED FULL TREE 16\n");
  printf ("nodes  = %d\n",Node_k::num_nodes());
  printf ("levels = %d\n",tree16->levels());
  image = tree16->create_image(n,line_width);
  write_image(image,n,n,"tree16-1.hdf5");
  delete image;

  // Optimize the tree

  printf ("\n");

  Memory::set_active(true);
  tree16->optimize();
  Memory::print();
  Memory::set_active(false);
  printf ("Bytes / node = %g\n",(float)Memory::bytes(0) / Node_k::num_nodes());

  printf ("OPTIMIZED FULL TREE 16\n");
  printf ("nodes  = %d\n",Node_k::num_nodes());
  printf ("levels = %d\n",tree16->levels());
  image = tree16->create_image(n,line_width);
  write_image(image,n,n,"tree16-2.hdf5");
  delete image;

  delete tree16;

  // Create a new tree and refine to the image with non-full nodes

  // 4

  tree16 = new Tree_k(4);

  Memory::set_active(true);
  tree16->refine(level_array,n0,n1,max_level,full_nodes=false);
  Memory::print();
  Memory::set_active(false);
  printf ("Bytes / node = %g\n",(float)Memory::bytes(0) / Node_k::num_nodes());

  // Determine image size
  n = cell_size + 2*line_width;
  for (int i=0; i<tree16->levels(); i++) {
    n = 2*n - line_width;
  }

  printf ("\n");
  image = tree16->create_image(n,line_width);
  printf ("INITIAL NON-FULL TREE 16\n");
  printf ("nodes  = %d\n",Node_k::num_nodes());
  printf ("levels = %d\n",tree16->levels());
  write_image(image,n,n,"tree16-3.hdf5");
  delete image;

  // Balance the tree

  printf ("\n");

  Memory::set_active(true);
  tree16->balance(full_nodes);
  Memory::print();
  Memory::set_active(false);
  printf ("Bytes / node = %g\n",(float)Memory::bytes(0) / Node_k::num_nodes());

  printf ("BALANCED NON-FULL TREE 16\n");
  printf ("nodes  = %d\n",Node_k::num_nodes());
  printf ("levels = %d\n",tree16->levels());
  image = tree16->create_image(n,line_width);
  write_image(image,n,n,"tree16-4.hdf5");
  delete image;

  // Optimize the tree

  printf ("\n");

  Memory::set_active(true);
  tree16->optimize();
  Memory::print();
  Memory::set_active(false);
  printf ("Bytes / node = %g\n",(float)Memory::bytes(0) / Node_k::num_nodes());

  printf ("OPTIMIZED NON-FULL TREE 16\n");
  printf ("nodes  = %d\n",Node_k::num_nodes());
  printf ("levels = %d\n",tree16->levels());
  image = tree16->create_image(n,line_width);
  write_image(image,n,n,"tree16-5.hdf5");
  delete image;
  
  delete tree16;



}
