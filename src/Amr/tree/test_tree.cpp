#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>

#include <string>

#include "cello.h"

#include "array.hpp"
#include "disk.hpp"
#include "tree4.h"
#include "node4.h"

const bool debug = false;
const int  cell_size = 1;
const int  gray_threshold = 127;
const bool full_nodes = false;
const int  max_level = 10;

//#include "dot_UR.h"  // point near center in upper right
//#include "dot_UL.h"  // point near center in upper left
//#include "dot_DL.h"  // point near center in lower left
//#include "dot_DR.h"  // point near center in lower right

//#include "egret-bg.h"  // point near center in lower right
//#include "peace3.h"
//#include "sdsc.h"
//#include "sdsc-logo.h"
//#include "sdsc-logo-512.h"
// #include "sdsc-logo-2048.h"

//#include "norman.h"

#include "image.h"

// read in the gimp-generated image data into a bitmap array

bool * create_bitmap(int * n0, int * n1)
{
  int size = (width > height) ? width: height;

  bool * bitmap_array = new bool [size*size];
  for (int i=0; i<size*size; i++) bitmap_array[i] = false;

  *n0 = size;
  *n1 = size;

  int pixel[3];
  char * data = header_data;

  int ix0 = (size - width) / 2;
  int iy0 = (size - height) / 2;

  for (int iy=0; iy<height; iy++) {
    for (int ix=0; ix<width; ix++) {
      HEADER_PIXEL(data,pixel);
      if (debug) printf ("%d %d  %d %d %d\n",ix,iy,pixel[0],pixel[1],pixel[2]);
      int i = (iy+iy0) + size*(ix+ix0);
      bitmap_array[i] = ((pixel[0] + pixel[1] + pixel[2]) > gray_threshold*3);
    }
  }
  return bitmap_array;
}


main()
{
  int n0,n1;

  // read in the gimp image into bitmap

  bool * bitmap = create_bitmap(&n0,&n1);

  // Create a tree and refine to the image

  Tree4 tree;

  tree.refine(bitmap,n0,n1,max_level,full_nodes);

  printf ("nodes  = %d\n",Node4::num_nodes());
  printf ("levels = %d\n",tree.levels());
  printf ("sizeof (Tree4) = %d\n",int(sizeof(Tree4)));

  // Generate an array containing an "image" of the tree

  // Determine image size
  int n = cell_size + 2;
  for (int i=0; i<tree.levels(); i++) {
    n = 2*n - 1;
  }
  float * image = tree.create_image(n);
  printf ("Image size: %d x %d\n",n,n);

  // Write the image to disk

  Hdf5 hdf5;
  hdf5.file_open("tree.hdf5","w");
  ArraySerial tree_array (image,n,n,1);
  hdf5.dataset_open ("tree_image",tree_array);
  hdf5.write(tree_array);
  hdf5.dataset_close ();
  hdf5.file_close();

  delete image;

  // Normalize the tree

  tree.normalize(full_nodes);

  printf ("nodes  = %d\n",Node4::num_nodes());

  // Generate an array containing an "image" of the tree

  // Determine image size
  n = cell_size + 2;
  for (int i=0; i<tree.levels(); i++) {
    n = 2*n - 1;
  }

  // Generate the image
  image = tree.create_image(n);

  // Write the image to disk

  hdf5.file_open("tree-normal.hdf5","w");
  ArraySerial tree_array_normal (image,n,n,1);
  hdf5.dataset_open ("tree_image",tree_array_normal);
  hdf5.write(tree_array_normal);
  hdf5.dataset_close ();
  hdf5.file_close();

  delete image;

  // Optimize the tree

  tree.optimize(full_nodes);

  printf ("nodes  = %d\n",Node4::num_nodes());

  // Generate an array containing an "image" of the tree

  // Determine image size
  n = cell_size + 2;
  for (int i=0; i<tree.levels(); i++) {
    n = 2*n - 1;
  }

  // Generate the image
  image = tree.create_image(n);

  // Write the image to disk

  hdf5.file_open("tree-optimal.hdf5","w");
  ArraySerial tree_array_optimal (image,n,n,1);
  hdf5.dataset_open ("tree_image",tree_array_optimal);
  hdf5.write(tree_array_optimal);
  hdf5.dataset_close ();
  hdf5.file_close();

  delete image;
  
}
