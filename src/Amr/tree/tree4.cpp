
#include <stdio.h>
#include "cello.h"
#include "tree4.h"
#include "node4.h"

Tree4::Tree4()
  : root_(new Node4()),
    levels_(0)
{
}

// Refine down to array
void Tree4::refine(const bool * mask_array, int nd0, int nd1, int max_level)
{
  levels_ = root_->refine(mask_array,nd0,nd1,0,nd0,0,nd1,0,max_level);
  printf ("%d\n",levels_);
}

// Remove level-jumps 
void Tree4::normalize()
{
  // Repeatedly normalize
  int pass = 0;
  do {
    printf ("Pass %d\n",pass++);
  } while (root_->normalize_pass());
  root_->normalize_pass();
  root_->normalize_pass();
  root_->normalize_pass();
}

/// Create an hdf5 file of tree, assuming given source bitmap size
/// 
float * Tree4::create_image (int n)
{
  float * image = new float [n*n];
  
  root_->fill_image(image,n,n,0,n-1,0,n-1,0,levels_);
  return image;
}

