
#include <stdio.h>
#include "cello.h"
#include "tree4.h"
#include "node4.h"

const bool debug = true;

Tree4::Tree4()
  : root_(new Node4()),
    levels_(0)
{
}

// Refine down to array
void Tree4::refine
(
 const bool * mask_array, 
 int nd0, int nd1, 
 int max_level,
 bool is_full
 )
{
  levels_ = root_->refine(mask_array,nd0,nd1,0,nd0,0,nd1,0,max_level,is_full);
  if (debug) printf ("%d\n",levels_);
}

// Remove level-jumps 
void Tree4::normalize(bool is_full)
{
  // Repeatedly normalize
  int pass = 0;
  bool tree_changed;
  do {
    tree_changed = false;
    root_->normalize_pass(tree_changed,is_full);
    if (debug) printf ("Normalize pass %d\n",pass++,tree_changed);
  } while (tree_changed);
}

// Replace uniformly-refined patch with single node
void Tree4::optimize(bool is_full)
{
  // Repeatedly optimize
  int pass = 0;
  bool tree_changed;
  do {
    tree_changed = false;
    root_->optimize_pass(tree_changed,is_full);
    if (debug) printf ("Optimize pass %d\n",pass++,tree_changed);
  } while (tree_changed);
}

/// Create an hdf5 file of tree, assuming given source bitmap size
/// 
float * Tree4::create_image (int n)
{
  float * image = new float [n*n];
  
  root_->fill_image(image,n,n,0,n-1,0,n-1,0,levels_);
  return image;
}

