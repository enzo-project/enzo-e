
#include <stdio.h>
#include "cello.h"
#include "node16.h"
#include "tree16.h"

const bool debug = false;

Tree16::Tree16()
  : levels_(0),
    root_(new Node16())
    
{
}

// Refine down to array
void Tree16::refine
(
 const int * level_array, 
 int nd0, int nd1, 
 int max_level,
 bool is_full
 )
{
  levels_ = root_->refine(level_array,nd0,nd1,0,nd0,0,nd1,0,max_level,is_full);
  if (debug) printf ("%d\n",levels_);
}

// Remove level-jumps 
void Tree16::normalize(bool is_full)
{
  // Repeatedly normalize
  int pass = 0;
  bool tree_changed;
  do {
    tree_changed = false;
    root_->normalize_pass(tree_changed,is_full);
    if (debug) printf ("Normalize pass %d\n",pass);
    pass++;
  } while (tree_changed);
  printf ("passes = %d\n",pass);
}

// Replace uniformly-refined patch with single node
void Tree16::optimize(bool is_full)
{
  // Repeatedly optimize
  int pass = 0;
  bool tree_changed;
  do {
    tree_changed = false;
    root_->optimize_pass(tree_changed,is_full);
    if (debug) printf ("Optimize pass %d\n",pass);
    pass++;
  } while (tree_changed);
  printf ("passes = %d\n",pass);
}

/// Create an hdf5 file of tree, assuming given source bitmap size
/// 
float * Tree16::create_image (int n,int line_width)
{
  printf ("n = %d\n",n);
  float * image = new float [n*n];
  
  root_->fill_image(image,n,n,0,n-1,0,n-1,0,levels_,line_width);
  return image;
}

