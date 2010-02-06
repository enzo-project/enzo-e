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
 * SYNOPSIS:
 *
 *    
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
#include "cello.h"
#include "amr_node4.hpp"
#include "amr_tree4.hpp"

const bool debug = false;

Tree4::Tree4()
  : levels_(0),
    root_(new Node4())

{
}

// Refine down to array
void Tree4::refine
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
void Tree4::balance(bool is_full)
{
  // Repeatedly balance
  int pass = 0;
  bool tree_changed;
  do {
    tree_changed = false;
    root_->balance_pass(tree_changed,is_full);
    if (debug) printf ("Balance pass %d\n",pass);
    pass++;
  } while (tree_changed);
  printf ("passes = %d\n",pass);
}

// Replace uniformly-refined patch with single node
void Tree4::optimize()
{
  // Repeatedly optimize
  int pass = 0;
  bool tree_changed;
  do {
    tree_changed = false;
    root_->optimize_pass(tree_changed);
    if (debug) printf ("Optimize pass %d\n",pass);
    pass++;
  } while (tree_changed);
  printf ("passes = %d\n",pass);
}

/// Create an hdf5 file of tree, assuming given source bitmap size
/// 
float * Tree4::create_image (int n,int line_width)
{
  float * image = new float [n*n];
  
  root_->fill_image(image,n,n,0,n-1,0,n-1,0,levels_,line_width);
  return image;
}

