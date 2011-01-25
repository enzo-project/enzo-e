// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_tree16.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-09-18
/// @brief    Implementation of the Tree16 class

#include "cello.hpp"

#include "mesh_node16.hpp"
#include "mesh_tree16.hpp"

const bool debug = false;

Tree16::Tree16()
  ///   
  : levels_(0),
    root_(new Node16())
{
}

void Tree16::refine
(
 const int * level_array, 
 int nd0, int nd1, 
 int max_level,
 bool is_full
 )
/// @param    level_array Array of levels to refine to
/// @param    nd0       x-dimension of level_array[]
/// @param    nd1       y-dimension of level_array[]
/// @param    max_level Maximum refinement level
/// @param    is_full   Whether nodes always have a full complement of children
{
  levels_ = root_->refine(level_array,nd0,nd1,0,nd0,0,nd1,0,max_level,is_full);
  if (debug) printf ("%d\n",levels_);
}

void Tree16::balance(bool is_full)
/// @param    is_full   Whether nodes always have a full complement of children
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

void Tree16::optimize()
///
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

float * Tree16::create_image (int n,int line_width)
/// @param    n         Size of the image along each dimension
/// @param    line_width Width of lines bounding nodes
{
  float * image = new float [n*n];
  
  root_->fill_image(image,n,n,0,n-1,0,n-1,0,levels_,line_width);
  return image;
}

