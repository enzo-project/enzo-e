//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      amr_tree2k.cpp
 * @brief     
 * @author    James Bordner (jobordner@ucsd.edu)
 * @date      
 * @bug       
 * @note      
 *
 * DESCRIPTION 
 * 
 *    
 *
 * PACKAGES
 *
 *    
 * 
 * INCLUDES
 *  
 *    
 *
 * PUBLIC FUNCTIONS
 *  
 *    
 *
 * PRIVATE FUCTIONS
 *  
 *    
 *
 * $Id$
 *
 *********************************************************************
 */
 

#include <stdio.h>
#include "cello.h"
#include "amr_node2k.hpp"
#include "amr_tree2k.hpp"

const bool debug = false;

Tree2K::Tree2K(int r)
  : r_(r),
    levels_(0),
    root_(new Node2K(r))
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Create a Tree2K object
 *
 *********************************************************************
 */
   
{
}

void Tree2K::refine
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Refine down to array
 *
 *********************************************************************
 */
(
 const int * level_array, 
 int ndx, int ndy,
 int max_level,
 bool is_full
 )
{
  levels_ = root_->refine
    (level_array,ndx,ndy,0,ndx,0,ndy,0,max_level,is_full);
  if (debug) printf ("%d\n",levels_);
}

void Tree2K::balance(bool is_full)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Remove level-jumps 
 *
 *********************************************************************
 */
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

void Tree2K::optimize()
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Replace uniformly-refined patch with single node
 *
 *********************************************************************
 */
{
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

float * Tree2K::create_image (int n,int line_width)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Create an hdf5 file of tree, assuming given source bitmap size
 *
 *********************************************************************
 */
{
  printf ("n = %d\n",n);
  float * image = new float [n*n];
  
  root_->fill_image(image,n,n,0,n-1,0,n-1,0,levels_,line_width);
  return image;
}

