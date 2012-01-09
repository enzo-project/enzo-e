// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Tree.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    
///
/// 

#include "mesh.hpp"

//----------------------------------------------------------------------

Tree::Tree(int d, int r) throw ()
  : d_(d), r_(r), root_(new Node), num_nodes_(1)
{
}

//----------------------------------------------------------------------

Tree::~Tree() throw ()
{
  delete root_;
}

//----------------------------------------------------------------------

Tree::Tree(const Tree & Tree) throw ()
/// @param     Tree  Object being copied
{
  INCOMPLETE("Tree::Tree(Tree)");
}

//----------------------------------------------------------------------

Tree & Tree::operator= (const Tree & Tree) throw ()
/// @param     Tree  Source object of the assignment
/// @return    The target assigned object
{
  INCOMPLETE("Tree::operator=");
  return *this;
}

//======================================================================

