// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_node4.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     
/// @brief    Node class for 2^2-trees Tree4

#include <stdio.h>

#include "cello.h"

#include "error.hpp"

#include "mesh_node4.hpp"

//----------------------------------------------------------------------

const bool debug = false;

//----------------------------------------------------------------------

Node4::Node4(int level_adjust) 
    /// @param    level_adjust difference: actual mesh level - tree level
  : level_adjust_(level_adjust)
{ 
  child_[0] = NULL;
  child_[1] = NULL;
  child_[2] = NULL;
  child_[3] = NULL;

  neighbor_[0] = NULL;
  neighbor_[1] = NULL;
  neighbor_[2] = NULL;
  neighbor_[3] = NULL;

  parent_ = NULL;
}

// Delete the node and all descendents

//----------------------------------------------------------------------

Node4::~Node4()
///
{ 
  // recursively delete children

  if (child_[UL]) delete child_[UL];
  if (child_[DL]) delete child_[DL];
  if (child_[UR]) delete child_[UR];
  if (child_[DR]) delete child_[DR];

  // update neighbor's neighbors

  if (neighbor_[R]) neighbor_[R]->neighbor_[L] = NULL;
  if (neighbor_[U]) neighbor_[U]->neighbor_[D] = NULL;
  if (neighbor_[L]) neighbor_[L]->neighbor_[R] = NULL;
  if (neighbor_[D]) neighbor_[D]->neighbor_[U] = NULL;

  neighbor_[R] = NULL;
  neighbor_[U] = NULL;
  neighbor_[L] = NULL;
  neighbor_[D] = NULL;

  // Update parent's child

  if (parent_) {
    if (parent_->child_[UL] == this) parent_->child_[UL] = NULL;
    if (parent_->child_[DL] == this) parent_->child_[DL] = NULL;
    if (parent_->child_[UR] == this) parent_->child_[UR] = NULL;
    if (parent_->child_[DR] == this) parent_->child_[DR] = NULL;
  }

  parent_ = NULL;
}

//----------------------------------------------------------------------

Node4::Node4(const Node4 & node4) throw()
{
  INCOMPLETE_MESSAGE("Node4::Node4","");
}

//----------------------------------------------------------------------

Node4 & Node4::operator= (const Node4 & node4) throw()
{
  INCOMPLETE_MESSAGE("Node4::operator =","");
  return *this;
}


//----------------------------------------------------------------------

int Node4::num_nodes()
{
  int node_count = 0;
  if (child(UL)) node_count += child(UL)->num_nodes();
  if (child(DL)) node_count += child(DL)->num_nodes();
  if (child(UR)) node_count += child(UR)->num_nodes();
  if (child(DR)) node_count += child(DR)->num_nodes();
  return (1 + node_count);
}

//----------------------------------------------------------------------

inline Node4 * Node4::child (corner_type corner) 
/// @param    corner     Corner [UL|DL|UR|DR] of the child node to return
{ 
  return child_[corner]; 
}

//----------------------------------------------------------------------

inline Node4 * Node4::neighbor (face_type face) 
/// @param    face      Face 0 <= (face = [XYZ][MP]) < 6
{ 
  return neighbor_[face]; 
}

//----------------------------------------------------------------------

inline Node4 * Node4::cousin (face_type face, corner_type corner)
/// @param    face      Face 0 <= (face = [XYZ][MP]) < 6
/// @param    corner     Corner [UL|DL|UR|DR] of the child node to return
{ 
  if (neighbor_[face] && neighbor_[face]->child_[corner]) {
    return neighbor_[face]->child_[corner];
  } else {
    return NULL;
  }
}

//----------------------------------------------------------------------

inline Node4 * Node4::parent ()
///
{ 
  return parent_; 
}

//----------------------------------------------------------------------

// Set two nodes to be neighbors.  Friend function since nodes may be NULL
inline void make_neighbors 
(
 Node4 * node_1,
 Node4 * node_2,
 face_type face_1
 )
/// @param    node_1    First neighbor node pointer 
/// @param    node_2    Second neighbor node pointer
/// @param    face_1    Face 0 <= face_1 < 6 of node_1 that is adjacent to node_2
{
  if (node_1 != NULL) node_1->neighbor_[face_1] = node_2;
  if (node_2 != NULL) node_2->neighbor_[(face_1+2)%4] = node_1;
}

//----------------------------------------------------------------------
int Node4::refine 
(
 const int * level_array, 
 int nd0,  int nd1,
 int low0, int up0,  
 int low1, int up1,
 int level, 
 int max_level,
 bool full_nodes
 )
/// @param    level_array Array of levels to refine to
/// @param    nd0       x-dimension of level_array[]
/// @param    nd1       y-dimension of level_array[]
/// @param    low0      Lowest x-index of array for this node
/// @param    up0       Upper bound on x-index of array for this node
/// @param    low1      Lowest y-index of array for this node
/// @param    up1       Upper bound on y-index of array for this node
/// @param    level     Level of this node
/// @param    max_level Maximum refinement level
/// @param    full_nodes Whether nodes always have a full complement of children
{
  int depth = 0;

  if ( level < max_level && low0 < up0-1 && low1 < up1-1 ) {

    // determine whether to refine the node


    int mid0 = (up0 + low0)/2;
    int mid1 = (up1 + low1)/2;

    int depth_child[4] = {0};

    if (full_nodes) {

      // Refine if any bits in the level_array are in this node

      bool refine_node = false;
      for (int i1=low1; i1<up1 && !refine_node; i1++) {
	for (int i0=low0; i0<up0 && !refine_node; i0++) {
	  if (level_array[i0 + nd0 * i1] >= level) refine_node = true;
	}
      }

      // refine the node if needed

      if (refine_node) {

	create_children_();

	update_children_();

	depth_child[UL] = child_[UL]->refine 
	  (level_array,nd0,nd1,low0,mid0,low1,mid1,level+1,max_level,full_nodes);
	depth_child[DL] = child_[DL]->refine 
	  (level_array,nd0,nd1,mid0,up0,low1,mid1,level+1,max_level,full_nodes);
	depth_child[UR] = child_[UR]->refine 
	  (level_array,nd0,nd1,low0,mid0,mid1,up1,level+1,max_level,full_nodes);
	depth_child[DR] = child_[DR]->refine 
	  (level_array,nd0,nd1,mid0,up0,mid1,up1,level+1,max_level,full_nodes);

      }

    } else {

      // Refine each child separately if any bits in the level_array
      // are in in them

      bool refine_child[4] = {false};

      for (int i1=low1; i1<mid1; i1++) {
	for (int i0=low0; i0<mid0; i0++) {
	  if (level_array[i0 + nd0 * i1] >= level) refine_child[UL] = true;
	}
      }
      for (int i1=low1; i1<mid1; i1++) {
	for (int i0=mid0; i0<up0; i0++) {
	  if (level_array[i0 + nd0 * i1] >= level) refine_child[DL] = true;
	}
      }
      for (int i1=mid1; i1<up1; i1++) {
	for (int i0=low0; i0<mid0; i0++) {
	  if (level_array[i0 + nd0 * i1] >= level) refine_child[UR] = true;
	}
      }
      for (int i1=mid1; i1<up1; i1++) {
	for (int i0=mid0; i0<up0; i0++) {
	  if (level_array[i0 + nd0 * i1] >= level) refine_child[DR] = true;
	}
      }

      // refine each child if needed

      if (refine_child[UL]) {
	create_child_(UL);
	update_child_(UL);
	depth_child[UL] = child_[UL]->refine 
	  (level_array,nd0,nd1,low0,mid0,low1,mid1,level+1,max_level,full_nodes);
      }

      if (refine_child[DL]) {
	create_child_(DL);
	update_child_(DL);
	depth_child[DL] = child_[DL]->refine 
	  (level_array,nd0,nd1,mid0,up0,low1,mid1,level+1,max_level,full_nodes);
      }

      if (refine_child[UR]) {
	create_child_(UR);
	update_child_(UR);
	depth_child[UR] = child_[UR]->refine 
	  (level_array,nd0,nd1,low0,mid0,mid1,up1,level+1,max_level,full_nodes);
      }

      if (refine_child[DR]) {
	create_child_(DR);
	update_child_(DR);
	depth_child[DR] = child_[DR]->refine 
	  (level_array,nd0,nd1,mid0,up0,mid1,up1,level+1,max_level,full_nodes);
      }

    }

    depth = (depth_child[UL] > depth) ? depth_child[UL] : depth;
    depth = (depth_child[DL] > depth) ? depth_child[DL] : depth;
    depth = (depth_child[UR] > depth) ? depth_child[UR] : depth;
    depth = (depth_child[DR] > depth) ? depth_child[DR] : depth;

    ++depth;

  } // if not at bottom of recursion

  return depth;
}

//----------------------------------------------------------------------

void Node4::create_children_()
///
{
  create_child_(UL);
  create_child_(DL);
  create_child_(UR);
  create_child_(DR);
}

//----------------------------------------------------------------------

void Node4::update_children_()
///
{
  update_child_ (UL);
  update_child_ (DL);
  update_child_ (UR);
  update_child_ (DR);
}

//----------------------------------------------------------------------

void Node4::create_child_(corner_type corner)
/// @param    corner     Corner [UL|DL|UR|DR] of the child node to return
{
  child_[corner] = new Node4();
}

//----------------------------------------------------------------------

void Node4::update_child_ (corner_type corner)
/// @param    corner     Corner [UL|DL|UR|DR] of the child node to return
{
  if (child(corner)) {

    child(corner)->parent_ = this;

    if (corner == UL) {

      make_neighbors (child (UL),child (UR),R);
      make_neighbors (child (UL),child (DL),D);
      make_neighbors (child (UL),cousin(L,UR),L);
      make_neighbors (child (UL),cousin(U,DL),U);
      
    } else if (corner == DL) {

      make_neighbors (child (DL),child (UL),U);
      make_neighbors (child (DL),child (DR),R);
      make_neighbors (child (DL),cousin(D,UL),D);
      make_neighbors (child (DL),cousin(L,DR),L);
      
    } else if (corner == UR) {

      make_neighbors (child (UR),child (DR),D);
      make_neighbors (child (UR),child (UL),L);
      make_neighbors (child (UR),cousin(U,DR),U);
      make_neighbors (child (UR),cousin(R,UL),R);
      
    } else if (corner == DR) {

      make_neighbors (child (DR),child (UR),U);
      make_neighbors (child (DR),child (DL),L);
      make_neighbors (child (DR),cousin(D,UR),D);
      make_neighbors (child (DR),cousin(R,DL),R);
      
    }
  }
}

// Perform a pass of trying to remove level-jumps 
//----------------------------------------------------------------------

void Node4::balance_pass(bool & refined_tree, bool full_nodes)
{

  if (full_nodes) {

    // is a leaf
    if (! any_children()) {

      bool refine_node = false;

      // Check for adjacent nodes two levels finer

      corner_type any = UL;

      refine_node =
	( cousin(R,UL) && ( cousin(R,UL)->child(any) ) ) ||
	( cousin(D,UL) && ( cousin(D,UL)->child(any) ) ) ||
	( cousin(R,DL) && ( cousin(R,DL)->child(any) ) ) ||
	( cousin(U,DL) && ( cousin(U,DL)->child(any) ) ) ||
	( cousin(L,UR) && ( cousin(L,UR)->child(any) ) ) ||
	( cousin(D,UR) && ( cousin(D,UR)->child(any) ) ) ||
	( cousin(U,DR) && ( cousin(U,DR)->child(any) ) ) ||
	( cousin(L,DR) && ( cousin(L,DR)->child(any) ) );

      if (refine_node) {

	refined_tree = true;

	create_children_();
	update_children_();

	child(UL)->balance_pass(refined_tree,full_nodes);
	child(DL)->balance_pass(refined_tree,full_nodes);
	child(UR)->balance_pass(refined_tree,full_nodes);
	child(DR)->balance_pass(refined_tree,full_nodes);

      }

    } else {
      child(UL)->balance_pass(refined_tree,full_nodes);
      child(DL)->balance_pass(refined_tree,full_nodes);
      child(UR)->balance_pass(refined_tree,full_nodes);
      child(DR)->balance_pass(refined_tree,full_nodes);
    }

  } else { 

    // not full node refinement


    bool refine_child[4] = {false};

    if (! child(UL)) {

      refine_child[UL] = 
	(cousin(U,DL) && 
	 ( cousin(U,DL)->child(DL) ||
	   cousin(U,DL)->child(DR))) ||
	(cousin(L,UR) && 
	 ( cousin(L,UR)->child(UR) ||
	   cousin(L,UR)->child(DR))) ||
	(child(UR) && 
	 ( child(UR)->child(UL) ||
	   child(UR)->child(DL))) ||
	(child(DL) && 
	 ( child(DL)->child(UL) ||
	   child(DL)->child(UR)));

      if (refine_child[UL]) {
	create_child_(UL); 
	update_child_(UL);
	child(UL)->balance_pass(refined_tree,full_nodes);
	refined_tree = true;
      }
    } else {
      child(UL)->balance_pass(refined_tree,full_nodes);
    }
	 
    if ( ! child(DL) ) {

      refine_child[DL] =
	(cousin(L,DR) && 
	 ( cousin(L,DR)->child(UR) ||
	   cousin(L,DR)->child(DR))) ||
	(cousin(D,UL) && 
	 ( cousin(D,UL)->child(UL) ||
	   cousin(D,UL)->child(UR))) ||
	(child(DR) && 
	 ( child(DR)->child(UL) ||
	   child(DR)->child(DL))) ||
	(child(UL) && 
	 ( child(UL)->child(DL) ||
	   child(UL)->child(DR)));

      if (refine_child[DL]) {
	create_child_(DL); 
	update_child_(DL);
	child(DL)->balance_pass(refined_tree,full_nodes);
	refined_tree = true;
      }
    } else {
      child(DL)->balance_pass(refined_tree,full_nodes);
    }

    if ( ! child(UR) ) {

      refine_child[UR] = 
	(cousin(R,UL) && 
	 ( cousin(R,UL)->child(UL) ||
	   cousin(R,UL)->child(DL))) ||
	(cousin(U,DR) && 
	 ( cousin(U,DR)->child(DL) ||
	   cousin(U,DR)->child(DR))) ||
	(child(UL) && 
	 ( child(UL)->child(UR) ||
	   child(UL)->child(DR))) ||
	(child(DR) && 
	 ( child(DR)->child(UL) ||
	   child(DR)->child(UR)));

      if (refine_child[UR]) {
	create_child_(UR); 
	update_child_(UR);
	child(UR)->balance_pass(refined_tree,full_nodes);
	refined_tree = true;
      }
    } else {
      child(UR)->balance_pass(refined_tree,full_nodes);
    }
       
    if ( ! child(DR) ) {
      refine_child[DR] = 
	(cousin(R,DL) && 
	 ( cousin(R,DL)->child(UL) ||
	   cousin(R,DL)->child(DL))) ||
	(cousin(D,UR) && 
	 ( cousin(D,UR)->child(UL) ||
	   cousin(D,UR)->child(UR))) ||
	(child(DL) && 
	 ( child(DL)->child(UR) ||
	   child(DL)->child(DR))) ||
	(child(UR) && 
	 ( child(UR)->child(DL) ||
	   child(UR)->child(DR)));

      if (refine_child[DR]) {
	create_child_(DR); 
	update_child_(DR);
	child(DR)->balance_pass(refined_tree,full_nodes);
	refined_tree = true;
      }
    } else {
      child(DR)->balance_pass(refined_tree,full_nodes);
    }
  }
}

// Perform a pass of trying to optimize uniformly-refined nodes
//----------------------------------------------------------------------

void Node4::optimize_pass(bool & refined_tree)
{

  if (child_[0] && ! child_[0]->any_children() &&
      child_[1] && ! child_[1]->any_children() &&
      child_[2] && ! child_[2]->any_children() &&
      child_[3] && ! child_[3]->any_children() &&
      child_[0]->level_adjust_ == child_[1]->level_adjust_ &&
      child_[1]->level_adjust_ == child_[2]->level_adjust_ &&
      child_[2]->level_adjust_ == child_[3]->level_adjust_ ) {

    // adjust effective resolution

    level_adjust_ += 1 + child_[0]->level_adjust_; 

    delete child_[0];
    delete child_[1];
    delete child_[2];
    delete child_[3];

    refined_tree = true;

  } else {
    if (child_[0]) child_[0]->optimize_pass(refined_tree);
    if (child_[1]) child_[1]->optimize_pass(refined_tree);
    if (child_[2]) child_[2]->optimize_pass(refined_tree);
    if (child_[3]) child_[3]->optimize_pass(refined_tree);
  }

}

// Fill the image region with values
//----------------------------------------------------------------------

void Node4::fill_image
(
 float * image,
 int nd0,  int nd1,
 int low0, int high0,  
 int low1, int high1,
 int level,
 int num_levels,
 int line_width
 )
{
  int i0,i1,i;

  level += level_adjust_;

  // Fill interior

  for (i1=low1; i1<=high1; i1++) {
    for (i0=low0; i0<=high0; i0++) {
      i = i0 + nd0 * i1;
      image[i] = 2*num_levels - level;
    }
  }

  // Draw border
  for (i0=low0; i0<=high0; i0++) {
    i1 = low1;
    for (int k=0; k < line_width; k++) {
      image[i0 + nd0*(i1+k)] = 0;
    }
    i1 = high1;
    for (int k=0; k < line_width; k++) {
      image[i0 + nd0*(i1+k)] = 0;
    }
  }

  for (i1=low1; i1<=high1; i1++) {
    i0 = low0;
    for (int k=0; k < line_width; k++) {
      image[(i0+k) + nd0*i1] = 0;
    }
    i0 = high0;
    for (int k=0; k < line_width; k++) {
      image[(i0+k) + nd0*i1] = 0;
    }
  }
    

  // Recurse

  int mid0 = (high0 + low0)/2;
  int mid1 = (high1 + low1)/2;

  if (child_[UL]) {
    child_[UL]->fill_image 
      (image,nd0,nd1,low0,mid0, low1,mid1, level + 1, num_levels,line_width);
  }
  if (child_[DL]) {
    child_[DL]->fill_image 
      (image,nd0,nd1,mid0,high0,low1,mid1, level + 1, num_levels,line_width);
  }
  if (child_[UR]) {
    child_[UR]->fill_image 
      (image,nd0,nd1,low0,mid0, mid1,high1,level + 1, num_levels,line_width);
  }
  if (child_[DR]) {
    child_[DR]->fill_image 
      (image,nd0,nd1,mid0,high0,mid1,high1,level + 1, num_levels,line_width);
  }
}

