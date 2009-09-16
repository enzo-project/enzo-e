#include <stdio.h>
#include "cello.h"
#include "node4.h"
#include <assert.h>

const bool debug = false;

Node4::Node4(int level_adjust) 

  : level_adjust_(level_adjust)

{ 
  ++Node4::num_nodes_;

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

Node4::~Node4() 
{ 
  --Node4::num_nodes_;

  // recursively delete children

  if (child_[UL]) delete child_[UL];
  if (child_[DL]) delete child_[DL];
  if (child_[UR]) delete child_[UR];
  if (child_[DR]) delete child_[DR];

  // update neighbor's neighbor pointers

  if (neighbor_[R]) neighbor_[R]->neighbor_[L] = NULL;
  if (neighbor_[U]) neighbor_[U]->neighbor_[D] = NULL;
  if (neighbor_[L]) neighbor_[L]->neighbor_[R] = NULL;
  if (neighbor_[D]) neighbor_[D]->neighbor_[U] = NULL;

  neighbor_[R] = NULL;
  neighbor_[U] = NULL;
  neighbor_[L] = NULL;
  neighbor_[D] = NULL;

  // Update parent's child

  if (parent_->child_[UL] == this) parent_->child_[UL] = NULL;
  if (parent_->child_[DL] == this) parent_->child_[DL] = NULL;
  if (parent_->child_[UR] == this) parent_->child_[UR] = NULL;
  if (parent_->child_[DR] == this) parent_->child_[DR] = NULL;

  parent_ = NULL;
}

inline Node4 * Node4::child (corner_type corner) 
{ 
  return child_[corner]; 
}

inline Node4 * Node4::neighbor (face_type face) 
{ 
  return neighbor_[face]; 
}

inline Node4 * Node4::cousin (face_type face, corner_type corner) 
{ 
  if (neighbor_[face] && neighbor_[face]->child_[corner]) {
    return neighbor_[face]->child_[corner];
  } else {
    return NULL;
  }
}

inline Node4 * Node4::parent () 
{ 
  return parent_; 
}

// Set two nodes to be neighbors.  Friend function since nodes may be NULL
inline void make_neighbors 
(
 Node4 * node_1, face_type face_1,
 Node4 * node_2, face_type face_2
 )
{
  if (node_1 != NULL) node_1->neighbor_[face_1] = node_2;
  if (node_2 != NULL) node_2->neighbor_[face_2] = node_1;
}


// Create 4 empty child nodes

int Node4::refine 
(
 const bool * mask_array, 
 int nd0,  int nd1,
 int low0, int up0,  
 int low1, int up1,
 int level, 
 int max_level,
 bool full_nodes
 )
{
  int depth = 0;

  if ( level < max_level && low0 < up0-1 && low1 < up1-1 ) {

    // determine whether to refine the node


    int mid0 = (up0 + low0)/2;
    int mid1 = (up1 + low1)/2;

    int depth_child[4] = {0};

    if (full_nodes) {

      // Refine if any bits in the mask_array are in this node

      bool refine_node = false;
      for (int i1=low1; i1<up1 && !refine_node; i1++) {
	for (int i0=low0; i0<up0 && !refine_node; i0++) {
	  if (mask_array[i0 + nd0 * i1]) refine_node = true;
	}
      }

      // refine the node if needed

      if (refine_node) {

	create_children_();

	update_children_();

	depth_child[UL] = child_[UL]->refine 
	  (mask_array,nd0,nd1,low0,mid0,low1,mid1,level+1,max_level,full_nodes);
	depth_child[DL] = child_[DL]->refine 
	  (mask_array,nd0,nd1,mid0,up0,low1,mid1,level+1,max_level,full_nodes);
	depth_child[UR] = child_[UR]->refine 
	  (mask_array,nd0,nd1,low0,mid0,mid1,up1,level+1,max_level,full_nodes);
	depth_child[DR] = child_[DR]->refine 
	  (mask_array,nd0,nd1,mid0,up0,mid1,up1,level+1,max_level,full_nodes);

      }

    } else {

      // Refine each child separately if any bits in the mask_array
      // are in in them

      bool refine_child[4] = {false};

      for (int i1=low1; i1<mid1 && !refine_child[UL]; i1++) {
	for (int i0=low0; i0<mid0 && !refine_child[UL]; i0++) {
	  if (mask_array[i0 + nd0 * i1]) refine_child[UL] = true;
	}
      }
      for (int i1=low1; i1<mid1 && !refine_child[DL]; i1++) {
	for (int i0=mid0; i0<up0 && !refine_child[DL]; i0++) {
	  if (mask_array[i0 + nd0 * i1]) refine_child[DL] = true;
	}
      }
      for (int i1=mid1; i1<up1 && !refine_child[UR]; i1++) {
	for (int i0=low0; i0<mid0 && !refine_child[UR]; i0++) {
	  if (mask_array[i0 + nd0 * i1]) refine_child[UR] = true;
	}
      }
      for (int i1=mid1; i1<up1 && !refine_child[DR]; i1++) {
	for (int i0=mid0; i0<up0 && !refine_child[DR]; i0++) {
	  if (mask_array[i0 + nd0 * i1]) refine_child[DR] = true;
	}
      }

      // refine each child if needed

      if (refine_child[UL]) {
	create_child_(UL);
	update_child_(UL);
	depth_child[UL] = child_[UL]->refine 
	  (mask_array,nd0,nd1,low0,mid0,low1,mid1,level+1,max_level,full_nodes);
      }

      if (refine_child[DL]) {
	create_child_(DL);
	update_child_(DL);
	depth_child[DL] = child_[DL]->refine 
	  (mask_array,nd0,nd1,mid0,up0,low1,mid1,level+1,max_level,full_nodes);
      }

      if (refine_child[UR]) {
	create_child_(UR);
	update_child_(UR);
	depth_child[UR] = child_[UR]->refine 
	  (mask_array,nd0,nd1,low0,mid0,mid1,up1,level+1,max_level,full_nodes);
      }

      if (refine_child[DR]) {
	create_child_(DR);
	update_child_(DR);
	depth_child[DR] = child_[DR]->refine 
	  (mask_array,nd0,nd1,mid0,up0,mid1,up1,level+1,max_level,full_nodes);
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

void Node4::create_children_()
{
  create_child_(UL);
  create_child_(DL);
  create_child_(UR);
  create_child_(DR);
}

void Node4::update_children_()
{
  update_child_ (UL);
  update_child_ (DL);
  update_child_ (UR);
  update_child_ (DR);
}

void Node4::create_child_(corner_type corner)
{
  child_[corner] = new Node4();
}

void Node4::update_child_ (corner_type corner)
{
  if (child(corner)) {

    child(corner)->parent_ = this;

    if (corner == UL) {

      make_neighbors (child (UL),R,child (UR),  L);
      make_neighbors (child (UL),D,child (DL),  U);
      make_neighbors (child (UL),L,cousin(L,UR),R);
      make_neighbors (child (UL),U,cousin(U,DL),D);

    } else if (corner == DL) {

      make_neighbors (child (DL),U,child (UL),  D);
      make_neighbors (child (DL),R,child (DR),  L);
      make_neighbors (child (DL),D,cousin(D,UL),U);
      make_neighbors (child (DL),L,cousin(L,DR),R);

    } else if (corner == UR) {

      make_neighbors (child (UR),D,child (DR),  U);
      make_neighbors (child (UR),L,child (UL),  R);
      make_neighbors (child (UR),U,cousin(U,DR),D);
      make_neighbors (child (UR),R,cousin(R,UL),L);

    } else if (corner == DR) {

      make_neighbors (child (DR),U,child (UR),  D);
      make_neighbors (child (DR),L,child (DL),  R);
      make_neighbors (child (DR),D,cousin(D,UR),U);
      make_neighbors (child (DR),R,cousin(R,DL),L);

    }
  }
}

// Perform a pass of trying to remove level-jumps 
void Node4::normalize_pass(bool & refined_tree, bool full_nodes)
{

  bool refine_UL = false;
  bool refine_DL = false;
  bool refine_UR = false;
  bool refine_DR = false;

  if (is_leaf()) {

    refine_DR = refine_DR ||
      (neighbor_[R] && 
       neighbor_[R]->child_[UL] && 
       (neighbor_[R]->child_[UL]->child_[UL] ||
	neighbor_[R]->child_[UL]->child_[DL]));
    refine_DR = refine_DR ||
      (neighbor_[D] && 
       neighbor_[D]->child_[UL] && 
       (neighbor_[D]->child_[UL]->child_[UL] ||
	neighbor_[D]->child_[UL]->child_[UR]));

    refine_UR = refine_UR ||
      (neighbor_[R] && 
       neighbor_[R]->child_[DL] && 
       (neighbor_[R]->child_[DL]->child_[UL] ||
	neighbor_[R]->child_[DL]->child_[DL]));
    refine_UR = refine_UR ||
      (neighbor_[U] && 
       neighbor_[U]->child_[DL] && 
       (neighbor_[U]->child_[DL]->child_[DL] ||
	neighbor_[U]->child_[DL]->child_[DR]));

    refine_DL = refine_DL ||
      (neighbor_[L] && 
       neighbor_[L]->child_[UR] && 
       (neighbor_[L]->child_[UR]->child_[UR] ||
	neighbor_[L]->child_[UR]->child_[DR]));
    refine_DL = refine_DL || 
      (neighbor_[D] && 
       neighbor_[D]->child_[UR] && 
       (neighbor_[D]->child_[UR]->child_[UL] ||
	neighbor_[D]->child_[UR]->child_[UR]));
       
    refine_UL = refine_UL ||
      (neighbor_[U] && 
       neighbor_[U]->child_[DR] && 
       (neighbor_[U]->child_[DR]->child_[DL] ||
	neighbor_[U]->child_[DR]->child_[DR]));
    refine_UL = refine_UL || 
      (neighbor_[L] && 
       neighbor_[L]->child_[DR] && 
       (neighbor_[L]->child_[DR]->child_[UR] ||
	neighbor_[L]->child_[DR]->child_[DR]));
       

    if (refine_UL || refine_DL || refine_UR || refine_DR) {

      refined_tree = true;

      if (full_nodes) {

	// refine the node completely

	create_children_();

	update_children_();

	child_[UL]->normalize_pass(refined_tree,full_nodes);
	child_[DL]->normalize_pass(refined_tree,full_nodes);
	child_[UR]->normalize_pass(refined_tree,full_nodes);
	child_[DR]->normalize_pass(refined_tree,full_nodes);

      } else {

	// refine only specified children

	if (refine_UL) {
	  create_child_(UL); 
	  update_child_(UL);
	  child(UL)->normalize_pass(refined_tree,full_nodes);
	}
	if (refine_DL) {
	  create_child_(DL); 
	  update_child_(DL);
	  child(DL)->normalize_pass(refined_tree,full_nodes);
	}
	if (refine_UR) {
	  create_child_(UR); 
	  update_child_(UR);
	  child(UR)->normalize_pass(refined_tree,full_nodes);
	}
	if (refine_DR) {
	  create_child_(DR); 
	  update_child_(DR);
	  child(DR)->normalize_pass(refined_tree,full_nodes);
	}
      }
    }

  } else {
    
    if (full_nodes) {
      child_[UL]->normalize_pass(refined_tree,full_nodes);
      child_[DL]->normalize_pass(refined_tree,full_nodes);
      child_[UR]->normalize_pass(refined_tree,full_nodes);
      child_[DR]->normalize_pass(refined_tree,full_nodes);
    } else {
      if (child(UL)) child(UL)->normalize_pass(refined_tree,full_nodes);
      if (child(DL)) child(DL)->normalize_pass(refined_tree,full_nodes);
      if (child(UR)) child(UR)->normalize_pass(refined_tree,full_nodes);
      if (child(DR)) child(DR)->normalize_pass(refined_tree,full_nodes);
    }

  }

}

// Perform a pass of trying to optimize uniformly-refined nodes
void Node4::optimize_pass(bool & refined_tree, bool full_nodes)
{
  int any = 0;
  if (child_[0] && ! child_[0]->child_[any] &&
      child_[1] && ! child_[1]->child_[any] &&
      child_[2] && ! child_[2]->child_[any] &&
      child_[3] && ! child_[3]->child_[any] &&
      child_[0]->level_adjust_ == child_[1]->level_adjust_ &&
      child_[1]->level_adjust_ == child_[2]->level_adjust_ &&
      child_[2]->level_adjust_ == child_[3]->level_adjust_ ) {

    // adjust effective resolution

    level_adjust_ += 1 + child_[any]->level_adjust_; 

    delete child_[0];
    delete child_[1];
    delete child_[2];
    delete child_[3];

    refined_tree = true;

  } else {
    if (child_[0]) child_[0]->optimize_pass(refined_tree,full_nodes);
    if (child_[1]) child_[1]->optimize_pass(refined_tree,full_nodes);
    if (child_[2]) child_[2]->optimize_pass(refined_tree,full_nodes);
    if (child_[3]) child_[3]->optimize_pass(refined_tree,full_nodes);
  }

}

// Fill the image region with values
void Node4::fill_image
(
 float * image,
 int nd0,  int nd1,
 int low0, int high0,  
 int low1, int high1,
 int level,
 int num_levels
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
    image[i0 + nd0*i1] = 0;
    i1 = high1;
    image[i0 + nd0*i1] = 0;
  }

  for (i1=low1; i1<=high1; i1++) {
    i0 = low0;
    image[i0 + nd0*i1] = 0;
    i0 = high0;
    image[i0 + nd0*i1] = 0;
  }
    

  // Recurse

  int mid0 = (high0 + low0)/2;
  int mid1 = (high1 + low1)/2;

  if (child_[UL]) {
    child_[UL]->fill_image 
      (image,nd0,nd1,low0,mid0, low1,mid1, level + 1, num_levels);
  }
  if (child_[DL]) {
    child_[DL]->fill_image 
      (image,nd0,nd1,mid0,high0,low1,mid1, level + 1, num_levels);
  }
  if (child_[UR]) {
    child_[UR]->fill_image 
      (image,nd0,nd1,low0,mid0, mid1,high1,level + 1, num_levels);
  }
  if (child_[DR]) {
    child_[DR]->fill_image 
      (image,nd0,nd1,mid0,high0,mid1,high1,level + 1, num_levels);
  }
}

int Node4::num_nodes_ = 0;

