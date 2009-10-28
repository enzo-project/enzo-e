#include <stdio.h>
#include "cello.h"
#include "amr_node_k_2.hpp"
#include <assert.h>

Node_K::Node_K(int k, int level_adjust) 
  : k_(k),
    level_adjust_(level_adjust)

{ 
  ++Node_K::num_nodes_;

  neighbor_ = new Node_K * [k_];

  for (int ix=0; ix<k_; ix++) {
    neighbor_[ix] = NULL;
  }

  child_    = new Node_K * [k_*k_];

  for (int i=0; i<k_*k_; i++) {
      child_[i] = NULL;
  }

  parent_ = NULL;
}

@@@

// Delete the node and all descendents

Node_K::~Node_K() 
{ 
  --Node_K::num_nodes_;

  // recursively delete children

  for (int ix=0; ix<k_; ix++) {
    for (int iy=0; iy<k_; iy++) {
      if (child_[ix][iy]) delete child_[ix][iy];
    }
  }

  // update neighbor's neighbors

  if (neighbor_[R]) neighbor_[R]->neighbor_[L] = NULL;
  if (neighbor_[U]) neighbor_[U]->neighbor_[D] = NULL;
  if (neighbor_[L]) neighbor_[L]->neighbor_[R] = NULL;
  if (neighbor_[D]) neighbor_[D]->neighbor_[U] = NULL;

  neighbor_[R] = NULL;
  neighbor_[U] = NULL;
  neighbor_[L] = NULL;
  neighbor_[D] = NULL;

  // Update parent's children

  if (parent_) {
    for (int ix=0; ix<k_; ix++) {
      for (int iy=0; iy<k_; iy++) {
	if (parent_->child_[ix][iy] == this) parent_->child_[ix][iy] = NULL;
      }
    }
  }

  parent_ = NULL;
}

inline Node_K * Node_K::child (int ix, int iy) 
{ 
  return child_[ix][iy]; 
}

inline Node_K * Node_K::neighbor (face_type face) 
{ 
  return neighbor_[face]; 
}

inline Node_K * Node_K::cousin (face_type face, int ix, int iy) 
{ 
  if (neighbor_[face] && neighbor_[face]->child_[ix][iy]) {
    return neighbor_[face]->child_[ix][iy];
  } else {
    return NULL;
  }
}

inline Node_K * Node_K::parent () 
{ 
  return parent_; 
}

// Set two nodes to be neighbors.  Friend function since nodes may be NULL
inline void make_neighbors 
(
 Node_K * node_1, 
 Node_K * node_2, 
 face_type face_1
 )
{
  if (node_1 != NULL) node_1->neighbor_[face_1] = node_2;
  if (node_2 != NULL) node_2->neighbor_[(face_1+2)%k_] = node_1;
}


// Create 4 empty child nodes

int Node_K::refine 
(
 const int * level_array, 
 int ndx,  int ndy,
 int lowx, int upx,  
 int lowy, int upy,
 int level, 
 int max_level,
 bool full_nodes
 )
{
  int depth = 0;
  if ( level < max_level && lowx < upx-1 && lowy < upy-1 ) {

    // determine whether to refine the node

    int dx = (upx - lowx)/k_;
    int dy = (upy - lowy)/k_;
    int ix2[3] = {lowx, lowx+dx, upx};
    int iy2[3] = {lowy, lowy+dy, upy};

    int depth_child[k_][k_] = {{0}};

    if (full_nodes) {

      // Refine if any bits in the level_array are in this node

      bool refine_node = false;
      for (int iy=lowy; iy<upy && !refine_node; iy++) {
	for (int ix=lowx; ix<upx && !refine_node; ix++) {
	  if (level_array[ix + ndx * iy] >= level) refine_node = true;
	}
      }

      // refine the node if needed

      if (refine_node) {

	create_children_();

	update_children_();

	for (int ix=0; ix<k_; ix++) {
	  for (int iy=0; iy<k_; iy++) {
	    depth_child[ix][iy] = child_[ix][iy]->refine 
	      (level_array,ndx,ndy,ix2[ix],ix2[ix+1],iy2[iy],iy2[iy+1],
	       level+2,max_level,full_nodes);
	  }
	}
      }

    } else {

      // Refine each child separately if any bits in the level_array
      // are in in them

      bool refine_child[k_][k_] = {{false}};

      // loop over children

      for (int ix=0; ix<k_; ix++) {
	for (int iy=0; iy<k_; iy++) {

	  // loop over values in the level array

	  for (int ly=iy2[iy]; ly<iy2[iy+1]; ly++) {
	    for (int lx=ix2[ix]; lx<ix2[ix+1]; lx++) {

	      if (level_array[lx + ndx * ly] >= level) {
		refine_child[ix][iy] = true;

	      }
	    }
	  }

	}
      }

      // refine each child if needed

      for (int ix=0; ix<k_; ix++) {
	for (int iy=0; iy<k_; iy++) {
	  if (refine_child[ix][iy]) {
	    create_child_(ix,iy);
	    update_child_(ix,iy);
	    depth_child[ix][iy] = child_[ix][iy]->refine 
	      (level_array,ndx,ndy,ix2[ix],ix2[ix+1],iy2[iy],iy2[iy+1],
	       level+2,max_level,full_nodes);
	  }
	}
      }

    }

    for (int ix=0; ix<k_; ix++) {
      for (int iy=0; iy<k_; iy++) {
	depth = (depth_child[ix][iy] > depth) ? depth_child[ix][iy] : depth;
      }
    }
    depth += 2;

  } // if not at bottom of recursion

  return depth;
}

void Node_K::create_children_()
{
  for (int ix=0; ix<k_; ix++) {
    for (int iy=0; iy<k_; iy++) {
      create_child_(ix,iy);
    }
  }
}

void Node_K::update_children_()
{
  for (int ix=0; ix<k_; ix++) {
    for (int iy=0; iy<k_; iy++) {
      update_child_(ix,iy);
    }
  }
}

void Node_K::create_child_(int ix, int iy)
{
  child_[ix][iy] = new Node_K();
}

void Node_K::update_child_ (int ix, int iy)
{
  if (child(ix,iy)) {

    child(ix,iy)->parent_ = this;

    // Right neighbor

    if (ix < 3) {
      make_neighbors (child (ix,iy),child (ix+1,iy),R);
    } else {
      make_neighbors (child (ix,iy),cousin (R,0,iy),R);
    }

    // Left neighbor

    if (ix > 0) {
      make_neighbors (child (ix,iy),child (ix-1,iy),L);
    } else {
      make_neighbors (child (ix,iy),cousin (L,3,iy),L);
    }

    // Up neighbor

    if (iy < 3) {
      make_neighbors (child (ix,iy),child (ix,iy+1),U);
    } else {
      make_neighbors (child (ix,iy),cousin (U,ix,0),U);
    }

    // Down neighbor

    if (iy > 0) {
      make_neighbors (child (ix,iy),child (ix,iy-1),D);
    } else {
      make_neighbors (child (ix,iy),cousin (D,ix,3),D);
    }
  }
}

// Perform a pass of trying to remove level-jumps 
void Node_K::balance_pass(bool & refined_tree, bool full_nodes)
{
  if (full_nodes) {

    // is a leaf
    if (! any_children()) {

      bool refine_node = false;

      // Check for adjacent nodes two levels finer

      for (int i=0; i<k_; i++) {
	
	refine_node = refine_node ||
	  ( cousin(R,0,i) && ( cousin(R,0,i)->child(0,0) ) ) ||
	  ( cousin(L,3,i) && ( cousin(L,3,i)->child(0,0) ) ) ||
	  ( cousin(U,i,0) && ( cousin(U,i,0)->child(0,0) ) ) ||
	  ( cousin(D,i,3) && ( cousin(D,i,3)->child(0,0) ) );
      }

      if (refine_node) {

	refined_tree = true;

	create_children_();
	update_children_();

	for (int ix=0; ix<k_; ix++) {
	  for (int iy=0; iy<k_; iy++) {
	    if (child(ix,iy)) {
	      child(ix,iy)->balance_pass(refined_tree,full_nodes);
	    }
	  }
	}
      }
    }

  } else { 

    // not full node refinement
    
    if (! all_children ()) {

      bool refine_child[k_][k_] = {{false}};

      for (int ix=0; ix<k_; ix++) {
	for (int iy=0; iy<k_; iy++) {
	  if (! child(ix,iy)) {
	    // right neighbor
	    if (ix < 3 && child(ix+1,iy)) {
	      for (int k=0; k<k_; k++) {
		refine_child[ix][iy] = refine_child[ix][iy] ||
		  child(ix+1,iy)->child(0,k);
	      }
	    } else if (ix == 3 && cousin(R,0,iy)) {
	      for (int k=0; k<k_; k++) {
		refine_child[ix][iy] = refine_child[ix][iy] ||
		  cousin(R,0,iy)->child(0,k);
	      }
	    }
	    // left neighbor
	    if (ix > 0 && child(ix-1,iy)) {
	      for (int k=0; k<k_; k++) {
		refine_child[ix][iy] = refine_child[ix][iy] ||
		  child(ix-1,iy)->child(3,k);
	      }
	    } else if (ix == 0 && cousin(L,3,iy)) {
	      for (int k=0; k<k_; k++) {
		refine_child[ix][iy] = refine_child[ix][iy] ||
		  cousin(L,3,iy)->child(3,k);
	      }
	    }
	    // up neighbor
	    if (iy < 3 && child(ix,iy+1)) {
	      for (int k=0; k<k_; k++) {
		refine_child[ix][iy] = refine_child[ix][iy] ||
		  child(ix,iy+1)->child(k,0);
	      }
	    } else if (iy == 3 && cousin(U,ix,0)) {
	      for (int k=0; k<k_; k++) {
		refine_child[ix][iy] = refine_child[ix][iy] ||
		  cousin(U,ix,0)->child(k,0);
	      }
	    }
	    // down neighbor
	    if (iy > 0 && child(ix,iy-1)) {
	      for (int k=0; k<k_; k++) {
		refine_child[ix][iy] = refine_child[ix][iy] ||
		  child(ix,iy-1)->child(k,3);
	      }
	    } else if (iy == 0 && cousin(D,ix,3)) {
	      for (int k=0; k<k_; k++) {
		refine_child[ix][iy] = refine_child[ix][iy] ||
		  cousin(D,ix,3)->child(k,3);
	      }
	    }
	  } // k
	} // iy
      } // ix

      for (int ix=0; ix<k_; ix++) {
	for (int iy=0; iy<k_; iy++) {
	  if (refine_child[ix][iy]) {
	    create_child_(ix,iy); 
	    update_child_(ix,iy);
	    child(ix,iy)->balance_pass(refined_tree,full_nodes);
	    refined_tree = true;
	  }
	}
      }
    }
  }

  // not a leaf: recurse

  for (int ix=0; ix<k_; ix++) {
    for (int iy=0; iy<k_; iy++) {
      if (child(ix,iy)) {
	child(ix,iy)->balance_pass(refined_tree,full_nodes);
      }
    }
  }

    //  }
}
  // Perform a pass of trying to optimize uniformly-refined nodes
void Node_K::optimize_pass(bool & refined_tree)
{

  bool single_children = true;
  bool same_level = true;

  for (int ix=0; ix<k_; ix++) {
    for (int iy=0; iy<k_; iy++) {
      single_children = single_children && 
	child_[ix][iy] && (! child_[ix][iy]->any_children());
    }
  }

  if (single_children) {

    // assert: all children exist

    int level = child(0,0)->level_adjust_;
    for (int ix=0; ix<k_; ix++) {
      for (int iy=0; iy<k_; iy++) {
	same_level = same_level &&
	  (child(ix,iy)->level_adjust_ == level);
      }
    }
  }

  if (single_children && same_level) {

    // adjust effective resolution

    level_adjust_ += 2 + child_[0][0]->level_adjust_; 

    for (int ix=0; ix<k_; ix++) {
      for (int iy=0; iy<k_; iy++) {
	delete child_[ix][iy];
	child_[ix][iy]=NULL;
      }
    }

    refined_tree = true;

  } else {
    
    // recurse

    for (int ix=0; ix<k_; ix++) {
      for (int iy=0; iy<k_; iy++) {
	if (child_[ix][iy]) {
	  child_[ix][iy]->optimize_pass(refined_tree);
	}
      }
    }
  }

}

// Fill the image region with values
void Node_K::fill_image
(
 float * image,
 int ndx,  int ndy,
 int lowx, int upx,  
 int lowy, int upy,
 int level,
 int num_levels,
 int line_width
 )
{
  int ix,iy,i;

  level += level_adjust_;

  // Fill interior

  for (iy=lowy; iy<=upy; iy++) {
    for (ix=lowx; ix<=upx; ix++) {
      i = ix + ndx * iy;
      image[i] = 2*num_levels - level; 
    }
  }

  // Draw border
  for (ix=lowx; ix<=upx; ix++) {
    iy = lowy;
    for (int k=0; k < line_width; k++) {
      image[ix + ndx*(iy+k)] = 0;
    }
    iy = upy;
    for (int k=0; k < line_width; k++) {
      image[ix + ndx*(iy+k)] = 0;
    }
  }

  for (iy=lowy; iy<=upy; iy++) {
    ix = lowx;
    for (int k=0; k < line_width; k++) {
      image[(ix+k) + ndx*iy] = 0;
    }
    ix = upx;
    for (int k=0; k < line_width; k++) {
      image[(ix+k) + ndx*iy] = 0;
    }
  }
    

  // Recurse


  int dx = (upx - lowx)/k_;
  int dy = (upy - lowy)/k_;
  int ix2[5] = {lowx, lowx+dx, upx};
  int iy2[5] = {lowy, lowy+dy, upy};

  for (int ix=0; ix<k_; ix++) {
    for (int iy=0; iy<k_; iy++) {
      if (child_[ix][iy]) {
	child_[ix][iy]->fill_image 
	  (image,ndx,ndy,ix2[ix],ix2[ix+1], iy2[iy],iy2[iy+1], level + 1, num_levels,line_width);
      }
    }
  }
}

int Node_K::num_nodes_ = 0;

