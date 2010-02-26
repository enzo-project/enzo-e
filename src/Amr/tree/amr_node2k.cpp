/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @note      
 *
 *--------------------------------------------------------------------
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
#include <assert.h>

#include "cello.h"

#include "amr_node2k.hpp"

//----------------------------------------------------------------------

Node2K::Node2K(int k, int level_adjust) 
  : k_(k),
    child_(0),
    neighbor_(0),
    parent_(0),
    level_adjust_(level_adjust)

{ 
  ++Node2K::num_nodes_;

  allocate_neighbors_ ();
  
//   allocate_children_ ();

  parent_ = NULL;
}

//----------------------------------------------------------------------

// Delete the node and all descendents

Node2K::~Node2K() 
{ 
  --Node2K::num_nodes_;

  deallocate_children_();

  // update neighbor's neighbors

  for (int i=0; i<num_faces_(); i++) {
    int io = opposite_face_(face_type(i));
    if (neighbor_[i]) neighbor_[i]->neighbor_[io] = NULL;
    neighbor_[i] = NULL;
  }

  deallocate_neighbors_();

  // Update parent's children

  if (parent_) {
    for (int i=0; i<num_children_(); i++) {
      if (parent_->child_[i] == this) parent_->child_[i] = NULL;
    }
  }

  parent_ = NULL;
}

//----------------------------------------------------------------------

inline Node2K * Node2K::child (int ix, int iy) 
{
  if (child_==NULL) return NULL;
  return child_[index_(ix,iy)];
}

//----------------------------------------------------------------------

inline Node2K * Node2K::neighbor (face_type face) 
{ 
  return neighbor_[face]; 
}

//----------------------------------------------------------------------

inline Node2K * Node2K::cousin (face_type face, int ix, int iy) 
{ 
  if (neighbor_[face] && neighbor_[face]->child(ix,iy)) {
    return neighbor_[face]->child(ix,iy);
  } else {
    return NULL;
  }
}

//----------------------------------------------------------------------

inline Node2K * Node2K::parent () 
{ 
  return parent_; 
}

// Set two nodes to be neighbors.  Friend function since nodes may be NULL
inline void make_neighbors 
(
 Node2K * node_1, 
 Node2K * node_2, 
 face_type face_1
 )
{
  if (node_1 != NULL) node_1->neighbor_[face_1] = node_2;
  if (node_2 != NULL) {
    face_type face_2 = face_type(node_2->opposite_face_(face_1));
    node_2->neighbor_[face_2] = node_1;
  }
}


//----------------------------------------------------------------------

// Create empty child nodes

int Node2K::refine 
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
  int increment = level_increment_();

  if ( level < max_level && 
       lowx < upx-1 && 
       lowy < upy-1 ) {

    // determine whether to refine the node

    int dx = (upx - lowx)/k_;
    int dy = (upy - lowy)/k_;

    int * ixk = new int [k_+1];
    int * iyk = new int [k_+1];

    for (int i=0; i<=k_; i++) {
      ixk[i] = lowx + i*dx;
      iyk[i] = lowy + i*dy;
    }

    int * depth_child = new int [num_children_()];

    for (int i=0; i<num_children_(); i++) {
      depth_child[i] = 0;
    }

    if (full_nodes) {

      // Refine if any bits in the level_array are in this node

      bool refine_node = false;

      for (int iy=lowy; iy<upy && !refine_node; iy++) {
	for (int ix=lowx; ix<upx && !refine_node; ix++) {
	  if (level_array[ix + ndx*iy] >= level) {
	    refine_node = true;
	  }
	}
      }

      // refine the node if needed

      if (refine_node) {

	create_children_();

	update_children_();

	for (int iy=0; iy<k_; iy++) {
	  for (int ix=0; ix<k_; ix++) {
	    depth_child[index_(ix,iy)] = child(ix,iy)->refine 
	      (level_array,
	       ndx,ndy,
	       ixk[ix],ixk[ix+1],
	       iyk[iy],iyk[iy+1],
	       level + increment,
	       max_level,full_nodes);
	  }
	}
      }
      
    } else {

      // Refine each child separately if any bits in the level_array
      // are in in them

      bool * refine_child = new bool [num_children_()];

      for (int i=0; i<num_children_(); i++) {
	refine_child[i] = false;
      }
      

      // loop over children

      for (int iy=0; iy<k_; iy++) {
	for (int ix=0; ix<k_; ix++) {

	  int i = index_(ix,iy);

	  // loop over values in the level array

	  for (int ly=iyk[iy]; ly<iyk[iy+1] && ! refine_child[i]; ly++) {
	    for (int lx=ixk[ix]; lx<ixk[ix+1] && ! refine_child[i]; lx++) {

	      int l = lx + ndx * ly;
	      if (level_array[l] >= level) {
		refine_child[i] = true;
	      }

	    }
	  }

	}
      }

      // refine each child if needed

      for (int iy=0; iy<k_; iy++) {
	for (int ix=0; ix<k_; ix++) {

	  int i = index_(ix,iy);

	  if (refine_child[i]) {
	    create_child_(ix,iy);
	    update_child_(ix,iy);
	    depth_child[i] = child(ix,iy)->refine 
	      (level_array,
	       ndx,ndy,
	       ixk[ix],ixk[ix+1],
	       iyk[iy],iyk[iy+1],
	       level + increment,
	       max_level,full_nodes);
	  }

	}
      }

      delete [] refine_child;
    }

    // determine depth as depth of deepest child + level increment

    for (int iy=0; iy<k_; iy++) {
      for (int ix=0; ix<k_; ix++) {
	int i = index_(ix,iy);
	depth = (depth_child[i] > depth) ? depth_child[i] : depth;
      }
    }

    delete [] depth_child;

    depth += increment;

  } // if not at bottom of recursion

  return depth;
}

//----------------------------------------------------------------------

void Node2K::create_children_()
{
  for (int iy=0; iy<k_; iy++) {
    for (int ix=0; ix<k_; ix++) {
      create_child_(ix,iy);
    }
  }
}

//----------------------------------------------------------------------

void Node2K::update_children_()
{
  for (int iy=0; iy<k_; iy++) {
    for (int ix=0; ix<k_; ix++) {
      update_child_(ix,iy);
    }
  }
}

//----------------------------------------------------------------------

void Node2K::create_child_(int ix, int iy)
{
  if (child_ == NULL) allocate_children_();
  child_[index_(ix,iy)] = new Node2K(k_);
}

//----------------------------------------------------------------------

void Node2K::update_child_ (int ix, int iy)
{
  if (child(ix,iy)) {

    child(ix,iy)->parent_ = this;

    int nx = k_;
    int ny = k_;

    // XM-face neighbors

    if (ix > 0) {
      make_neighbors (child (ix,iy), child (ix-1,iy),XM);
    } else {
      make_neighbors (child (ix,iy), cousin (XM,nx-1,iy),XM);
    }

    // XP-face neighbors

    if (ix < nx-1) {
      make_neighbors (child (ix,iy), child (ix+1,iy), XP);
    } else {
      make_neighbors (child (ix,iy), cousin (XP,0,iy), XP);
    }

    // YM-face neighbor

    if (iy > 0) {
      make_neighbors (child (ix,iy), child (ix,iy-1),YM);
    } else {
      make_neighbors (child (ix,iy), cousin (YM,ix,ny-1),YM);
    }

    // YP-face neighbor

    if (iy < ny-1) {
      make_neighbors (child (ix,iy), child (ix,iy+1),YP);
    } else {
      make_neighbors (child (ix,iy), cousin (YP,ix,0),YP);
    }

  }
}

//----------------------------------------------------------------------

// Perform a pass of trying to remove level-jumps 

void Node2K::balance_pass(bool & refined_tree, bool full_nodes)
{
  int nx = k_;
  int ny = k_;

  if (full_nodes) {

    // if is a leaf

    if (! any_children()) {

      // Check for adjacent nodes two levels finer

      bool refine_node = false;

      // X faces

      for (int iy=0; iy<ny; iy++) {
	refine_node = refine_node ||
	  (cousin(XP,   0,iy) && 
	   cousin(XP,   0,iy)->any_children() ) ||
	  (cousin(XM,nx-1,iy) && 
	   cousin(XM,nx-1,iy)->any_children() );
      }

      // Y faces

      for (int ix=0; ix<nx; ix++) {
	refine_node = refine_node ||
	  (cousin(YP,ix,   0) && 
	   cousin(YP,ix,   0)->any_children() ) ||
	  (cousin(YM,ix,ny-1) && 
	   cousin(YM,ix,ny-1)->any_children() );
      }

      if (refine_node) {

	refined_tree = true;

	create_children_();
	update_children_();

	for (int iy=0; iy<ny; iy++) {
	  for (int ix=0; ix<nx; ix++) {
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

      bool * refine_child = new bool [num_children_()];

      for (int i=0; i<num_children_(); i++) {
	refine_child[i] = false;
      }

      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {

	  int i = index_(ix,iy);

	  bool r = refine_child[i];

	  if (! child(ix,iy)) {

	    // XM-face neighbor

	    if (ix > 0 && child(ix-1,iy)) {
	      for (int ky=0; ky<ny; ky++) {
		r = r || child(ix-1,iy)->child(nx-1,ky);
	      }
	    } else if (ix == 0 && cousin(XM,nx-1,iy)) {
	      for (int ky=0; ky<ny; ky++) {
		r = r || cousin(XM,nx-1,iy)->child(nx-1,ky);
	      }
	    }

	    // XP-face neighbor

	    if (ix < nx-1 && child(ix+1,iy)) {
	      for (int ky=0; ky<ny; ky++) {
		r = r || child(ix+1,iy)->child(0,ky);
	      }
	    } else if (ix == nx-1 && cousin(XP,0,iy)) {
	      for (int ky=0; ky<ny; ky++) {
		r = r || cousin(XP,0,iy)->child(0,ky);
	      }
	    }

	    // YM-face neighbor

	    if (iy > 0 && child(ix,iy-1)) {
	      for (int kx=0; kx<nx; kx++) {
		r = r || child(ix,iy-1)->child(kx,ny-1);
	      }
	    } else if (iy == 0 && cousin(YM,ix,ny-1)) {
	      for (int kx=0; kx<nx; kx++) {
		r = r || cousin(YM,ix,ny-1)->child(kx,ny-1);
	      }
	    }

	    // YP-face neighbor

	    if (iy < ny-1 && child(ix,iy+1)) {
	      for (int kx=0; kx<nx; kx++) {
		r = r || child(ix,iy+1)->child(kx,0);
	      }
	    } else if (iy == ny-1 && cousin(YP,ix,0)) {
	      for (int kx=0; kx<nx; kx++) {
		r = r || cousin(YP,ix,0)->child(kx,0);
	      }
	    }

	  } // if ! child

	  refine_child[i] = r;

	} // ix
      } // iy

      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  if (refine_child[index_(ix,iy)]) {
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

  for (int iy=0; iy<ny; iy++) {
    for (int ix=0; ix<nx; ix++) {
      if (child(ix,iy)) {
	child(ix,iy)->balance_pass(refined_tree,full_nodes);
      }
    }
  }
}


//----------------------------------------------------------------------

// Perform a pass of trying to optimize uniformly-refined nodes

void Node2K::optimize_pass(bool & refined_tree)
{

  bool single_children = true;
  bool same_level = true;

  int nx = k_;
  int ny = k_;

  for (int iy=0; iy<ny; iy++) {
    for (int ix=0; ix<nx; ix++) {
      single_children = single_children && 
	child(ix,iy) && (! child(ix,iy)->any_children());
    }
  }

  if (single_children) {

    // assert: all children exist

    int level = child(0,0)->level_adjust_;
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	same_level = same_level &&
	  (child(ix,iy)->level_adjust_ == level);
      }
    }
  }

  if (single_children && same_level) {

    // adjust effective resolution

    int increment = level_increment_();

    level_adjust_ += increment + child(0,0)->level_adjust_; 

    deallocate_children_();

    refined_tree = true;

  } else {
    
    // recurse

    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	if (child(ix,iy)) {
	  child(ix,iy)->optimize_pass(refined_tree);
	}
      }
    }
  }

}

//----------------------------------------------------------------------

// Fill the image region with values

void Node2K::fill_image
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

  int * ixk = new int [k_+1];
  int * iyk = new int [k_+1];
  for (int i=0; i<=k_; i++) {
    ixk[i] = lowx + i*dx;
    iyk[i] = lowy + i*dy;
  }

  for (int ix=0; ix<k_; ix++) {
    for (int iy=0; iy<k_; iy++) {
      if (child(ix,iy)) {
	child(ix,iy)->fill_image 
	  (image,ndx,ndy,ixk[ix],ixk[ix+1], iyk[iy],iyk[iy+1], level + 1, num_levels,line_width);
      }
    }
  }
}

//----------------------------------------------------------------------

// Fill the image region with values

void Node2K::geomview
(
 FILE * fpr,
 double lowx, double upx,  
 double lowy, double upy,
 double lowz, double upz,
 bool full )
{

  if (full) {
    fprintf (fpr,"VECT\n");
    fprintf (fpr,"6 18 2\n");
    fprintf (fpr,"1 1 8 3 3 2\n");
    fprintf (fpr,"1 0 1 0 0 0\n");

    // Print points at domain boundaries to provide geomview with bounding box

    fprintf (fpr,"0 0 0\n");
    fprintf (fpr,"1 1 1\n");

  }
  fprintf (fpr,"%g %g %g\n",lowx,lowy,lowz);
  fprintf (fpr,"%g %g %g\n",upx,lowy,lowz);
  fprintf (fpr,"%g %g %g\n",upx,upy,lowz);
  fprintf (fpr,"%g %g %g\n",lowx,upy,lowz);
  fprintf (fpr,"%g %g %g\n",lowx,lowy,lowz);
  fprintf (fpr,"%g %g %g\n",lowx,lowy,upz);
  fprintf (fpr,"%g %g %g\n",upx,lowy,upz);
  fprintf (fpr,"%g %g %g\n",upx,lowy,lowz);
  fprintf (fpr,"%g %g %g\n",lowx,upy,lowz);
  fprintf (fpr,"%g %g %g\n",lowx,upy,upz);
  fprintf (fpr,"%g %g %g\n",lowx,lowy,upz);
  fprintf (fpr,"%g %g %g\n",upx,upy,lowz);
  fprintf (fpr,"%g %g %g\n",upx,upy,upz);
  fprintf (fpr,"%g %g %g\n",upx,lowy,upz);
  fprintf (fpr,"%g %g %g\n",lowx,upy,upz);
  fprintf (fpr,"%g %g %g\n",upx,upy,upz);

  if (full) {
    fprintf (fpr,"1 1 1 1\n");
    fprintf (fpr,"1 1 1 0\n");
  }

  double * xk = new double [k_+1];
  double * yk = new double [k_+1];
  double hx = (upx-lowx) / k_;
  double hy = (upy-lowy) / k_;

  for (int i=0; i<k_+1; i++) {
    xk[i] = lowx + hx*i;
    yk[i] = lowy + hy*i;
  }
  for (int iy=0; iy<k_; iy++) {
    for (int ix=0; ix<k_; ix++) {
      if (child(ix,iy)) {
	child(ix,iy)->geomview
	  (fpr,
	   xk[ix],xk[ix+1], 
	   yk[iy],yk[iy+1],
	   lowz+0.25,upz+0.25,false);
      }
    }
  }
}

//----------------------------------------------------------------------

void Node2K::allocate_neighbors_ ()
{
  neighbor_ = new Node2K * [num_faces_()];
  for (int i=0; i<num_faces_(); i++) {
    neighbor_[i] = NULL;
  }
}

//----------------------------------------------------------------------

void Node2K::deallocate_neighbors_ ()
{
  delete [] neighbor_;
  neighbor_ = NULL;
}

//----------------------------------------------------------------------

void Node2K::allocate_children_ ()
{
  child_    = new Node2K * [num_children_()];
  for (int i=0; i<num_children_(); i++) {
    child_[i] = NULL;
  }
}

//----------------------------------------------------------------------

void Node2K::deallocate_children_ ()
{
  if (child_) {
    for (int i=0; i<num_children_(); i++) {
      delete child_[i];
    }

    delete [] child_;
    child_ = NULL;
  }
}



int Node2K::num_nodes_ = 0;

