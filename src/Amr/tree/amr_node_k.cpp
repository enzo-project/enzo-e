#include <stdio.h>
#include "cello.h"
#include "amr_node_k.hpp"
#include <assert.h>

Node_k::Node_k(int k, int level_adjust) 
  : k_(k),
    level_adjust_(level_adjust)

{ 
  ++Node_k::num_nodes_;

  neighbor_ = new Node_k * [num_faces_()];
  for (int i=0; i<num_faces_(); i++) {
    neighbor_[i] = NULL;
  }

  
  child_    = new Node_k * [num_children_()];
  for (int i=0; i<num_children_(); i++) {
      child_[i] = NULL;
  }

  parent_ = NULL;
}

// Delete the node and all descendents

Node_k::~Node_k() 
{ 
  --Node_k::num_nodes_;

  for (int i=0; i<num_children_(); i++) {
    delete child_[i];
  }

  delete [] child_;

  // update neighbor's neighbors

  for (int i=0; i<num_faces_(); i++) {
    int io = opposite_face_(face_type(i));
    if (neighbor_[i]) neighbor_[i]->neighbor_[io] = NULL;
    neighbor_[i] = NULL;
  }

  delete [] neighbor_;

  // Update parent's children

  if (parent_) {
    for (int i=0; i<num_children_(); i++) {
      if (parent_->child_[i] == this) parent_->child_[i] = NULL;
    }
  }

  parent_ = NULL;
}

inline Node_k * Node_k::child (int ix, int iy, int iz) 
{ 
  return child_[index_(ix,iy,iz)]; 
}

inline Node_k * Node_k::neighbor (face_type face) 
{ 
  return neighbor_[face]; 
}

inline Node_k * Node_k::cousin (face_type face, int ix, int iy, int iz) 
{ 
  if (neighbor_[face] && neighbor_[face]->child(ix,iy,iz)) {
    return neighbor_[face]->child(ix,iy,iz);
  } else {
    return NULL;
  }
}

inline Node_k * Node_k::parent () 
{ 
  return parent_; 
}

// Set two nodes to be neighbors.  Friend function since nodes may be NULL
inline void make_neighbors 
(
 Node_k * node_1, 
 Node_k * node_2, 
 face_type face_1
 )
{
  face_type face_2 = face_type(node_2->opposite_face_(face_1));
  if (node_1 != NULL) node_1->neighbor_[face_1] = node_2;
  if (node_2 != NULL) node_2->neighbor_[face_2] = node_1;
}


// Create 4 empty child nodes

int Node_k::refine 
(
 const int * level_array, 
 int ndx,  int ndy, int ndz,
 int lowx, int upx,  
 int lowy, int upy,
 int lowz, int upz,
 int level, 
 int max_level,
 bool full_nodes
 )
{
  int depth = 0;
  if ( level < max_level && 
       lowx < upx-1 && 
       lowy < upy-1 &&
       lowz < upz-1 
       ) {

    // determine whether to refine the node

    int dx = (upx - lowx)/k_;
    int dy = (upy - lowy)/k_;
    int dz = (upz - lowz)/k_;

    int * ixk = new int [k_+1];
    int * iyk = new int [k_+1];
    int * izk = new int [k_+1];

    for (int i=0; i<=k_; i++) {
      ixk[i] = lowx + i*dx;
      iyk[i] = lowy + i*dy;
      izk[i] = lowz + i*dz;
    }

    int * depth_child = new int [num_children_()];

    for (int i=0; i<num_children_(); i++) {
      depth_child[i] = 0;
    }

    if (full_nodes) {

      // Refine if any bits in the level_array are in this node

      bool refine_node = false;

      for (int iz=lowz; iz<upz && !refine_node; iz++) {
	for (int iy=lowy; iy<upy && !refine_node; iy++) {
	  for (int ix=lowx; ix<upx && !refine_node; ix++) {
	    if (level_array[index_(ix,iy,iz)] >= level) refine_node = true;
	  }
	}
      }

      // refine the node if needed

      if (refine_node) {

	create_children_();

	update_children_();

	int increment = level_increment_();

	for (int iz=0; iz<k_; iz++) {
	  for (int iy=0; iy<k_; iy++) {
	    for (int ix=0; ix<k_; ix++) {
	      depth_child[index_(ix,iy,iz)] = child(ix,iy,iz)->refine 
		(level_array,
		 ndx,ndy,ndz,
		 ixk[ix],ixk[ix+1],
		 iyk[iy],iyk[iy+1],
		 izk[iz],izk[iz+1],
		 level + increment,
		 max_level,full_nodes);
	    }
	  }
	}
      }
      
    } else {

      // Refine each child separately if any bits in the level_array
      // are in in them

      bool * refine_child = new bool [num_children_()];
      for (int i=0; i<num_children_(); i++) refine_child[i] = false;
      

      // loop over children

      for (int iz=0; iz<k_; iz++) {
	for (int iy=0; iy<k_; iy++) {
	  for (int ix=0; ix<k_; ix++) {

	    // loop over values in the level array

	    for (int lz=izk[iz]; lz<izk[iz+1]; lz++) {
	      for (int ly=iyk[iy]; ly<iyk[iy+1]; ly++) {
		for (int lx=ixk[ix]; lx<ixk[ix+1]; lx++) {

		  int l = lx + ndx * (ly + ndy * lz);
		  if (level_array[l] >= level) {
		    refine_child[index_(ix,iy,iz)] = true;
		  }

		}
	      }
	    }

	  }
	}
      }

      // refine each child if needed

      int increment = level_increment_();

      for (int iz=0; iz<k_; iz++) {
	for (int iy=0; iy<k_; iy++) {
	  for (int ix=0; ix<k_; ix++) {

	    if (refine_child[index_(ix,iy,iz)]) {
	      create_child_(ix,iy,iz);
	      update_child_(ix,iy,iz);
	      depth_child[index_(ix,iy,iz)] = child(ix,iy,iz)->refine 
		(level_array,
		 ndx,ndy,ndz,
		 ixk[ix],ixk[ix+1],
		 iyk[iy],iyk[iy+1],
		 izk[iz],izk[iz+1],
		 level + increment,
		 max_level,full_nodes);
	    }

	  }
	}
      }
    }

    // determine depth as depth of deepest child + level increment

    for (int iz=0; iz<k_; iz++) {
      for (int iy=0; iy<k_; iy++) {
	for (int ix=0; ix<k_; ix++) {
	  depth = (depth_child[index_(ix,iy,iz)] > depth) ? 
	    depth_child[index_(ix,iy,iz)] : depth;
	}
      }
    }

    depth += level_increment_();

  } // if not at bottom of recursion

  return depth;
}

void Node_k::create_children_()
{
  for (int iz=0; iz<k_; iz++) {
    for (int iy=0; iy<k_; iy++) {
      for (int ix=0; ix<k_; ix++) {
	create_child_(ix,iy,iz);
      }
    }
  }
}

void Node_k::update_children_()
{
  for (int iz=0; iz<k_; iz++) {
    for (int iy=0; iy<k_; iy++) {
      for (int ix=0; ix<k_; ix++) {
	update_child_(ix,iy,iz);
      }
    }
  }
}

void Node_k::create_child_(int ix, int iy, int iz)
{
  child_[index_(ix,iy,iz)] = new Node_k(k_);
}

void Node_k::update_child_ (int ix, int iy, int iz)
{
  if (child(ix,iy,iz)) {

    child(ix,iy,iz)->parent_ = this;

    int nx = k_;
    int ny = k_;
    int nz = d_ >= 3 ? k_ : 1;

    // XM-face neighbors

    if (ix > 0) {
      make_neighbors (child (ix,iy,iz), child (ix-1,iy,iz),XM);
    } else {
      make_neighbors (child (ix,iy,iz), cousin (XM,nx-1,iy,iz),XM);
    }

    // XP-face neighbors

    if (ix < nx-1) {
      make_neighbors (child (ix,iy,iz), child (ix+1,iy,iz), XP);
    } else {
      make_neighbors (child (ix,iy,iz), cousin (XP,0,iy,iz), XP);
    }

    // YM-face neighbor

    if (iy > 0) {
      make_neighbors (child (ix,iy,iz), child (ix,iy-1,iz),YM);
    } else {
      make_neighbors (child (ix,iy,iz), cousin (YM,ix,ny-1,iz),YM);
    }

    // YP-face neighbor

    if (iy < ny-1) {
      make_neighbors (child (ix,iy,iz), child (ix,iy+1,iz),YP);
    } else {
      make_neighbors (child (ix,iy,iz), cousin (YP,ix,0,iz),YP);
    }

    if (d_ > 2) {

      // ZM-face neighbor

      if (iz > 0) {
	make_neighbors (child (ix,iy,iz), child (ix,iy,iz-1),ZM);
      } else {
	make_neighbors (child (ix,iy,iz), cousin (ZM,ix,ny,iz-1),ZM);
      }

      // ZP-face neighbor

      if (iz < nz-1) {
	make_neighbors (child (ix,iy,iz), child (ix,iy,iz+1),ZP);
      } else {
	make_neighbors (child (ix,iy,iz), cousin (ZP,ix,iy,0),ZP);
      }

    }

  }
}

// Perform a pass of trying to remove level-jumps 
void Node_k::balance_pass(bool & refined_tree, bool full_nodes)
{
  int nx = k_;
  int ny = k_;
  int nz = d_ >= 3 ? k_ : 1;

  if (full_nodes) {

    // if is a leaf

    if (! any_children()) {

      // Check for adjacent nodes two levels finer

      bool refine_node = false;

      // X faces

      for (int iz=0; iz<nz; iz++) {
	for (int iy=0; iy<ny; iy++) {
	  refine_node = refine_node ||
	    (cousin(XP,0,iy,iz) && 
	     cousin(XP,0,iy,iz)->any_children() ) ||
	    (cousin(XM,nx-1,iy,iz) && 
	     cousin(XM,nx-1,iy,iz)->any_children() );
	}
      }

      // Y faces

      for (int iz=0; iz<nz; iz++) {
	for (int ix=0; ix<nx; ix++) {
	  refine_node = refine_node ||
	    (cousin(YP,ix,0,iz) && 
	     cousin(YP,ix,0,iz)->any_children() ) ||
	    (cousin(YM,ix,ny-1,iz) && 
	     cousin(YM,ix,ny-1,iz)->any_children() );
	}
      }

      // Z faces

      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  refine_node = refine_node ||
	    (cousin(ZP,ix,iy,0) && 
	     cousin(ZP,ix,iy,0)->any_children() ) ||
	    (cousin(ZM,ix,iy,nz-1) && 
	     cousin(ZM,ix,iy,nz-1)->any_children() );
	}
      }

      if (refine_node) {

	refined_tree = true;

	create_children_();
	update_children_();

	for (int iz=0; iz<nz; iz++) {
	  for (int iy=0; iy<ny; iy++) {
	    for (int ix=0; ix<nx; ix++) {
	      if (child(ix,iy,iz)) {
		child(ix,iy,iz)->balance_pass(refined_tree,full_nodes);
	      }
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
	refine_child[i] = 0;
      }

      for (int iz=0; iz<nz; iz++) {
	for (int iy=0; iy<ny; iy++) {
	  for (int ix=0; ix<nx; ix++) {

	    int i = index_(ix,iy,iz);

	    if (! child(ix,iy,iz)) {

	      // XM-face neighbor

	      if (ix > 0 && child(ix-1,iy,iz)) {
		for (int kz=0; kz<nz; kz++) {
		  for (int ky=0; ky<ny; ky++) {
		    refine_child[i] = refine_child[i] ||
		      child(ix-1,iy,iz)->child(nx-1,ky,kz);
		  }
		}
	      } else if (ix == 0 && cousin(XM,nx-1,iy,iz)) {
		for (int kz=0; kz<nz; kz++) {
		  for (int ky=0; ky<ny; ky++) {
		    refine_child[i] = refine_child[i] ||
		      cousin(XM,nx-1,iy,iz)->child(nx-1,ky,kz);
		  }
		}
	      }

	      // XP-face neighbor

	      if (ix < nx-1 && child(ix+1,iy,iz)) {
		for (int kz=0; kz<nz; kz++) {
		  for (int ky=0; ky<ny; ky++) {
		    int i = index_(ix,iy,iz);
		    refine_child[i] = refine_child[i] ||
		      child(ix+1,iy,iz)->child(0,ky,kz);
		  }
		}
	      } else if (ix == nx-1 && cousin(XP,0,iy,iz)) {
		for (int kz=0; kz<nz; kz++) {
		  for (int ky=0; ky<ny; ky++) {
		    int i = index_(ix,iy,iz);
		    refine_child[i] = refine_child[i] ||
		      cousin(XP,0,iy,iz)->child(0,ky,kz);
		  }
		}
	      }

	      // YM-face neighbor

	      if (iy > 0 && child(ix,iy-1,iz)) {
		for (int kz=0; kz<nz; kz++) {
		  for (int kx=0; kx<nx; kx++) {
		    refine_child[i] = refine_child[i] ||
		      child(ix,iy-1,iz)->child(kx,ny-1,kz);
		  }
		}
	      } else if (iy == 0 && cousin(YM,ix,ny-1,iz)) {
		for (int kz=0; kz<nz; kz++) {
		  for (int kx=0; kx<nx; kx++) {
		    refine_child[i] = refine_child[i] ||
		      cousin(XM,ix,ny-1,iz)->child(kx,ny-1,kz);
		  }
		}
	      }

	      // YP-face neighbor

	      if (iy < ny-1 && child(ix,iy+1,iz)) {
		for (int kz=0; kz<nz; kz++) {
		  for (int kx=0; kx<nx; kx++) {
		    refine_child[i] = refine_child[i] ||
		      child(ix,iy+1,iz)->child(kx,0,kz);
		  }
		}

	      } else if (iy == ny-1 && cousin(YP,ix,0,iz)) {
		for (int kz=0; kz<nz; kz++) {
		  for (int kx=0; kx<nx; kx++) {
		    refine_child[i] = refine_child[i] ||
		      cousin(YP,ix,0,iz)->child(kx,0,kz);
		  }
		}
	      }

	      if (d_ >= 3) {

		// ZM-face neighbor

		if (iz > 0 && child(ix,iy,iz-1)) {
		  for (int ky=0; ky<ny; ky++) {
		    for (int kx=0; kx<nx; kx++) {
		      refine_child[i] = refine_child[i] ||
			child(ix,iy,iz-1)->child(kx,ky,nz-1);
		    }
		  }
		} else if (iz == 0 && cousin(ZM,ix,iy,nz-1)) {
		  for (int ky=0; ky<ny; ky++) {
		    for (int kx=0; kx<nx; kx++) {
		      refine_child[i] = refine_child[i] ||
			cousin(ZM,ix,iy,nz-1)->child(kx,ky,nz-1);
		    }
		  }
		}

		// ZP-face neighbor

		if (iz < nz-1 && child(ix,iy,iz+1)) {
		  for (int ky=0; ky<ny; ky++) {
		    for (int kx=0; kx<nx; kx++) {
		      refine_child[i] = refine_child[i] ||
			child(ix,iy,iz+1)->child(kx,ky,0);
		    }
		  }

		} else if (iz == nz-1 && cousin(ZP,ix,iy,0)) {
		  for (int ky=0; ky<ny; ky++) {
		    for (int kx=0; kx<nx; kx++) {
		      refine_child[i] = refine_child[i] ||
			cousin(ZP,ix,iy,0)->child(kx,ky,0);
		    }
		  }
		}

	      }
	    } // if ! child

	  } // ix
	} // iy
      } // iz

      for (int iz=0; iz<nz; iz++) {
	for (int iy=0; iy<ny; iy++) {
	  for (int ix=0; ix<nx; ix++) {
	    if (refine_child[index_(ix,iy,iz)]) {
	      create_child_(ix,iy,iz); 
	      update_child_(ix,iy,iz);
	      child(ix,iy,iz)->balance_pass(refined_tree,full_nodes);
	      refined_tree = true;
	    }
	  }
	}
      }
    }
  }

  // not a leaf: recurse

  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	if (child(ix,iy,iz)) {
	  child(ix,iy,iz)->balance_pass(refined_tree,full_nodes);
	}
      }
    }
  }
}


// Perform a pass of trying to optimize uniformly-refined nodes
void Node_k::optimize_pass(bool & refined_tree)
{

  bool single_children = true;
  bool same_level = true;

  int nx = k_;
  int ny = k_;
  int nz = d_ >= 3 ? k_ : 1;

  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	single_children = single_children && 
	  child(ix,iy,iz) && (! child(ix,iy,iz)->any_children());
      }
    }
  }

  if (single_children) {

    // assert: all children exist

    int level = child(0,0)->level_adjust_;
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  same_level = same_level &&
	    (child(ix,iy,iz)->level_adjust_ == level);
	}
      }
    }
  }

  if (single_children && same_level) {

    // adjust effective resolution

    int increment = level_increment_();

    level_adjust_ += increment + child(0,0,0)->level_adjust_; 

    for (int i=0; i<num_children_(); i++) {
      delete child_[i];
      child_[i] = NULL;
    }

    refined_tree = true;

  } else {
    
    // recurse

    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  if (child(ix,iy,iz)) {
	    child(ix,iy,iz)->optimize_pass(refined_tree);
	  }
	}
      }
    }
  }

}

// Fill the image region with values
void Node_k::fill_image
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

int Node_k::num_nodes_ = 0;

