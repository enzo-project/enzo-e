// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Node3K.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-29
/// @brief    Implementation of the Node3K class

#include "mesh.hpp"

//----------------------------------------------------------------------

Node3K::Node3K(int k, int level_adjust) 
  : k_(k),
    child_(0),
    neighbor_(0),
    parent_(0),
    level_adjust_(level_adjust),
    data_(0)
    /// @param    k            refinement factor   
    /// @param    level_adjust difference: actual mesh level - tree level
{ 
  allocate_neighbors_();

  //  allocate_children_();

  parent_ = NULL;
}

//----------------------------------------------------------------------

Node3K::~Node3K()
///
{ 
  deallocate_children_();

  // update neighbor's neighbors

  for (int i=0; i<num_faces_(); i++) {
    int io = opposite_face_(face_enum(i));
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

Node3K::Node3K(const Node3K & node3k) throw()
{
  INCOMPLETE_MESSAGE("Node3K::Node3K","");
}

//----------------------------------------------------------------------

Node3K & Node3K::operator= (const Node3K & node3k) throw()
{
  INCOMPLETE_MESSAGE("Node3K::operator =","");
  return *this;
}

//----------------------------------------------------------------------

int Node3K::num_nodes()
{
  int node_count = 0;
  for (int iz=0; iz<k_; iz++) {
    for (int iy=0; iy<k_; iy++) {
      for (int ix=0; ix<k_; ix++) {
	if (child(ix,iy,iz)) {
	  node_count += child(ix,iy,iz)->num_nodes();
	}
      }
    }
  }

  return (1 + node_count);
}

//----------------------------------------------------------------------

inline Node3K * Node3K::child (int ix, int iy, int iz) 
/// @param    ix        Index 0 <= ix < k of cell in grid block
/// @param    iy        Index 0 <= iy < k of cell in grid block
/// @param    iz        Index 0 <= iz < k of cell in grid block
{ 
  return (child_ == NULL) ? NULL : child_[index_(ix,iy,iz)];
}

//----------------------------------------------------------------------

inline Node3K * Node3K::neighbor (face_enum face) 
/// @param    face      Face 0 <= (face = [XYZ][MP]) < 6
{ 
  return neighbor_[face]; 
}

//----------------------------------------------------------------------

inline Node3K * Node3K::cousin (face_enum face, int ix, int iy, int iz) 
/// @param    face      Face 0 <= (face = [XYZ][MP]) < 6
/// @param    ix        Index 0 <= ix < k of cell in grid block
/// @param    iy        Index 0 <= iy < k of cell in grid block
/// @param    iz        Index 0 <= iz < k of cell in grid block
{ 
  if (neighbor_[face] && neighbor_[face]->child(ix,iy,iz)) {
    return neighbor_[face]->child(ix,iy,iz);
  } else {
    return NULL;
  }
}

//----------------------------------------------------------------------

inline Node3K * Node3K::parent () 
///
{ 
  return parent_; 
}

inline void make_neighbors 
(
 Node3K * node_1, 
 Node3K * node_2, 
 face_enum face_1
 )
/// @param    node_1    First neighbor node pointer 
/// @param    node_2    Second neighbor node pointer
/// @param    face_1    Face 0 <= face_1 < 6 of node_1 that is adjacent to node_2
{
  if (node_1 != NULL) node_1->neighbor_[face_1] = node_2;
  if (node_2 != NULL) {
    face_enum face_2 = face_enum(node_2->opposite_face_(face_1));
    node_2->neighbor_[face_2] = node_1;
  }
}

//----------------------------------------------------------------------

int Node3K::refine 
(
 const int * level_array, 
 int ndx,  int ndy, int ndz,
 int nxm, int nxp,  
 int nym, int nyp,
 int nzm, int nzp,
 int level, 
 int max_level,
 bool full_nodes
 )
/// @param    level_array Array of levels to refine to
/// @param    ndx       x-dimension of level_array[]
/// @param    ndy       y-dimension of level_array[]
/// @param    ndz       z-dimension of level_array[]
/// @param    nxm       Lowest x-index of array for this node
/// @param    nxp       Upper bound on x-index of array for this node
/// @param    nym       Lowest y-index of array for this node
/// @param    nyp       Upper bound on y-index of array for this node
/// @param    nzm       Lowest z-index of array for this node
/// @param    nzp       Upper bound on z-index of array for this node
/// @param    level     Level of this node
/// @param    max_level Maximum refinement level
/// @param    full_nodes Whether nodes always have a full complement of children
{

  int depth = 0;
  int increment = level_increment_();

  if ( level < max_level && 
       nxm < nxp-1 && 
       nym < nyp-1 &&
       nzm < nzp-1  ) {

    // determine whether to refine the node

    int dx = (nxp - nxm)/k_;
    int dy = (nyp - nym)/k_;
    int dz = (nzp - nzm)/k_;

    int * ixk = new int [k_+1];
    int * iyk = new int [k_+1];
    int * izk = new int [k_+1];

    for (int i=0; i<=k_; i++) {
      ixk[i] = nxm + i*dx;
      iyk[i] = nym + i*dy;
      izk[i] = nzm + i*dz;
    }

    int * depth_child = new int [num_children_()];

    for (int i=0; i<num_children_(); i++) {
      depth_child[i] = 0;
    }

    if (full_nodes) {

      // Refine if any bits in the level_array are in this node

      bool refine_node = false;

      for (int iz=nzm; iz<nzp && !refine_node; iz++) {
	for (int iy=nym; iy<nyp && !refine_node; iy++) {
	  for (int ix=nxm; ix<nxp && !refine_node; ix++) {
	    if (level_array[ix + ndx*(iy + ndy*iz)] >= level) {
	      refine_node = true;
	    }
	  }
	}
      }

      // refine the node if needed

      if (refine_node) {

	create_children_();

	update_children_();

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

      for (int i=0; i<num_children_(); i++) {
	refine_child[i] = false;
      }
      

      // loop over children

      for (int iz=0; iz<k_; iz++) {
	for (int iy=0; iy<k_; iy++) {
	  for (int ix=0; ix<k_; ix++) {

	    int i = index_(ix,iy,iz);
	  
	    // loop over values in the level array

	    for (int lz=izk[iz]; lz<izk[iz+1] && ! refine_child[i]; lz++) {
	      for (int ly=iyk[iy]; ly<iyk[iy+1] && ! refine_child[i]; ly++) {
		for (int lx=ixk[ix]; lx<ixk[ix+1] && ! refine_child[i]; lx++) {

		  int l = lx + ndx * (ly + ndy * lz);
		  if (level_array[l] >= level) {
		    refine_child[i] = true;
		  }

		}
	      }
	    }

	  }
	}
      }

      // refine each child if needed

      for (int iz=0; iz<k_; iz++) {
	for (int iy=0; iy<k_; iy++) {
	  for (int ix=0; ix<k_; ix++) {

	    int i = index_(ix,iy,iz);
	  
	    if (refine_child[i]) {
	      create_child_(ix,iy,iz);
	      update_child_(ix,iy,iz);
	      depth_child[i] = child(ix,iy,iz)->refine 
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

      delete [] refine_child;
    }

    // determine depth as depth of deepest child + level increment

    for (int iz=0; iz<k_; iz++) {
      for (int iy=0; iy<k_; iy++) {
	for (int ix=0; ix<k_; ix++) {
	  int i = index_(ix,iy,iz);
	  depth = (depth_child[i] > depth) ? depth_child[i] : depth;
	}
      }
    }

    delete [] depth_child;

    depth += increment;

  } // if not at bottom of recursion

  return depth;
}

//----------------------------------------------------------------------

void Node3K::create_children_()
///
{
  for (int iz=0; iz<k_; iz++) {
    for (int iy=0; iy<k_; iy++) {
      for (int ix=0; ix<k_; ix++) {
	create_child_(ix,iy,iz);
      }
    }
  }
}

//----------------------------------------------------------------------

void Node3K::update_children_()
///
{
  for (int iz=0; iz<k_; iz++) {
    for (int iy=0; iy<k_; iy++) {
      for (int ix=0; ix<k_; ix++) {
	update_child_(ix,iy,iz);
      }
    }
  }
}

//----------------------------------------------------------------------

void Node3K::create_child_(int ix, int iy, int iz)
{
  if (child_ == NULL) allocate_children_();
  child_[index_(ix,iy,iz)] = new Node3K(k_);
}

//----------------------------------------------------------------------

void Node3K::update_child_ (int ix, int iy, int iz)
/// @param    ix        Index 0 <= ix < k of cell in grid block
/// @param    iy        Index 0 <= iy < k of cell in grid block
/// @param    iz        Index 0 <= iz < k of cell in grid block
{
  if (child(ix,iy,iz)) {

    child(ix,iy,iz)->parent_ = this;

    int nx = k_;
    int ny = k_;
    int nz = k_;

    // lower_x-face neighbors

    if (ix > 0) {
      make_neighbors (child (ix,iy,iz), child (ix-1,iy,iz),face_lower_x);
    } else {
      make_neighbors (child (ix,iy,iz), cousin (face_lower_x,nx-1,iy,iz),face_lower_x);
    }

    // upper_x-face neighbors

    if (ix < nx-1) {
      make_neighbors (child (ix,iy,iz), child (ix+1,iy,iz), face_upper_x);
    } else {
      make_neighbors (child (ix,iy,iz), cousin (face_upper_x,0,iy,iz), face_upper_x);
    }

    // lower_y-face neighbor

    if (iy > 0) {
      make_neighbors (child (ix,iy,iz), child (ix,iy-1,iz),face_lower_y);
    } else {
      make_neighbors (child (ix,iy,iz), cousin (face_lower_y,ix,ny-1,iz),face_lower_y);
    }

    // upper_y-face neighbor

    if (iy < ny-1) {
      make_neighbors (child (ix,iy,iz), child (ix,iy+1,iz),face_upper_y);
    } else {
      make_neighbors (child (ix,iy,iz), cousin (face_upper_y,ix,0,iz),face_upper_y);
    }

    // lower_z-face neighbor

    if (iz > 0) {
      make_neighbors (child (ix,iy,iz), child (ix,iy,iz-1),face_lower_z);
    } else {
      make_neighbors (child (ix,iy,iz), cousin (face_lower_z,ix,iy,nz-1),face_lower_z);
    }

    // upper_z-face neighbor

    if (iz < nz-1) {
      make_neighbors (child (ix,iy,iz), child (ix,iy,iz+1),face_upper_z);
    } else {
      make_neighbors (child (ix,iy,iz), cousin (face_upper_z,ix,iy,0),face_upper_z);
    }
  }
}

//----------------------------------------------------------------------

// Perform a pass of trying to remove level-jumps 

void Node3K::balance_pass(bool & refined_tree, bool full_nodes)
/// @param    refined_tree Whether tree has been refined
/// @param    full_nodes   Whether nodes always have a full complement of children
{
  int nx = k_;
  int ny = k_;
  int nz = k_;

  if (full_nodes) {

    // if is a leaf

    if (! any_children()) {

      // Check for adjacent nodes two levels finer

      bool refine_node = false;

      // X faces

      for (int iz=0; iz<nz; iz++) {
	for (int iy=0; iy<ny; iy++) {
	  refine_node = refine_node ||
	    (cousin(face_upper_x,   0,iy,iz) && 
	     cousin(face_upper_x,   0,iy,iz)->any_children() ) ||
	    (cousin(face_lower_x,nx-1,iy,iz) && 
	     cousin(face_lower_x,nx-1,iy,iz)->any_children() );
	}
      }

      // Y faces

      for (int iz=0; iz<nz; iz++) {
	for (int ix=0; ix<nx; ix++) {
	  refine_node = refine_node ||
	    (cousin(face_upper_y,ix,   0,iz) && 
	     cousin(face_upper_y,ix,   0,iz)->any_children() ) ||
	    (cousin(face_lower_y,ix,ny-1,iz) && 
	     cousin(face_lower_y,ix,ny-1,iz)->any_children() );
	}
      }

      // Z faces

      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  refine_node = refine_node ||
	    (cousin(face_upper_z,ix,iy,   0) && 
	     cousin(face_upper_z,ix,iy,   0)->any_children() ) ||
	    (cousin(face_lower_z,ix,iy,nz-1) && 
	     cousin(face_lower_z,ix,iy,nz-1)->any_children() );
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
	refine_child[i] = false;
      }

      for (int iz=0; iz<nz; iz++) {
	for (int iy=0; iy<ny; iy++) {
	  for (int ix=0; ix<nx; ix++) {

	    int i = index_(ix,iy,iz);

	    bool r = refine_child[i];

	    if (! child(ix,iy,iz)) {

	      // lower_x-face neighbor

	      if (ix > 0 && child(ix-1,iy,iz)) {
		for (int kz=0; kz<nz; kz++) {
		  for (int ky=0; ky<ny; ky++) {
		    r = r || child(ix-1,iy,iz)->child(nx-1,ky,kz);
		  }
		}
	      } else if (ix == 0 && cousin(face_lower_x,nx-1,iy,iz)) {
		for (int kz=0; kz<nz; kz++) {
		  for (int ky=0; ky<ny; ky++) {
		    r = r || cousin(face_lower_x,nx-1,iy,iz)->child(nx-1,ky,kz);
		  }
		}
	      }

	      // upper_x-face neighbor

	      if (ix < nx-1 && child(ix+1,iy,iz)) {
		for (int kz=0; kz<nz; kz++) {
		  for (int ky=0; ky<ny; ky++) {
		    r = r || child(ix+1,iy,iz)->child(0,ky,kz);
		  }
		}
	      } else if (ix == nx-1 && cousin(face_upper_x,0,iy,iz)) {
		for (int kz=0; kz<nz; kz++) {
		  for (int ky=0; ky<ny; ky++) {
		    r = r || cousin(face_upper_x,0,iy,iz)->child(0,ky,kz);
		  }
		}
	      }

	      // lower_y-face neighbor

	      if (iy > 0 && child(ix,iy-1,iz)) {
		for (int kz=0; kz<nz; kz++) {
		  for (int kx=0; kx<nx; kx++) {
		    r = r || child(ix,iy-1,iz)->child(kx,ny-1,kz);
		  }
		}
	      } else if (iy == 0 && cousin(face_lower_y,ix,ny-1,iz)) {
		for (int kz=0; kz<nz; kz++) {
		  for (int kx=0; kx<nx; kx++) {
		    r = r || cousin(face_lower_y,ix,ny-1,iz)->child(kx,ny-1,kz);
		  }
		}
	      }

	      // upper_y-face neighbor

	      if (iy < ny-1 && child(ix,iy+1,iz)) {
		for (int kz=0; kz<nz; kz++) {
		  for (int kx=0; kx<nx; kx++) {
		    r = r || child(ix,iy+1,iz)->child(kx,0,kz);
		  }
		}

	      } else if (iy == ny-1 && cousin(face_upper_y,ix,0,iz)) {
		for (int kz=0; kz<nz; kz++) {
		  for (int kx=0; kx<nx; kx++) {
		    r = r || cousin(face_upper_y,ix,0,iz)->child(kx,0,kz);
		  }
		}
	      }

	      // lower_z-face neighbor

	      if (iz > 0 && child(ix,iy,iz-1)) {
		for (int ky=0; ky<ny; ky++) {
		  for (int kx=0; kx<nx; kx++) {
		    r = r || child(ix,iy,iz-1)->child(kx,ky,nz-1);
		  }
		}
	      } else if (iz == 0 && cousin(face_lower_z,ix,iy,nz-1)) {
		for (int ky=0; ky<ny; ky++) {
		  for (int kx=0; kx<nx; kx++) {
		    r = r || cousin(face_lower_z,ix,iy,nz-1)->child(kx,ky,nz-1);
		  }
		}
	      }

	      // upper_z-face neighbor

	      if (iz < nz-1 && child(ix,iy,iz+1)) {
		for (int ky=0; ky<ny; ky++) {
		  for (int kx=0; kx<nx; kx++) {
		    r = r || child(ix,iy,iz+1)->child(kx,ky,0);
		  }
		}
	      } else if (iz == nz-1 && cousin(face_upper_z,ix,iy,0)) {
		for (int ky=0; ky<ny; ky++) {
		  for (int kx=0; kx<nx; kx++) {
		    r = r || cousin(face_upper_z,ix,iy,0)->child(kx,ky,0);
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

      delete [] refine_child;

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


//----------------------------------------------------------------------

// Perform a pass of trying to optimize uniformly-refined nodes

void Node3K::optimize_pass(bool & refined_tree)
/// @param    refined_tree Whether tree has been refined
{

  bool single_children = true;
  bool same_level = true;

  int nx = k_;
  int ny = k_;
  int nz = k_;

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

    int level = child(0,0,0)->level_adjust_;
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

    deallocate_children_();

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

//----------------------------------------------------------------------

// Fill the image region with values

void Node3K::fill_image
(
 float * image,
 int ndx,  int ndy, int ndz,
 int nxm, int nxp,  
 int nym, int nyp,
 int nzm, int nzp,
 int level,
 int num_levels,
 int line_width,
 int ia
 )
/// @param    image     Array of colormap indices
/// @param    ndx       x-dimension of image[]
/// @param    ndy       y-dimension of image[]
/// @param    ndz       z-dimension of image[]
/// @param    nxm       Lowest x-index of image[] for this node
/// @param    nxp       Upper bound on x-index of image[] for this node
/// @param    nym       Lowest y-index of image[] for this node
/// @param    nyp       Upper bound on y-index of image[] for this node
/// @param    nzm       Lowest z-index of image[] for this node
/// @param    nzp       Upper bound on z-index of image[] for this node
/// @param    level     Level of this node
/// @param    num_levels Total number of levels
/// @param    line_width Width of lines bounding nodes
/// @param    ia        Axis along which center slice is taken
{

  level += level_adjust_;

  int jx = (ia + 1) % 3;
  int jy = (ia + 2) % 3;

  int nm3[3] = {nxm,nym,nzm};
  int np3[3] = {nxp,nyp,nzp};
  int n3[3] = {ndx,ndy,ndz};
  // Exit if we're not in the middle slice orthogonal to axis ia
  if ( ! (nm3[ia] <= n3[ia]/2 && n3[ia]/2 <= np3[ia])) return;

  int n = n3[jx];

  // Fill interior

  int ix,iy;

  for (iy=nm3[jy]; iy<=np3[jy]; iy++) {
    for (ix=nm3[jx]; ix<=np3[jx]; ix++) {
      int i = ix + iy*n;
      image[i] = 2*num_levels - level; 
    }
  }

  // Draw border

  for (int k=0; k < line_width; k++) {

    for (ix=nm3[jx]; ix<=np3[jx]; ix++) {
      image[ix + (nm3[jy]+k)*n] = 0;
      image[ix + (np3[jy]+k)*n] = 0;
    }

    for (iy=nm3[jy]; iy<=np3[jy]; iy++) {
      image[(nm3[jx]+k) + iy*n] = 0;
      image[(np3[jx]+k) + iy*n] = 0;
    }

  }
  
  // Recurse
  
  int dx = (nxp - nxm)/k_;
  int dy = (nyp - nym)/k_;
  int dz = (nzp - nzm)/k_;

  int * ixk = new int [k_+1];
  int * iyk = new int [k_+1];
  int * izk = new int [k_+1];

  for (int i=0; i<=k_; i++) {
    ixk[i] = nxm + i*dx;
    iyk[i] = nym + i*dy;
    izk[i] = nzm + i*dz;
  }

  for (int iz=0; iz<k_; iz++) {
    for (int iy=0; iy<k_; iy++) {
      for (int ix=0; ix<k_; ix++) {
	if (child(ix,iy,iz)) {
	  child(ix,iy,iz)->fill_image 
	    (image,
	     ndx,ndy,ndz,
	     ixk[ix],ixk[ix+1], 
	     iyk[iy],iyk[iy+1],
	     izk[iz],izk[iz+1], 
	     level + 1, num_levels,line_width,ia);
	}
      }
    }
  }

  delete [] ixk;
  delete [] iyk;
  delete [] izk;
}


//----------------------------------------------------------------------

// Fill the image region with values

void Node3K::geomview
(
 FILE * fpr,
 double nxm, double nxp,  
 double nym, double nyp,
 double nzm, double nzp,
 bool full )
/// @param    fpr       Open file pointer for geomview file
/// @param    nxm       Lowest x-index of image[] for this node
/// @param    nxp       Upper bound on x-index of image[] for this node
/// @param    nym       Lowest y-index of image[] for this node
/// @param    nyp       Upper bound on y-index of image[] for this node
/// @param    nzm       Lowest z-index of image[] for this node
/// @param    nzp       Upper bound on z-index of image[] for this node
/// @param    full      Whether this is a standalone geomview node file or part of a tree
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
  fprintf (fpr,"%g %g %g\n",nxm,nym,nzm);
  fprintf (fpr,"%g %g %g\n",nxp,nym,nzm);
  fprintf (fpr,"%g %g %g\n",nxp,nyp,nzm);
  fprintf (fpr,"%g %g %g\n",nxm,nyp,nzm);
  fprintf (fpr,"%g %g %g\n",nxm,nym,nzm);
  fprintf (fpr,"%g %g %g\n",nxm,nym,nzp);
  fprintf (fpr,"%g %g %g\n",nxp,nym,nzp);
  fprintf (fpr,"%g %g %g\n",nxp,nym,nzm);
  fprintf (fpr,"%g %g %g\n",nxm,nyp,nzm);
  fprintf (fpr,"%g %g %g\n",nxm,nyp,nzp);
  fprintf (fpr,"%g %g %g\n",nxm,nym,nzp);
  fprintf (fpr,"%g %g %g\n",nxp,nyp,nzm);
  fprintf (fpr,"%g %g %g\n",nxp,nyp,nzp);
  fprintf (fpr,"%g %g %g\n",nxp,nym,nzp);
  fprintf (fpr,"%g %g %g\n",nxm,nyp,nzp);
  fprintf (fpr,"%g %g %g\n",nxp,nyp,nzp);

  if (full) {
    fprintf (fpr,"1 1 1 1\n");
    fprintf (fpr,"1 1 1 0\n");
  }

  double * xk = new double [k_+1];
  double * yk = new double [k_+1];
  double * zk = new double [k_+1];
  double hx = (nxp-nxm) / k_;
  double hy = (nyp-nym) / k_;
  double hz = (nzp-nzm) / k_;

  for (int i=0; i<k_+1; i++) {
    xk[i] = nxm + hx*i;
    yk[i] = nym + hy*i;
    zk[i] = nzm + hz*i;
  }
  for (int iz=0; iz<k_; iz++) {
    for (int iy=0; iy<k_; iy++) {
      for (int ix=0; ix<k_; ix++) {
	if (child(ix,iy,iz)) {
	  child(ix,iy,iz)->geomview
	    (fpr,
	     xk[ix],xk[ix+1], 
	     yk[iy],yk[iy+1],
	     zk[iz],zk[iz+1],false);
	}
      }
    }
  }
}

//----------------------------------------------------------------------

void Node3K::allocate_neighbors_ ()
///
{
  neighbor_ = new Node3K * [num_faces_()];
  for (int i=0; i<num_faces_(); i++) {
    neighbor_[i] = NULL;
  }
}

//----------------------------------------------------------------------

void Node3K::deallocate_neighbors_ ()
///
{
  delete [] neighbor_;
  neighbor_ = NULL;
}


//----------------------------------------------------------------------

void Node3K::allocate_children_ ()
///
{
  child_    = new Node3K * [num_children_()];
  for (int i=0; i<num_children_(); i++) {
    child_[i] = NULL;
  }
}


//----------------------------------------------------------------------

void Node3K::deallocate_children_ ()
///
{
  if (child_) {
    for (int i=0; i<num_children_(); i++) {
      delete child_[i];
    }

    delete [] child_;
    child_ = NULL;
  }
}
