// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Node2K.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-27
/// @brief    Implementation of the Node2K class 

#include "mesh_tree.hpp"

//----------------------------------------------------------------------

Node2K::Node2K(int k, int level_adjust) 
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

//   allocate_children_();

  parent_ = NULL;
}

//----------------------------------------------------------------------

Node2K::~Node2K() 
///
{ 
  deallocate_children_();

  // update neighbor's neighbors

  for (int i=0; i<num_faces_(); i++) {
    int io = opposite_face_(face_axis_enum(i));
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

Node2K::Node2K(const Node2K & node2k) throw()
{
  INCOMPLETE("Node2K::Node2K");
}

//----------------------------------------------------------------------

Node2K & Node2K::operator= (const Node2K & node2k) throw()
{
  INCOMPLETE("Node2K::operator =");
  return *this;
}

//----------------------------------------------------------------------

int Node2K::num_nodes()
{
  int node_count = 0;
  for (int iy=0; iy<k_; iy++) {
    for (int ix=0; ix<k_; ix++) {
      if (child(ix,iy)) {
	node_count += child(ix,iy)->num_nodes();
      }
    }
  }

  return (1 + node_count);
}

//----------------------------------------------------------------------

inline Node2K * Node2K::child (int ix, int iy) 
/// @param    ix        Index 0 <= ix < k of cell in grid block
/// @param    iy        Index 0 <= iy < k of cell in grid block
{
  return (child_ == NULL) ? NULL : child_[index_(ix,iy)];
}

//----------------------------------------------------------------------

inline Node2K * Node2K::neighbor (face_axis_enum face_axis) 
/// @param    face_axis      Face_Axis 0 <= (face_axis = [XYZ][MP]) < 6
{ 
  return neighbor_[face_axis]; 
}

//----------------------------------------------------------------------

inline Node2K * Node2K::cousin (face_axis_enum face_axis, int ix, int iy) 
/// @param    face_axis      Face_Axis 0 <= (face_axis = [XYZ][MP]) < 6
/// @param    ix        Index 0 <= ix < k of cell in grid block
/// @param    iy        Index 0 <= iy < k of cell in grid block
{ 
  if (neighbor_[face_axis] && neighbor_[face_axis]->child(ix,iy)) {
    return neighbor_[face_axis]->child(ix,iy);
  } else {
    return NULL;
  }
}

//----------------------------------------------------------------------

inline Node2K * Node2K::parent () 
///
{ 
  return parent_; 
}

inline void make_neighbors 
(
 Node2K * node_1, 
 Node2K * node_2, 
 face_axis_enum face_axis_1
 )
/// @param    node_1    First neighbor node pointer 
/// @param    node_2    Second neighbor node pointer
/// @param    face_axis_1    Face_Axis 0 <= face_axis_1 < 4 of node_1 that is adjacent to node_2
{
  if (node_1 != NULL) node_1->neighbor_[face_axis_1] = node_2;
  if (node_2 != NULL) {
    face_axis_enum face_axis_2 = face_axis_enum(node_2->opposite_face_(face_axis_1));
    node_2->neighbor_[face_axis_2] = node_1;
  }
}

//----------------------------------------------------------------------

int Node2K::refine 
(
 const int * level_array, 
 int ndx,  int ndy,
 int nxm, int nxp,  
 int nym, int nyp,
 int level, 
 int max_level,
 bool full_nodes
 )
/// @param    level_array Array of levels to refine to
/// @param    ndx       x-dimension of level_array[]
/// @param    ndy       y-dimension of level_array[]
/// @param    nxm       Lowest x-index of array for this node
/// @param    nxp       Upper bound on x-index of array for this node
/// @param    nym       Lowest y-index of array for this node
/// @param    nyp       Upper bound on y-index of array for this node
/// @param    level     Level of this node
/// @param    max_level Maximum refinement level
/// @param    full_nodes Whether nodes always have a full complement of children
{

  int depth = 0;
  int increment = level_increment_();

  if ( level < max_level && 
       nxm < nxp-1 && 
       nym < nyp-1 ) {

    // determine whether to refine the node

    int dx = (nxp - nxm)/k_;
    int dy = (nyp - nym)/k_;

    int * ixk = new int [k_+1];
    int * iyk = new int [k_+1];

    for (int i=0; i<=k_; i++) {
      ixk[i] = nxm + i*dx;
      iyk[i] = nym + i*dy;
    }

    int * depth_child = new int [num_children_()];

    for (int i=0; i<num_children_(); i++) {
      depth_child[i] = 0;
    }

    if (full_nodes) {

      // Refine if any bits in the level_array are in this node

      bool refine_node = false;

      for (int iy=nym; iy<nyp && !refine_node; iy++) {
	for (int ix=nxm; ix<nxp && !refine_node; ix++) {
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
    delete [] iyk;
    delete [] ixk;

    depth += increment;

  } // if not at bottom of recursion

  return depth;
}

//----------------------------------------------------------------------

void Node2K::create_children_()
///
{
  for (int iy=0; iy<k_; iy++) {
    for (int ix=0; ix<k_; ix++) {
      create_child_(ix,iy);
    }
  }
}

//----------------------------------------------------------------------

void Node2K::update_children_()
///
{
  for (int iy=0; iy<k_; iy++) {
    for (int ix=0; ix<k_; ix++) {
      update_child_(ix,iy);
    }
  }
}

//----------------------------------------------------------------------

void Node2K::create_child_(int ix, int iy)
/// @param    ix        Index 0 <= ix < k of cell in grid block
/// @param    iy        Index 0 <= iy < k of cell in grid block
{
  if (child_ == NULL) allocate_children_();
  child_[index_(ix,iy)] = new Node2K(k_);
}

//----------------------------------------------------------------------

void Node2K::update_child_ (int ix, int iy)
/// @param    ix        Index 0 <= ix < k of cell in grid block
/// @param    iy        Index 0 <= iy < k of cell in grid block
{
  if (child(ix,iy)) {

    child(ix,iy)->parent_ = this;

    int nx = k_;
    int ny = k_;

    // XM-face neighbors

    if (ix > 0) {
      make_neighbors (child (ix,iy), child (ix-1,iy),face_lower_axis_x);
    } else {
      make_neighbors (child (ix,iy), cousin (face_lower_axis_x,nx-1,iy),face_lower_axis_x);
    }

    // upper_x-face neighbors

    if (ix < nx-1) {
      make_neighbors (child (ix,iy), child (ix+1,iy), face_upper_axis_x);
    } else {
      make_neighbors (child (ix,iy), cousin (face_upper_axis_x,0,iy), face_upper_axis_x);
    }

    // lower_y-face neighbor

    if (iy > 0) {
      make_neighbors (child (ix,iy), child (ix,iy-1),face_lower_axis_y);
    } else {
      make_neighbors (child (ix,iy), cousin (face_lower_axis_y,ix,ny-1),face_lower_axis_y);
    }

    // upper_y-face neighbor

    if (iy < ny-1) {
      make_neighbors (child (ix,iy), child (ix,iy+1),face_upper_axis_y);
    } else {
      make_neighbors (child (ix,iy), cousin (face_upper_axis_y,ix,0),face_upper_axis_y);
    }

  }
}

//----------------------------------------------------------------------

// Perform a pass of trying to remove level-jumps 

void Node2K::balance_pass(bool & refined_tree, bool full_nodes)
/// @param    refined_tree Whether tree has been refined
/// @param    full_nodes   Whether nodes always have a full complement of children
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
	  (cousin(face_upper_axis_x,   0,iy) && 
	   cousin(face_upper_axis_x,   0,iy)->any_children() ) ||
	  (cousin(face_lower_axis_x,nx-1,iy) && 
	   cousin(face_lower_axis_x,nx-1,iy)->any_children() );
      }

      // Y faces

      for (int ix=0; ix<nx; ix++) {
	refine_node = refine_node ||
	  (cousin(face_upper_axis_y,ix,   0) && 
	   cousin(face_upper_axis_y,ix,   0)->any_children() ) ||
	  (cousin(face_lower_axis_y,ix,ny-1) && 
	   cousin(face_lower_axis_y,ix,ny-1)->any_children() );
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

	    // lower_x-face neighbor

	    if (ix > 0 && child(ix-1,iy)) {
	      for (int ky=0; ky<ny; ky++) {
		r = r || child(ix-1,iy)->child(nx-1,ky);
	      }
	    } else if (ix == 0 && cousin(face_lower_axis_x,nx-1,iy)) {
	      for (int ky=0; ky<ny; ky++) {
		r = r || cousin(face_lower_axis_x,nx-1,iy)->child(nx-1,ky);
	      }
	    }

	    // upper_x-face neighbor

	    if (ix < nx-1 && child(ix+1,iy)) {
	      for (int ky=0; ky<ny; ky++) {
		r = r || child(ix+1,iy)->child(0,ky);
	      }
	    } else if (ix == nx-1 && cousin(face_upper_axis_x,0,iy)) {
	      for (int ky=0; ky<ny; ky++) {
		r = r || cousin(face_upper_axis_x,0,iy)->child(0,ky);
	      }
	    }

	    // lower_y-face neighbor

	    if (iy > 0 && child(ix,iy-1)) {
	      for (int kx=0; kx<nx; kx++) {
		r = r || child(ix,iy-1)->child(kx,ny-1);
	      }
	    } else if (iy == 0 && cousin(face_lower_axis_y,ix,ny-1)) {
	      for (int kx=0; kx<nx; kx++) {
		r = r || cousin(face_lower_axis_y,ix,ny-1)->child(kx,ny-1);
	      }
	    }

	    // upper_y-face neighbor

	    if (iy < ny-1 && child(ix,iy+1)) {
	      for (int kx=0; kx<nx; kx++) {
		r = r || child(ix,iy+1)->child(kx,0);
	      }
	    } else if (iy == ny-1 && cousin(face_upper_axis_y,ix,0)) {
	      for (int kx=0; kx<nx; kx++) {
		r = r || cousin(face_upper_axis_y,ix,0)->child(kx,0);
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
/// @param    refined_tree Whether tree has been refined
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
 int nxm, int nxp,  
 int nym, int nyp,
 int level,
 int num_levels,
 int line_width
 )
/// @param    image     Array of colormap indices
/// @param    ndx       x-dimension of image[]
/// @param    ndy       y-dimension of image[]
/// @param    nxm      Lowest x-index of image[] for this node
/// @param    nxp       Upper bound on x-index of image[] for this node
/// @param    nym      Lowest y-index of image[] for this node
/// @param    nyp       Upper bound on y-index of image[] for this node
/// @param    level     Level of this node
/// @param    num_levels Total number of levels
/// @param    line_width Width of lines bounding nodes
{
  int ix,iy,i;

  level += level_adjust_;

  // Fill interior

  for (iy=nym; iy<=nyp; iy++) {
    for (ix=nxm; ix<=nxp; ix++) {
      i = ix + ndx * iy;
      image[i] = 2*num_levels - level; 
    }
  }

  // Draw border
  for (ix=nxm; ix<=nxp; ix++) {
    iy = nym;
    for (int k=0; k < line_width; k++) {
      image[ix + ndx*(iy+k)] = 0;
    }
    iy = nyp;
    for (int k=0; k < line_width; k++) {
      image[ix + ndx*(iy+k)] = 0;
    }
  }

  for (iy=nym; iy<=nyp; iy++) {
    ix = nxm;
    for (int k=0; k < line_width; k++) {
      image[(ix+k) + ndx*iy] = 0;
    }
    ix = nxp;
    for (int k=0; k < line_width; k++) {
      image[(ix+k) + ndx*iy] = 0;
    }
  }

  // Recurse

  int dx = (nxp - nxm)/k_;
  int dy = (nyp - nym)/k_;

  int * ixk = new int [k_+1];
  int * iyk = new int [k_+1];
  for (int i=0; i<=k_; i++) {
    ixk[i] = nxm + i*dx;
    iyk[i] = nym + i*dy;
  }

  for (int ix=0; ix<k_; ix++) {
    for (int iy=0; iy<k_; iy++) {
      if (child(ix,iy)) {
	child(ix,iy)->fill_image 
	  (image,ndx,ndy,ixk[ix],ixk[ix+1], iyk[iy],iyk[iy+1], level + 1, num_levels,line_width);
      }
    }
  }
  delete [] ixk;
  delete [] iyk;
}

//----------------------------------------------------------------------

void Node2K::geomview
(
 FILE * fpr,
 double nxm, double nxp,  
 double nym, double nyp,
 double nzm, double nzp,
 bool full )
/// @param    fpr       File pointer of geomview file opened for output
/// @param    nxm       Lowest x-index of this node
/// @param    nxp       Upper bound on x-index of this node (warning: different meanings)
/// @param    nym       Lowest y-index of this node (warning: different meanings)
/// @param    nyp       Upper bound on y-index of this node (warning: different meanings)
/// @param    nzm       Lowest z-index of this node (warning: different meanings)
/// @param    nzp       Upper bound on z-index of this node (warning: different meanings)
/// @param    full      Whether a refined node is always fully refined
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
  double hx = (nxp-nxm) / k_;
  double hy = (nyp-nym) / k_;

  for (int i=0; i<k_+1; i++) {
    xk[i] = nxm + hx*i;
    yk[i] = nym + hy*i;
  }
  for (int iy=0; iy<k_; iy++) {
    for (int ix=0; ix<k_; ix++) {
      if (child(ix,iy)) {
	child(ix,iy)->geomview
	  (fpr,
	   xk[ix],xk[ix+1], 
	   yk[iy],yk[iy+1],
	   nzm+0.25,nzp+0.25,false);
      }
    }
  }
  delete [] xk;
  delete [] yk;
}

//----------------------------------------------------------------------

void Node2K::allocate_neighbors_ ()
///
{
  neighbor_ = new Node2K * [num_faces_()];
  for (int i=0; i<num_faces_(); i++) {
    neighbor_[i] = NULL;
  }
}

//----------------------------------------------------------------------

void Node2K::deallocate_neighbors_ ()
///
{
  delete [] neighbor_;
  neighbor_ = 0;
}

//----------------------------------------------------------------------

void Node2K::allocate_children_ ()
///
{
  child_    = new Node2K * [num_children_()];
  for (int i=0; i<num_children_(); i++) {
    child_[i] = NULL;
  }
}

//----------------------------------------------------------------------

void Node2K::deallocate_children_ ()
///
{
  if (child_) {
    for (int i=0; i<num_children_(); i++) {
      delete child_[i];
      child_[i] = 0;
    }

    delete [] child_;
    child_ = 0;
  }
}

