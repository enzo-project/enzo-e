// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_Layout.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Implementation of the Layout class

//----------------------------------------------------------------------

#include "error.hpp"
#include "parallel.hpp"

//----------------------------------------------------------------------

Layout::Layout() throw()
{
  for (int i=0; i<3; i++) {
    np_[i] = 1;
    nt_[i] = 1;
    nd_[i] = 1;
  }
}
    
//----------------------------------------------------------------------

void Layout::set_process_blocks(int p0, int p1, int p2) throw()
{
  np_[0] = p0;
  np_[1] = p1;
  np_[2] = p2; 
}


//----------------------------------------------------------------------

void Layout::set_thread_blocks(int t0, int t1, int t2) throw()
{ 
  nt_[0] = t0;
  nt_[1] = t1;
  nt_[2] = t2; 
}

//----------------------------------------------------------------------

void Layout::set_data_blocks(int d0, int d1, int d2) throw()
{ 
  nd_[0] = d0;
  nd_[1] = d1;
  nd_[2] = d2; 
}

//----------------------------------------------------------------------

int Layout::process_count () throw()
{ 
  return np_[0]*np_[1]*np_[2]; 
}

//----------------------------------------------------------------------

int Layout::thread_count () throw()
{
  return  nt_[0]*nt_[1]*nt_[2]; 
}

//----------------------------------------------------------------------

int Layout::data_blocks_per_process () throw()
{ 
  return thread_count() * data_blocks_per_thread(); 
}

//----------------------------------------------------------------------

int Layout::data_blocks_per_thread () throw()
{
  return nd_[0]*nd_[1]* nd_[2]; 
}

//----------------------------------------------------------------------

bool Layout::neighbor_is_internal (int ip, int it, int id,
				   axis_type axis, int face)
{
  int k = neighbor_project_(ip,it,id,axis,face);
  return (0 <= k && k < nd_[axis] * nt_[axis] * np_[axis]);
}

//----------------------------------------------------------------------

int Layout::neighbor_process (int ip, int it, int id,
			      axis_type axis, int face)  throw()
{

  // Project indices along axis of interest and return updated block index
  int i = neighbor_project_(ip,it,id,axis,face);

  // convert blocks to process blocks, and return relative change
  return index_1_to_z(i,nd_[axis],nt_[axis],np_[axis]) - ip;

}

//----------------------------------------------------------------------

int Layout::neighbor_thread (int ip, int it, int id,
			     axis_type axis, int face) throw()
{

  // Project indices along axis of interest and return updated block index
  int i = neighbor_project_(ip,it,id,axis,face);

  // convert blocks to thread block, and return relative change
  return index_1_to_y(i,nd_[axis],nt_[axis],np_[axis]) - it;

}

//----------------------------------------------------------------------

void Layout::box_extent (int ip, int it, int id,
			 double lower_extent[3],    
			 double upper_extent[3])
{

  int block_index[3];
  block_indices_(ip,it,id,block_index);

  lower_extent[0] = 1.0 * block_index[0] / (nd_[0]*nt_[0]*np_[0]);
  lower_extent[1] = 1.0 * block_index[1] / (nd_[1]*nt_[1]*np_[1]);
  lower_extent[2] = 1.0 * block_index[2] / (nd_[2]*nt_[2]*np_[2]);

  upper_extent[0] = 1.0 * (block_index[0] + 1) / (nd_[0]*nt_[0]*np_[0]);
  upper_extent[1] = 1.0 * (block_index[1] + 1) / (nd_[1]*nt_[1]*np_[1]);
  upper_extent[2] = 1.0 * (block_index[2] + 1) / (nd_[2]*nt_[2]*np_[2]);

}

//----------------------------------------------------------------------

void Layout::array_indices (int ip, int it, int id,
			    int nx, int ny, int nz,
			    int lower_index[3],
			    int upper_index[3]) throw ()
{

  int block_index[3];
  block_indices_(ip,it,id,block_index);

  lower_index[0] = block_index[0] * nx /(nd_[0]*nt_[0]*np_[0]);
  lower_index[1] = block_index[1] * ny /(nd_[1]*nt_[1]*np_[1]);
  lower_index[2] = block_index[2] * nz /(nd_[2]*nt_[2]*np_[2]);

  upper_index[0] = (block_index[0]+1) * nx / (nd_[0]*nt_[0]*np_[0]);
  upper_index[1] = (block_index[1]+1) * ny / (nd_[1]*nt_[1]*np_[1]);
  upper_index[2] = (block_index[2]+1) * nz / (nd_[2]*nt_[2]*np_[2]);

}

//----------------------------------------------------------------------

void Layout::set_periodic (axis_type axis, bool periodic)
{
  periodic_[axis] = periodic; 
}

//----------------------------------------------------------------------

bool Layout::is_periodic (axis_type axis)
{
  return periodic_[axis]; 
}

//----------------------------------------------------------------------

int Layout::neighbor_project_(int ip, int it, int id, axis_type axis, int face)
{
  int ipa,ita,ida;
  switch (axis) {
  case axis_x:
    ipa = index_1_to_x(ip,np_[0],np_[1],np_[2]);
    ita = index_1_to_x(it,nt_[0],nt_[1],nt_[2]);
    ida = index_1_to_x(id,nd_[0],nd_[1],nd_[2]);
    break;
  case axis_y:
    ipa = index_1_to_y(ip,np_[0],np_[1],np_[2]);
    ita = index_1_to_y(it,nt_[0],nt_[1],nt_[2]);
    ida = index_1_to_y(id,nd_[0],nd_[1],nd_[2]);
    break;
  case axis_z:
    ipa = index_1_to_z(ip,np_[0],np_[1],np_[2]);
    ita = index_1_to_z(it,nt_[0],nt_[1],nt_[2]);
    ida = index_1_to_z(id,nd_[0],nd_[1],nd_[2]);
    break;
  default:
    break;
  }

  // Move "face" blocks (could be < 0) along axis
  ida += face;

  // Convert to total data blocks
  int i = index_3_to_1(ida,ita,ipa,nd_[axis],nt_[axis],np_[axis]);

  // Wrap if periodic
  if (periodic_[axis]) {
    int n = nd_[axis]*nt_[axis]*np_[axis];
    i = (i + n) % n;
  }
  return i;
}

//----------------------------------------------------------------------

void Layout::block_indices_ 
(int ip, int it, int id, int block_index[3]) throw ()
 
{
  // Project block id's along each axis
  int ip3[3],it3[3],id3[3];

  index_1_to_3(ip,ip3[0],ip3[1],ip3[2],np_[0],np_[1],np_[2]);
  index_1_to_3(it,it3[0],it3[1],it3[2],nt_[0],nt_[1],nt_[2]);
  index_1_to_3(id,id3[0],id3[1],id3[2],nd_[0],nd_[1],nd_[2]);

  // Convert from (ip,it,id) relative blocks to absolute data blocks

  block_index[0] = index_3_to_1(id3[0],it3[0],ip3[0],nd_[0],nt_[0],np_[0]);
  block_index[1] = index_3_to_1(id3[1],it3[1],ip3[1],nd_[1],nt_[1],np_[1]);
  block_index[2] = index_3_to_1(id3[2],it3[2],ip3[2],nd_[2],nt_[2],np_[2]);

}
